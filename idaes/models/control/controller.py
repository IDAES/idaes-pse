#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
PID controller model module
"""

__author__ = ["John Eslick", "Jinliang Ma"]

import enum
import pyomo.environ as pyo
import pyomo.dae as pyodae

from idaes.core import UnitModelBlockData, declare_process_block_class
from pyomo.common.config import ConfigValue, In
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_bound


class ControllerType(enum.Enum):
    """Controller types."""

    P = 1
    PI = 2
    PD = 3
    PID = 4


class ControllerMVBoundType(enum.Enum):
    """Manipulated value bound type.

    NONE: No bound on manipulated value output.
    SMOOTH_BOUND: Use a smoothed version of mv = min(max(mv_unbound, ub), lb)
    LOGISTIC: Use a logistic function too keep mv between the bounds
    """

    NONE = 1
    SMOOTH_BOUND = 2
    LOGISTIC = 3


@declare_process_block_class(
    "PIDController",
    doc="PID controller model block.  To use this the model must be dynamic.",
)
class PIDControllerData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "process_var",
        ConfigValue(
            default=None,
            description="Process variable to be controlled",
            doc=(
                "A Pyomo Var, Expression, or Reference for the measured"
                " process variable. Should be indexed by time. Slices are okay."
            ),
        ),
    )
    CONFIG.declare(
        "manipulated_var",
        ConfigValue(
            default=None,
            description="Manipulated variable",
            doc=(
                "A Pyomo Var, Reference to a Var, or something that can be"
                " used to construct a Reference to a Var for the controlled"
                " variable. The final Reference should be indexed by time."
                " Slices are okay."
            ),
        ),
    )
    CONFIG.declare(
        "mv_bound_type",
        ConfigValue(
            default=ControllerMVBoundType.NONE,
            domain=In(
                [
                    ControllerMVBoundType.NONE,
                    ControllerMVBoundType.SMOOTH_BOUND,
                    ControllerMVBoundType.LOGISTIC,
                ]
            ),
            description="Type of bounds to apply to the manipulated variable (mv)).",
            doc=(
                """Type of bounds to apply to the manipulated variable output. If,
bounds are applied, the model parameters **mv_lb** and **mv_ub** set the bounds.
The **default** is ControllerMVBoundType.NONE. See the controller documentation
for details on the mathematical formulation. The options are:
**ControllerMVBoundType.NONE** no bounds, **ControllerMVBoundType.SMOOTH_BOUND**
smoothed mv = min(max(mv_unbound, ub), lb), and **ControllerMVBoundType.LOGISTIC**
logistic function to enforce bounds.
"""
            ),
        ),
    )
    CONFIG.declare(
        "calculate_initial_integral",
        ConfigValue(
            default=True,
            domain=bool,
            description="Calculate the initial integral term value if True",
            doc="Calculate the initial integral term value if True",
        ),
    )
    CONFIG.declare(
        "type",
        ConfigValue(
            default=ControllerType.PI,
            domain=In(
                [
                    ControllerType.P,
                    ControllerType.PI,
                    ControllerType.PD,
                    ControllerType.PID,
                ]
            ),
            description="Control type",
            doc="""Controller type. The **default** = ControllerType.PI and the
options are: **ControllerType.P** Proportional, **ControllerType.PI**
proportional and integral, **ControllerType.PD** proportional and derivative, and
**ControllerType.PID** proportional, integral, and derivative

""",
        ),
    )

    def build(self):
        """
        Build the PID block
        """
        super().build()

        # Check for required config
        if self.config.process_var is None or self.config.manipulated_var is None:
            raise ConfigurationError("Controller config requires 'pv' and 'mv'")
        self.process_var = pyo.Reference(self.config.process_var)
        self.manipulated_var = pyo.Reference(self.config.manipulated_var)

        # Shorter pointers to time set information
        time_set = self.flowsheet().time
        time_units = self.flowsheet().time_units
        time_0 = time_set.first()

        # Type Check
        if not issubclass(self.process_var[time_0].ctype, (pyo.Var, pyo.Expression)):
            raise TypeError(
                f"process_var must reference a Var or Expression not {self.process_var[time_0].ctype}"
            )
        if not issubclass(self.process_var[time_0].ctype, pyo.Var):
            raise TypeError(
                f"manipulated_var must reference a Var not {self.process_var[time_0].ctype}"
            )

        # Get the appropriate units for various contoller varaibles
        mv_units = pyo.units.get_units(self.manipulated_var[time_0])
        pv_units = pyo.units.get_units(self.process_var[time_0])
        if mv_units is None:
            mv_units = pyo.units.dimensionless
        if pv_units is None:
            pv_units = pyo.units.dimensionless
        gain_p_units = mv_units / pv_units

        # Parameters
        self.mv_lb = pyo.Param(
            mutable=True,
            initialize=0.05,
            doc="Controller output lower bound",
            units=mv_units,
        )
        self.mv_ub = pyo.Param(
            mutable=True,
            initialize=1,
            doc="Controller output upper bound",
            units=mv_units,
        )
        self.smooth_eps = pyo.Param(
            mutable=True,
            initialize=1e-4,
            doc="Smoothing parameter for controller output limits when the bound"
            " type is SMOOTH_BOUND",
            units=mv_units,
        )
        self.logistic_bound_k = pyo.Param(
            mutable=True,
            initialize=4,
            doc="Smoothing parameter for controller output limits when the bound"
            " type is LOGISTIC",
            units=mv_units,
        )

        # Variable for basic controller settings may change with time.
        self.setpoint = pyo.Var(
            time_set, initialize=0.5, doc="Setpoint", units=pv_units
        )
        self.gain_p = pyo.Var(
            time_set,
            initialize=0.1,
            doc="Gain for proportional part",
            units=gain_p_units,
        )
        if self.config.type in [ControllerType.PI, ControllerType.PID]:
            self.gain_i = pyo.Var(
                time_set,
                initialize=0.1,
                doc="Gain for integral part",
                units=gain_p_units / time_units,
            )
        if self.config.type in [ControllerType.PD, ControllerType.PID]:
            self.gain_d = pyo.Var(
                time_set,
                initialize=0.01,
                doc="Gain for derivative part",
                units=gain_p_units * time_units,
            )
        self.mv_ref = pyo.Var(
            time_set,
            initialize=0.5,
            doc="Controller bias",
            units=mv_units,
        )

        # Error expression or variable (variable required for derivative term)
        if self.config.type in [ControllerType.PD, ControllerType.PID]:

            self.error = pyo.Var(
                time_set, initialize=0, doc="Error variable", units=pv_units
            )

            @self.Constraint(time_set, doc="Error constraint")
            def error_eqn(b, t):
                return b.error[t] == b.setpoint[t] - b.process_var[t]

            self.derivative_of_error = pyodae.DerivativeVar(
                self.error,
                wrt=self.flowsheet().time,
                initialize=0,
                units=pv_units / time_units,
            )

        else:

            @self.Expression(time_set, doc="Error expression")
            def error(b, t):
                return b.setpoint[t] - b.process_var[t]

        # integral term written de_i(t)/dt = e(t)
        if self.config.type in [ControllerType.PI, ControllerType.PID]:
            self.integral_of_error = pyo.Var(
                time_set,
                initialize=0,
                doc="Integral term calculated from de_i(t)/dt = e(t)",
                units=pv_units * time_units,
            )
            self.integral_of_error_dot = pyodae.DerivativeVar(
                self.integral_of_error,
                wrt=time_set,
                initialize=0,
                units=pv_units,
                doc="de_i(t)/dt",
            )

            @self.Constraint(time_set, doc="de_i(t)/dt = e(t)")
            def error_from_integral_eqn(b, t):
                return b.error[t] == b.integral_of_error_dot[t]

            if self.config.calculate_initial_integral:
                t0 = time_set.first()

                @self.Constraint(doc="Calculate initial e_i based on output")
                def initial_integral_error_eqn(b):
                    if self.config.type == ControllerType.PI:
                        return (
                            self.integral_of_error[t0]
                            == (
                                b.manipulated_var[t0]
                                - b.mv_ref[t0]
                                - b.gain_p[t0] * b.error[t0]
                            )
                            / b.gain_i[t0]
                        )
                    return (
                        self.integral_of_error[t0]
                        == (
                            b.manipulated_var[t0]
                            - b.mv_ref[t0]
                            - b.gain_p[t0] * b.error[t0]
                            - b.gain_d[t0] * b.derivative_of_error[t0]
                        )
                        / b.gain_i[t0]
                    )

        @self.Expression(time_set, doc="Unbounded output for manipulated variable")
        def mv_unbounded(b, t):
            if self.config.type == ControllerType.PID:
                return (
                    b.mv_ref[t]
                    + b.gain_p[t] * b.error[t]
                    + b.gain_i[t] * b.integral_of_error[t]
                    + b.gain_d[t] * b.derivative_of_error[t]
                )
            elif self.config.type == ControllerType.PI:
                return (
                    b.mv_ref[t]
                    + b.gain_p[t] * b.error[t]
                    + b.gain_i[t] * b.integral_of_error[t]
                )
            elif self.config.type == ControllerType.PD:
                return (
                    b.mv_ref[t]
                    + b.gain_p[t] * b.error[t]
                    + b.gain_d[t] * b.derivative_of_error[t]
                )
            elif self.config.type == ControllerType.P:
                return b.mv_ref[t] + b.gain_p[t] * b.error[t]
            else:
                raise ConfigurationError(
                    f"{self.config.type} is not a valid PID controller type"
                )

        @self.Constraint(time_set, doc="Bounded output of manipulated variable")
        def mv_eqn(b, t):
            if self.config.mv_bound_type == ControllerMVBoundType.SMOOTH_BOUND:
                return b.manipulated_var[t] == smooth_bound(
                    b.mv_unbounded[t], lb=b.mv_lb, ub=b.mv_ub, eps=b.smooth_eps
                )
            elif self.config.mv_bound_type == ControllerMVBoundType.LOGISTIC:
                return (
                    (b.manipulated_var[t] - b.mv_lb)
                    * (
                        1
                        + pyo.exp(
                            -b.logistic_bound_k
                            / (b.mv_ub - b.mv_lb)
                            * (b.mv_unbounded[t] - (b.mv_lb + b.mv_ub) / 2)
                        )
                    )
                ) == b.mv_ub - b.mv_lb
            return b.manipulated_var[t] == b.mv_unbounded[t]

        # deactivate the time 0 mv_eqn instead of skip, should be fine since
        # first time step always exists.
        self.mv_eqn[time_set.first()].deactivate()
