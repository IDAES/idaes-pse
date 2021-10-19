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
PID controller block, Revised to use pyomo.DAE
to calculate integral and differential parts
"""

__author__ = ["John Eslick", "Jinliang Ma"]

import enum
import pyomo.environ as pyo
import pyomo.dae as pyodae

from idaes.core import UnitModelBlockData, declare_process_block_class
from pyomo.common.config import ConfigValue, In
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_bound
import idaes.logger as idaeslog
from idaes.core.util import get_solver


class ControllerType(enum.Enum):
    """Controller types."""

    P = 1
    PI = 2
    PD = 3
    PID = 4


@declare_process_block_class(
    "PIDController",
    doc="PID controller model block.  To use this the model must be dynamic.",
)
class PIDControllerData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare(
        "pv",
        ConfigValue(
            default=None,
            description="Process variable to be controlled",
            doc=(
                "A Pyomo Var, Expression, or Reference for the measured"
                " process variable. Should be indexed by time, slices are okay."
            ),
        ),
    )
    CONFIG.declare(
        "mv",
        ConfigValue(
            default=None,
            description="Manipulated process variable",
            doc=(
                "A Pyomo Var, Expression, or Reference for the controlled"
                " process variable. Should be indexed by time, slices are okay."
            ),
        ),
    )
    CONFIG.declare(
        "bounded_output",
        ConfigValue(
            default=False,
            description="Flag to bound manipulated variable",
            doc=(
                "If True, the output for the manipulated variable is bounded,"
                " between the 'mv_lb' and 'mv_ub' parameters. (Default=False)"
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
            doc="""Controller type options including
- ControllerType.P: Proportional only
- ControllerType.PI: Proportional and integral only
- ControllerType.PD: Proportional and derivative only
- ControllerType.PID: Proportional, integral and derivative
Default is ControllerType.PI""",
        ),
    )

    def build(self):
        """
        Build the PID block
        """
        super().build()

        # Check for required config
        if self.config.pv is None or self.config.mv is None:
            raise ConfigurationError("Controller config requires 'pv' and 'mv'")

        # Make local references to the measured process varaible (pv) and the
        # manipulated variable (mv)
        self.pv = pyo.Reference(self.config.pv)
        self.mv = pyo.Reference(self.config.mv)

        # Shorter pointers to time set information
        time_set = self.flowsheet().time
        time_units = self.flowsheet().time_units

        # Get the appropriate units for various contoller varaibles
        mv_units = pyo.units.get_units(self.mv[time_set.first()])
        pv_units = pyo.units.get_units(self.pv[time_set.first()])
        if mv_units is None:
            mv_units = pyo.units.dimensionless
        if pv_units is None:
            pv_units = pyo.units.dimensionless
        if time_units is None:
            time_units = pyo.units.dimensionless
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
            doc="Smoothing parameter for controller output limits",
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
                return b.error[t] == b.setpoint[t] - b.pv[t]

            self.derivative_of_error = pyodae.DerivativeVar(
                self.error,
                wrt=self.flowsheet().time,
                initialize=0,
                units=pv_units / time_units,
            )

        else:

            @self.Expression(time_set, doc="Error expression")
            def error(b, t):
                return b.setpoint[t] - b.pv[t]

        # integral term written de_i(t)/dt = e(t)
        if self.config.type in [ControllerType.PI, ControllerType.PID]:
            self.integral_of_error = pyo.Var(
                time_set, initialize=0, doc="Integral term", units=pv_units * time_units
            )
            self.integral_of_error_dot = pyodae.DerivativeVar(
                self.integral_of_error, wrt=time_set, initialize=0, units=pv_units
            )

            @self.Constraint(time_set, doc="Error calculated by derivative of integral")
            def error_from_integral_eqn(b, t):
                if t == time_set.first():
                    if self.config.calculate_initial_integral:
                        if self.config.type == ControllerType.PI:
                            return (
                                self.integral_of_error[t]
                                == (b.mv[t] - b.mv_ref[t] - b.gain_p[t] * b.error[t])
                                / b.gain_i[t]
                            )
                        return (
                            self.integral_of_error[t]
                            == (
                                b.mv[t]
                                - b.mv_ref[t]
                                - b.gain_p[t] * b.error[t]
                                - b.gain_d[t] * b.derivative_of_error[t]
                            )
                            / b.gain_i[t]
                        )
                    return pyo.Constraint.Skip
                return b.error[t] == b.integral_of_error_dot[t]

        @self.Expression(time_set, doc="Unbounded output for manimulated variable")
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
            if t == time_set.first():
                return pyo.Constraint.Skip
            else:
                if self.config.bounded_output is True:
                    l = b.mv_lb
                    h = b.mv_ub
                    mv = b.mv_unbounded[t]
                    e = b.smooth_eps
                    return b.mv[t] == smooth_bound(mv, lb=l, ub=h, eps=e)
                else:
                    return b.mv[t] == b.mv_unbounded[t]
