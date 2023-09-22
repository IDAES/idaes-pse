#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
PID controller model module
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

__author__ = ["John Eslick", "Jinliang Ma"]

import enum

import pyomo.environ as pyo
import pyomo.dae as pyodae
from pyomo.common.config import ConfigValue, In, Bool

from idaes.core import UnitModelBlockData, declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_bound
from idaes.core.util import scaling as iscale


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
    LOGISTIC: Use a logistic function to keep mv between the bounds
    """

    NONE = 1
    SMOOTH_BOUND = 2
    LOGISTIC = 3


class ControllerAntiwindupType(enum.Enum):
    """Antiwindup type.

    NONE: No antiwindup scheme.
    CONDITIONAL_INTEGRATION: Error integrates only when the MV is not at a bound. Note that this scheme is "dumb"
    because it does not distinguish between the "correct" bound and "incorrect" bound being active for a particular
    value of integrated error. Note that switching between integrating and not integrating can cause the DAE solver to
    slow down.
    BACK_CALCULATION: Takes the difference between the actual MV output and the unbounded MV output, multiplies it by
    a gain, then subtracts that product from the error integral term. Benefits from having only one non-smooth term (the
    MV bounding equation) and supposedly can give better control performance if tuned well. The downside is an additional
    tuning parameter is needed.
    """

    NONE = 1
    CONDITIONAL_INTEGRATION = 2
    BACK_CALCULATION = 3


def smooth_heaviside(x, k):
    return 1 / (1 + pyo.exp(-2 * k * x))


@declare_process_block_class(
    "PIDController",
    doc="PID controller model block.  To use this the model must be dynamic.",
)
class PIDControllerData(UnitModelBlockData):
    """
    PID controller class.
    """

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
            domain=Bool,
            description="Calculate the initial integral term value if True",
            doc="Calculate the initial integral term value if True",
        ),
    )
    CONFIG.declare(
        "controller_type",
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
    CONFIG.declare(
        "antiwindup_type",
        ConfigValue(
            default=ControllerAntiwindupType.NONE,
            domain=In(
                [
                    ControllerAntiwindupType.NONE,
                    ControllerAntiwindupType.CONDITIONAL_INTEGRATION,
                    ControllerAntiwindupType.BACK_CALCULATION,
                ]
            ),
            description="Type of antiwindup technique to use.",
            doc=(
                """Type of antiwindup technique to use. Options are **ControllerAntiwindupType.NONE**,
                **ControllerAntiwindupType.CONDITIONAL_INTEGRATION**, and **ControllerAntiwindupType.BACK_CALCULATION**.
                See the controller documentation for details on the mathematical formulation.
                """
            ),
        ),
    )
    CONFIG.declare(
        "derivative_on_error",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Whether basing derivative action on process var or error",
            doc="""Naive implementations of derivative action can cause large spikes in
                control when the setpoint is changed. One solution is to use the (negative)
                derivative of the process variable to calculate derivative action instead
                of using the derivative of setpoint error. If **True**, use the derivative of
                setpoint error to calculate derivative action. If **False** (default), use the
                (negative) derivative of the process variable instead.
                """,
        ),
    )

    def build(self):
        """
        Build the PID block
        """
        super().build()

        if self.config.dynamic is False:
            raise ConfigurationError(
                "PIDControllers work only with dynamic flowsheets."
            )

        # Check for required config
        if self.config.process_var is None or self.config.manipulated_var is None:
            raise ConfigurationError(
                "Controller config requires specifying process_var and manipulated_var"
            )
        self.process_var = pyo.Reference(self.config.process_var)
        self.manipulated_var = pyo.Reference(self.config.manipulated_var)

        # Shorter pointers to time set information
        time_set = self.flowsheet().time
        time_units = self.flowsheet().time_units
        if time_units is None:
            time_units = pyo.units.dimensionless
        t0 = time_set.first()

        # Type Check
        if not issubclass(self.process_var[t0].ctype, (pyo.Var, pyo.Expression)):
            raise TypeError(
                f"process_var must reference a Var or Expression not {self.process_var[t0].ctype}"
            )
        if not issubclass(self.manipulated_var[t0].ctype, pyo.Var):
            raise TypeError(
                f"manipulated_var must reference a Var not {self.manipulated_var[t0].ctype}"
            )

        if not self.config.antiwindup_type == ControllerAntiwindupType.NONE:
            if not self.config.controller_type in [
                ControllerType.PI,
                ControllerType.PID,
            ]:
                raise ConfigurationError(
                    "User specified antiwindup method for controller without integral action."
                )
            if self.config.mv_bound_type == ControllerMVBoundType.NONE:
                raise ConfigurationError(
                    "User specified antiwindup method for unbounded MV."
                )

        # Get the appropriate units for various controller variables
        mv_units = pyo.units.get_units(self.manipulated_var[t0])
        pv_units = pyo.units.get_units(self.process_var[t0])
        if mv_units is None:
            mv_units = pyo.units.dimensionless
        if pv_units is None:
            pv_units = pyo.units.dimensionless
        gain_p_units = mv_units / pv_units

        if not self.config.mv_bound_type == ControllerMVBoundType.NONE:
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
        if self.config.mv_bound_type == ControllerMVBoundType.SMOOTH_BOUND:
            self.smooth_eps = pyo.Param(
                mutable=True,
                initialize=1e-4,
                doc="Smoothing parameter for controller output limits when the bound"
                " type is SMOOTH_BOUND",
                units=mv_units,
            )
        elif self.config.mv_bound_type == ControllerMVBoundType.LOGISTIC:
            self.logistic_bound_k = pyo.Param(
                mutable=True,
                initialize=4,
                doc="Smoothing parameter for controller output limits when the bound"
                " type is LOGISTIC",
                units=pyo.units.dimensionless,
            )
        if (
            self.config.antiwindup_type
            == ControllerAntiwindupType.CONDITIONAL_INTEGRATION
        ):
            self.conditional_integration_k = pyo.Param(
                mutable=True,
                initialize=200,
                doc="Parameter governing steepness of transition between integrating and not integrating."
                "A larger value means a steeper transition.",
                units=pyo.units.dimensionless,
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
        if self.config.controller_type in [ControllerType.PI, ControllerType.PID]:
            self.gain_i = pyo.Var(
                time_set,
                initialize=0.1,
                doc="Gain for integral part",
                units=gain_p_units / time_units,
            )
        if self.config.controller_type in [ControllerType.PD, ControllerType.PID]:
            self.gain_d = pyo.Var(
                time_set,
                initialize=0.01,
                doc="Gain for derivative part",
                units=gain_p_units * time_units,
            )

        if self.config.antiwindup_type == ControllerAntiwindupType.BACK_CALCULATION:
            self.gain_b = pyo.Var(
                time_set,
                initialize=0.1,
                doc="Gain for back calculation antiwindup",
                units=1 / time_units,
            )

        self.mv_ref = pyo.Var(
            time_set,
            initialize=0.5,
            doc="Controller bias",
            units=mv_units,
        )

        # Error expression or variable (variable required for derivative term)
        if (
            self.config.controller_type in [ControllerType.PD, ControllerType.PID]
            and self.config.derivative_on_error
        ):
            self.error = pyo.Var(
                time_set, initialize=0, doc="Error variable", units=pv_units
            )

            @self.Constraint(time_set, doc="Error constraint")
            def error_eqn(b, t):
                return b.error[t] == b.setpoint[t] - b.process_var[t]

            self.derivative_term = pyodae.DerivativeVar(
                self.error,
                wrt=self.flowsheet().time,
                initialize=0,
                units=pv_units / time_units,
            )

        else:

            @self.Expression(time_set, doc="Error expression")
            def error(b, t):
                return b.setpoint[t] - b.process_var[t]

            if (
                self.config.controller_type in [ControllerType.PD, ControllerType.PID]
                and not self.config.derivative_on_error
            ):
                # Need to create a Var because process_var might be an Expression
                self.negative_pv = pyo.Var(
                    time_set,
                    initialize=0,
                    doc="Negative of process variable",
                    units=pv_units,
                )

                @self.Constraint(time_set, doc="Negative process variable equation")
                def negative_pv_eqn(b, t):
                    return b.negative_pv[t] == -b.process_var[t]

                self.derivative_term = pyodae.DerivativeVar(
                    self.negative_pv,
                    wrt=self.flowsheet().time,
                    initialize=0,
                    units=pv_units / time_units,
                )

        # integral term written de_i(t)/dt = e(t)
        if self.config.controller_type in [ControllerType.PI, ControllerType.PID]:
            self.mv_integral_component = pyo.Var(
                time_set,
                initialize=0,
                doc="Integral contribution to control action",
                units=mv_units,
            )
            self.mv_integral_component_dot = pyodae.DerivativeVar(
                self.mv_integral_component,
                wrt=time_set,
                initialize=0,
                units=mv_units / time_units,
                doc="Rate of change of integral contribution to control action",
            )

            if self.config.calculate_initial_integral:

                @self.Constraint(doc="Calculate initial e_i based on output")
                def initial_integral_error_eqn(b):
                    if self.config.controller_type == ControllerType.PI:
                        return b.mv_integral_component[t0] == (
                            b.manipulated_var[t0]
                            - b.mv_ref[t0]
                            - b.gain_p[t0] * b.error[t0]
                        )
                    return b.mv_integral_component[t0] == (
                        b.manipulated_var[t0]
                        - b.mv_ref[t0]
                        - b.gain_p[t0] * b.error[t0]
                        - b.gain_d[t0] * b.derivative_term[t0]
                    )

        @self.Expression(time_set, doc="Unbounded output for manipulated variable")
        def mv_unbounded(b, t):
            if self.config.controller_type == ControllerType.PID:
                return (
                    b.mv_ref[t]
                    + b.gain_p[t] * b.error[t]
                    + b.mv_integral_component[t]
                    + b.gain_d[t] * b.derivative_term[t]
                )
            elif self.config.controller_type == ControllerType.PI:
                return (
                    b.mv_ref[t] + b.gain_p[t] * b.error[t] + b.mv_integral_component[t]
                )
            elif self.config.controller_type == ControllerType.PD:
                return (
                    b.mv_ref[t]
                    + b.gain_p[t] * b.error[t]
                    + b.gain_d[t] * b.derivative_term[t]
                )
            elif self.config.controller_type == ControllerType.P:
                return b.mv_ref[t] + b.gain_p[t] * b.error[t]
            else:
                raise ConfigurationError(
                    f"{self.config.controller_type} is not a valid PID controller type"
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
                            * (b.mv_unbounded[t] - (b.mv_lb + b.mv_ub) / 2)
                            / (b.mv_ub - b.mv_lb)
                        )
                    )
                ) == b.mv_ub - b.mv_lb
            return b.manipulated_var[t] == b.mv_unbounded[t]

        # deactivate the time 0 mv_eqn instead of skip, should be fine since
        # first time step always exists.
        if self.config.calculate_initial_integral:
            self.mv_eqn[t0].deactivate()

        if self.config.controller_type in [ControllerType.PI, ControllerType.PID]:

            @self.Constraint(time_set, doc="de_i(t)/dt = e(t)")
            def mv_integration_eqn(b, t):
                if (
                    self.config.antiwindup_type
                    == ControllerAntiwindupType.CONDITIONAL_INTEGRATION
                ):
                    # This expression is not sensitive to whether the "right" or "wrong" bound is active for a given
                    # expression of error.
                    return b.mv_integral_component_dot[t] == b.gain_i[t] * b.error[
                        t
                    ] * (
                        smooth_heaviside(
                            (b.mv_unbounded[t] - b.mv_lb) / (b.mv_ub - b.mv_lb),
                            b.conditional_integration_k,
                        )
                        # 1
                        - smooth_heaviside(
                            (b.mv_unbounded[t] - b.mv_ub) / (b.mv_ub - b.mv_lb),
                            b.conditional_integration_k,
                        )
                    )
                elif (
                    self.config.antiwindup_type
                    == ControllerAntiwindupType.BACK_CALCULATION
                ):
                    return b.mv_integral_component_dot[t] == b.gain_i[t] * b.error[
                        t
                    ] + b.gain_b[t] * (b.manipulated_var[t] - b.mv_unbounded[t])
                else:
                    return b.mv_integral_component_dot[t] == b.gain_i[t] * b.error[t]

            self.mv_integration_eqn[t0].deactivate()

    def calculate_scaling_factors(self):
        super().calculate_scaling_factors()
        gsf = iscale.get_scaling_factor
        ssf = lambda c, v: iscale.set_scaling_factor(c, v, overwrite=False)
        cst = lambda c, v: iscale.constraint_scaling_transform(c, v, overwrite=False)

        # orig_pv = self.config.process_var
        # orig_mv = self.config.manipulated_var
        time_set = self.flowsheet().time
        t0 = time_set.first()

        sf_pv = iscale.get_scaling_factor(self.config.process_var[t0])
        if sf_pv is None:
            sf_pv = iscale.get_scaling_factor(self.process_var[t0], default=1)
        sf_mv = iscale.get_scaling_factor(self.config.manipulated_var[t0])
        if sf_mv is None:
            sf_mv = iscale.get_scaling_factor(self.manipulated_var[t0], default=1)

        # Don't like calling scaling laterally like this, but we need scaling factors for the pv and mv
        # Except this causes a StackOverflow with flowsheet-level PVs or MVs---put this on ice for now
        # if sf_pv is None:
        #     try:
        #         iscale.calculate_scaling_factors(self.config.process_var.parent_block())
        #     except RecursionError:
        #         raise ConfigurationError(
        #             f"Circular scaling dependency detected in Controller {self.name}. The only way this should be "
        #             "able to happen is if a loop of controllers exists manipulating each others setpoints without "
        #             "terminating in an actual process variable."
        #         )
        # if sf_mv is None:
        #     try:
        #         iscale.calculate_scaling_factors(self.config.manipulated_var.parent_block())
        #     except RecursionError:
        #         raise ConfigurationError(
        #             f"Circular scaling dependency detected in Controller {self.name}. The only way this should be "
        #             "able to happen is if a loop of controllers exists manipulating each others setpoints without "
        #             "terminating in an actual process variable."
        #         )

        if self.config.calculate_initial_integral:
            sf_mv = gsf(self.manipulated_var[t0], default=1, warning=True)
            cst(self.initial_integral_error_eqn, sf_mv)

        for t in time_set:
            sf_pv = gsf(self.process_var[t], default=1, warning=True)
            sf_mv = gsf(self.manipulated_var[t], default=1, warning=True)

            ssf(self.setpoint[t], sf_pv)
            ssf(self.mv_ref[t], sf_mv)
            cst(self.mv_eqn[t], sf_mv)

            if self.config.controller_type in [ControllerType.PD, ControllerType.PID]:
                if self.config.derivative_on_error:
                    ssf(self.error[t], sf_pv)
                    cst(self.error_eqn[t], sf_pv)
                else:
                    ssf(self.negative_pv[t], sf_pv)
                    cst(self.negative_pv_eqn[t], sf_pv)

            if self.config.controller_type in [ControllerType.PI, ControllerType.PID]:
                ssf(self.mv_integral_component[t], sf_mv)

                cst(self.mv_integration_eqn[t], sf_pv)
