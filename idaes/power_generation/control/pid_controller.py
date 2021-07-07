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

__author__ = "John Eslick and Jinliang Ma"


from pyomo.environ import Var, Param, Constraint, Reference,\
    exp, log
from pyomo.dae import DerivativeVar

from idaes.core import UnitModelBlockData, declare_process_block_class
from pyomo.common.config import ConfigValue, In
from idaes.core.util.exceptions import ConfigurationError


@declare_process_block_class("PIDController", doc="""This is a PID controller
                             block. Based on UnitModelBlockData to configure
                             the dynamic flag based on parent""")
class PIDControllerData(UnitModelBlockData):
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("pv", ConfigValue(
        default=None,
        description="Process variable to be controlled",
        doc="A Pyomo Var, Expression, or Reference for the measured"
            " process variable. Should be indexed by time."))
    CONFIG.declare("mv", ConfigValue(
        default=None,
        description="Manipulated process variable",
        doc="A Pyomo Var, Expression, or Reference for the controlled"
            " process variable. Should be indexed by time."))
    CONFIG.declare("bounded_output", ConfigValue(
        default=False,
        description="Flag to bound manipulated variable",
        doc="""Indicating if the output for the manipulated variable is bounded
Default: False.
If True, user need to set the lower and upper bound parameters"""))
    CONFIG.declare("type", ConfigValue(
        default="PI",
        domain=In(['P', 'PI', 'PD', 'PID']),
        description="Control type",
        doc="""Controller type options including
- P: Proportional only
- PI: Proportional and integral only
- PD: Proportional and derivative only
- PID: Proportional, integral and derivative
Default is PI"""))

    def build(self):
        """
        Build the PID block
        """
        super().build()

        # Do nothing if steady-state
        if self.config.dynamic is True:
            # Check for required config
            if self.config.pv is None:
                raise ConfigurationError("Controller configuration"
                                         " requires 'pv'")
            if self.config.mv is None:
                raise ConfigurationError("Controller configuration"
                                         " requires 'mv'")
            # Shorter pointers to time set information
            time_set = self.flowsheet().time
            self.pv = Reference(self.config.pv)
            self.mv = Reference(self.config.mv)

            # Parameters
            self.mv_lb = Param(mutable=True,
                               initialize=0.05,
                               doc="Controller output lower bound")
            self.mv_ub = Param(mutable=True,
                               initialize=1,
                               doc="Controller output upper bound")

            # Variable for basic controller settings may change with time.
            self.setpoint = Var(time_set,
                                initialize=0.5,
                                doc="Setpoint")
            self.gain_p = Var(time_set,
                              initialize=0.1,
                              doc="Gain for proportional part")
            if self.config.type == 'PI' or self.config.type == 'PID':
                self.gain_i = Var(time_set,
                                  initialize=0.1,
                                  doc="Gain for integral part")
            if self.config.type == 'PD' or self.config.type == 'PID':
                self.gain_d = Var(time_set,
                                  initialize=0.01,
                                  doc="Gain for derivative part")
            self.mv_ref = Var(initialize=0.5,
                              doc="bias value of manipulated variable")

            if self.config.type == 'P' or self.config.type == 'PI':
                @self.Expression(time_set,
                                 doc="Error expression")
                def error(b, t):
                    return b.setpoint[t] - b.pv[t]
            else:
                self.error = Var(time_set, initialize=0, doc="Error variable")
                @self.Constraint(time_set, doc="Error variable")
                def error_eqn(b, t):
                    return b.error[t] == b.setpoint[t] - b.pv[t]

            if self.config.type == 'PI' or self.config.type == 'PID':
                self.integral_of_error = Var(time_set,
                                             initialize=0,
                                             doc="Integral term")
                self.error_from_integral = DerivativeVar(
                    self.integral_of_error,
                    wrt=self.flowsheet().time,
                    initialize=0)

                @self.Constraint(time_set,
                                 doc="Error calculated by"
                                 " derivative of integral")
                def error_from_integral_eqn(b, t):
                    return b.error[t] == b.error_from_integral[t]

            if self.config.type == 'PID' or self.config.type == 'PD':
                self.derivative_of_error = DerivativeVar(
                    self.error,
                    wrt=self.flowsheet().time,
                    initialize=0)

            @self.Expression(time_set, doc="Proportional output")
            def mv_p_only(b, t):
                return b.gain_p[t]*b.error[t]

            @self.Expression(time_set,
                             doc="Proportional output and reference")
            def mv_p_only_with_ref(b, t):
                return b.gain_p[t]*b.error[t]+b.mv_ref

            if self.config.type == 'PI' or self.config.type == 'PID':
                @self.Expression(time_set,
                                 doc="Integral output")
                def mv_i_only(b, t):
                    return b.gain_i[t] * b.integral_of_error[t]

            if self.config.type == 'PD' or self.config.type == 'PID':
                @self.Expression(time_set,
                                 doc="Derivative output")
                def mv_d_only(b, t):
                    return b.gain_d[t] * b.derivative_of_error[t]

            @self.Expression(time_set,
                             doc="Unbounded output for manimulated variable")
            def mv_unbounded(b, t):
                if self.config.type == 'PID':
                    return (b.mv_ref + b.gain_p[t]*b.error[t]
                            + b.gain_i[t]*b.integral_of_error[t]
                            + b.gain_d[t]*b.derivative_of_error[t])
                elif self.config.type == 'PI':
                    return (b.mv_ref + b.gain_p[t]*b.error[t]
                            + b.gain_i[t]*b.integral_of_error[t])
                elif self.config.type == 'PD':
                    return (b.mv_ref + b.gain_p[t]*b.error[t]
                            + b.gain_d[t]*b.derivative_of_error[t])
                else:
                    return b.mv_ref + b.gain_p[t]*b.error[t]

            @self.Constraint(time_set,
                             doc="Bounded output of manipulated variable")
            def mv_eqn(b, t):
                if t == b.flowsheet().time.first():
                    return Constraint.Skip
                else:
                    if self.config.bounded_output is True:
                        return (b.mv[t]-b.mv_lb) * \
                            (1+exp(-4/(b.mv_ub-b.mv_lb)
                                   * (b.mv_unbounded[t]
                                      - (b.mv_lb
                                         + b.mv_ub)/2))) == b.mv_ub-b.mv_lb
                    else:
                        return b.mv[t] == b.mv_unbounded[t]
            if self.config.bounded_output is True:
                if self.config.type == 'PI' or self.config.type == 'PID':
                    @self.Expression(time_set,
                                     doc="Integral error"
                                     " at error 0 and mv_ref")
                    def integral_of_error_ref(b, t):
                        return ((b.mv_lb+b.mv_ub)/2
                                - b.mv_ref-log((b.mv_ub-b.mv_lb)
                                               / (b.mv_ref-b.mv_lb)-1)
                                / 4*(b.mv_ub-b.mv_lb))/b.gain_i[t]

                    @self.Expression(time_set,
                                     doc="Integral error at error 0 and"
                                     " output value at current mv")
                    def integral_of_error_mv(b, t):
                        return ((b.mv_lb+b.mv_ub)/2-b.mv_ref
                                - log((b.mv_ub-b.mv_lb)/(b.mv[t]-b.mv_lb)-1)
                                / 4*(b.mv_ub-b.mv_lb))/b.gain_i[t]
