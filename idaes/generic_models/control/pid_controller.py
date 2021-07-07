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
PID controller block
"""

__author__ = "John Eslick"

from enum import Enum

import pyomo.environ as pyo

from idaes.core import ProcessBlockData, declare_process_block_class
from pyomo.common.config import ConfigValue, In
from idaes.core.util.math import smooth_max, smooth_min
from idaes.core.util.exceptions import ConfigurationError
from pyomo.dae import ContinuousSet


class PIDForm(Enum):
    """Enum for the pid ``pid_form`` option.
    Either standard or velocity form."""
    standard = 1
    velocity = 2


@declare_process_block_class("PIDBlock", doc="""
This is a PID controller block. The PID Controller block must be added
after the DAE transformation.""")
class PIDBlockData(ProcessBlockData):
    CONFIG = ProcessBlockData.CONFIG()
    CONFIG.declare("pv", ConfigValue(
        default=None,
        description="Measured process variable",
        doc="A Pyomo Var, Expression, or Reference for the measured"
            " process variable. Should be indexed by time."))
    CONFIG.declare("output", ConfigValue(
        default=None,
        description="Controlled process variable",
        doc="A Pyomo Var, Expression, or Reference for the controlled"
            " process variable. Should be indexed by time."))
    CONFIG.declare("upper", ConfigValue(
        default=1.0,
        domain=float,
        description="Output upper limit",
        doc="The upper limit for the controller output, default=1"))
    CONFIG.declare("lower", ConfigValue(
        default=0.0,
        domain=float,
        description="Output lower limit",
        doc="The lower limit for the controller output, default=0"))
    CONFIG.declare("calculate_initial_integral", ConfigValue(
        default=True,
        domain=bool,
        description="Calculate the initial integral term value if true, "
                    " otherwise provide a variable err_i0, which can be fixed",
        doc="Calculate the initial integral term value if true, otherwise"
            " provide a variable err_i0, which can be fixed, default=True"))
    CONFIG.declare("pid_form", ConfigValue(
        default=PIDForm.velocity,
        domain=In(PIDForm),
        description="Velocity or standard form",
        doc="Velocity or standard form"))
    # TODO<jce> options for P, PI, and PD, you can currently do PI by setting
    #           the derivative time to 0, this class should handle PI and PID
    #           controllers. Proportional, only controllers are sufficiently
    #           different that another class should be implemented.
    # TODO<jce> Anti-windup the integral term can keep accumulating error when
    #           the controller output is at a bound. This can cause trouble,
    #           and ways to deal with it should be implemented
    # TODO<jce> Implement way to better deal with the integral term for setpoint
    #           changes (see bumpless).  I need to look into the more, but this
    #           would basically use the calculation like the one already used
    #           for the first time point to calculate integral error to keep the
    #           controller output from suddenly jumping in response to a set
    #           point change or transition from manual to automatic control.

    def _build_standard(self, time_set, t0):
        # Want to fix the output variable at the first time step to make
        # solving easier. This calculates the initial integral error to line up
        # with the initial output value, keeps the controller from initially
        # jumping.
        if self.config.calculate_initial_integral:
            @self.Expression(doc="Initial integral error")
            def err_i0(b):
                return (b.time_i[t0]*(b.output[t0] - b.gain[t0]*b.pterm[t0] -
                                      b.gain[t0]*b.time_d[t0]*b.err_d[t0]) /
                        b.gain[t0])
        # integral error
        @self.Expression(time_set, doc="Integral error")
        def err_i(b, t_end):
            return b.err_i0 + sum((b.iterm[t] + b.iterm[time_set.prev(t)]) *
                                  (t - time_set.prev(t))/2.0
                                  for t in time_set if t <= t_end and t > t0)
        # Calculate the unconstrained controller output
        @self.Expression(time_set, doc="Unconstrained controller output")
        def unconstrained_output(b, t):
            return b.gain[t]*(
                b.pterm[t] +
                1.0/b.time_i[t]*b.err_i[t] +
                b.time_d[t]*b.err_d[t]
            )

        @self.Expression(doc="Initial integral error at the end")
        def err_i_end(b):
            return b.err_i[time_set.last()]

    def _build_velocity(self, time_set, t0):
        if self.config.calculate_initial_integral:
            @self.Expression(doc="Initial integral error")
            def err_i0(b):
                return (b.time_i[t0]*(b.output[t0] - b.gain[t0]*b.pterm[t0] -
                                      b.gain[t0]*b.time_d[t0]*b.err_d[t0]) /
                        b.gain[t0])

        # Calculate the unconstrained controller output
        @self.Expression(time_set, doc="Unconstrained controller output")
        def unconstrained_output(b, t):
            if t == t0:
                # do the standard first step so I have a previous time
                # for the rest of the velocity form
                return b.gain[t]*(
                    b.pterm[t] +
                    1.0/b.time_i[t]*b.err_i0 +
                    b.time_d[t]*b.err_d[t]
                )
            tb = time_set.prev(t)  # time back a step
            return self.output[tb] + self.gain[t]*(
                b.pterm[t] - b.pterm[tb] +
                (t - tb)/b.time_i[t]*(b.err[t] + b.err[tb])/2 +
                b.time_d[t]*(b.err_d[t] - b.err_d[tb])
            )

        @self.Expression(doc="Initial integral error at the end")
        def err_i_end(b):
            tl = time_set.last()
            return (b.time_i[tl]*(b.output[tl] - b.gain[tl]*b.pterm[tl] -
                                  b.gain[tl]*b.time_d[tl]*b.err_d[tl]) /
                    b.gain[tl])

    def build(self):
        """
        Build the PID block
        """
        if isinstance(self.flowsheet().time, ContinuousSet):
            # time may not be a continuous set if you have a steady state model
            # in the steady state model case obviously the controller should
            # not be active, but you can still add it.
            if 'scheme' not in self.flowsheet().time.get_discretization_info():
                # if you have a dynamic model, must do time discretization
                # before adding the PID model
                raise RuntimeError(
                    "PIDBlock must be added after time discretization")

        super().build()  # do the ProcessBlockData voodoo for config
        # Check for required config
        if self.config.pv is None:
            raise ConfigurationError("Controller configuration requires 'pv'")
        if self.config.output is None:
            raise ConfigurationError(
                "Controller configuration requires 'output'")

        # Shorter pointers to time set information
        time_set = self.flowsheet().time
        t0 = time_set.first()

        # Get units of time domain, PV and output
        t_units = self.flowsheet().time_units
        pv_units = self.config.pv.get_units()
        out_units = self.config.output.get_units()
        gain_units = out_units/pv_units if pv_units is not None else None
        err_d_units = pv_units/t_units if pv_units is not None else None
        err_i_units = pv_units*t_units if pv_units is not None else None

        # Variable for basic controller settings may change with time.
        self.setpoint = pyo.Var(time_set,
                                doc="Setpoint",
                                units=pv_units)
        self.gain = pyo.Var(time_set,
                            doc="Controller gain",
                            units=gain_units)
        self.time_i = pyo.Var(time_set, doc="Integral time", units=t_units)
        self.time_d = pyo.Var(time_set, doc="Derivative time", units=t_units)

        # Make the initial derivative term a variable so you can set it. This
        # should let you carry on from the end of another time period
        self.err_d0 = pyo.Var(doc="Initial derivative term",
                              initialize=0,
                              units=err_d_units)
        self.err_d0.fix()

        if not self.config.calculate_initial_integral:
            self.err_i0 = pyo.Var(doc="Initial integral term",
                                  initialize=0,
                                  units=err_i_units)
            self.err_i0.fix()

        # Make references to the output and measured variables
        self.pv = pyo.Reference(self.config.pv)  # No duplicate
        self.output = pyo.Reference(self.config.output)  # No duplicate

        # Create an expression for error from setpoint
        @self.Expression(time_set, doc="Setpoint error")
        def err(b, t):
            return self.setpoint[t] - self.pv[t]

        # Use expressions to allow the some future configuration
        @self.Expression(time_set)
        def pterm(b, t):
            return -self.pv[t]

        @self.Expression(time_set)
        def dterm(b, t):
            return -self.pv[t]

        @self.Expression(time_set)
        def iterm(b, t):
            return self.err[t]
        # Output limits parameter
        self.limits = pyo.Param(["l", "h"],
                                mutable=True,
                                doc="controller output limits",
                                initialize={"l": self.config.lower,
                                            "h": self.config.upper},
                                units=out_units)

        # Smooth min and max are used to limit output, smoothing parameter here
        self.smooth_eps = pyo.Param(
            mutable=True,
            initialize=1e-4,
            doc="Smoothing parameter for controller output limits",
            units=out_units)

        # This is ugly, but want integral and derivative error as expressions,
        # nice implementation with variables is harder to initialize and solve
        @self.Expression(time_set, doc="Derivative error.")
        def err_d(b, t):
            if t == t0:
                return self.err_d0
            else:
                return ((b.dterm[t] - b.dterm[time_set.prev(t)]) /
                        (t - time_set.prev(t)))

        if self.config.pid_form == PIDForm.standard:
            self._build_standard(time_set, t0)
        else:
            self._build_velocity(time_set, t0)

        # Add the controller output constraint and limit it with smooth min/max
        e = self.smooth_eps
        h = self.limits["h"]
        l = self.limits["l"]
        @self.Constraint(time_set, doc="Controller output constraint")
        def output_constraint(b, t):
            if t == t0:
                return pyo.Constraint.Skip
            else:
                return self.output[t] ==\
                    smooth_min(
                        smooth_max(self.unconstrained_output[t], l, e), h, e)
