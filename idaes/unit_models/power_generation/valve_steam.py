##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
This provides valve models for steam and liquid water.  These are for
steam cycle control valves and the turbine throttle valves.
"""


from __future__ import division

__Author__ = "John Eslick"

import logging
_log = logging.getLogger(__name__)

from pyomo.common.config import In, ConfigValue
from pyomo.environ import Var, Expression, SolverFactory, value, Constraint
from pyomo.opt import TerminationCondition

from idaes.core import declare_process_block_class
from idaes.unit_models.pressure_changer import PressureChangerData
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.ui.report import degrees_of_freedom
from .valve_steam_config import _define_config, ValveFunctionType

def _linear_rule(b, t):
    return b.valve_opening[t]

def _quick_open_rule(b, t):
    return sqrt(b.valve_opening[t])

def _equal_percentage_rule(b, t):
    return b.alpha**(b.valve_opening[t] - 1)

def _liquid_pressure_flow_rule(b, t):
    """
    For liquid F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    Cv = b.Cv
    fun = b.valve_function[t]
    return 1e-3*F**2 == 1e-3*Cv**2*(Pi - Po)*fun**2

def _vapor_pressure_flow_rule(b, t):
    """
    For vapor F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    Cv = b.Cv
    fun = b.valve_function[t]
    return 1e-6*F**2 == 1e-6*Cv**2*(Pi**2 - Po**2)*fun**2


@declare_process_block_class("SteamValve", doc="Basic steam valve models")
class SteamValveData(PressureChangerData):
    # Same settings as the default pressure changer, but force to expander with
    # isentropic efficiency
    CONFIG = PressureChangerData.CONFIG()
    _define_config(CONFIG)

    def build(self):
        super().build()

        self.valve_opening = Var(self.time_ref, initialize=1,
            doc="Fraction open for valve from 0 to 1")
        self.Cv = Var(initialize=0.01, doc="Valve flow coefficent, for vapor "
            "[mol/s/Pa] for liquid [mol/s/Pa^0.5]")
        self.Cv.fix()
        self.valve_opening.fix()

        # set up the valve function rule.  I'm not sure these matter too much
        # for us, but the options are easy enough to provide.
        if self.config.valve_function == ValveFunctionType.linear:
            rule = _linear_rule
        elif self.config.valve_function == ValveFunctionType.quick_opening:
            rule = _quick_open_rule
        elif self.config.valve_function == ValveFunctionType.equal_percentage:
            self.alpha = Var(initialize=1, doc="Valve function parameter")
            self.alpha.fix()
            rule = equal_percentage_rule
        else:
            rule = self.config.valve_function_rule

        self.valve_function = Expression(self.time_ref, rule=rule,
                doc="Valve function expression")

        if self.config.phase == "Liq":
            rule = _liquid_pressure_flow_rule
        else:
            rule = _vapor_pressure_flow_rule

        self.pressure_flow_equation = Constraint(self.time_ref, rule=rule)

    def initialize(self, state_args={}, outlvl=0, solver='ipopt',
        optarg={'tol': 1e-6, 'max_iter':30}):
        """
        Initialize the turbine stage model.  This deactivates the
        specialized constraints, then does the isentropic turbine initialization,
        then reactivates the constraints and solves.

        Args:
            state_args (dict): Initial state for property initialization
            outlvl (int): Amount of output (0 to 3) 0 is lowest
            solver (str): Solver to use for initialization
            optarg (dict): Solver arguments dictionary
        """
        stee = True if outlvl >= 3 else False
        # sp is what to save to make sure state after init is same as the start
        #   saves value, fixed, and active state, doesn't load originally free
        #   values, this makes sure original problem spec is same but initializes
        #   the values of free vars
        sp = StoreSpec.value_isfixed_isactive(only_fixed=True)
        istate = to_json(self, return_dict=True, wts=sp)

        self.deltaP[:].unfix()
        self.ratioP[:].unfix()

        # fix inlet and free outlet
        for t in self.time_ref:
            for k, v in self.inlet.vars.items():
                v[t].fix()
            for k, v in self.outlet.vars.items():
                v[t].unfix()
            # to calculate outlet pressure
            Pout = self.outlet.pressure[t]
            Pin = self.inlet.pressure[t]
            if self.deltaP[t].value is not None:
                prdp = value((self.deltaP[t] - Pin)/Pin)
            else:
                prdp = -100 # crazy number to say don't use deltaP as guess
            if value(Pout/Pin) > 1 or value(Pout/Pin) < 0.0:
                if value(self.ratioP[t]) <= 1 and value(self.ratioP[t]) >= 0:
                    Pout.value = value(Pin*self.ratioP[t])
                elif prdp <= 1 and prdp >= 0:
                    Pout.value = value(prdp*Pin)
                else:
                    Pout.value = value(Pin*0.95)
            self.deltaP[t] = value(Pout - Pin)
            self.ratioP[t] = value(Pout/Pin)

        # Make sure the initialization problem has no degrees of freedom
        # This shouldn't happen here unless there is a bug in this
        dof = degrees_of_freedom(self)
        print(dof)
        try:
            assert(dof == 0)
        except:
            _log.exception("degrees_of_freedom = {}".format(dof))
            raise

        # one bad thing about reusing this is that the log messages aren't
        # really compatible with being nested inside another initialization
        super().initialize(state_args=state_args,
            outlvl=outlvl, solver=solver, optarg=optarg)

        # reload original spec
        from_json(self, sd=istate, wts=sp)
