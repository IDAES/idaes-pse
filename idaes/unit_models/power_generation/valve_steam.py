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
from pyomo.environ import Var, Expression, SolverFactory, value, sqrt
from pyomo.opt import TerminationCondition

from idaes.core import declare_process_block_class
from idaes.unit_models.pressure_changer import PressureChangerData
from idaes.core.util import from_json, to_json, StoreSpec
from idaes.ui.report import degrees_of_freedom
from .valve_steam_config import _define_config, ValveFunctionType

def _linear_rule(b, t):
    return b.x[t]

def _quick_open_rule(b, t):
    return sqrt(b.x[t])

def _equal_percentage_rule(b, t):
    return b.alpha**(b.x[t] - 1)

def _liquid_pressure_flow_rule(b, t):
    """
    For liquid F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    Cv = self.Cv
    fun = self.valve_function[t]
    return F**2 == Cv**2*(Pi - Po)*fun**2

def _vapor_pressure_flow_rule(b, t):
    """
    For vapor F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    Po = b.control_volume.properties_out[t].pressure
    Pi = b.control_volume.properties_in[t].pressure
    F = b.control_volume.properties_in[t].flow_mol
    Cv = self.Cv
    fun = self.valve_function[t]
    return F**2 == Cv**2*(Pi**2 - Po**2)*fun**2


@declare_process_block_class("SteamValve", doc="Basic steam valve models")
class SteamValveData(PressureChangerData):
    # Same settings as the default pressure changer, but force to expander with
    # isentropic efficiency
    CONFIG = PressureChangerData.CONFIG()
    _define_config(CONFIG)

    def build(self):
        super(TurbineStageData, self).build()

        self.valve_opening = Var(self.time_ref, initialize=1,
            doc="Fraction open for valve from 0 to 1")
        self.Cv = Var(initialize=1, doc="Valve flow coefficent, for vapor "
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
