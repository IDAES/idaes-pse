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

from pyutilib.enum import Enum

ValveFunctionType = Enum("linear", "quick_opening", "equal_percentage")

@declare_process_block_class("SteamValve", doc="Basic steam valve models")
class SteamValveData(PressureChangerData):
    # Same settings as the default pressure changer, but force to expander with
    # isentropic efficiency
    CONFIG = PressureChangerData.CONFIG()
    CONFIG.compressor = False
    CONFIG.get('compressor')._default = False
    CONFIG.get('compressor')._domain = In([False])
    CONFIG.thermodynamic_assumption = 'adiabatic'
    CONFIG.get('thermodynamic_assumption')._default = 'adiabatic'
    CONFIG.get('thermodynamic_assumption')._domain = In(['adiabatic'])
    CONFIG.declare("valve_function", ConfigValue(
        default=ValveFunctionType.linear,
        domain=In(ValveFunctionType),
        description="Valve function type",
        doc="""The type of valve function
**default** - ValveFunctionType.linear
**Valid values:** {
ValveFunctionType.linear,
ValveFunctionType.quick_opening
ValveFunctionType.equal_percentage}"""))
    def build(self):
        super(TurbineStageData, self).build()
