##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################

"""
Define configuration block for the SteamValve model.
"""


__author__ = "John Eslick"

from enum import Enum

from pyomo.common.config import ConfigBlock, ConfigValue, In, ConfigList

from idaes.unit_models.pressure_changer import ThermodynamicAssumption


class ValveFunctionType(Enum):
    linear = 1
    quick_opening = 2
    equal_percentage = 3
    custom = 4


def _define_config(config):
    config.compressor = False
    config.get('compressor')._default = False
    config.get('compressor')._domain = In([False])
    config.thermodynamic_assumption = ThermodynamicAssumption.adiabatic
    config.get('thermodynamic_assumption')._default = \
        ThermodynamicAssumption.adiabatic
    config.get('thermodynamic_assumption')._domain = \
        In([ThermodynamicAssumption.adiabatic])
    config.declare("valve_function", ConfigValue(
        default=ValveFunctionType.linear,
        domain=In(ValveFunctionType),
        description="Valve function type, if custom provide an expression rule",
        doc="""The type of valve function, if custom provide an expression rule
with the valve_function_rule argument.
**default** - ValveFunctionType.linear
**Valid values** - {
ValveFunctionType.linear,
ValveFunctionType.quick_opening,
ValveFunctionType.equal_percentage,
ValveFunctionType.custom}"""))
    config.declare("valve_function_rule", ConfigValue(
        default=None,
        description="This is a rule that returns a time indexed valve function expression.",
        doc="""This is a rule that returns a time indexed valve function expression.
This is required only if valve_function==ValveFunctionType.custom"""))
    config.declare("phase", ConfigValue(
        default="Vap",
        domain=In(("Vap", "Liq")),
        description='Expected phase of fluid in valve in {"Liq", "Vap"}'))
