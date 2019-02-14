"""
Define configuration block for the SteamValve model.
"""


__author__ = "John Eslick"

from pyomo.common.config import ConfigBlock, ConfigValue, In, ConfigList
from pyutilib.enum import Enum

ValveFunctionType = Enum("linear", "quick_opening", "equal_percentage", "custom")

def _define_config(config):
    config.compressor = False
    config.get('compressor')._default = False
    config.get('compressor')._domain = In([False])
    config.thermodynamic_assumption = 'adiabatic'
    config.get('thermodynamic_assumption')._default = 'adiabatic'
    config.get('thermodynamic_assumption')._domain = In(['adiabatic'])
    config.declare("valve_function", ConfigValue(
        default=ValveFunctionType.linear,
        domain=In(ValveFunctionType),
        description="Valve function type, if custom provide an expression rule",
        doc="""The type of valve function, if custom provide an expression rule
with the valve_function_rule argument.
**default** - ValveFunctionType.linear
**Valid values - ** {
ValveFunctionType.linear,
ValveFunctionType.quick_opening,
ValveFunctionType.equal_percentage,
ValveFunctionType.custom}"""))
    config.declare("valve_function_rule", ConfigValue(
        default=None,
        description="Rule with time index returns valve function expression",
        doc="This is a rule with time index that returns an expression for the "
            "this is required only if valve_function==ValveFunctionType.custom"))
    config.declare("phase", ConfigValue(
        default="Vap",
        domain=In(("Vap", "Liq")),
        description='Expected phase of fluid in valve in {"Liq", "Vap"}'))
