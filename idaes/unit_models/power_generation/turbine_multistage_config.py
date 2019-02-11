__author__ = "John Eslick"

from pyutilib.enum import Enum
from pyomo.common.config import ConfigBlock, ConfigValue, In, ConfigList
from idaes.core import (EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block

_ReheatType = Enum('none', 'single', 'double')

def _define_turbine_mutlistage_config(config):
    config.declare("dynamic", ConfigValue(
        domain=In([True, False]),
        default=False,
        description="Dynamic model flag",
        doc="Indicates whether the model is dynamic."))
    config.declare("has_holdup", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
    config.declare("has_phase_equilibrium", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Calculate phase equilibrium in mixed stream",
        doc="""Argument indicating whether phase equilibrium should be
calculated for the resulting mixed stream,
**default** - False.
**Valid values:** {
**True** - calculate phase equilibrium in mixed stream,
**False** - do not calculate equilibrium in mixed stream.}"""))
    config.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    config.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))
    config.declare("num_hp", ConfigValue(
        default=2,
        domain=int,
        description="Number of high pressure stages not including inlet stage",
        doc="Number of high pressure stages not including inlet stage"))
    config.declare("num_ip", ConfigValue(
        default=10,
        domain=int,
        description="Number of intermediate pressure stages",
        doc="Number of intermediate pressure stages"))
    config.declare("num_lp", ConfigValue(
        default=5,
        domain=int,
        description="Number of low pressure stages not including outlet stage",
        doc="Number of low pressure stages not including outlet stage"))
    config.declare("hp_split_locations", ConfigList(
        default=[],
        domain=int,
        description="Locations of splitters in HP section",
        doc="A list of index locations of splitters in the HP section. The "
            "indexes indicate after which stage to include splitters.  0 is "
            "between the inlet stage and the first regular HP stage."))
    config.declare("ip_split_locations", ConfigList(
        default=[],
        domain=int,
        description="Locations of splitters in IP section",
        doc="A list of index locations of splitters in the IP section. The "
            "indexes indicate after which stage to include splitters."))
    config.declare("lp_split_locations", ConfigList(
        default=[],
        domain=int,
        description="Locations of splitter in LP section",
        doc="A list of index locations of splitters in the LP section. The "
            "indexes indicate after which stage to include splitters."))
    config.declare("hp_disconnect", ConfigList(
        default=[],
        domain=int,
        description="HP Turbine stages to not connect to next with an arc.",
        doc="HP Turbine stages to not connect to next with an arc. This is "
            "usually used to insert addtional units between stages on a "
            "flowsheet, such as a reheater"))
    config.declare("ip_disconnect", ConfigList(
        default=[],
        domain=int,
        description="IP Turbine stages to not connect to next with an arc.",
        doc="IP Turbine stages to not connect to next with an arc. This is "
            "usually used to insert addtional units between stages on a "
            "flowsheet, such as a reheater"))
    config.declare("lp_disconnect", ConfigList(
        default=[],
        domain=int,
        description="LP Turbine stages to not connect to next with an arc.",
        doc="LP Turbine stages to not connect to next with an arc. This is "
            "usually used to insert addtional units between stages on a "
            "flowsheet, such as a reheater"))
    config.declare("hp_split_num_outlets", ConfigValue(
        default={},
        domain=dict,
        description="Dict, hp split index: number of splitter outlets"))
    config.declare("ip_split_num_outlets", ConfigValue(
        default={},
        domain=dict,
        description="Dict, ip split index: number of splitter outlets"))
    config.declare("lp_split_num_outlets", ConfigValue(
        default={},
        domain=dict,
        description="Dict, lp split index: number of splitter outlets"))
