__author__ = "John Eslick"

from pyutilib.enum import Enum
from pyomo.common.config import ConfigBlock, ConfigValue, In
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
    config.declare("hp_split_loc")
