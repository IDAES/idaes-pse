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
Function to make 0D feedwater heater model config blocks.
"""

__author__ = "John Eslick"

from pyomo.common.config import ConfigValue, In, ConfigBlock
from idaes.unit_models import MomentumMixingType
from idaes.unit_models.heat_exchanger import HeatExchangerData
from idaes.core import useDefault
from idaes.core.util.config import is_physical_parameter_block

def _define_feedwater_heater_0D_config(config):
    config.declare("has_drain_mixer", ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Add a mixer to the inlet of the condensing section",
            doc="""Add a mixer to the inlet of the condensing section to add
water from the drain of another feedwaterheater to the steam, if True"""))
    config.declare("has_desuperheat", ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Add a mixer desuperheat section to the heat exchanger",
            doc="Add a mixer desuperheat section to the heat exchanger"))
    config.declare("has_drain_cooling", ConfigValue(
            default=True,
            domain=In([True, False]),
            description="Add a section after condensing section cool condensate.",
            doc="Add a section after condensing section to cool condensate."))
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
    config.declare("condense", HeatExchangerData.CONFIG())
    config.declare("desuperheat", HeatExchangerData.CONFIG())
    config.declare("cooling", HeatExchangerData.CONFIG())
