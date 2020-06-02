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
Basic heater/cooler models
"""

__author__ = "John Eslick"

# Import Pyomo libraries
from pyomo.environ import Reference
from pyomo.common.config import ConfigValue, In

# Import IDAES cores
from idaes.core import declare_process_block_class
from idaes.generic_models.unit_models.balance import BalanceBlockData
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def _make_heater_config_block(config):
    """
    Declare configuration options for HeaterData block.  In addtion to the
    balance block config options the heat and work transfer options are fixed, so
    there is heat transfer is on and work transfer is off. Having the options
    allows the a user to easily deternin wheather the haet and work transfer
    terms exist in a model.
    """
    config.declare("has_work_transfer", ConfigValue(
            default=False,
            domain=In([False]),
            description="Heater does not have work transfer term.",
        )
    )
    config.declare("has_heat_transfer", ConfigValue(
            default=True,
            domain=In([True]),
            description="Heater has heat transfer term.",
        )
    )


@declare_process_block_class("Heater", doc="Simple 0D heater/cooler model.")
class HeaterData(BalanceBlockData):
    """
    Simple 0D heater unit.
    Unit model to add or remove heat from a material.
    """
    CONFIG = BalanceBlockData.CONFIG()
    _make_heater_config_block(CONFIG)

    def build(self):
        """Building model

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super().build() # Add control volume, balance equations, and ports
        self.heat_duty = Reference(self.control_volume.heat)


    def _get_performance_contents(self, time_point=0):
        return {"vars":{"Heat Duty":self.heat_duty[time_point]}}
