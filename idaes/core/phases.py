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
IDAES Phase objects

Created on Tue Feb 18 10:54:52 2020

@author: alee
"""
from pyomo.common.config import ConfigBlock, ConfigValue

from .process_base import (declare_process_block_class,
                           ProcessBlockData)


@declare_process_block_class("Phase")
class PhaseData(ProcessBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare("component_list", ConfigValue(
            default=None,
            domain=list,
            description="List of components in phase",
            doc="List of components which are present in phase. This is used "
            "to construct the phase-component Set for the property package."))

    def build(self):
        super(PhaseData, self).build()

    # For the base Phase class, determine phase type based on component name
    # Derived classes will overload these and return the correct type
    # This will handle backwards compatability for old-style property packages
    def is_liquid_phase(self):
        if "Liq" in self.name:
            return True
        else:
            return False

    def is_solid_phase(self):
        if "Liq" in self.name:
            return True
        else:
            return False

    def is_vapor_phase(self):
        if "Liq" in self.name:
            return True
        else:
            return False


@declare_process_block_class("LiquidPhase")
class LiquidPhaseData(PhaseData):
    def is_liquid_phase(self):
        return True

    def is_solid_phase(self):
        return False

    def is_vapor_phase(self):
        return False


@declare_process_block_class("SolidPhase")
class SolidPhaseData(PhaseData):
    def is_liquid_phase(self):
        return False

    def is_solid_phase(self):
        return True

    def is_vapor_phase(self):
        return False


@declare_process_block_class("VaporPhase")
class VaporPhaseData(PhaseData):
    def is_liquid_phase(self):
        return False

    def is_solid_phase(self):
        return False

    def is_vapor_phase(self):
        return True
