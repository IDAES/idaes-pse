#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
IDAES Phase objects

Created on Tue Feb 18 10:54:52 2020

@author: alee
"""
from enum import Enum

from pyomo.environ import Set
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.core.base.process_base import declare_process_block_class, ProcessBlockData


# Enumerate recognised Phase types
class PhaseType(Enum):
    undefined = 0
    liquidPhase = 1
    vaporPhase = 2
    solidPhase = 3
    aqueousPhase = 4


# TODO: Document EoS options and parameter_Data
@declare_process_block_class("Phase")
class PhaseData(ProcessBlockData):
    CONFIG = ConfigBlock()
    CONFIG.declare(
        "component_list",
        ConfigValue(
            default=None,
            domain=list,
            description="List of components in phase",
            doc="List of components which are present in phase. This is used "
            "to construct the phase-component Set for the property package.",
        ),
    )
    CONFIG.declare(
        "equation_of_state",
        ConfigValue(
            default=None,
            description="Equation of state for phase",
            doc="""A valid Python class with the necessary methods for
                constructing the desired equation of state (or similar
                model).""",
        ),
    )
    CONFIG.declare(
        "equation_of_state_options",
        ConfigValue(
            default=None,
            description="Options for equation of state",
            doc="""A dict or ConfigBlock of options to be used when setting
                up equation of state for phase.""",
        ),
    )
    CONFIG.declare(
        "parameter_data",
        ConfigValue(
            default={},
            domain=dict,
            description="Dict containing initialization data for parameters",
        ),
    )
    CONFIG.declare(
        "_phase_list_exists",
        ConfigValue(
            default=False,
            doc="Internal config argument indicating whether phase_list "
            "needs to be populated.",
        ),
    )

    CONFIG.declare(
        "therm_cond_phase",
        ConfigValue(description="Method to calculate thermal conductivity phase"),
    )
    CONFIG.declare(
        "surf_tens_phase",
        ConfigValue(description="Method to calculate surface tension of phase"),
    )
    CONFIG.declare(
        "visc_d_phase",
        ConfigValue(description="Method to calculate dynamic viscosity of phase"),
    )

    def build(self):
        super(PhaseData, self).build()

        # If the phase_list does not exist, add a reference to the new Phase
        # The IF is mostly for backwards compatability, to allow for old-style
        # property packages where the phase_list already exists but we need to
        # add new Phase objects
        if not self.config._phase_list_exists:
            self.__add_to_phase_list()

    # For the base Phase class, determine phase type based on component name
    # Derived classes will overload these and return the correct type
    # This will handle backwards compatability for old-style property packages
    def is_liquid_phase(self):
        if "Liq" in self.name:
            return True
        else:
            return False

    def is_solid_phase(self):
        if "Sol" in self.name:
            return True
        else:
            return False

    def is_vapor_phase(self):
        if "Vap" in self.name:
            return True
        else:
            return False

    def is_aqueous_phase(self):
        # Returns bool indicating if this phase involve electrolytes
        return False

    def __add_to_phase_list(self):
        """
        Method to add reference to new Phase in phase_list
        """
        parent = self.parent_block()
        try:
            phase_list = getattr(parent, "phase_list")
            phase_list.add(self.local_name)
        except AttributeError:
            # Parent does not have a phase_list yet, so create one
            parent.phase_list = Set(initialize=[self.local_name], ordered=True)


@declare_process_block_class("LiquidPhase", block_class=Phase)
class LiquidPhaseData(PhaseData):
    def is_liquid_phase(self):
        return True

    def is_solid_phase(self):
        return False

    def is_vapor_phase(self):
        return False


@declare_process_block_class("SolidPhase", block_class=Phase)
class SolidPhaseData(PhaseData):
    def is_liquid_phase(self):
        return False

    def is_solid_phase(self):
        return True

    def is_vapor_phase(self):
        return False


@declare_process_block_class("VaporPhase", block_class=Phase)
class VaporPhaseData(PhaseData):
    def is_liquid_phase(self):
        return False

    def is_solid_phase(self):
        return False

    def is_vapor_phase(self):
        return True


@declare_process_block_class("AqueousPhase", block_class=LiquidPhase)
class AqueousPhaseData(LiquidPhaseData):
    # Special phase type for liquid phases involving electrolytes
    # This is used to determine if we need to do the more complex component
    # list determinations
    def is_aqueous_phase(self):
        return True


# List of all Phase types to use for validation
__all_phases__ = [Phase, LiquidPhase, SolidPhase, VaporPhase, AqueousPhase]
