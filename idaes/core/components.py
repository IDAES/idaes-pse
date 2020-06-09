##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
IDAES Component objects

@author: alee
"""
from pyomo.environ import Set, Param, Var
from pyomo.common.config import ConfigBlock, ConfigValue

from .process_base import (declare_process_block_class,
                           ProcessBlockData)
from .phases import PhaseType as PT
from .util.config import list_of_phase_types
from .util.exceptions import ConfigurationError


@declare_process_block_class("Component")
class ComponentData(ProcessBlockData):
    CONFIG = ConfigBlock()

    CONFIG.declare("valid_phase_types", ConfigValue(
            domain=list_of_phase_types,
            doc="List of valid PhaseTypes (Enums) for this Component."))

    CONFIG.declare("elemental_composition", ConfigValue(
            domain=dict,
            description="Elemental composition of component",
            doc="Dict containing elemental composition in the form element "
                ": stoichiometry"))

    CONFIG.declare("henry_component", ConfigValue(
            domain=dict,
            description="Phases in which component follows Henry's Law",
            doc="Dict indicating phases in which component follows Herny's "
                "Law (keys) with values indicating form of law."))

    CONFIG.declare("dens_mol_liq_comp", ConfigValue(
        description="Method to use to calculate liquid phase molar density"))
    CONFIG.declare("enth_mol_liq_comp", ConfigValue(
        description="Method to calculate liquid component molar enthalpies"))
    CONFIG.declare("enth_mol_ig_comp", ConfigValue(
        description="Method to calculate ideal gas component molar enthalpies"
        ))
    CONFIG.declare("entr_mol_liq_comp", ConfigValue(
        description="Method to calculate liquid component molar entropies"))
    CONFIG.declare("entr_mol_ig_comp", ConfigValue(
        description="Method to calculate ideal gas component molar entropies"))
    CONFIG.declare("pressure_sat_comp", ConfigValue(
        description="Method to use to calculate saturation pressure"))

    CONFIG.declare("phase_equilibrium_form", ConfigValue(
        domain=dict,
        description="Form of phase equilibrium constraints for component"))

    CONFIG.declare("parameter_data", ConfigValue(
        default={},
        domain=dict,
        description="Dict containing initialization data for parameters"))

    CONFIG.declare("_component_list_exists", ConfigValue(
            default=False,
            doc="Internal config argument indicating whether component_list "
            "needs to be populated."))

    def build(self):
        super(ComponentData, self).build()

        # If the component_list does not exist, add reference to new Component
        # The IF is mostly for backwards compatability, to allow for old-style
        # property packages where the component_list already exists but we
        # need to add new Component objects
        if not self.config._component_list_exists:
            self.__add_to_component_list()

        # Create Param for molecular weight if provided
        if "mw" in self.config.parameter_data:
            self.mw = Param(initialize=self.config.parameter_data["mw"])

        # Create Vars for common parameters
        for p in ["pressure_crit", "temperature_crit", "omega"]:
            if p in self.config.parameter_data:
                self.add_component(p, Var(
                    initialize=self.config.parameter_data[p]))

    def is_solute(self):
        raise TypeError(
            "{} Generic Component objects do not support is_solute() method. "
            "Use a Solvent or Solute Component instead."
            .format(self.name))

    def is_solvent(self):
        raise TypeError(
            "{} Generic Component objects do not support is_solvent() method. "
            "Use a Solvent or Solute Component instead."
            .format(self.name))

    def __add_to_component_list(self):
        """
        Method to add reference to new Component in component_list
        """
        parent = self.parent_block()
        try:
            comp_list = getattr(parent, "component_list")
            comp_list.add(self.local_name)
        except AttributeError:
            # Parent does not have a component_list yet, so create one
            parent.component_list = Set(initialize=[self.local_name],
                                        ordered=True)

    def _is_phase_valid(self, phase):
        # If no valid phases assigned, assume all are valid
        if self.config.valid_phase_types is None:
            return True

        # Check for behaviour of phase, and see if that is a valid behaviour
        # for component.
        if (phase.is_liquid_phase() and
                PT.liquidPhase in self.config.valid_phase_types):
            return True
        elif (phase.is_vapor_phase() and
                PT.vaporPhase in self.config.valid_phase_types):
            return True
        elif (phase.is_solid_phase() and
                PT.solidPhase in self.config.valid_phase_types):
            return True
        else:
            return False


# TODO : What about LLE systems where a species is a solvent in one liquid
# phase, but a solute in another?
@declare_process_block_class("Solute")
class SoluteData(ComponentData):
    """
    Component type for species which should be considered as solutes in
    LiquidPhases.
    """

    def is_solute(self):
        return True

    def is_solvent(self):
        return False


# TODO : What about LLE systems where a species is a solvent in one liquid
# phase, but a solute in another?
@declare_process_block_class("Solvent")
class SolventData(ComponentData):
    """
    Component type for species which should be considered as solvents in
    LiquidPhases.
    """

    def is_solute(self):
        return False

    def is_solvent(self):
        return True


@declare_process_block_class("Ion")
class IonData(SoluteData):
    """
    Component type for ionic species. These can exist only in LiquidPhases,
    and are always solutes.
    """
    CONFIG = SoluteData.CONFIG()

    # Remove valid_phase_types argument, as ions are liquid phase only
    CONFIG.__delitem__("valid_phase_types")

    CONFIG.declare("charge", ConfigValue(
            domain=int,
            doc="Charge of ionic species."))

    def _is_phase_valid(self, phase):
        return phase.is_liquid_phase()


@declare_process_block_class("Anion")
class AnionData(IonData):
    """
    Component type for anionic species. These can exist only in LiquidPhases,
    and are always solutes.
    """
    CONFIG = IonData.CONFIG()

    def build(self):
        super().build()

        # Validate charge config argument
        if self.config.charge is None:
            raise ConfigurationError(
                "{} was not provided with a value for charge."
                .format(self.name))
        elif self.config.charge >= 0:
            raise ConfigurationError(
                "{} received invalid value for charge configuration argument."
                " Anions must have a negative charge.".format(self.name))


@declare_process_block_class("Cation")
class CationData(IonData):
    """
    Component type for cationic species. These can exist only in LiquidPhases,
    and are always solutes.
    """
    CONFIG = IonData.CONFIG()

    def build(self):
        super().build()

        # Validate charge config argument
        if self.config.charge is None:
            raise ConfigurationError(
                "{} was not provided with a value for charge."
                .format(self.name))
        elif self.config.charge <= 0:
            raise ConfigurationError(
                "{} received invalid value for charge configuration argument."
                " Cations must have a positive charge.".format(self.name))


__all_components__ = [Component, Solute, Solvent, Ion, Anion, Cation]
