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
Base class for control volumes
"""

from __future__ import division

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Constraint, Param, Reals, Var
from pyomo.dae import DerivativeVar

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        ControlVolumeBlockData,
                        FlowDirection,
                        MaterialFlowBasis)
from idaes.core.util.exceptions import (BalanceTypeNotSupportedError,
                                        ConfigurationError,
                                        PropertyNotSupportedError,
                                        PropertyPackageError)
from idaes.core.util.misc import add_object_reference

__author__ = "Andrew Lee"


_log = logging.getLogger(__name__)

# TODO : Custom terms in material balances, other types of material balances
# TODO : Improve flexibility for get_material_flow_terms and associated


@declare_process_block_class("ControlVolume0DBlock", doc="""
    ControlVolume0DBlock is a specialized Pyomo block for IDAES non-discretized
    control volume blocks, and contains instances of ControlVolume0DBlockData.

    ControlVolume0DBlock should be used for any control volume with a defined volume
    and distinct inlets and outlets which does not require spatial
    discretization. This encompases most basic unit models used in process
    modeling.""")
class ControlVolume0DBlockData(ControlVolumeBlockData):
    """
    0-Dimensional (Non-Discretised) ControlVolume Class

    This class forms the core of all non-discretized IDAES models. It provides
    methods to build property and reaction blocks, and add mass, energy and
    momentum balances. The form of the terms used in these constraints is
    specified in the chosen property package.
    """
    def build(self):
        """
        Build method for ControlVolume0DBlock blocks.

        Returns:
            None
        """
        # Call build method from base class
        super(ControlVolume0DBlockData, self).build()

    def add_geometry(self):
        """
        Method to create volume Var in ControlVolume.

        Args:
            None

        Returns:
            None
        """
        l_units = self.config.property_package.get_metadata().default_units[
                                                                      "length"]
        self.volume = Var(self.time_ref, initialize=1.0,
                          doc='Holdup Volume [{}^3]'.format(l_units))

    def add_state_blocks(self,
                         information_flow=FlowDirection.forward,
                         has_phase_equilibrium=None):
        """
        This method constructs the inlet and outlet state blocks for the
        control volume.

        Args:
            information_flow: a FlowDirection Enum indicating whether
                               information flows from inlet-to-outlet or
                               outlet-to-inlet
            has_phase_equilibrium: indicates whether equilibrium calculations
                                    will be required in state blocks
            package_arguments: dict-like object of arguments to be passed to
                                state blocks as construction arguments
        Returns:
            None
        """
        if has_phase_equilibrium is None:
            raise ConfigurationError(
                    "{} add_state_blocks method was not provided with a "
                    "has_phase_equilibrium_argument.".format(self.name))
        elif has_phase_equilibrium not in [True, False]:
            raise ConfigurationError(
                    "{} add_state_blocks method was provided with an invalid "
                    "has_phase_equilibrium_argument. Must be True or False"
                    .format(self.name))

        tmp_dict = dict(**self.config.property_package_args)
        tmp_dict["has_phase_equilibrium"] = has_phase_equilibrium
        tmp_dict["parameters"] = self.config.property_package

        if information_flow == FlowDirection.forward:
            tmp_dict["defined_state"] = True
        elif information_flow == FlowDirection.backward:
            tmp_dict["defined_state"] = False
        else:
            raise ConfigurationError(
                    '{} invalid value for information_flow argument. '
                    'Valid values are FlowDirection.forward and '
                    'FlowDirection.backward'.format(self.name))

        try:
            self.properties_in = (
                    self.config.property_package.state_block_class(
                            self.time_ref,
                            doc="Material properties at inlet",
                            default=tmp_dict))

            # Reverse defined_state
            tmp_dict["defined_state"] = not tmp_dict["defined_state"]

            self.properties_out = (
                    self.config.property_package.state_block_class(
                            self.time_ref,
                            doc="Material properties at outlet",
                            default=tmp_dict))
        except AttributeError:
            raise PropertyPackageError(
                    "{} physical property package has not implemented the "
                    "state_block_class attribute. Please contact the "
                    "developer of the physical property package."
                    .format(self.name))

    def add_reaction_blocks(self, has_equilibrium=None):
        """
        This method constructs the reaction block for the control volume.

        Args:
            has_equilibrium: indicates whether equilibrium calculations
                              will be required in reaction block
            package_arguments: dict-like object of arguments to be passed to
                                reaction block as construction arguments

        Returns:
            None
        """
        if has_equilibrium is None:
            raise ConfigurationError(
                    "{} add_reaction_blocks method was not provided with a "
                    "has_equilibrium_argument.".format(self.name))
        elif has_equilibrium not in [True, False]:
            raise ConfigurationError(
                    "{} add_reaction_blocks method was provided with an "
                    "invalid has_equilibrium_argument. Must be True or False"
                    .format(self.name))

        tmp_dict = dict(**self.config.reaction_package_args)
        tmp_dict["state_block"] = self.properties_out
        tmp_dict["has_equilibrium"] = has_equilibrium
        tmp_dict["parameters"] = self.config.reaction_package

        try:
            self.reactions = (
                    self.config.reaction_package.reaction_block_class(
                            self.time_ref,
                            doc="Reaction properties in control volume",
                            default=tmp_dict))
        except AttributeError:
            raise PropertyPackageError(
                    "{} reaction property package has not implemented the "
                    "reaction_block_class attribute. Please contact the "
                    "developer of the reaction property package."
                    .format(self.name))

    def add_phase_component_balances(self,
                                     has_rate_reactions=False,
                                     has_equilibrium_reactions=False,
                                     has_phase_equilibrium=False,
                                     has_mass_transfer=False,
                                     custom_molar_term=None,
                                     custom_mass_term=None):
        """
        This method constructs a set of 0D material balances indexed by time,
        phase and component.

        Args:
            has_rate_reactions: whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions: whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium: whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer: whether generic mass transfer terms should be
                    included in material balances
            custom_molar_term: a Pyomo Expression representing custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, phase list and
                    component list
            custom_mass_term: a Pyomo Expression representing custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, phase list and
                    component list

        Returns:
            Constraint object representing material balances
        """
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        # Check that reaction block exists if required
        if has_rate_reactions or has_equilibrium_reactions:
            try:
                rblock = self.reactions
            except AttributeError:
                raise ConfigurationError(
                        "{} does not contain a Reaction Block, but material "
                        "balances have been set to contain reaction terms. "
                        "Please construct a reaction block before adding "
                        "balance equations.".format(self.name))

        if has_equilibrium_reactions:
            # Check that reaction block is set to calculate equilibrium
            for t in self.time_ref:
                if self.reactions[t].config.has_equilibrium is False:
                    raise ConfigurationError(
                            "{} material balance was set to include "
                            "equilibrium reactions, however the associated "
                            "ReactionBlock was not set to include equilibrium "
                            "constraints (has_equilibrium_reactions=False). "
                            "Please correct your configuration arguments."
                            .format(self.name))

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            for t in self.time_ref:
                if not self.properties_out[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the associated outlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                if not self.properties_in[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the associated inlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))

        # Get units from property package
        units = {}
        for u in ['length', 'holdup', 'amount', 'time']:
            try:
                units[u] = \
                   self.config.property_package.get_metadata().default_units[u]
            except KeyError:
                units[u] = '-'

        # Test for components that must exist prior to calling this method
        if has_holdup:
            if not hasattr(self, "volume"):
                raise ConfigurationError(
                        "{} control volume must have volume defined to have "
                        "holdup and/or rate reaction terms. Please call the "
                        "add_geometry method before adding balance equations."
                        .format(self.name))

        # Material holdup and accumulation
        if has_holdup:
            self.material_holdup = Var(self.time_ref,
                                       self.phase_list_ref,
                                       self.component_list_ref,
                                       domain=Reals,
                                       doc="Material holdup in unit [{}]"
                                           .format(units['holdup']))
        if dynamic:
            self.material_accumulation = DerivativeVar(
                    self.material_holdup,
                    wrt=self.time_ref,
                    doc="Material accumulation in unit [{}/{}]"
                        .format(units['holdup'], units['time']))

        # Get phase component list(s)
        phase_component_list = self._get_phase_comp_list()

        # Create material balance terms as required
        # Kinetic reaction generation
        if has_rate_reactions:
            try:
                add_object_reference(
                        self,
                        "rate_reaction_idx_ref",
                        self.config.reaction_package.rate_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of rate "
                    "reactions (rate_reaction_idx), thus does not support "
                    "rate-based reactions.".format(self.name))
            self.rate_reaction_generation = Var(
                        self.time_ref,
                        self.phase_list_ref,
                        self.component_list_ref,
                        domain=Reals,
                        doc="Amount of component generated in "
                            "unit by kinetic reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Equilibrium reaction generation
        if has_equilibrium_reactions:
            try:
                add_object_reference(
                    self,
                    "equilibrium_reaction_idx_ref",
                    self.config.reaction_package.equilibrium_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of "
                    "equilibrium reactions (equilibrium_reaction_idx), thus "
                    "does not support equilibrium-based reactions."
                    .format(self.name))
            self.equilibrium_reaction_generation = Var(
                        self.time_ref,
                        self.phase_list_ref,
                        self.component_list_ref,
                        domain=Reals,
                        doc="Amount of component generated in unit "
                            "by equilibrium reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Phase equilibrium generation
        if has_phase_equilibrium:
            try:
                add_object_reference(
                    self,
                    "phase_equilibrium_idx_ref",
                    self.config.property_package.phase_equilibrium_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Property package does not contain a list of phase "
                    "equilibrium reactions (phase_equilibrium_idx), thus does "
                    "not support phase equilibrium.".format(self.name))
            self.phase_equilibrium_generation = Var(
                        self.time_ref,
                        self.phase_equilibrium_idx_ref,
                        domain=Reals,
                        doc="Amount of generation in unit by phase "
                            "equilibria [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Material transfer term
        if has_mass_transfer:
            self.mass_transfer_term = Var(
                        self.time_ref,
                        self.phase_list_ref,
                        self.component_list_ref,
                        domain=Reals,
                        doc="Component material transfer into unit [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Create rules to substitute material balance terms
        # Accumulation term
        def accumulation_term(b, t, p, j):
            return b.material_accumulation[t, p, j] if dynamic else 0

        def kinetic_term(b, t, p, j):
            return (b.rate_reaction_generation[t, p, j] if has_rate_reactions
                    else 0)

        def equilibrium_term(b, t, p, j):
            return (b.equilibrium_reaction_generation[t, p, j]
                    if has_equilibrium_reactions else 0)

        def phase_equilibrium_term(b, t, p, j):
            if has_phase_equilibrium:
                sd = {}
                sblock = self.properties_out[t]
                for r in b.phase_equilibrium_idx_ref:
                    if sblock.phase_equilibrium_list[r][0] == j:
                        if sblock.phase_equilibrium_list[r][1][0] == p:
                            sd[r] = 1
                        elif sblock.phase_equilibrium_list[r][1][1] == p:
                            sd[r] = -1
                        else:
                            sd[r] = 0
                    else:
                        sd[r] = 0

                return sum(b.phase_equilibrium_generation[t, r]*sd[r]
                           for r in b.phase_equilibrium_idx_ref)
            else:
                return 0

        def transfer_term(b, t, p, j):
            return (b.mass_transfer_term[t, p, j] if has_mass_transfer else 0)

        def user_term_mol(b, t, p, j):
            if custom_molar_term is not None:
                flow_basis = b.properties_out[t].get_material_flow_basis()
                if flow_basis == MaterialFlowBasis.molar:
                    return custom_molar_term(t, p, j)
                elif flow_basis == MaterialFlowBasis.mass:
                    try:
                        return (custom_molar_term(t, p, j) *
                                b.properties_out[t].mw[j])
                    except AttributeError:
                        raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances."
                                .format(self.name))
                else:
                    raise ConfigurationError(
                            "{} contained a custom_molar_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name))
            else:
                return 0

        def user_term_mass(b, t, p, j):
            if custom_mass_term is not None:
                flow_basis = b.properties_out[t].get_material_flow_basis()
                if flow_basis == MaterialFlowBasis.mass:
                    return custom_mass_term(t, p, j)
                elif flow_basis == MaterialFlowBasis.molar:
                    try:
                        return (custom_mass_term(t, p, j) /
                                b.properties_out[t].mw[j])
                    except AttributeError:
                        raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances."
                                .format(self.name))
                else:
                    raise ConfigurationError(
                            "{} contained a custom_mass_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name))
            else:
                return 0

        # Add component balances
        @self.Constraint(self.time_ref,
                         self.phase_list_ref,
                         self.component_list_ref,
                         doc="Material balances")
        def material_balances(b, t, p, j):
            if j in phase_component_list[p]:
                return accumulation_term(b, t, p, j) == (
                        b.properties_in[t].get_material_flow_terms(p, j) -
                        b.properties_out[t].get_material_flow_terms(p, j) +
                        kinetic_term(b, t, p, j) +
                        equilibrium_term(b, t, p, j) +
                        phase_equilibrium_term(b, t, p, j) +
                        transfer_term(b, t, p, j) +
                        user_term_mol(b, t, p, j) + user_term_mass(b, t, p, j))
            else:
                return Constraint.Skip

        # TODO: Need to set material_holdup = 0 for non-present component-phase
        # pairs. Not ideal, but needed to close DoF. Is there a better way?

        # Material Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time_ref,
                             self.phase_list_ref,
                             self.component_list_ref,
                             doc="Material holdup calculations")
            def material_holdup_calculation(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.material_holdup[t, p, j] == (
                          b.volume[t]*self.phase_fraction[t, p] *
                          b.properties_out[t].get_material_density_terms(p, j))
                else:
                    return b.material_holdup[t, p, j] == 0

        if has_rate_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.rate_reaction_extent = Var(
                    self.time_ref,
                    self.rate_reaction_idx_ref,
                    domain=Reals,
                    doc="Extent of kinetic reactions[{}/{}]"
                        .format(units['holdup'], units['time']))

            @self.Constraint(self.time_ref,
                             self.phase_list_ref,
                             self.component_list_ref,
                             doc="Kinetic reaction stoichiometry constraint")
            def rate_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.rate_reaction_generation[t, p, j] == (
                        sum(rblock[t].rate_reaction_stoichiometry[r, p, j] *
                            b.rate_reaction_extent[t, r]
                            for r in b.rate_reaction_idx_ref))
                else:
                    return Constraint.Skip

        if has_equilibrium_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.equilibrium_reaction_extent = Var(
                            self.time_ref,
                            self.equilibrium_reaction_idx_ref,
                            domain=Reals,
                            doc="Extent of equilibrium reactions[{}/{}]"
                                .format(units['holdup'], units['time']))

            @self.Constraint(self.time_ref,
                             self.phase_list_ref,
                             self.component_list_ref,
                             doc="Equilibrium reaction stoichiometry")
            def equilibrium_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.equilibrium_reaction_generation[t, p, j] == (
                            sum(rblock[t].
                                equilibrium_reaction_stoichiometry[r, p, j] *
                                b.equilibrium_reaction_extent[t, r]
                                for r in b.equilibrium_reaction_idx_ref))
                else:
                    return Constraint.Skip

        return self.material_balances

    def add_total_component_balances(self,
                                     has_rate_reactions=False,
                                     has_equilibrium_reactions=False,
                                     has_phase_equilibrium=False,
                                     has_mass_transfer=False,
                                     custom_molar_term=None,
                                     custom_mass_term=None):
        """
        This method constructs a set of 0D material balances indexed by time
        and component.

        Args:
            has_rate_reactions - whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions - whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium - whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer - whether generic mass transfer terms should be
                    included in material balances
            custom_molar_term - a Pyomo Expression representing custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, phase list and
                    component list
            custom_mass_term - a Pyomo Expression representing custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, phase list and
                    component list

        Returns:
            Constraint object representing material balances
        """
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        # Check that reaction block exists if required
        if has_rate_reactions or has_equilibrium_reactions:
            try:
                rblock = self.reactions
            except AttributeError:
                raise ConfigurationError(
                        "{} does not contain a Reaction Block, but material "
                        "balances have been set to contain reaction terms. "
                        "Please construct a reaction block before adding "
                        "balance equations.".format(self.name))

        if has_equilibrium_reactions:
            # Check that reaction block is set to calculate equilibrium
            for t in self.time_ref:
                if self.reactions[t].config.has_equilibrium is False:
                    raise ConfigurationError(
                            "{} material balance was set to include "
                            "equilibrium reactions, however the associated "
                            "ReactionBlock was not set to include equilibrium "
                            "constraints (has_equilibrium_reactions=False). "
                            "Please correct your configuration arguments."
                            .format(self.name))

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            for t in self.time_ref:
                if not self.properties_out[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the associated outlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                if not self.properties_in[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the associated inlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))

        # Get units from property package
        units = {}
        for u in ['length', 'holdup', 'amount', 'time']:
            try:
                units[u] = \
                   self.config.property_package.get_metadata().default_units[u]
            except KeyError:
                units[u] = '-'

        # Test for components that must exist prior to calling this method
        if has_holdup:
            if not hasattr(self, "volume"):
                raise ConfigurationError(
                        "{} control volume must have volume defined to have "
                        "holdup and/or rate reaction terms. Please call the "
                        "add_geometry method before adding balance equations."
                        .format(self.name))

        # Material holdup and accumulation
        if has_holdup:
            self.material_holdup = Var(self.time_ref,
                                       self.phase_list_ref,
                                       self.component_list_ref,
                                       domain=Reals,
                                       doc="Material holdup in unit [{}]"
                                           .format(units['holdup']))
        if dynamic:
            self.material_accumulation = DerivativeVar(
                    self.material_holdup,
                    wrt=self.time_ref,
                    doc="Material accumulation in unit [{}/{}]"
                        .format(units['holdup'], units['time']))

        # Get phase component list(s)
        phase_component_list = self._get_phase_comp_list()

        # Create material balance terms as required
        # Kinetic reaction generation
        if has_rate_reactions:
            try:
                add_object_reference(
                        self,
                        "rate_reaction_idx_ref",
                        self.config.reaction_package.rate_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of rate "
                    "reactions (rate_reaction_idx), thus does not support "
                    "rate-based reactions.".format(self.name))
            self.rate_reaction_generation = Var(
                        self.time_ref,
                        self.phase_list_ref,
                        self.component_list_ref,
                        domain=Reals,
                        doc="Amount of component generated in "
                            "unit by kinetic reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Equilibrium reaction generation
        if has_equilibrium_reactions:
            try:
                add_object_reference(
                    self,
                    "equilibrium_reaction_idx_ref",
                    self.config.reaction_package.equilibrium_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of "
                    "equilibrium reactions (equilibrium_reaction_idx), thus "
                    "does not support equilibrium-based reactions."
                    .format(self.name))
            self.equilibrium_reaction_generation = Var(
                        self.time_ref,
                        self.phase_list_ref,
                        self.component_list_ref,
                        domain=Reals,
                        doc="Amount of component generated in unit "
                            "by equilibrium reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Material transfer term
        if has_mass_transfer:
            self.mass_transfer_term = Var(
                        self.time_ref,
                        self.phase_list_ref,
                        self.component_list_ref,
                        domain=Reals,
                        doc="Component material transfer into unit [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Create rules to substitute material balance terms
        # Accumulation term
        def accumulation_term(b, t, p, j):
            return b.material_accumulation[t, p, j] if dynamic else 0

        def kinetic_term(b, t, p, j):
            return (b.rate_reaction_generation[t, p, j] if has_rate_reactions
                    else 0)

        def equilibrium_term(b, t, p, j):
            return (b.equilibrium_reaction_generation[t, p, j]
                    if has_equilibrium_reactions else 0)

        def transfer_term(b, t, p, j):
            return (b.mass_transfer_term[t, p, j] if has_mass_transfer else 0)

        def user_term_mol(b, t, j):
            if custom_molar_term is not None:
                flow_basis = b.properties_out[t].get_material_flow_basis()
                if flow_basis == MaterialFlowBasis.molar:
                    return custom_molar_term(t, j)
                elif flow_basis == MaterialFlowBasis.mass:
                    try:
                        return (custom_molar_term(t, j) *
                                b.properties_out[t].mw[j])
                    except AttributeError:
                        raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances."
                                .format(self.name))
                else:
                    raise ConfigurationError(
                            "{} contained a custom_molar_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name))
            else:
                return 0

        def user_term_mass(b, t, j):
            if custom_mass_term is not None:
                flow_basis = b.properties_out[t].get_material_flow_basis()
                if flow_basis == MaterialFlowBasis.mass:
                    return custom_mass_term(t, j)
                elif flow_basis == MaterialFlowBasis.molar:
                    try:
                        return (custom_mass_term(t, j) /
                                b.properties_out[t].mw[j])
                    except AttributeError:
                        raise PropertyNotSupportedError(
                                "{} property package does not support "
                                "molecular weight (mw), which is required for "
                                "using custom terms in material balances."
                                .format(self.name))
                else:
                    raise ConfigurationError(
                            "{} contained a custom_mass_term argument, but "
                            "the property package used an undefined basis "
                            "(MaterialFlowBasis.other). Custom terms can "
                            "only be used when the property package declares "
                            "a molar or mass flow basis.".format(self.name))
            else:
                return 0

        # Add component balances
        @self.Constraint(self.time_ref,
                         self.component_list_ref,
                         doc="Material balances")
        def material_balances(b, t, j):
            cplist = []
            for p in self.phase_list_ref:
                if j in phase_component_list[p]:
                    cplist.append(p)
            return (
                sum(accumulation_term(b, t, p, j) for p in cplist) ==
                sum(b.properties_in[t].get_material_flow_terms(p, j)
                    for p in cplist) -
                sum(b.properties_out[t].get_material_flow_terms(p, j)
                    for p in cplist) +
                sum(kinetic_term(b, t, p, j) for p in cplist) +
                sum(equilibrium_term(b, t, p, j) for p in cplist) +
                sum(transfer_term(b, t, p, j) for p in cplist) +
                user_term_mol(b, t, j) + user_term_mass(b, t, j))

        # TODO: Need to set material_holdup = 0 for non-present component-phase
        # pairs. Not ideal, but needed to close DoF. Is there a better way?

        # Material Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time_ref,
                             self.phase_list_ref,
                             self.component_list_ref,
                             doc="Material holdup calculations")
            def material_holdup_calculation(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.material_holdup[t, p, j] == (
                          b.volume[t]*self.phase_fraction[t, p] *
                          b.properties_out[t].get_material_density_terms(p, j))
                else:
                    return b.material_holdup[t, p, j] == 0

        if has_rate_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.rate_reaction_extent = Var(
                    self.time_ref,
                    self.rate_reaction_idx_ref,
                    domain=Reals,
                    doc="Extent of kinetic reactions[{}/{}]"
                        .format(units['holdup'], units['time']))

            @self.Constraint(self.time_ref,
                             self.phase_list_ref,
                             self.component_list_ref,
                             doc="Kinetic reaction stoichiometry constraint")
            def rate_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.rate_reaction_generation[t, p, j] == (
                        sum(rblock[t].rate_reaction_stoichiometry[r, p, j] *
                            b.rate_reaction_extent[t, r]
                            for r in b.rate_reaction_idx_ref))
                else:
                    return Constraint.Skip

        if has_equilibrium_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.equilibrium_reaction_extent = Var(
                            self.time_ref,
                            self.equilibrium_reaction_idx_ref,
                            domain=Reals,
                            doc="Extent of equilibrium reactions[{}/{}]"
                                .format(units['holdup'], units['time']))

            @self.Constraint(self.time_ref,
                             self.phase_list_ref,
                             self.component_list_ref,
                             doc="Equilibrium reaction stoichiometry")
            def equilibrium_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.equilibrium_reaction_generation[t, p, j] == (
                            sum(rblock[t].
                                equilibrium_reaction_stoichiometry[r, p, j] *
                                b.equilibrium_reaction_extent[t, r]
                                for r in b.equilibrium_reaction_idx_ref))
                else:
                    return Constraint.Skip

        return self.material_balances

    def add_total_element_balances(self,
                                   has_rate_reactions=False,
                                   has_equilibrium_reactions=False,
                                   has_phase_equilibrium=False,
                                   has_mass_transfer=False,
                                   custom_elemental_term=None):
        """
        This method constructs a set of 0D element balances indexed by time.

        Args:
            has_rate_reactions - whether default generation terms for rate
                    reactions should be included in material balances
            has_equilibrium_reactions - whether generation terms should for
                    chemical equilibrium reactions should be included in
                    material balances
            has_phase_equilibrium - whether generation terms should for phase
                    equilibrium behaviour should be included in material
                    balances
            has_mass_transfer - whether generic mass transfer terms should be
                    included in material balances
            custom_elemental_term - a Pyomo Expression representing custom
                    terms to be included in material balances on a molar
                    elemental basis. Expression must be indexed by time and
                    element list

        Returns:
            Constraint object representing material balances
        """
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        if has_rate_reactions:
            raise ConfigurationError(
                    "{} add_total_element_balances method as provided with "
                    "argument has_rate_reactions = True. Total element "
                    "balances do not support rate based reactions, "
                    "please correct your configuration arguments"
                    .format(self.name))

        # Check that property package supports element balances
        try:
            add_object_reference(self,
                                 "element_list_ref",
                                 self.config.property_package.element_list)
        except AttributeError:
            raise PropertyNotSupportedError(
                    "{} property package provided does not contain a list of "
                    "elements (element_list), and thus does not support "
                    "elemental material balances. Please choose another type "
                    "of material balance or a property pakcage which supports "
                    "elemental balances.")

        # Check that reaction block exists if required
        if has_equilibrium_reactions:
            try:
                rblock = self.reactions
            except AttributeError:
                raise ConfigurationError(
                        "{} does not contain a Reaction Block, but material "
                        "balances have been set to contain reaction terms. "
                        "Please construct a reaction block before adding "
                        "balance equations.".format(self.name))

        if has_equilibrium_reactions:
            # Check that reaction block is set to calculate equilibrium
            for t in self.time_ref:
                if self.reactions[t].config.has_equilibrium is False:
                    raise ConfigurationError(
                            "{} material balance was set to include "
                            "equilibrium reactions, however the associated "
                            "ReactionBlock was not set to include equilibrium "
                            "constraints (has_equilibrium_reactions=False). "
                            "Please correct your configuration arguments."
                            .format(self.name))
                try:
                    add_object_reference(
                        self,
                        "equilibrium_reaction_idx_ref",
                        self.config.reaction_package.equilibrium_reaction_idx)
                except AttributeError:
                    raise PropertyNotSupportedError(
                        "{} Reaction package does not contain a list of "
                        "equilibrium reactions (equilibrium_reaction_idx), "
                        "thus does not support equilibrium-based reactions."
                        .format(self.name))

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            for t in self.time_ref:
                if not self.properties_out[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the associated outlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                if not self.properties_in[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the associated inlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                try:
                    add_object_reference(
                        self,
                        "phase_equilibrium_idx_ref",
                        self.config.property_package.phase_equilibrium_idx)
                except AttributeError:
                    raise PropertyNotSupportedError(
                        "{} Property package does not contain a list of phase "
                        "equilibrium reactions (phase_equilibrium_idx), thus "
                        "does not support phase equilibrium."
                        .format(self.name))

        # Test for components that must exist prior to calling this method
        if has_holdup:
            if not hasattr(self, "volume"):
                raise ConfigurationError(
                        "{} control volume must have volume defined to have "
                        "holdup terms. Please call the "
                        "add_geometry method before adding balance equations."
                        .format(self.name))

        # Get units from property package
        units = {}
        for u in ['amount', 'time']:
            try:
                units[u] = \
                   self.config.property_package.get_metadata().default_units[u]
            except KeyError:
                units[u] = '-'

        # Add Material Balance terms
        if has_holdup:
            self.element_holdup = Var(
                    self.time_ref,
                    self.element_list_ref,
                    domain=Reals,
                    doc="Elemental holdup in unit [{}]"
                        .format(units['amount']))

        if dynamic:
            self.element_accumulation = DerivativeVar(
                    self.element_holdup,
                    wrt=self.time_ref,
                    doc="Elemental accumulation in unit [{}/{}]"
                        .format(units['amount'], units['time']))

        @self.Expression(self.time_ref,
                         self.phase_list_ref,
                         self.element_list_ref,
                         doc="Inlet elemental flow terms [{}/{}]"
                             .format(units['amount'], units['time']))
        def elemental_flow_in(b, t, p, e):
            return sum(b.properties_in[t].get_material_flow_terms(p, j) *
                       b.properties_out[t].config.parameters.element_comp[j][e]
                       for j in b.component_list_ref)

        @self.Expression(self.time_ref,
                         self.phase_list_ref,
                         self.element_list_ref,
                         doc="Outlet elemental flow terms [{}/{}]"
                             .format(units['amount'], units['time']))
        def elemental_flow_out(b, t, p, e):
            return sum(b.properties_out[t].get_material_flow_terms(p, j) *
                       b.properties_out[t].config.parameters.element_comp[j][e]
                       for j in b.component_list_ref)

        # Create material balance terms as needed
        if has_mass_transfer:
            self.elemental_mass_transfer_term = Var(
                            self.time_ref,
                            self.element_list_ref,
                            domain=Reals,
                            doc="Element material transfer into unit [{}/{}]"
                            .format(units['amount'], units['time']))

        # Create rules to substitute material balance terms
        # Accumulation term
        def accumulation_term(b, t, e):
            return b.element_accumulation[t, e] if dynamic else 0

        # Mass transfer term
        def transfer_term(b, t, e):
            return (b.elemental_mass_transfer_term[t, e]
                    if has_mass_transfer else 0)

        # Custom term
        def user_term(t, e):
            if custom_elemental_term is not None:
                return custom_elemental_term(t, e)
            else:
                return 0

        # Element balances
        @self.Constraint(self.time_ref,
                         self.element_list_ref,
                         doc="Elemental material balances")
        def element_balances(b, t, e):
            return accumulation_term(b, t, e) == (
                        sum(b.elemental_flow_in[t, p, e]
                            for p in b.phase_list_ref) -
                        sum(b.elemental_flow_out[t, p, e]
                            for p in b.phase_list_ref) +
                        transfer_term(b, t, e) +
                        user_term(t, e))

        # Elemental Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time_ref,
                             self.element_list_ref,
                             doc="Elemental holdup calculation")
            def elemental_holdup_calculation(b, t, e):
                return b.element_holdup[t, e] == (
                    b.volume[t] *
                    sum(b.phase_fraction[t, p] *
                        b.properties_out[t].get_material_density_terms(p, j) *
                        b.properties_out[t]
                        .config.parameters.element_comp[j][e]
                        for p in b.phase_list_ref
                        for j in b.component_list_ref))

        return self.element_balances

    def add_total_material_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_total_material_balances (yet)."
                .format(self.name))

    def add_total_enthalpy_balances(self,
                                    has_heat_of_reaction=False,
                                    has_heat_transfer=False,
                                    has_work_transfer=False,
                                    custom_term=None):
        """
        This method constructs a set of 0D enthalpy balances indexed by time
        and phase.

        Args:
            has_heat_of_reaction - whether terms for heat of reaction should
                    be included in enthalpy balance
            has_heat_transfer - whether terms for heat transfer should be
                    included in enthalpy balances
            has_work_transfer - whether terms for work transfer should be
                    included in enthalpy balances
            custom_term - a Pyomo Expression representing custom terms to
                    be included in enthalpy balances.
                    Expression must be indexed by time and phase list

        Returns:
            Constraint object representing enthalpy balances
        """
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        # Test for components that must exist prior to calling this method
        if has_holdup:
            if not hasattr(self, "volume"):
                raise ConfigurationError(
                        "{} control volume must have volume defined to have "
                        "holdup terms. Please call the "
                        "add_geometry method before adding balance equations."
                        .format(self.name))
        if has_heat_of_reaction:
            if not (hasattr(self, "rate_reaction_extent") or
                    hasattr(self, "equilibrium_reaction_extent")):
                raise ConfigurationError(
                        "{} extent of reaction terms must exist in order to "
                        "have heat of reaction terms. Please ensure that "
                        "add_material_balance (or equivalent) is called before"
                        " adding energy balances.".format(self.name))

        # Get units from property package
        units = {}
        for u in ['energy', 'time']:
            try:
                units[u] = \
                   self.config.property_package.get_metadata().default_units[u]
            except KeyError:
                units[u] = '-'

        # Create variables
        if has_holdup:
            self.enthalpy_holdup = Var(
                        self.time_ref,
                        self.phase_list_ref,
                        domain=Reals,
                        doc="Enthalpy holdup in unit [{}]"
                        .format(units['energy']))

        if dynamic is True:
            self.enthalpy_accumulation = DerivativeVar(
                        self.enthalpy_holdup,
                        wrt=self.time_ref,
                        doc="Enthaly holdup in unit [{}/{}]"
                        .format(units['energy'], units['time']))

        # Create scaling factor
        self.scaling_factor_energy = Param(
                        default=1e-6,
                        mutable=True,
                        doc='Energy balance scaling parameter')

        # Create energy balance terms as needed
        # Heat transfer term
        if has_heat_transfer:
            self.heat = Var(self.time_ref,
                            domain=Reals,
                            initialize=0.0,
                            doc="Heat transfered in unit [{}/{}]"
                                .format(units['energy'], units['time']))

        # Work transfer
        if has_work_transfer:
            self.work = Var(self.time_ref,
                            domain=Reals,
                            initialize=0.0,
                            doc="Work transfered in unit [{}/{}]"
                                .format(units['energy'], units['time']))

        # Heat of Reaction
        if has_heat_of_reaction:
            @self.Expression(self.time_ref,
                             doc="Heat of reaction term [{}/{}]"
                                 .format(units['energy'], units['time']))
            def heat_of_reaction(b, t):
                if hasattr(self, "rate_reaction_extent"):
                    rate_heat = -sum(b.rate_reaction_extent[t, r] *
                                    b.reactions[t].dh_rxn[r]
                                    for r in self.rate_reaction_idx_ref)
                else:
                    rate_heat = 0

                if hasattr(self, "equilibrium_reaction_extent"):
                    equil_heat = -sum(
                            b.equilibrium_reaction_extent[t, e] *
                            b.reactions[t].dh_rxn[e]
                            for e in self.equilibrium_reaction_idx_ref)
                else:
                    equil_heat = 0

                return rate_heat + equil_heat

        # Create rules to substitute energy balance terms
        # Accumulation term
        def accumulation_term(b, t, p):
            return b.enthalpy_accumulation[t, p] if dynamic else 0

        def heat_term(b, t):
            return b.heat[t] if has_heat_transfer else 0

        def work_term(b, t):
            return b.work[t] if has_work_transfer else 0

        def rxn_heat_term(b, t):
            return b.heat_of_reaction[t] if has_heat_of_reaction else 0

        # Custom term
        def user_term(t):
            if custom_term is not None:
                return custom_term(t)
            else:
                return 0

        # Energy balance equation
        @self.Constraint(self.time_ref, doc="Energy balances")
        def enthalpy_balances(b, t):
            return (sum(accumulation_term(b, t, p) for p in b.phase_list_ref) *
                    b.scaling_factor_energy) == (
                        sum(b.properties_in[t].get_enthalpy_flow_terms(p)
                            for p in b.phase_list_ref) *
                        b.scaling_factor_energy -
                        sum(self.properties_out[t].get_enthalpy_flow_terms(p)
                            for p in b.phase_list_ref) *
                        b.scaling_factor_energy +
                        heat_term(b, t)*b.scaling_factor_energy +
                        work_term(b, t)*b.scaling_factor_energy +
                        rxn_heat_term(b, t)*b.scaling_factor_energy +
                        user_term(t)*b.scaling_factor_energy)

        # Energy Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time_ref,
                             self.phase_list_ref,
                             doc="Enthalpy holdup constraint")
            def enthalpy_holdup_calculation(b, t, p):
                return b.enthalpy_holdup[t, p] == (
                            b.volume[t]*self.phase_fraction[t, p] *
                            b.properties_out[t].get_enthalpy_density_terms(p))

        return self.enthalpy_balances

    def add_phase_enthalpy_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_phase_enthalpy_balances."
                .format(self.name))

    def add_phase_energy_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_phase_energy_balances."
                .format(self.name))

    def add_total_energy_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_total_energy_balances."
                .format(self.name))

    def add_total_pressure_balances(self,
                                    has_pressure_change=False,
                                    custom_term=None):
        """
        This method constructs a set of 0D pressure balances indexed by time.

        Args:
            has_pressure_change - whether terms for pressure change should be
                    included in enthalpy balances
            custom_term - a Pyomo Expression representing custom terms to
                    be included in pressure balances.
                    Expression must be indexed by time

        Returns:
            Constraint object representing pressure balances
        """
        # Get dynamic and holdup flags from config block
        dynamic = self.config.dynamic
        has_holdup = self.config.has_holdup

        if dynamic:
            _log.info("{} add_total_pressure_balances was provided with "
                      "argument dynamic = True. Total pressure balances do "
                      "not support dynamic terms (yet), and this argument "
                      "will be ignored.".format(self.name))

        if has_holdup:
            _log.info("{} add_total_pressure_balances was provided with "
                      "argument has_holdup = True. Total pressure balances do "
                      "not support holdup terms (yet), and this argument "
                      "will be ignored.".format(self.name))

        # Get units from property package
        try:
            p_units = (self.config.property_package.get_metadata().
                       default_units['pressure'])
        except KeyError:
            p_units = '-'

        # Add Momentum Balance Variables as necessary
        if has_pressure_change:
            self.deltaP = Var(self.time_ref,
                              domain=Reals,
                              doc="Pressure difference across unit [{}]"
                                  .format(p_units))

        # Create rules to substitute energy balance terms
        # Pressure change term
        def deltaP_term(b, t):
            return b.deltaP[t] if has_pressure_change else 0

        # Custom term
        def user_term(t):
            if custom_term is not None:
                return custom_term(t)
            else:
                return 0

        # Create scaling factor
        self.scaling_factor_pressure = Param(
                    default=1e-4,
                    mutable=True,
                    doc='Momentum balance scaling parameter')

        # Momentum balance equation
        @self.Constraint(self.time_ref, doc='Momentum balance')
        def pressure_balance(b, t):
            return 0 == (b.properties_in[t].pressure *
                         b.scaling_factor_pressure -
                         b.properties_out[t].pressure *
                         b.scaling_factor_pressure +
                         deltaP_term(b, t)*b.scaling_factor_pressure +
                         user_term(t)*b.scaling_factor_pressure)

        return self.pressure_balance

    def add_phase_pressure_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_phase_pressure_balances."
                .format(self.name))

    def add_phase_momentum_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_phase_momentum_balances."
                .format(self.name))

    def add_total_momentum_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_total_momentum_balances."
                .format(self.name))

    def model_check(blk):
        """
        This method executes the model_check methods on the associated state
        blocks (if they exist). This method is generally called by a unit model
        as part of the unit's model_check method.

        Args:
            None

        Returns:
            None
        """
        # Try property block model check
        for t in blk.time_ref:
            try:
                blk.properties_in[t].model_check()
            except AttributeError:
                _log.warning('{} ControlVolume inlet property block has no '
                             'model checks. To correct this, add a model_check'
                             ' method to the associated StateBlock class.'
                             .format(blk.name))
            try:
                blk.properties_out[t].model_check()
            except AttributeError:
                _log.warning('{} ControlVolume outlet property block has no '
                             'model checks. To correct this, add a '
                             'model_check method to the associated '
                             'StateBlock class.'.format(blk.name))

            try:
                blk.reactions[t].model_check()
            except AttributeError:
                _log.warning('{} ControlVolume outlet reaction block has no '
                             'model check. To correct this, add a '
                             'model_check method to the associated '
                             'ReactionBlock class.'.format(blk.name))

    def initialize(blk, state_args=None, outlvl=0, optarg=None,
                   solver='ipopt', hold_state=True):
        '''
        Initialisation routine for 0D control volume (default solver ipopt)

        Keyword Arguments:
            state_args : a dict of arguments to be passed to the property
                         package(s) to provide an initial state for
                         initialization (see documentation of the specific
                         property package) (default = {}).
            outlvl : sets output level of initialisation routine. **Valid
                     values:** **0** - no output (default), **1** - return
                     solver state for each step in routine, **2** - include
                     solver output infomation (tee=True)
            optarg : solver options dictionary object (default=None)
            solver : str indicating whcih solver to use during
                     initialization (default = 'ipopt')
            hold_state : flag indicating whether the initialization routine
                     should unfix any state variables fixed during
                     initialization, **default** - True. **Valid values:**
                     **True** - states variables are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the release_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        '''
        # Get inlet state if not provided
        if state_args is None:
            state_args = {}
            state_dict = \
                blk.properties_in[blk.time_ref.first()].define_port_members()

            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # Initialize state blocks
        flags = blk.properties_in.initialize(outlvl=outlvl-1,
                                             optarg=optarg,
                                             solver=solver,
                                             hold_state=hold_state,
                                             **state_args)

        blk.properties_out.initialize(outlvl=outlvl-1,
                                      optarg=optarg,
                                      solver=solver,
                                      hold_state=False,
                                      **state_args)

        try:
            blk.reactions.initialize(outlvl=outlvl-1,
                                     optarg=optarg,
                                     solver=solver)
        except AttributeError:
            pass

        if outlvl > 0:
            _log.info('{} Initialisation Complete'.format(blk.name))

        return flags

    def release_state(blk, flags, outlvl=0):
        '''
        Method to release state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of logging

        Returns:
            None
        '''
        blk.properties_in.release_state(flags, outlvl=outlvl-1)

    def _add_phase_fractions(self):
        """
        This method constructs the phase_fraction variables for the control
        volume, and the associated constraint on the sum of phase_fractions
        == 1. For systems with only one phase, phase_fraction is created as a
        Pyomo Expression with a value of 1.

        Args:
            None

        Returns:
            None
        """
        if len(self.phase_list_ref) > 1:
            self.phase_fraction = Var(
                            self.time_ref,
                            self.phase_list_ref,
                            initialize=1/len(self.phase_list_ref),
                            doc='Volume fraction of holdup by phase')

            @self.Constraint(self.time_ref,
                             doc='Sum of phase fractions == 1')
            def sum_of_phase_fractions(self, t):
                return 1 == sum(self.phase_fraction[t, p]
                                for p in self.phase_list_ref)
        else:
            @self.Expression(self.time_ref,
                             self.phase_list_ref,
                             doc='Volume fraction of holdup by phase')
            def phase_fraction(self, t, p):
                return 1
