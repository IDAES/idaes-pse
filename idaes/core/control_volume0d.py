##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
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
                        ControlVolumeBase,
                        FlowDirection,
                        useDefault)
from idaes.core.util.exceptions import (BalanceTypeNotSupportedError,
                                        ConfigurationError,
                                        PropertyNotSupportedError)

__author__ = "Andrew Lee"


_log = logging.getLogger(__name__)

# TODO : Custom terms in material balances, other types of material balances

@declare_process_block_class("ControlVolume0D", doc="""
    ControlVolume0D is a specialized Pyomo block for IDAES non-discretized
    control volume blocks, and contains instances of ControlVolume0dData.

    ControlVolume0D should be used for any control volume with a defined volume
    and distinct inlets and outlets which does not require spatial
    discretization. This encompases most basic unit models used in process
    modeling.""")
class ControlVolume0dData(ControlVolumeBase):
    """
    0-Dimensional (Non-Discretised) ControlVolume Class

    This class forms the core of all non-discretized IDAES models. It provides
    methods to build property and reaction blocks, and add mass, energy and
    momentum balances. The form of the terms used in these constraints is
    specified in the chosen property package.
    """

    def build(self):
        """
        Build method for ControlVolume0D blocks.

        Args:
            None

        Returns:
            None
        """
        # Call build method from base class
        super(ControlVolume0dData, self).build()

    # TODO : add autobuild method

    def add_geometry(self):
        """
        Method to create volume Var in ControlVolume.

        Args:
            None

        Returns:
            None
        """
        l_units = \
            self.config.property_package.get_metadata().default_units["length"]
        self.volume = Var(self.time, initialize=1.0,
                          doc='Holdup Volume [{}^3]'.format(l_units))

    def add_state_blocks(self,
                         information_flow=FlowDirection.forward,
                         has_phase_equilibrium=False,
                         package_arguments={}):
        """
        This method constructs the inlet and outlet state blocks for the
        control volume.

        Args:
            information_flow - a FlowDirection Enum indicating whether
                               information flows from inlet-to-outlet or
                               outlet-to-inlet
            has_phase_equilibrium - indicates whether equilibrium calculations
                                    will be required in state blocks
            package_arguments - dict-like object of arguments to be passed to
                                state blocks as construction arguments

        Returns:
            None
        """
        tmp_dict = package_arguments
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

        self.properties_in = self._property_module.StateBlock(
                self.time,
                doc="Material properties at inlet",
                default=tmp_dict)

        # Reverse defined_state
        tmp_dict["defined_state"] = not tmp_dict["defined_state"]

        self.properties_out = self._property_module.StateBlock(
                self.time,
                doc="Material properties at outlet",
                default=tmp_dict)

    def add_reaction_blocks(self,
                            has_equilibrium=False,
                            package_arguments={}):
        """
        This method constructs the reaction block for the control volume.

        Args:
            has_equilibrium - indicates whether equilibrium calculations
                              will be required in reaction block
            package_arguments - dict-like object of arguments to be passed to
                                reaction block as construction arguments

        Returns:
            None
        """
        tmp_dict = package_arguments
        tmp_dict["state_block"] = self.properties_out
        tmp_dict["has_equilibrium"] = has_equilibrium
        tmp_dict["parameters"] = self.config.reaction_package

        self.reactions = self._reaction_module.ReactionBlock(
                self.time,
                doc="Reaction properties in control volume",
                default=tmp_dict)

    def add_phase_component_balances(self,
                                     dynamic=useDefault,
                                     has_holdup=False,
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
            dynamic - argument indicating whether material balances should
                    include temporal derivative terms. If not provided,
                    will use the dynamic flag of the control volume block
            has_holdup - whether material holdup terms should be included in
                    material balances. Must be True if dynamic = True
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
            custom_molar_term - a Pyomo Expression reresenting custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, phase list and
                    component list
            custom_mass_term - a Pyomo Expression reresenting custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, phase list and
                    component list

        Returns:
            Constraint object representing material balances
        """
        # Validate arguments
        dynamic, has_holdup = self._validate_add_balance_arguments(
                                            dynamic=dynamic,
                                            has_holdup=has_holdup)

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
            for t in self.time:
                if self.reactions[t].config.has_equilibrium is False:
                    raise ConfigurationError(
                            "{} material balance was set to include "
                            "equilibrium reactions, however the assoicated "
                            "ReactionBlock was not set to include equilibrium "
                            "constraints (has_equilibrium_reactions=False). "
                            "Please correct your configuration arguments."
                            .format(self.name))

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            for t in self.time:
                if not self.properties_out[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the assoicated outlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                if not self.properties_in[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the assoicated inlet "
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
        if has_holdup or has_rate_reactions:
            if not hasattr(self, "volume"):
                raise ConfigurationError(
                        "{} control volume must have volume defined to have "
                        "holdup and/or rate reaction terms. Please call the "
                        "add_geometry method before adding balance equations."
                        .format(self.name))

        # Material holdup and accumulation
        if has_holdup:
            self.material_holdup = Var(self.time,
                                       self.phase_list,
                                       self.component_list,
                                       domain=Reals,
                                       doc="Material holdup in unit [{}]"
                                           .format(units['holdup']))
        if dynamic:
            self.material_accumulation = DerivativeVar(
                    self.material_holdup,
                    wrt=self.time,
                    doc="Material accumulation in unit [{}/{}]"
                        .format(units['holdup'], units['time']))

        # Get phase component list(s)
        phase_component_list = self._get_phase_comp_list()

        # Create material balance terms as required
        # Kinetic reaction generation
        if has_rate_reactions:
            try:
                # TODO : replace with Reference
                object.__setattr__(
                        self,
                        "rate_reaction_idx",
                        self.config.reaction_package.rate_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of rate "
                    "reactions (rate_reaction_idx), thus does not support "
                    "rate-based reactions.".format(self.name))
            self.rate_reaction_generation = Var(
                        self.time,
                        self.phase_list,
                        self.component_list,
                        domain=Reals,
                        doc="Amount of component generated in "
                            "unit by kinetic reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Equilibrium reaction generation
        if has_equilibrium_reactions:
            try:
                # TODO : replace with Reference
                object.__setattr__(
                    self,
                    "equilibrium_reaction_idx",
                    self.config.reaction_package.equilibrium_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of "
                    "equilibrium reactions (equilibrium_reaction_idx), thus "
                    "does not support equilibrium-based reactions."
                    .format(self.name))
            self.equilibrium_reaction_generation = Var(
                        self.time,
                        self.phase_list,
                        self.component_list,
                        domain=Reals,
                        doc="Amount of component generated in unit "
                            "by equilibrium reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Phase equilibrium generation
        if has_phase_equilibrium:
            try:
                # TODO : replace with Reference
                object.__setattr__(
                    self,
                    "phase_equilibrium_idx",
                    self.config.property_package.phase_equilibrium_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Property package does not contain a list of phase "
                    "equilibrium reactions (phase_equilibrium_idx), thus does "
                    "not support phase equilibrium.".format(self.name))
            self.phase_equilibrium_generation = Var(
                        self.time,
                        self.phase_equilibrium_idx,
                        domain=Reals,
                        doc="Amount of generation in unit by phase "
                            "equilibria [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Material transfer term
        if has_mass_transfer:
            self.mass_transfer_term = Var(
                        self.time,
                        self.phase_list,
                        self.component_list,
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
                for r in b.phase_equilibrium_idx:
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
                           for r in b.phase_equilibrium_idx)

        def transfer_term(b, t, p, j):
            return (b.mass_transfer_term[t, p, j] if has_mass_transfer else 0)

        # TODO : Add custom terms

        # Add component balances
        @self.Constraint(self.time,
                         self.phase_list,
                         self.component_list,
                         doc="Material balances")
        def material_balances(b, t, p, j):
            if j in phase_component_list[p]:
                return accumulation_term(b, t, p, j) == (
                        b.properties_in[t].get_material_flow_terms(p, j) -
                        b.properties_out[t].get_material_flow_terms(p, j) +
                        kinetic_term(b, t, p, j) +
                        equilibrium_term(b, t, p, j) +
                        phase_equilibrium_term(b, t, p, j) +
                        transfer_term(b, t, p, j))
            else:
                return Constraint.Skip

        # TODO: Need to set material_holdup = 0 for non-present component-phase
        # pairs. Not ideal, but needed to close DoF. Is there a better way?

        # Material Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time,
                             self.phase_list,
                             self.component_list,
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
                    self.time,
                    self.rate_reaction_idx,
                    domain=Reals,
                    doc="Extent of kinetic reactions[{}/{}]"
                        .format(units['holdup'], units['time']))

            @self.Constraint(self.time,
                             self.phase_list,
                             self.component_list,
                             doc="Kinetic reaction stoichiometry constraint")
            def rate_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.rate_reaction_generation[t, p, j] == (
                        sum(rblock[t].rate_reaction_stoichiometry[r, p, j] *
                            b.rate_reaction_extent[t, r]
                            for r in b.rate_reaction_idx))
                else:
                    return Constraint.Skip

            try:
                @self.Constraint(self.time,
                                 self.rate_reaction_idx,
                                 doc="Kinetic reaction extents constraint")
                def rate_reaction_extents_constraint(b, t, r):
                    return b.rate_reaction_extent[t, r] == (
                            rblock[t].reaction_rate[r]*b.volume[t])
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a reaction_rate "
                    "variable, thus does not support kinetic equilibrium."
                    .format(self.name))

        if has_equilibrium_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.equilibrium_reaction_extent = Var(
                            self.time,
                            self.equilibrium_reaction_idx,
                            domain=Reals,
                            doc="Extent of equilibrium reactions[{}/{}]"
                                .format(units['holdup'], units['time']))

            @self.Constraint(self.time,
                             self.phase_list,
                             self.component_list,
                             doc="Equilibrium reaction stoichiometry")
            def equilibrium_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.equilibrium_reaction_generation[t, p, j] == (
                            sum(rblock[t].
                                equilibrium_reaction_stoichiometry[r, p, j] *
                                b.equilibrium_reaction_extent[t, r]
                                for r in b.equilibrium_reaction_idx))
                else:
                    return Constraint.Skip

        return self.material_balances

    def add_total_component_balances(self,
                                     dynamic=useDefault,
                                     has_holdup=False,
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
            dynamic - argument indicating whether material balances should
                    include temporal derivative terms. If not provided,
                    will use the dynamic flag of the control volume block
            has_holdup - whether material holdup terms should be included in
                    material balances. Must be True if dynamic = True
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
            custom_molar_term - a Pyomo Expression reresenting custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, phase list and
                    component list
            custom_mass_term - a Pyomo Expression reresenting custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, phase list and
                    component list

        Returns:
            Constraint object representing material balances
        """
        # Validate arguments
        dynamic, has_holdup = self._validate_add_balance_arguments(
                                            dynamic=dynamic,
                                            has_holdup=has_holdup)

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
            for t in self.time:
                if self.reactions[t].config.has_equilibrium is False:
                    raise ConfigurationError(
                            "{} material balance was set to include "
                            "equilibrium reactions, however the assoicated "
                            "ReactionBlock was not set to include equilibrium "
                            "constraints (has_equilibrium_reactions=False). "
                            "Please correct your configuration arguments."
                            .format(self.name))

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            for t in self.time:
                if not self.properties_out[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the assoicated outlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                if not self.properties_in[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the assoicated inlet "
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
        if has_holdup or has_rate_reactions:
            if not hasattr(self, "volume"):
                raise ConfigurationError(
                        "{} control volume must have volume defined to have "
                        "holdup and/or rate reaction terms. Please call the "
                        "add_geometry method before adding balance equations."
                        .format(self.name))

        # Material holdup and accumulation
        if has_holdup:
            self.material_holdup = Var(self.time,
                                       self.phase_list,
                                       self.component_list,
                                       domain=Reals,
                                       doc="Material holdup in unit [{}]"
                                           .format(units['holdup']))
        if dynamic:
            self.material_accumulation = DerivativeVar(
                    self.material_holdup,
                    wrt=self.time,
                    doc="Material accumulation in unit [{}/{}]"
                        .format(units['holdup'], units['time']))

        # Get phase component list(s)
        phase_component_list = self._get_phase_comp_list()

        # Create material balance terms as required
        # Kinetic reaction generation
        if has_rate_reactions:
            try:
                # TODO : replace with Reference
                object.__setattr__(
                        self,
                        "rate_reaction_idx",
                        self.config.reaction_package.rate_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of rate "
                    "reactions (rate_reaction_idx), thus does not support "
                    "rate-based reactions.".format(self.name))
            self.rate_reaction_generation = Var(
                        self.time,
                        self.phase_list,
                        self.component_list,
                        domain=Reals,
                        doc="Amount of component generated in "
                            "unit by kinetic reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Equilibrium reaction generation
        if has_equilibrium_reactions:
            try:
                # TODO : replace with Reference
                object.__setattr__(
                    self,
                    "equilibrium_reaction_idx",
                    self.config.reaction_package.equilibrium_reaction_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a list of "
                    "equilibrium reactions (equilibrium_reaction_idx), thus "
                    "does not support equilibrium-based reactions."
                    .format(self.name))
            self.equilibrium_reaction_generation = Var(
                        self.time,
                        self.phase_list,
                        self.component_list,
                        domain=Reals,
                        doc="Amount of component generated in unit "
                            "by equilibrium reactions [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Phase equilibrium generation
        if has_phase_equilibrium:
            try:
                # TODO : replace with Reference
                object.__setattr__(
                    self,
                    "phase_equilibrium_idx",
                    self.config.property_package.phase_equilibrium_idx)
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Property package does not contain a list of phase "
                    "equilibrium reactions (phase_equilibrium_idx), thus does "
                    "not support phase equilibrium.".format(self.name))
            self.phase_equilibrium_generation = Var(
                        self.time,
                        self.phase_equilibrium_idx,
                        domain=Reals,
                        doc="Amount of generation in unit by phase "
                            "equilibria [{}/{}]"
                            .format(units['holdup'], units['time']))

        # Material transfer term
        if has_mass_transfer:
            self.mass_transfer_term = Var(
                        self.time,
                        self.phase_list,
                        self.component_list,
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
                for r in b.phase_equilibrium_idx:
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
                           for r in b.phase_equilibrium_idx)

        def transfer_term(b, t, p, j):
            return (b.mass_transfer_term[t, p, j] if has_mass_transfer else 0)

        # TODO : Add custom terms

        # Add component balances
        @self.Constraint(self.time,
                         self.component_list,
                         doc="Material balances")
        def material_balances(b, t, j):
            cplist = []
            for p in self.phase_list:
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
                sum(transfer_term(b, t, p, j) for p in cplist))

        # TODO: Need to set material_holdup = 0 for non-present component-phase
        # pairs. Not ideal, but needed to close DoF. Is there a better way?

        # Material Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time,
                             self.phase_list,
                             self.component_list,
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
                    self.time,
                    self.rate_reaction_idx,
                    domain=Reals,
                    doc="Extent of kinetic reactions[{}/{}]"
                        .format(units['holdup'], units['time']))

            @self.Constraint(self.time,
                             self.phase_list,
                             self.component_list,
                             doc="Kinetic reaction stoichiometry constraint")
            def rate_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.rate_reaction_generation[t, p, j] == (
                        sum(rblock[t].rate_reaction_stoichiometry[r, p, j] *
                            b.rate_reaction_extent[t, r]
                            for r in b.rate_reaction_idx))
                else:
                    return Constraint.Skip

            try:
                @self.Constraint(self.time,
                                 self.rate_reaction_idx,
                                 doc="Kinetic reaction extents constraint")
                def rate_reaction_extents_constraint(b, t, r):
                    return b.rate_reaction_extent[t, r] == (
                            rblock[t].reaction_rate[r]*b.volume[t])
            except AttributeError:
                raise PropertyNotSupportedError(
                    "{} Reaction package does not contain a reaction_rate "
                    "variable, thus does not support kinetic equilibrium."
                    .format(self.name))

        if has_equilibrium_reactions:
            # Add extents of reaction and stoichiometric constraints
            self.equilibrium_reaction_extent = Var(
                            self.time,
                            self.equilibrium_reaction_idx,
                            domain=Reals,
                            doc="Extent of equilibrium reactions[{}/{}]"
                                .format(units['holdup'], units['time']))

            @self.Constraint(self.time,
                             self.phase_list,
                             self.component_list,
                             doc="Equilibrium reaction stoichiometry")
            def equilibrium_reaction_stoichiometry_constraint(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.equilibrium_reaction_generation[t, p, j] == (
                            sum(rblock[t].
                                equilibrium_reaction_stoichiometry[r, p, j] *
                                b.equilibrium_reaction_extent[t, r]
                                for r in b.equilibrium_reaction_idx))
                else:
                    return Constraint.Skip

        return self.material_balances

    def add_total_element_balances(self,
                                   dynamic=useDefault,
                                   has_holdup=False,
                                   has_rate_reactions=False,
                                   has_equilibrium_reactions=False,
                                   has_phase_equilibrium=False,
                                   has_mass_transfer=False,
                                   custom_molar_term=None,
                                   custom_mass_term=None):
        """
        This method constructs a set of 0D element balances indexed by time.

        Args:
            dynamic - argument indicating whether material balances should
                    include temporal derivative terms. If not provided,
                    will use the dynamic flag of the control volume block
            has_holdup - whether material holdup terms should be included in
                    material balances. Must be True if dynamic = True
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
            custom_molar_term - a Pyomo Expression reresenting custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, phase list and
                    component list
            custom_mass_term - a Pyomo Expression reresenting custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, phase list and
                    component list

        Returns:
            Constraint object representing material balances
        """

        # Validate arguments
        dynamic, has_holdup = self._validate_add_balance_arguments(
                                            dynamic=dynamic,
                                            has_holdup=has_holdup)

        if has_rate_reactions:
            raise ConfigurationError(
                    "{} add_total_element_balances method as provided with "
                    "argument has_rate_reactions = True. Total element "
                    "balances do not support rate based reactions, "
                    "please correct your configuration arguments"
                    .format(self.name))

        # Check that property package supports element balances
        try:
            object.__setattr__(self,
                               "element_list",
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
            for t in self.time:
                if self.reactions[t].config.has_equilibrium is False:
                    raise ConfigurationError(
                            "{} material balance was set to include "
                            "equilibrium reactions, however the assoicated "
                            "ReactionBlock was not set to include equilibrium "
                            "constraints (has_equilibrium_reactions=False). "
                            "Please correct your configuration arguments."
                            .format(self.name))
                try:
                # TODO : replace with Reference
                    object.__setattr__(
                        self,
                        "equilibrium_reaction_idx",
                        self.config.reaction_package.equilibrium_reaction_idx)
                except AttributeError:
                    raise PropertyNotSupportedError(
                        "{} Reaction package does not contain a list of "
                        "equilibrium reactions (equilibrium_reaction_idx), "
                        "thus does not support equilibrium-based reactions."
                        .format(self.name))

        if has_phase_equilibrium:
            # Check that state blocks are set to calculate equilibrium
            for t in self.time:
                if not self.properties_out[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the assoicated outlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                if not self.properties_in[t].config.has_phase_equilibrium:
                    raise ConfigurationError(
                            "{} material balance was set to include phase "
                            "equilibrium, however the assoicated inlet "
                            "StateBlock was not set to include equilibrium "
                            "constraints (has_phase_equilibrium=False). Please"
                            " correct your configuration arguments."
                            .format(self.name))
                try:
                    # TODO : replace with Reference
                    object.__setattr__(
                        self,
                        "phase_equilibrium_idx",
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
                    self.time,
                    self.element_list,
                    domain=Reals,
                    doc="Elemental holdup in unit [{}]"
                        .format(units['amount']))

        if dynamic:
            self.element_accumulation = DerivativeVar(
                    self.element_holdup,
                    wrt=self.time,
                    doc="Elemental accumulation in unit [{}/{}]"
                        .format(units['amount'], units['time']))

        @self.Expression(self.time,
                         self.phase_list,
                         self.element_list,
                         doc="Inlet elemental flow terms [{}/{}]"
                             .format(units['amount'], units['time']))
        def elemental_flow_in(b, t, p, e):
            return sum(b.properties_in[t].get_material_flow_terms(p, j) *
                       b.properties_out[t].element_comp[j][e]
                       for j in b.component_list)

        @self.Expression(self.time,
                         self.phase_list,
                         self.element_list,
                         doc="Outlet elemental flow terms [{}/{}]"
                             .format(units['amount'], units['time']))
        def elemental_flow_out(b, t, p, e):
            return sum(b.properties_out[t].get_material_flow_terms(p, j) *
                       b.properties_out[t].element_comp[j][e]
                       for j in b.component_list)

        # Create material balance terms as needed
        if has_mass_transfer:
            self.elemental_mass_transfer_term = Var(
                            self.time,
                            self.element_list,
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

        # Element balances
        @self.Constraint(self.time,
                         self.element_list,
                         doc="Elemental material balances")
        def element_balances(b, t, e):
            return accumulation_term(b, t, e) == (
                        sum(b.elemental_flow_in[t, p, e]
                            for p in b.phase_list) -
                        sum(b.elemental_flow_out[t, p, e]
                            for p in b.phase_list) +
                        transfer_term(b, t, e))

        # Elemental Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time,
                             self.element_list,
                             doc="Elemental holdup calcuation")
            def elemental_holdup_calculation(b, t, e):
                return b.element_holdup[t, e] == (
                    b.volume[t] *
                    sum(b.phase_fraction[t, p] *
                        b.properties_out[t].get_material_density_terms(p, j) *
                        b.properties_out[t].element_comp[j][e]
                        for p in b.phase_list
                        for j in b.component_list))

        return self.element_balances

    def add_total_material_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_total_material_balances (yet)."
                .format(self.name))

    def add_phase_enthalpy_balances(self,
                                    dynamic=useDefault,
                                    has_holdup=False,
                                    has_heat_transfer=False,
                                    has_work_transfer=False,
                                    custom_term=None):
        """
        This method constructs a set of 0D enthalpy balances indexed by time
        and phase.

        Args:
            dynamic - argument indicating whether enthalpy balances should
                    include temporal derivative terms. If not provided,
                    will use the dynamic flag of the control volume block
            has_holdup - whether enthalpy holdup terms should be included in
                    material balances. Must be True if dynamic = True
            has_heat_transfer - whether terms for heat transfer should be
                    included in enthalpy balances
            has_work_transfer - whether terms for work transfer should be
                    included in enthalpy balances
            custom_term - a Pyomo Expression reresenting custom terms to
                    be included in enthalpy balances.
                    Expression must be indexed by time and phase list

        Returns:
            Constraint object representing enthalpy balances
        """
        # Validate arguments
        dynamic, has_holdup = self._validate_add_balance_arguments(
                                            dynamic=dynamic,
                                            has_holdup=has_holdup)

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
        for u in ['energy', 'time']:
            try:
                units[u] = \
                   self.config.property_package.get_metadata().default_units[u]
            except KeyError:
                units[u] = '-'

        # Create scaling factor
        self.scaling_factor_energy = Param(
                        default=1e-6,
                        mutable=True,
                        doc='Energy balance scaling parameter')

        # Create energy balance terms as needed
        # Heat transfer term
        if has_heat_transfer:
            self.heat = Var(self.time,
                            domain=Reals,
                            initialize=0.0,
                            doc="Heat transfered in unit [{}/{}]"
                                .format(units['energy'], units['time']))

        # Work transfer
        if has_work_transfer:
            self.work = Var(self.time,
                            domain=Reals,
                            initialize=0.0,
                            doc="Work transfered in unit [{}/{}]"
                                .format(units['energy'], units['time']))

        # Create rules to substitute energy balance terms
        # Accumulation term
        def accumulation_term(b, t, p):
            return b.energy_accumulation[t, p] if b.config.dynamic else 0

        def heat_term(b, t):
            return b.heat[t] if b.config.has_heat_transfer else 0

        def work_term(b, t):
            return b.work[t] if b.config.has_work_transfer else 0

        # Energy balance equation
        @self.Constraint(self.time, doc="Energy balances")
        def enthalpy_balance(b, t):
            return (sum(accumulation_term(b, t, p) for p in b.phase_list) *
                    b.scaling_factor_energy) == (
                        sum(b.properties_in[t].get_enthalpy_flow_terms(p)
                            for p in b.phase_list) *
                        b.scaling_factor_energy -
                        sum(self.properties_out[t].get_enthalpy_flow_terms(p)
                            for p in b.phase_list) *
                        b.scaling_factor_energy +
                        heat_term(b, t)*b.scaling_factor_energy +
                        work_term(b, t)*b.scaling_factor_energy)

        # Energy Holdup
        if has_holdup:
            if not hasattr(self, "phase_fraction"):
                self._add_phase_fractions()

            @self.Constraint(self.time,
                             self.phase_list,
                             doc="Enthalpy holdup constraint")
            def enthalpy_holdup_calculation(b, t, p):
                return b.enthalpy_holdup[t, p] == (
                            b.volume[t]*self.phase_fraction[t, p] *
                            b.properties_out[t].get_energy_density_terms(p))

        return self.enthalpy_balances

    def add_total_enthalpy_balances(self, *args, **kwargs):
        raise BalanceTypeNotSupportedError(
                "{} OD control volumes do not support "
                "add_total_enthalpy_balances."
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
                                    dynamic=useDefault,
                                    has_holdup=False,
                                    has_pressure_change=False,
                                    custom_term=None):
        """
        This method constructs a set of 0D pressure balances indexed by time.

        Args:
            dynamic - argument indicating whether enthalpy balances should
                    include temporal derivative terms. If not provided,
                    will use the dynamic flag of the control volume block
            has_holdup - whether enthalpy holdup terms should be included in
                    material balances. Must be True if dynamic = True
            has_pressure_change - whether terms for pressure change should be
                    included in enthalpy balances
            custom_term - a Pyomo Expression reresenting custom terms to
                    be included in pressure balances.
                    Expression must be indexed by time

        Returns:
            Constraint object representing pressure balances
        """
        # Validate arguments
        dynamic, has_holdup = self._validate_add_balance_arguments(
                                            dynamic=dynamic,
                                            has_holdup=has_holdup)

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
            self.deltaP = Var(self.time,
                              domain=Reals,
                              doc="Pressure difference across unit [{}]"
                                  .format(p_units))

        # Create rules to substitute energy balance terms
        # Pressure change term
        def deltaP_term(b, t):
            return b.deltaP[t] if b.config.has_pressure_change else 0

        # Create scaling factor
        self.scaling_factor_pressure = Param(
                    default=1e-4,
                    mutable=True,
                    doc='Momentum balance scaling parameter')

        # Momentum balance equation
        @self.Constraint(self.time, doc='Momentum balance')
        def pressure_balance(b, t):
            return 0 == (b.properties_in[t].pressure *
                         b.scaling_factor_pressure -
                         b.properties_out[t].pressure *
                         b.scaling_factor_pressure +
                         deltaP_term(b, t)*b.scaling_factor_pressure)

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
        This method exectues the model_check methods on the associated property
        blocks (if they exist). This method is generally called by a unit model
        as part of the unit's model_check method.

        Args:
            None

        Returns:
            None
        """
        # Try property block model check
        for t in blk.time:
            try:
                blk.properties_in[t].model_check()
            except AttributeError:
                _log.warning('{} Holdup inlet property block has no model '
                             'check. To correct this, add a model_check '
                             'method to the associated PropertyBlock class.'
                             .format(blk.name))
            try:
                blk.properties_out[t].model_check()
            except AttributeError:
                _log.warning('{} Holdup outlet property block has no '
                             'model check. To correct this, add a '
                             'model_check method to the associated '
                             'PropertyBlock class.'.format(blk.name))

    def initialize(blk, state_args=None, outlvl=0, optarg=None,
                   solver='ipopt', hold_state=True):
        '''
        Initialisation routine for holdup (default solver ipopt)

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
                     **True** - states varaibles are not unfixed, and a dict of
                     returned containing flags for which states were fixed
                     during initialization, **False** - state variables are
                     unfixed after initialization by calling the relase_state
                     method.

        Returns:
            If hold_states is True, returns a dict containing flags for which
            states were fixed during initialization.
        '''
        # Get inlet state if not provided
        if state_args is None:
            state_args = {}
            state_dict = blk.properties_in[0].declare_port_members()
            for k in state_dict.keys():
                if state_dict[k].is_indexed():
                    state_args[k] = {}
                    for m in state_dict[k].keys():
                        state_args[k][m] = state_dict[k][m].value
                else:
                    state_args[k] = state_dict[k].value

        # Initialize property blocks
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

        if outlvl > 0:
            _log.info('{} Initialisation Complete'.format(blk.name))

        return flags

    def release_state(blk, flags, outlvl=0):
        '''
        Method to relase state variables fixed during initialisation.

        Keyword Arguments:
            flags : dict containing information of which state variables
                    were fixed during initialization, and should now be
                    unfixed. This dict is returned by initialize if
                    hold_state = True.
            outlvl : sets output level of of logging

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
        if len(self.phase_list) > 1:
            self.phase_fraction = Var(
                            self.time,
                            self.phase_list,
                            initialize=1/len(self.phase_list),
                            doc='Volume fraction of holdup by phase')

            @self.Constraint(self.time,
                             doc='Sum of phase fractions == 1')
            def sum_of_phase_fractions(self, t):
                return 1 == sum(self.phase_fraction[t, p]
                                for p in self.phase_list)
        else:
            @self.Expression(self.time,
                             self.phase_list,
                             doc='Volume fraction of holdup by phase')
            def phase_fraction(self, t, p):
                return 1
