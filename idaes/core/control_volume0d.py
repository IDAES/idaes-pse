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

# Import Pyomo libraries
from pyomo.environ import Constraint, Reals, Var
from pyomo.dae import DerivativeVar

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        ControlVolumeBase,
                        FlowDirection,
                        useDefault)
from idaes.core.util.exceptions import (ConfigurationError,
                                        PropertyNotSupportedError)

__author__ = "Andrew Lee"


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
        l_units = self.config.property_package.get_metadata().default_units[
                                                                      "length"]
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
            has_equilibrium: indicates whether equilibrium calculations
                              will be required in reaction block
            package_arguments: dict-like object of arguments to be passed to
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
            dynamic: argument indicating whether material balances should
                    include temporal derivative terms. If not provided,
                    will use the dynamic flag of the control volume block
            has_holdup: whether material holdup terms should be included in
                    material balances. Must be True if dynamic = True
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
            custom_molar_term: a Pyomo Expression reresenting custom terms to
                    be included in material balances on a molar basis.
                    Expression must be indexed by time, phase list and
                    component list
            custom_mass_term: a Pyomo Expression reresenting custom terms to
                    be included in material balances on a mass basis.
                    Expression must be indexed by time, phase list and
                    component list

        Returns:
            Constraint object representing material balances
        """
        # Validate arguments
        self._validate_add_balance_arguments(dynamic=dynamic,
                                             has_holdup=has_holdup)

        # Get units from property package
        units = {}
        for u in ['length', 'holdup', 'amount', 'time']:
            try:
                units[u] = self.config.property_package.get_package_units()[u]
            except KeyError:
                units[u] = '-'

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

        # Create material balance terms as required
        # Kinetic reaction generation
        if has_rate_reactions:
            try:
                # TODO : replace with Reference
                object.__setattr__(
                        self,
                        "rate_reaction_idx",
                        self.config.reaction_parameters.rate_reaction_idx)
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
                    self.config.reaction_parameters.equilibrium_reaction_idx)
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
                    self.config.property_parameters.phase_equilibrium_idx)
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
                for r in b.phase_equilibrium_idx:
                    if self.phase_equilibrium_list[r][0] == j:
                        if self.phase_equilibrium_list[r][1][0] == p:
                            sd[r] = 1
                        elif self.phase_equilibrium_list[r][1][1] == p:
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
        def material_balance(b, t, p, j):
            if j in phase_component_list[p]:
                return accumulation_term(b, t, p, j) == (
                        b.properties_in[t].get_material_balance_term[p, j] -
                        b.properties_out[t].get_material_balance_term[p, j] +
                        kinetic_term(b, t, p, j) +
                        equilibrium_term(b, t, p, j) +
                        phase_equilibrium_term(b, t, p, j) +
                        transfer_term(b, t, p, j))
            else:
                return Constraint.Skip

        # TODO: Need to set material_holdup = 0 for non-present component-phase
        # pairs. Not ideal, but needed to close DoF. Is there a better way?

        # Material Holdup
        if self.config.include_holdup:
            @self.Constraint(self.time,
                             self.phase_list,
                             self.component_list,
                             doc="Material holdup calculations")
            def material_holdup_calculation(b, t, p, j):
                if j in phase_component_list[p]:
                    return b.material_holdup[t, p, j] == (
                           b.volume[t]*self.phase_fraction[t, p] *
                           b.properties_out[t].material_density_term[p, j])
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

            @self.Constraint(self.time,
                             self.reaction_idx,
                             doc="Kinetic reaction extents constraint")
            def rate_reaction_extents_constraint(b, t, r):
                return b.rate_reaction_extent[t, r] == (
                        rblock[t].reaction_rate[r]*b.volume[t])

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

#    def _make_phase_frac(self):
#        """
#        This method constructs the phase_fraction variables for the control
#        volume, and the associated constraint on the sum of phase_fractions
#        == 1. For systems with only one phase, phase_fraction is created as a
#        Pyomo Expression with a value of 1.
#
#        Args:
#            None
#
#        Returns:
#            None
#        """
#        if len(self.phase_list) > 1:
#            self.phase_fraction = Var(
#                            self.time,
#                            self.phase_list,
#                            initialize=1/len(self.phase_list),
#                            doc='Volume fraction of holdup by phase')
#
#            @self.Constraint(self.time,
#                             doc='Sum of phase fractions == 1')
#            def sum_of_phase_fractions(self, t):
#                return 1 == sum(self.phase_fraction[t, p]
#                                for p in self.phase_list)
#        else:
#            @self.Expression(self.time,
#                             self.phase_list,
#                             doc='Volume fraction of holdup by phase')
#            def phase_fraction(self, t, p):
#                return 1
