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
from idaes.core.util.exceptions import ConfigurationError

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

        Args:
            None

        Returns:
            None
        """
        # Call build method from base class
        super(ControlVolume0dData, self).build()

    # TODO : add autobuild method

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
        if information_flow == FlowDirection.forward:
            inlet_defined = True
        elif information_flow == FlowDirection.backward:
            inlet_defined = False
        else:
            raise ConfigurationError(
                    '{} invalid value for information_flow argument. '
                    'Valid values are FlowDirection.forward and '
                    'FlowDirection.backward'.format(self.name))

        self.properties_in = self.property_module.StateBlock(
                self.time,
                doc="Material properties at inlet",
                defined_state=inlet_defined,
                has_phase_equilibrium=has_phase_equilibrium,
                parameters=self.config.property_package,
                **package_arguments)

        self.properties_out = self.property_module.PropertyBlock(
                self.time,
                doc="Material properties at outlet",
                defined_state=not inlet_defined,
                has_phase_equilibrium=has_phase_equilibrium,
                parameters=self.config.property_package,
                **package_arguments)

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
        self.reactions = self.reaction_module.ReactionBlock(
                self.time,
                doc="Reaction properties in control volume",
                state_block=self.properties_out,
                has_equilibrium=has_equilibrium,
                parameters=self.config.reaction_package,
                **package_arguments)

    def _add_volume(self, length_units):
        """
        Method to create volume Var in ControlVolume.

        Args:
            length_units - string to use for units for length

        Returns:
            None
        """
        self.volume = Var(self.time, initialize=1.0,
                          doc='Holdup Volume [{}^3]'.format(length_units))

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
        self._validate_add_balance_arguments(dynamic=dynamic,
                                             has_holdup=has_holdup)

        # Get units from property package
        units = {}
        for u in ['length', 'holdup', 'amount', 'time']:
            try:
                units[u] = self.config.property_package.get_package_units()[u]
            except KeyError:
                units[u] = '-'

        if self.config.has_holdup:
            if not hasattr(self, "volume"):
                self._add_volume(lenght_units=units["length"])

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

#        mbal_basis = self.properties_out[0].get_material_balance_term
#        def custom_molar_term(b, t, p, j):

        # Add component balances
        @self.Constraint(self.time,
                         self.phase_list,
                         self.component_list,
                         doc="Material balances")
        def material_balance(b, t, p, j):
            if j in phase_component_list[p]:
                return accumulation_term(b, t, p, j) == (
                        b.properties_in[t].material_balance_term[p, j] -
                        b.properties_out[t].material_balance_term[p, j] +
                        kinetic_term(b, t, p, j) +
                        equilibrium_term(b, t, p, j) +
                        phase_equilibrium_term(b, t, p, j) +
                        transfer_term(b, t, p, j))
            else:
                return Constraint.Skip

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
