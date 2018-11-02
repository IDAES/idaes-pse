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
from pyomo.environ import Var

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        ControlVolumeBase,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType,
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

#    def add_material_balances(balance_type=MaterialBalanceType.componentPhase,
#                              dynamic=useDefault,
#                              has_holdup=False,
#                              has_rate_reactions=False,
#                              has_equilibrium_reactions=False,
#                              has_phase_equilibrium=False,
#                              has_mass_transfer=False,
#                              custom_molar_term=None,
#                              custom_mass_term=None):
#        """
#        General method for adding material balances to a 0D control volume.
#        This method makes calls to specialised sub-methods for each type of
#        material balance.
#
#        Args:
#            balance_type - MaterialBalanceType Enum indicating which type of
#                    material balance should be constructed.
#            dynamic - argument indicating whether material balances should
#                    include temporal derivative terms. If not provided,
#                    will use the dynamic flag of the control volume block
#            has_holdup - whether material holdup terms should be included in
#                    material balances. Must be True if dynamic = True
#            has_rate_reactions - whether default generation terms for rate
#                    reactions should be included in material balances
#            has_equilibrium_reactions - whether generation terms should for
#                    chemical equilibrium reactions should be included in
#                    material balances
#            has_phase_equilibrium - whether generation terms should for phase
#                    equilibrium behaviour should be included in material
#                    balances
#            has_mass_transfer - whether generic mass transfer terms should be
#                    included in material balances
#            custom_molar_term - a Pyomo Expression reresenting custom terms to
#                    be included in material balances on a molar basis.
#                    Expression must be indexed by time, phase list and
#                    component list
#            custom_mass_term - a Pyomo Expression reresenting custom terms to
#                    be included in material balances on a mass basis.
#                    Expression must be indexed by time, phase list and
#                    component list
#
#        Returns:
#            None
#        """
#        pass

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
