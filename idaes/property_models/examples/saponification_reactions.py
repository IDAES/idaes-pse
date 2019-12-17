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
Example property package for the saponification of Ethyl Acetate with NaOH
Assumes dilute solutions with properties of H2O.
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           exp,
                           Param,
                           PositiveReals,
                           Set,
                           value,
                           Var)

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        ReactionParameterBlock,
                        ReactionBlockDataBase,
                        ReactionBlockBase)
from idaes.core.util.misc import add_object_reference

# Some more inforation about this module
__author__ = "Andrew Lee"


# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("SaponificationReactionParameterBlock")
class ReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(ReactionParameterData, self).build()

        self.reaction_block_class = ReactionBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Liq'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['H2O', 'NaOH',
                                              'EthylAcetate',
                                              'SodiumAcetate',
                                              'Ethanol'])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = {("R1", "Liq", "NaOH"): -1,
                                            ("R1", "Liq", "EthylAcetate"): -1,
                                            ("R1", "Liq", "SodiumAcetate"): 1,
                                            ("R1", "Liq", "Ethanol"): 1,
                                            ("R1", "Liq", "H2O"): 0}

        # Arrhenius Constant
        self.arrhenius = Param(default=3.132e6,
                               doc="Arrhenius constant")

        # Activation Energy
        self.energy_activation = Param(default=43000,
                                       doc="Activation energy [J/mol]")

        # Heat of Reaction
        dh_rxn_dict = {"R1": -49000}
        self.dh_rxn = Param(self.rate_reaction_idx,
                            initialize=dh_rxn_dict,
                            doc="Heat of reaction [J/mol]")

        # Gas Constant
        self.gas_const = Param(within=PositiveReals,
                               mutable=False,
                               default=8.314,
                               doc='Gas Constant [J/mol.K]')

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({
                'k_rxn': {'method': '_rate_constant', 'units': 'm^3/mol.s'},
                'reaction_rate': {'method': "_rxn_rate", 'units': 'mol/m^3.s'}
                })
        obj.add_default_units({'time': 's',
                               'length': 'm',
                               'mass': 'g',
                               'amount': 'mol',
                               'temperature': 'K',
                               'energy': 'J',
                               'holdup': 'mol'})


class _ReactionBlock(ReactionBlockBase):
    """
    This Class contains methods which should be applied to Reaction Blocks as a
    whole, rather than individual elements of indexed Reaction Blocks.
    """
    def initialize(blk, outlvl=0, **kwargs):
        '''
        Initialization routine for reaction package.

        Keyword Arguments:
            outlvl : sets output level of initialization routine

                     * 0 = no output (default)
                     * 1 = report after each step

        Returns:
            None
        '''
        if outlvl > 0:
            _log.info('{} Initialization Complete.'.format(blk.name))


@declare_process_block_class("ReactionBlock",
                             block_class=_ReactionBlock)
class ReactionBlockData(ReactionBlockDataBase):
    """
    An example reaction package for saponification of ethyl acetate
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(ReactionBlockData, self).build()

        # Create references to state vars
        # Concentration
        add_object_reference(self,
                             "conc_mol_comp_ref",
                             self.state_ref.conc_mol_comp)

        # Temperature
        add_object_reference(self,
                             "temperature_ref",
                             self.state_ref.temperature)

        # Heat of reaction - no _ref as this is the actual property
        add_object_reference(
                self,
                "dh_rxn",
                self.config.parameters.dh_rxn)

    # Rate constant method
    def _rate_constant(self):
        self.k_rxn = Var(initialize=1,
                         doc="Rate constant [m^3/mol.s]")

        try:
            self.arrhenius_eqn = Constraint(
                    expr=self.k_rxn == self._params.arrhenius *
                    exp(-self._params.energy_activation /
                        (self._params.gas_const*self.temperature_ref)))
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.k_rxn)
            self.del_component(self.arrhenius_eqn)
            raise

    # Rate of reaction method
    def _rxn_rate(self):
        self.reaction_rate = Var(self._params.rate_reaction_idx,
                                 initialize=0,
                                 doc="Rate of reaction [mol/m^3.s]")

        try:
            def rate_rule(b, r):
                return b.reaction_rate[r] == (
                            b.k_rxn*b.conc_mol_comp_ref["EthylAcetate"] *
                            b.conc_mol_comp_ref["NaOH"])
            self.rate_expression = Constraint(self._params.rate_reaction_idx,
                                              rule=rate_rule)
        except AttributeError:
            # If constraint fails, clean up so that DAE can try again later
            self.del_component(self.reaction_rate)
            self.del_component(self.rate_expression)
            raise

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar

    def model_check(blk):
        pass
