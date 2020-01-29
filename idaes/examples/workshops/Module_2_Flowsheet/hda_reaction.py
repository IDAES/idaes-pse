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
Property package for the hydrodealkylation of toluene to form benzene
"""

# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           exp,
                           Expression,
                           Param,
                           PositiveReals,
                           Set,
                           Var)

# Import IDAES cores
from idaes.core import (declare_process_block_class,
                        MaterialFlowBasis,
                        ReactionParameterBlock,
                        ReactionBlockDataBase,
                        ReactionBlockBase)
from idaes.core.util.misc import add_object_reference

# Set up logger
_log = logging.getLogger(__name__)


@declare_process_block_class("HDAReactionParameterBlock")
class HDAReactionParameterData(ReactionParameterBlock):
    """
    Property Parameter Block Class

    Contains parameters and indexing sets associated with properties for
    superheated steam.

    """
    def build(self):
        '''
        Callable method for Block construction.
        '''
        super(HDAReactionParameterData, self).build()

        self.reaction_block_class = HDAReactionBlock

        # List of valid phases in property package
        self.phase_list = Set(initialize=['Vap'])

        # Component list - a list of component identifiers
        self.component_list = Set(initialize=['benzene',
                                              'toluene',
                                              'hydrogen',
                                              'methane'])

        # Reaction Index
        self.rate_reaction_idx = Set(initialize=["R1"])

        # Reaction Stoichiometry
        self.rate_reaction_stoichiometry = {("R1", "Vap", "benzene"): 1,
                                            ("R1", "Vap", "toluene"): -1,
                                            ("R1", "Vap", "hydrogen"): -1,
                                            ("R1", "Vap", "methane"): 1,
                                            ("R1", "Liq", "benzene"): 0,
                                            ("R1", "Liq", "toluene"): 0,
                                            ("R1", "Liq", "hydrogen"): 0,
                                            ("R1", "Liq", "methane"): 0}

        # Heat of Reaction
        dh_rxn_dict = {"R1": -1.08e5}
        self.dh_rxn = Param(self.rate_reaction_idx,
                            initialize=dh_rxn_dict,
                            doc="Heat of reaction [J/mol]")

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


class ReactionBlock(ReactionBlockBase):
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


@declare_process_block_class("HDAReactionBlock",
                             block_class=ReactionBlock)
class HDAReactionBlockData(ReactionBlockDataBase):
    """
    An example reaction package for saponification of ethyl acetate
    """

    def build(self):
        """
        Callable method for Block construction
        """
        super(HDAReactionBlockData, self).build()

        # Heat of reaction - no _ref as this is the actual property
        add_object_reference(
                self,
                "dh_rxn",
                self.config.parameters.dh_rxn)

    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar
