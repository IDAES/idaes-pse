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
This module contains simple property packages for the purpose of testing.
"""

__author__ = "Robert Parker"


from pyomo.environ import (Set, SolverFactory, Var, Reals, Constraint, Param,
        exp)

from idaes.core import (declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlock,
                        StateBlockData,
                        ReactionParameterBlock,
                        ReactionBlockBase,
                        ReactionBlockDataBase,
                        MaterialFlowBasis,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MomentumBalanceType)


def get_default_solver():
    """
    Tries to set-up the default solver for testing, and returns None if not
    available
    """
    if SolverFactory('ipopt').available(exception_flag=False):
        solver = SolverFactory('ipopt')
        solver.options = {'tol': 1e-6,
                          'linear_solver': 'ma27'}
    else:
        solver = None

    return solver


@declare_process_block_class("AqueousEnzymeParameterBlock")
class ParameterData(PhysicalParameterBlock):
    """
    Parameter block for the aqueous enzyme reaction in the biochemical CSTR
    used by Christofides and Daoutidis, 1996, presented by Heineken et al, 1967
    """
    def build(self):
        super(ParameterData, self).build()

        # all components are in the aqueous phase
        self.phase_list = Set(initialize=['aq'])
        self.component_list = Set(initialize=['S', 'E', 'C', 'P'])

        self.state_block_class = AqueousEnzymeStateBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': 'min',
                               'length': 'm',
                               'amount': 'kmol',
                               'temperature': 'K',
                               'energy': 'kcal',
                               'holdup': 'kmol'})


class _AqueousEnzymeStateBlock(StateBlock):
    def initialize(blk):
        pass

@declare_process_block_class("AqueousEnzymeStateBlock",
                             block_class=_AqueousEnzymeStateBlock)
class AqueousEnzymeStateBlockData(StateBlockData):
    def build(self):
        super(AqueousEnzymeStateBlockData, self).build()

        self.conc_mol = Var(self._params.component_list,
                             domain=Reals,
                             doc='Component molar concentration [kmol/m^3]')
        
        self.flow_mol_comp = Var(self._params.component_list,
                                 domain=Reals, 
                                 doc='Molar component flow rate [kmol/min]')

        self.flow_rate = Var(domain=Reals,
                             doc='Volumetric flow rate out of reactor [m^3/min]')

        self.temperature = Var(initialize=303, domain=Reals,
                               doc='Temperature within reactor [K]')

        def flow_mol_comp_rule(b, j):
            return b.flow_mol_comp[j] == b.flow_rate*b.conc_mol[j]
        
        self.flow_mol_comp_eqn = Constraint(self._params.component_list,
                rule=flow_mol_comp_rule,
                doc='Outlet component molar flow rate equation')

    def get_material_density_terms(b, p, j):
        return b.conc_mol[j]

    def get_material_flow_terms(b, p, j):
        return b.flow_mol_comp[j]

    def get_material_flow_basis(b):
        return MaterialFlowBasis.molar

    def get_enthalpy_flow_terms(b, p):
        return b.flow_rate*b.temperature

    def get_energy_density_terms(b, p):
        return b.temperature

    def define_state_vars(b):
        return {'conc_mol': b.conc_mol,
                'flow_mol_comp': b.flow_mol_comp,
                'temperature': b.temperature,
                'flow_rate': b.flow_rate}

@declare_process_block_class('EnzymeReactionParameterBlock')
class EnzymeReactionParameterData(ReactionParameterBlock):
    '''
    Enzyme reaction:
    S + E <-> C -> P + E
    '''
    def build(self):
        super(EnzymeReactionParameterData, self).build()

        self.reaction_block_class = EnzymeReactionBlock

        self.rate_reaction_idx = Set(initialize=['R1', 'R2', 'R3'])
        self.rate_reaction_stoichiometry = {('R1', 'aq', 'S'): -1,
                                            ('R1', 'aq', 'E'): -1,
                                            ('R1', 'aq', 'C'): 1,
                                            ('R1', 'aq', 'P'): 0,
                                            ('R2', 'aq', 'S'): 1,
                                            ('R2', 'aq', 'E'): 1,
                                            ('R2', 'aq', 'C'): -1,
                                            ('R2', 'aq', 'P'): 0,
                                            ('R3', 'aq', 'S'): 0,
                                            ('R3', 'aq', 'E'): 1,
                                            ('R3', 'aq', 'C'): -1,
                                            ('R3', 'aq', 'P'): 1}

        self.act_energy = Param(self.rate_reaction_idx,
                initialize={'R1': 8.0e3,
                            'R2': 9.0e3,
                            'R3': 1.0e4},
                doc='Activation energy [kcal/kmol]')

        self.gas_const = Param(initialize=1.987, 
                doc='Gas constant R [kcal/kmol/K]')

        self.temperature_ref = Param(initialize=300.0, doc='Reference temperature')

        self.k_rxn = Param(self.rate_reaction_idx,
                initialize={'R1': 3.36e6,
                            'R2': 1.80e6,
                            'R3': 5.79e7},
                doc='Pre-exponential rate constant in Arrhenius expression')

    #    self.reaction_block_class = EnzymeReactionBlock

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units({'time': 'min',
                               'length': 'm',
                               'amount': 'kmol',
                               'energy': 'kcal'})

class _EnzymeReactionBlock(ReactionBlockBase):
    def initialize(blk):
        # initialize for reaction rates for each data object
        pass

@declare_process_block_class('EnzymeReactionBlock',
                             block_class=_EnzymeReactionBlock)
class EnzymeReactionBlockData(ReactionBlockDataBase):
    def build(self):
        super(EnzymeReactionBlockData, self).build()

        self.reaction_coef = Var(self._params.rate_reaction_idx,
                         domain=Reals, doc='Reaction rate coefficient')

        self.reaction_rate = Var(self._params.rate_reaction_idx,
                                 domain=Reals, 
                                 doc='Reaction rate [kmol/m^3/min]')

        self.dh_rxn = Param(self._params.rate_reaction_idx,
                domain=Reals, doc='Heat of reaction',
                initialize={'R1': 1e3/900/0.231,
                            'R2': 1e3/900/0.231,
                            'R3': 5e3/900/0.231})

        def reaction_rate_rule(b, r):
            if r == 'R1':
                return (b.reaction_rate[r] == 
                        b.reaction_coef[r]*
                        b.state_ref.conc_mol['S']*b.state_ref.conc_mol['E'])
            elif r == 'R2':
                return (b.reaction_rate[r] ==
                        b.reaction_coef[r]*
                        b.state_ref.conc_mol['C'])
            elif r == 'R3':
                return (b.reaction_rate[r] ==
                        b.reaction_coef[r]*
                        b.state_ref.conc_mol['C'])

        self.reaction_rate_eqn = Constraint(self._params.rate_reaction_idx,
                rule=reaction_rate_rule)

        def arrhenius_rule(b, r):
            return (b.reaction_coef[r] == b._params.k_rxn[r]*
                    exp(-b._params.act_energy[r]/b._params.gas_const/
                        b.state_ref.temperature))

        self.arrhenius_eqn = Constraint(self._params.rate_reaction_idx,
                rule=arrhenius_rule)
    
    def get_reaction_rate_basis(b):
        return MaterialFlowBasis.molar


