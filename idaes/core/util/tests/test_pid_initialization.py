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
Test for element-by-element initialization on PID controller. This is important
as PID controllers have additional time-linking variables beyond derivative
and differential variables.
"""

import pytest
from pyomo.environ import (Block, ConcreteModel, Constraint, Expression,
                           Set, SolverFactory, Var, value, Param, Reals,
                           TransformationFactory, TerminationCondition,
                           exp)
from pyomo.network import Arc, Port
from pyomo.dae import DerivativeVar

from idaes.core import (FlowsheetBlock, 
                        MaterialBalanceType, 
                        EnergyBalanceType,
                        MomentumBalanceType, 
                        declare_process_block_class,
                        PhysicalParameterBlock,
                        StateBlock,
                        StateBlockData,
                        ReactionParameterBlock,
                        ReactionBlockBase,
                        ReactionBlockDataBase,
                        MaterialFlowBasis)
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.generic_models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.generic_models.control import PIDBlock, PIDForm
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.initialization import (fix_state_vars,
                                            revert_state_vars,
                                            propagate_state,
                                            solve_indexed_blocks,
                                            initialize_by_time_element)
import idaes.logger as idaeslog

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


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


def make_model(horizon=6, ntfe=60, ntcp=2, inlet_E=11.91, inlet_S=12.92):

    time_set = [0, horizon]

    m = ConcreteModel(name='CSTR with level control')
    m.fs = FlowsheetBlock(default={'dynamic': True,
                                   'time_set': time_set})

    m.fs.properties = AqueousEnzymeParameterBlock()
    m.fs.reactions = EnzymeReactionParameterBlock(
            default={'property_package': m.fs.properties})
    m.fs.cstr = CSTR(default={'has_holdup': True,
                              'property_package': m.fs.properties,
                              'reaction_package': m.fs.reactions,
                              'material_balance_type': MaterialBalanceType.componentTotal,
                              'energy_balance_type': EnergyBalanceType.enthalpyTotal,
                              'momentum_balance_type': MomentumBalanceType.none,
                              'has_heat_of_reaction': True})

    m.fs.mixer = Mixer(default={
        'property_package': m.fs.properties,
        'material_balance_type': MaterialBalanceType.componentTotal,
        'momentum_mixing_type': MomentumMixingType.none,
        'num_inlets': 2,
        'inlet_list': ['S_inlet', 'E_inlet']})
    # Allegedly the proper energy balance is being used...

    # Add DerivativeVar for CSTR volume
    m.fs.cstr.control_volume.dVdt = DerivativeVar(
            m.fs.cstr.control_volume.volume,
            wrt=m.fs.time)

    # Time discretization
    disc = TransformationFactory('dae.collocation')
    disc.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU')

    m.fs.pid = PIDBlock(default={'pv': m.fs.cstr.volume,
                                 'output': m.fs.cstr.outlet.flow_rate,
                                 'upper': 5.0,
                                 'lower': 0.5,
                                 'calculate_initial_integral': True,
                                 # ^ Why would initial integral be calculated
                                 # to be nonzero?
                                 'pid_form': PIDForm.velocity})

    m.fs.pid.gain.fix(-1.0)
    m.fs.pid.time_i.fix(0.1)
    m.fs.pid.time_d.fix(0.0)
    m.fs.pid.setpoint.fix(1.0)

    # Fix initial condition for volume:
    m.fs.cstr.volume.unfix()
    m.fs.cstr.volume[m.fs.time.first()].fix(1.0)

    # Fix initial conditions for other variables:
    for p, j in m.fs.properties.phase_list*m.fs.properties.component_list:
        m.fs.cstr.control_volume.material_holdup[0, p, j].fix(0.001)
    # Note: Model does not solve when initial conditions are empty tank

    m.fs.mixer.E_inlet.conc_mol.fix(0)
    m.fs.mixer.S_inlet.conc_mol.fix(0)

    for t, j in m.fs.time*m.fs.properties.component_list:
        if j == 'E':
            m.fs.mixer.E_inlet.conc_mol[t, j].fix(inlet_E)
        elif j == 'S':
            m.fs.mixer.S_inlet.conc_mol[t, j].fix(inlet_S)

    m.fs.mixer.E_inlet.flow_rate.fix(0.1)
    m.fs.mixer.S_inlet.flow_rate.fix(2.1)

    # Specify a perturbation to substrate flow rate:
    for t in m.fs.time:
        if t < horizon/4:
            continue
        else:
            m.fs.mixer.S_inlet.flow_rate[t].fix(3.0)
    
    '''
    Not sure what the 'proper' way to enforce this balance is...
    Should have some sort of 'sum_flow_rate_eqn' constraint, but
    that doesn't really make sense for my aqueous property package
    with dilute components...
    '''
    @m.fs.mixer.Constraint(m.fs.time,
            doc='Solvent flow rate balance')
    def total_flow_balance(mx, t):
        return (mx.E_inlet.flow_rate[t] + mx.S_inlet.flow_rate[t]
                == mx.outlet.flow_rate[t])

    m.fs.mixer.E_inlet.temperature.fix(290)
    m.fs.mixer.S_inlet.temperature.fix(310)

    m.fs.inlet = Arc(source=m.fs.mixer.outlet, destination=m.fs.cstr.inlet)

    '''
    This constraint is in lieu of tracking the CSTR's level and allowing
    the outlet flow rate to be another degree of freedom.
    ^ Not sure best way to do this in IDAES.
    '''
#    @m.fs.cstr.Constraint(m.fs.time,
#        doc='Total flow rate balance')
#    def total_flow_balance(cstr, t):
#        return (cstr.inlet.flow_rate[t] == cstr.outlet.flow_rate[t])
    '''
    This constraint is omitted in the PID controlled case - outlet flow
    rate will be determined by controller
    '''
    @m.fs.cstr.Constraint(m.fs.time,
            doc='Total volume balance')
    def total_volume_balance(cstr, t):
        return (cstr.control_volume.dVdt[t] ==
                cstr.inlet.flow_rate[t] - cstr.outlet.flow_rate[t])

    # Fix "initial condition" for outlet flow rate, as here it cannot be
    # specified by the PID controller
    m.fs.cstr.outlet.flow_rate[m.fs.time.first()].fix(2.2)

    # Specify initial condition for energy
    m.fs.cstr.control_volume.energy_holdup[m.fs.time.first(), 'aq'].fix(300)

    # This generates constraints for my arc, but shouldn't I be able to
    # enforce that the variables in mixer.outlet are the same as cstr.inlet,
    # insteady of having two collections of 'the same' variables with an
    # equality constraint?
    TransformationFactory('network.expand_arcs').apply_to(m.fs)

    # Need to deactivate some arc equations because they over-specify.
    # Not sure how to avoid this...
    m.fs.inlet_expanded.flow_mol_comp_equality.deactivate()

    return m


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_initialize():
    '''Very rough test, just to make sure degrees of freedom are not violated.
    '''
    mod = make_model(horizon=2, ntfe=20, ntcp=1, inlet_E=11.91, inlet_S=12.92)
    assert degrees_of_freedom(mod) == 0
    initialize_by_time_element(mod.fs, mod.fs.time, solver=solver, 
            outlvl=idaeslog.DEBUG,
            fix_diff_only=False)
    assert degrees_of_freedom(mod) == 0


if __name__ == '__main__':
    test_initialize()
