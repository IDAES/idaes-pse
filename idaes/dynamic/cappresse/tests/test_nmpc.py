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
Test for Cappresse's module for NMPC.
"""

import pytest
from pyomo.environ import (Block, ConcreteModel,  Constraint, Expression,
                           Set, SolverFactory, Var, value, 
                           TransformationFactory, TerminationCondition)
from pyomo.network import Arc
from pyomo.kernel import ComponentSet

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.testing import (PhysicalParameterTestBlock,
        AqueousEnzymeParameterBlock, EnzymeReactionParameterBlock,
        EnzymeReactionBlock)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.generic_models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.dynamic.cappresse import nmpc
from idaes.dynamic.cappresse.nmpc import *
import idaes.logger as idaeslog
import pdb

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


def make_model(horizon=6, ntfe=60, ntcp=2, 
        inlet_E=11.91, inlet_S=12.92):
#    horizon = 6 # minutes
    time_set = [0, horizon]
#    ntfe = 60 # fe spacing: 6 seconds
#    ntcp = 2
#    sample_time = 0.5 # 5 fe per sample

    m = ConcreteModel(name='CSTR model for testing')
    m.fs = FlowsheetBlock(default={'dynamic': True,
                                   'time_set': time_set})

    m.fs.properties = AqueousEnzymeParameterBlock()
    m.fs.reactions = EnzymeReactionParameterBlock(
            default={'property_package': m.fs.properties})
    m.fs.cstr = CSTR(default={'property_package': m.fs.properties,
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

    # Time discretization
    disc = TransformationFactory('dae.collocation')
    disc.apply_to(m, wrt=m.fs.time, nfe=ntfe, ncp=ntcp, scheme='LAGRANGE-RADAU')

    # Fix geometry variables
    m.fs.cstr.volume.fix(1.0)

    # Fix initial conditions:
    for p, j in m.fs.properties.phase_list*m.fs.properties.component_list:
        m.fs.cstr.control_volume.material_holdup[0, p, j].fix(0)

    m.fs.mixer.E_inlet.conc_mol.fix(0)
    m.fs.mixer.S_inlet.conc_mol.fix(0)

    for t, j in m.fs.time*m.fs.properties.component_list:
        if j == 'E':
            m.fs.mixer.E_inlet.conc_mol[t, j].fix(inlet_E)
        elif j == 'S':
            m.fs.mixer.S_inlet.conc_mol[t, j].fix(inlet_S)

    m.fs.mixer.E_inlet.flow_rate.fix(0.1)
    m.fs.mixer.S_inlet.flow_rate.fix(2.1)
    
    # Not sure what the 'proper' way to enforce this balance is...
    # Should have some sort of 'sum_flow_rate_eqn' constraint, but
    # that doesn't really make sense for my aqueous property package
    # with dilute components...
    @m.fs.mixer.Constraint(m.fs.time,
            doc='Solvent flow rate balance')
    def total_flow_balance(mx, t):
        return (mx.E_inlet.flow_rate[t] + mx.S_inlet.flow_rate[t]
                == mx.outlet.flow_rate[t])

    m.fs.mixer.E_inlet.temperature.fix(290)
    m.fs.mixer.S_inlet.temperature.fix(310)

    m.fs.inlet = Arc(source=m.fs.mixer.outlet, destination=m.fs.cstr.inlet)

#    m.fs.cstr.outlet.flow_rate.fix(2.2)

    # This constraint is in lieu of tracking the CSTR's level and allowing
    # the outlet flow rate to be another degree of freedom.
    # ^ Not sure how to do this in IDAES.
    @m.fs.cstr.Constraint(m.fs.time,
        doc='Total flow rate balance')
    def total_flow_balance(cstr, t):
        return (cstr.inlet.flow_rate[t] == cstr.outlet.flow_rate[t])

    # This is an initial condition specification via an algebraic variable...
    m.fs.cstr.outlet.temperature[m.fs.time.first()].fix(300)

    # This generates constraints for my arc, but shouldn't I be able to
    # enforce that the variables in mixer.outlet are the same as cstr.inlet,
    # insteady of having two collections of 'the same' variables with an
    # equality constraint?
    TransformationFactory('network.expand_arcs').apply_to(m.fs)

    # Need to deactivate some arc equations because they over specify.
    # Not sure how to avoid this...
    m.fs.inlet_expanded.flow_mol_comp_equality.deactivate()

    return m


@pytest.fixture
def model():
    pass

# @ pytest something...?
def test_find_comp_in_block():
    m1 = ConcreteModel()

    @m1.Block([1,2,3])
    def b1(b):
        b.v = Var([1,2,3])

    m2 = ConcreteModel()

    @m2.Block([1,2,3])
    def b1(b):
        b.v = Var([1,2,3,4])

    @m2.Block([1,2,3])
    def b2(b):
        b.v = Var([1,2,3])

    v1 = m1.b1[1].v[1]

    assert find_comp_in_block(m2, m1, v1) is m2.b1[1].v[1]

    v2 = m2.b2[1].v[1]
    v3 = m2.b1[3].v[4]

    # These should result in Attribute/KeyErrors
    #find_comp_in_block(m1, m2, v2)
    #find_comp_in_block(m1, m2, v3)
    assert find_comp_in_block(m1, m2, v2, allow_miss=True) is None
    assert find_comp_in_block(m1, m2, v3, allow_miss=True) is None


def test_constructor():
    m_plant = make_model(horizon=6, ntfe=60, ntcp=2)
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2)
    sample_time = 0.5
    # Six samples per horizon, five elements per sample

    initial_plant_inputs = [m_plant.fs.mixer.S_inlet.flow_rate[0],
                            m_plant.fs.mixer.E_inlet.flow_rate[0]]

    nmpc = NMPCSim(m_plant.fs, m_controller.fs, initial_plant_inputs,
            solver=solver, outlvl=idaeslog.DEBUG, 
            sample_time=sample_time)
    # IPOPT output looks a little weird solving for initial conditions here...
    # has non-zero dual infeasibility, iteration 1 has a non-zero
    # regularization coefficient. (Would love to debug this with a
    # transparent NLP solver...)

    # Check that variables have been categorized properly
    p_mod = nmpc.p_mod
    init_input_set = ComponentSet(initial_plant_inputs)
    init_deriv_set = ComponentSet(
            [m_plant.fs.cstr.control_volume.material_accumulation[0, 'aq', 'S'],
             m_plant.fs.cstr.control_volume.material_accumulation[0, 'aq', 'E'],
             m_plant.fs.cstr.control_volume.material_accumulation[0, 'aq', 'C'],
             m_plant.fs.cstr.control_volume.material_accumulation[0, 'aq', 'P'],
             m_plant.fs.cstr.control_volume.energy_accumulation[0, 'aq']])
    init_diff_set = ComponentSet(
            [m_plant.fs.cstr.control_volume.material_holdup[0, 'aq', 'S'],
             m_plant.fs.cstr.control_volume.material_holdup[0, 'aq', 'E'],
             m_plant.fs.cstr.control_volume.material_holdup[0, 'aq', 'C'],
             m_plant.fs.cstr.control_volume.material_holdup[0, 'aq', 'P'],
             m_plant.fs.cstr.control_volume.energy_holdup[0, 'aq']])
    assert len(p_mod.input_vars) == 2
    assert p_mod.input_vars[0][0] in init_input_set
    assert p_mod.input_vars[1][0] in init_input_set
    assert len(p_mod.deriv_vars) == 5
    assert p_mod.deriv_vars[0][0] in init_deriv_set
    assert p_mod.deriv_vars[1][0] in init_deriv_set
    assert p_mod.deriv_vars[2][0] in init_deriv_set
    assert p_mod.deriv_vars[3][0] in init_deriv_set
    assert p_mod.deriv_vars[4][0] in init_deriv_set
    assert len(p_mod.diff_vars) == 5
    assert p_mod.diff_vars[0][0] in init_diff_set
    assert p_mod.diff_vars[1][0] in init_diff_set
    assert p_mod.diff_vars[2][0] in init_diff_set
    assert p_mod.diff_vars[3][0] in init_diff_set
    assert p_mod.diff_vars[4][0] in init_diff_set

    init_alg_set = ComponentSet([
             m_plant.fs.cstr.outlet.conc_mol[0, 'S'],
             m_plant.fs.cstr.outlet.conc_mol[0, 'E'],
             m_plant.fs.cstr.outlet.conc_mol[0, 'C'],
             m_plant.fs.cstr.outlet.conc_mol[0, 'P'],
             m_plant.fs.cstr.outlet.flow_mol_comp[0, 'S'],
             m_plant.fs.cstr.outlet.flow_mol_comp[0, 'E'],
             m_plant.fs.cstr.outlet.flow_mol_comp[0, 'C'],
             m_plant.fs.cstr.outlet.flow_mol_comp[0, 'P'],
             m_plant.fs.cstr.outlet.flow_rate[0],
             m_plant.fs.cstr.outlet.temperature[0],
             m_plant.fs.cstr.inlet.conc_mol[0, 'S'],
             m_plant.fs.cstr.inlet.conc_mol[0, 'E'],
             m_plant.fs.cstr.inlet.conc_mol[0, 'C'],
             m_plant.fs.cstr.inlet.conc_mol[0, 'P'],
             m_plant.fs.cstr.inlet.flow_mol_comp[0, 'S'],
             m_plant.fs.cstr.inlet.flow_mol_comp[0, 'E'],
             m_plant.fs.cstr.inlet.flow_mol_comp[0, 'C'],
             m_plant.fs.cstr.inlet.flow_mol_comp[0, 'P'],
             m_plant.fs.cstr.inlet.flow_rate[0],
             m_plant.fs.cstr.inlet.temperature[0],
             m_plant.fs.mixer.outlet.conc_mol[0, 'S'],
             m_plant.fs.mixer.outlet.conc_mol[0, 'E'],
             m_plant.fs.mixer.outlet.conc_mol[0, 'C'],
             m_plant.fs.mixer.outlet.conc_mol[0, 'P'],
             m_plant.fs.mixer.outlet.flow_mol_comp[0, 'S'],
             m_plant.fs.mixer.outlet.flow_mol_comp[0, 'E'],
             m_plant.fs.mixer.outlet.flow_mol_comp[0, 'C'],
             m_plant.fs.mixer.outlet.flow_mol_comp[0, 'P'],
             m_plant.fs.mixer.outlet.flow_rate[0],
             m_plant.fs.mixer.outlet.temperature[0],
             m_plant.fs.mixer.inlet.flow_mol_comp[0, 'S'],
             m_plant.fs.mixer.inlet.flow_mol_comp[0, 'E'],
             m_plant.fs.mixer.inlet.flow_mol_comp[0, 'C'],
             m_plant.fs.mixer.inlet.flow_mol_comp[0, 'P'],
             m_plant.fs.cstr.control_volume.reactions[0].reaction_coef['R1'],
             m_plant.fs.cstr.control_volume.reactions[0].reaction_coef['R2'],
             m_plant.fs.cstr.control_volume.reactions[0].reaction_coef['R3'],
             m_plant.fs.cstr.control_volume.reactions[0].reaction_rate['R1'],
             m_plant.fs.cstr.control_volume.reactions[0].reaction_rate['R2'],
             m_plant.fs.cstr.control_volume.reactions[0].reaction_rate['R3'],
             m_plant.fs.cstr.control_volume.rate_reaction_generation[0, 'R1'],
             m_plant.fs.cstr.control_volume.rate_reaction_generation[0, 'R2'],
             m_plant.fs.cstr.control_volume.rate_reaction_generation[0, 'R3'],
             m_plant.fs.cstr.control_volume.rate_reaction_extent[0, 'R1'],
             m_plant.fs.cstr.control_volume.rate_reaction_extent[0, 'R2'],
             m_plant.fs.cstr.control_volume.rate_reaction_extent[0, 'R3']
             ])
    assert len(p_mod.alg_vars) == 46
    for i in range(46):
        assert p_mod.alg_vars[i][0] in init_alg_set

    return nmpc


# TODO: dedicated test for validate_models


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_something():
    pass


if __name__ == '__main__':
    m = make_model()
#    assert degrees_of_freedom(m) == 0
#    # This is much slower than the CSTR without mixer. 
#    # Unclear to me what is causing the slowdown.
#    initialize_by_time_element(m.fs, m.fs.time, solver=solver)
#    assert degrees_of_freedom(m) == 0

    test_find_comp_in_block()

#    pdb.set_trace()

    nmpc = test_constructor()
    
    pdb.set_trace()
