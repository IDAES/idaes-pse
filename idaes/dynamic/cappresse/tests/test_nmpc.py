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
                           Set, SolverFactory, Var, value, Objective,
                           TransformationFactory, TerminationCondition)
from pyomo.network import Arc
from pyomo.kernel import ComponentSet
from pyomo.core.expr.visitor import identify_variables

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom, 
        activated_equalities_generator, unfixed_variables_generator)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.generic_models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.dynamic.cappresse import nmpc
from idaes.dynamic.cappresse.nmpc import *
import idaes.logger as idaeslog
from cstr_for_testing import make_model

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes'}
else:
    solver = None


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
    # This tests the same model constructed in the test_nmpc_constructor_1 file
    m_plant = make_model(horizon=6, ntfe=60, ntcp=2)
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
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

    assert hasattr(nmpc, 'p_mod')
    assert hasattr(nmpc, 'c_mod')
    assert hasattr(nmpc, 'sample_time')

    p_mod = nmpc.p_mod
    c_mod = nmpc.c_mod

    assert hasattr(p_mod, 'diff_vars')
    assert hasattr(p_mod, 'deriv_vars')
    assert hasattr(p_mod, 'alg_vars')
    assert hasattr(p_mod, 'input_vars')
    assert hasattr(p_mod, 'fixed_vars')
    assert hasattr(p_mod, 'scalar_vars')
    assert hasattr(p_mod, 'ic_vars')
    assert hasattr(p_mod, 'n_dv')
    assert hasattr(p_mod, 'n_iv')
    assert hasattr(p_mod, 'n_av')

    assert hasattr(c_mod, 'diff_vars')
    assert hasattr(c_mod, 'deriv_vars')
    assert hasattr(c_mod, 'alg_vars')
    assert hasattr(c_mod, 'input_vars')
    assert hasattr(c_mod, 'fixed_vars')
    assert hasattr(c_mod, 'scalar_vars')
    assert hasattr(c_mod, 'ic_vars')
    assert hasattr(c_mod, 'n_dv')
    assert hasattr(c_mod, 'n_iv')
    assert hasattr(c_mod, 'n_av')
    assert hasattr(c_mod, '_ncp')

    # Check that variables have been categorized properly
    ##################
    # In plant model #
    ##################
    assert p_mod is m_plant.fs
    init_input_set = ComponentSet(initial_plant_inputs)

    init_deriv_list = [p_mod.cstr.control_volume.energy_accumulation[0, 'aq']]
    init_diff_list = [p_mod.cstr.control_volume.energy_holdup[0, 'aq']]
    init_fixed_list = [p_mod.cstr.control_volume.volume[0],
                       p_mod.mixer.E_inlet.temperature[0],
                       p_mod.mixer.S_inlet.temperature[0]]

    init_ic_list = [p_mod.cstr.control_volume.energy_holdup[0, 'aq']]

    init_alg_list = [
        p_mod.cstr.outlet.flow_rate[0],
        p_mod.cstr.outlet.temperature[0],
        p_mod.cstr.inlet.flow_rate[0],
        p_mod.cstr.inlet.temperature[0],
        p_mod.mixer.outlet.flow_rate[0],
        p_mod.mixer.outlet.temperature[0]
        ]

    for j in p_mod.properties.component_list:
        init_deriv_list.append(
                p_mod.cstr.control_volume.material_accumulation[0, 'aq', j])
        init_diff_list.append(
                p_mod.cstr.control_volume.material_holdup[0, 'aq', j])
        
        init_fixed_list.append(p_mod.mixer.E_inlet.conc_mol[0, j])
        init_fixed_list.append(p_mod.mixer.S_inlet.conc_mol[0, j])

        init_ic_list.append(
                p_mod.cstr.control_volume.material_holdup[0, 'aq', j])

        init_alg_list.extend([
            p_mod.cstr.outlet.conc_mol[0, j],
            p_mod.cstr.outlet.flow_mol_comp[0, j],
            p_mod.cstr.inlet.conc_mol[0, j],
            p_mod.cstr.inlet.flow_mol_comp[0, j],
            p_mod.cstr.control_volume.rate_reaction_generation[0, 'aq', j],
            p_mod.mixer.outlet.conc_mol[0, j],
            p_mod.mixer.outlet.flow_mol_comp[0, j],
            p_mod.mixer.E_inlet.flow_mol_comp[0, j],
            p_mod.mixer.S_inlet.flow_mol_comp[0, j]
            ])

    for r in p_mod.reactions.rate_reaction_idx:
        init_alg_list.extend([
            p_mod.cstr.control_volume.reactions[0].reaction_coef[r],
            p_mod.cstr.control_volume.reactions[0].reaction_rate[r],
            p_mod.cstr.control_volume.rate_reaction_extent[0, r]
            ])

    init_deriv_set = ComponentSet(init_deriv_list)
    init_diff_set = ComponentSet(init_diff_list)
    init_fixed_set = ComponentSet(init_fixed_list)
    init_ic_set = ComponentSet(init_ic_list)
    init_alg_set = ComponentSet(init_alg_list)

    assert len(p_mod.input_vars) == len(init_input_set)
    for v in p_mod.input_vars:
        assert v[0] in init_input_set

    assert len(p_mod.deriv_vars) == len(init_deriv_set)
    for v in p_mod.deriv_vars:
        assert v[0] in init_deriv_set

    assert len(p_mod.diff_vars) == len(init_deriv_set)
    for v in p_mod.diff_vars:
        assert v[0] in init_diff_set

    assert len(p_mod.fixed_vars) == len(init_fixed_set)
    for v in p_mod.fixed_vars:
        assert v[0] in init_fixed_set

    assert len(p_mod.alg_vars) == len(init_alg_set)
    for v in p_mod.alg_vars:
        assert v[0] in init_alg_set

    assert len(p_mod.ic_vars) == len(init_ic_set)
    for v in p_mod.ic_vars:
        assert v[0] in init_ic_set

    assert len(p_mod.scalar_vars) == 0

    for var in p_mod.deriv_vars:
        assert len(var) == len(p_mod.time)
        assert var.index_set() is p_mod.time
    for var in p_mod.alg_vars:
        assert len(var) == len(p_mod.time)
        assert var.index_set() is p_mod.time

    #######################
    # In controller model #
    #######################
    assert c_mod is m_controller.fs
    init_controller_inputs = [c_mod.mixer.E_inlet.flow_rate[0],
                              c_mod.mixer.S_inlet.flow_rate[0]]
    init_input_set = ComponentSet(init_controller_inputs)

    init_deriv_list = [c_mod.cstr.control_volume.energy_accumulation[0, 'aq']]
    init_diff_list = [c_mod.cstr.control_volume.energy_holdup[0, 'aq']]
    init_fixed_list = [c_mod.cstr.control_volume.volume[0],
                       c_mod.mixer.E_inlet.temperature[0],
                       c_mod.mixer.S_inlet.temperature[0]]

    init_ic_list = [c_mod.cstr.control_volume.energy_holdup[0, 'aq']]

    init_alg_list = [
        c_mod.cstr.outlet.flow_rate[0],
        c_mod.cstr.outlet.temperature[0],
        c_mod.cstr.inlet.flow_rate[0],
        c_mod.cstr.inlet.temperature[0],
        c_mod.mixer.outlet.flow_rate[0],
        c_mod.mixer.outlet.temperature[0]
        ]

    for j in c_mod.properties.component_list:
        init_deriv_list.append(
                c_mod.cstr.control_volume.material_accumulation[0, 'aq', j])
        init_diff_list.append(
                c_mod.cstr.control_volume.material_holdup[0, 'aq', j])
        
        init_fixed_list.append(c_mod.mixer.E_inlet.conc_mol[0, j])
        init_fixed_list.append(c_mod.mixer.S_inlet.conc_mol[0, j])

        init_ic_list.append(
                c_mod.cstr.control_volume.material_holdup[0, 'aq', j])

        init_alg_list.extend([
            c_mod.cstr.outlet.conc_mol[0, j],
            c_mod.cstr.outlet.flow_mol_comp[0, j],
            c_mod.cstr.inlet.conc_mol[0, j],
            c_mod.cstr.inlet.flow_mol_comp[0, j],
            c_mod.cstr.control_volume.rate_reaction_generation[0, 'aq', j],
            c_mod.mixer.outlet.conc_mol[0, j],
            c_mod.mixer.outlet.flow_mol_comp[0, j],
            c_mod.mixer.E_inlet.flow_mol_comp[0, j],
            c_mod.mixer.S_inlet.flow_mol_comp[0, j]
            ])

    for r in c_mod.reactions.rate_reaction_idx:
        init_alg_list.extend([
            c_mod.cstr.control_volume.reactions[0].reaction_coef[r],
            c_mod.cstr.control_volume.reactions[0].reaction_rate[r],
            c_mod.cstr.control_volume.rate_reaction_extent[0, r]
            ])

    init_deriv_set = ComponentSet(init_deriv_list)
    init_diff_set = ComponentSet(init_diff_list)
    init_fixed_set = ComponentSet(init_fixed_list)
    init_alg_set = ComponentSet(init_alg_list)
    init_ic_set = ComponentSet(init_ic_list)

    assert len(c_mod.input_vars) == len(init_input_set)
    for v in c_mod.input_vars:
        assert v[0] in init_input_set
        
    assert len(c_mod.deriv_vars) == len(init_deriv_set)
    for v in c_mod.deriv_vars:
        assert v[0] in init_deriv_set

    assert len(c_mod.diff_vars) == len(init_diff_set)
    for v in c_mod.diff_vars:
        assert v[0] in init_diff_set

    assert len(c_mod.fixed_vars) == len(init_fixed_set)
    for v in c_mod.fixed_vars:
        assert v[0] in init_fixed_set

    assert len(c_mod.alg_vars) == len(init_alg_set)
    for v in c_mod.alg_vars:
        assert v[0] in init_alg_set

    assert len(c_mod.ic_vars) == len(init_ic_set)
    for v in c_mod.ic_vars:
        assert v[0] in init_ic_set

    assert len(c_mod.scalar_vars) == 0

    for var in c_mod.deriv_vars:
        assert len(var) == len(c_mod.time)
        assert var.index_set() is c_mod.time
    for var in c_mod.alg_vars:
        assert len(var) == len(c_mod.time)
        assert var.index_set() is c_mod.time

    return nmpc


def test_validate_setpoint(nmpc, m_steady):

    # set_point is a list of VarData, Value tuples
    # with VarDatas from the controller model,
    # at any point in time.
    # If multiple points in time are given, the last provided is used.
    # TODO: Should probably issue warning if a set point is overwritten in
    # this way...
    c_mod = nmpc.c_mod
    # This is just a random set point I chose, no specific physical
    # meaning, or any reason to think it will be feasible as is.
    set_point = [(c_mod.cstr.outlet.conc_mol[0, 'P'], 0.4),
                 (c_mod.cstr.outlet.conc_mol[0, 'S'], 0.0),
                 (c_mod.cstr.control_volume.energy_holdup[0, 'aq'], 300),
                 (c_mod.mixer.E_inlet.flow_rate[0], 0.1),
                 (c_mod.mixer.S_inlet.flow_rate[0], 2.0)]
    # Interestingly, this (st.st. set point) solve converges infeasible
    # if energy_holdup set point is not 300. (Needs higher weight?)

    weight_tolerance = 5e-7

    # Weight overwrite expects a list of VarData, value tuples
    # in the STEADY MODEL
    # User should /probably/ overwrite weights.
    weight_overwrite = [(m_steady.mixer.E_inlet.flow_rate[0], 20.0)]

    nmpc.validate_steady_setpoint(set_point, m_steady,
            outlvl=idaeslog.DEBUG, 
            solver=solver,
            weight_tolerance=weight_tolerance,
            weight_overwrite=weight_overwrite)

    for con in activated_equalities_generator(m_steady):
        assert value(con.body) - value(con.upper) < 1e-6
        assert value(con.lower) - value(con.body) < 1e-6

    # Validate that steady model has positive objective
    assert hasattr(m_steady, 'user_setpoint_objective')
    assert isinstance(m_steady.user_setpoint_objective, Objective)
    assert value(m_steady.user_setpoint_objective.expr) >= 0

    assert hasattr(m_steady, 'var_locator')

    real_diff_weights = [w for w in m_steady.diff_weights if w is not None]
    real_alg_weights = [w for w in m_steady.alg_weights if w is not None]
    real_input_weights = [w for w in m_steady.input_weights if w is not None]

    assert len(real_diff_weights) == 1
    assert len(real_alg_weights) == 2
    assert len(real_input_weights) == 2

    # Validate that correct bounds have been added
    for j in m_steady.properties.component_list:
        assert (m_steady.cstr.control_volume.material_holdup[0, 'aq', j].lb 
                == 0)

    assert m_steady.mixer.S_inlet.flow_rate[0].lb == 0.5
    assert m_steady.mixer.S_inlet.flow_rate[0].ub == 5
    assert m_steady.mixer.E_inlet.flow_rate[0].lb == 0.01

    # Validate that weights are correct
    C_P_location = m_steady.var_locator\
            [id(m_steady.cstr.outlet.conc_mol[0, 'P'])].location
    assert m_steady.alg_weights[C_P_location] == 1/(0.4-0.001)

    C_S_location = m_steady.var_locator\
            [id(m_steady.cstr.outlet.conc_mol[0, 'S'])].location
    assert m_steady.alg_weights[C_S_location] == 1/(0.001)

    energy_location = m_steady.var_locator\
            [id(m_steady.cstr.control_volume.energy_holdup[0, 'aq'])].location
    assert m_steady.diff_weights[energy_location] == 1/weight_tolerance

    E_inlet_location = m_steady.var_locator\
            [id(m_steady.mixer.E_inlet.flow_rate[0])].location
    assert m_steady.input_weights[E_inlet_location] == weight_overwrite[0][1]

    S_inlet_location = m_steady.var_locator\
            [id(m_steady.mixer.S_inlet.flow_rate[0])].location
    assert abs(m_steady.input_weights[S_inlet_location] - 1/(0.1)) < 1e-6

    # Validate correct objective function
    obj_expr = (m_steady.alg_weights[C_P_location]*
    (m_steady.cstr.outlet.conc_mol[0, 'P'] - m_steady.alg_sp[C_P_location])**2
              + m_steady.alg_weights[C_S_location]*
    (m_steady.cstr.outlet.conc_mol[0, 'S'] - m_steady.alg_sp[C_S_location])**2
              + m_steady.diff_weights[energy_location]*
    (m_steady.cstr.control_volume.energy_holdup[0, 'aq'] - 
                                        m_steady.diff_sp[energy_location])**2
              + m_steady.input_weights[E_inlet_location]*
    (m_steady.mixer.E_inlet.flow_rate[0] - 
                                    m_steady.input_sp[E_inlet_location])**2
              + m_steady.input_weights[S_inlet_location]*
    (m_steady.mixer.S_inlet.flow_rate[0] - 
                                    m_steady.input_sp[S_inlet_location])**2)
    
    assert (abs(value(m_steady.user_setpoint_objective.expr - obj_expr))
            /value(obj_expr)) < 1e-5

    # Validate that set point attributes have been correctly added to 
    # controller model
    assert hasattr(c_mod, 'diff_sp') 
    assert hasattr(c_mod, 'input_sp')
    for i, var in enumerate(m_steady.diff_vars):
        assert var[0].value == c_mod.diff_sp[i]
    for i, var in enumerate(m_steady.input_vars):
        assert var[0].value == c_mod.input_sp[i]


def test_construct_objective_weight_matrices(nmpc):
    
    c_mod = nmpc.c_mod
    dynamic_weight_tol = 5e-7
    dynamic_weight_overwrite = \
            [(nmpc.c_mod.cstr.control_volume.energy_holdup[0, 'aq'], 0.1)]

    nmpc.construct_objective_weight_matrices(nmpc.c_mod,
            overwrite=dynamic_weight_overwrite,
            tol=dynamic_weight_tol)

    # Validate that attributes were added
    assert hasattr(c_mod, 'diff_weights')
    assert hasattr(c_mod, 'input_weights')
    assert len(c_mod.diff_weights) == len(c_mod.diff_vars)
    assert len(c_mod.input_weights) == len(c_mod.input_weights)

    # Validate that weights are as expected
    for i, _slice in enumerate(c_mod.diff_vars):
        if not _slice[0].local_name.startswith('energy_holdup'):
            assert (c_mod.diff_weights[i] == 
                    1/abs(_slice[0].value - c_mod.diff_sp[i]))
        else:
            assert c_mod.diff_weights[i] == dynamic_weight_overwrite[0][1]


def test_add_objective_function(nmpc):

    c_mod = nmpc.c_mod
    time = c_mod.time
    nmpc.add_objective_function(c_mod,
            control_penalty_type='action',
            name='tracking_objective')

    # Validate that something called 'tracking_objective' has been added
    assert hasattr(c_mod, 'tracking_objective')
    assert value(c_mod.tracking_objective.expr) > 0

    obj_state_term = sum(sum(c_mod.diff_weights[i]*(var[t] - c_mod.diff_sp[i])**2
                        for i, var in enumerate(c_mod.diff_vars))
                        for t in time)

    obj_control_term = sum(sum(c_mod.input_weights[i]*(var[time[k]] - var[time[k-1]])**2
                        for i, var in enumerate(c_mod.input_vars))
                        for k in range(2, len(time)+1))

    obj_expr = obj_state_term + obj_control_term

    assert value(obj_expr) == value(c_mod.tracking_objective.expr)
    # Controller model has not been initialized yet, so value of
    # objective function may not be meaningful


def test_add_pwc_constraints(nmpc):
    sample_time = 0.5
    nmpc.add_pwc_constraints(sample_time=sample_time)

    c_mod = nmpc.c_mod

    assert nmpc.sample_time == sample_time
    assert nmpc.c_mod._samples_per_horizon == 6

    # Test that components were added
    assert hasattr(c_mod, '_pwc_input_0') 
    assert hasattr(c_mod, '_pwc_input_1') 

    # Test that constraints have the correct indexing set
    n_sample = int(c_mod.time.last()/sample_time)
    sample_points = [sample_time*i
            for i in range(n_sample+1)]

    for t in sample_points:
        assert t not in c_mod._pwc_input_0
        assert t not in c_mod._pwc_input_1

    # Rough test the the constraints are correct - contain the correct 
    # variables
    for i, t in enumerate(c_mod.time):
        if t not in sample_points:
            t_next = c_mod.time[i+2]
            var_in_0 = [id(v) for v in 
                    identify_variables(c_mod._pwc_input_0[t].expr)]
            var_in_1 = [id(v) for v in 
                    identify_variables(c_mod._pwc_input_1[t].expr)]
            assert len(var_in_0) == 2
            assert len(var_in_1) == 2
            assert id(c_mod.input_vars[0][t]) in var_in_0
            assert id(c_mod.input_vars[0][t_next]) in var_in_0
            assert id(c_mod.input_vars[1][t]) in  var_in_1
            assert id(c_mod.input_vars[1][t_next]) in var_in_1


def test_initialization_by_simulation(nmpc):

    nmpc.initialize_control_problem(strategy='from_simulation')

    c_mod = nmpc.c_mod
    time = c_mod.time

    # Validate that model has been correctly unfixed.
    # (At least at non-initial time)
    for _slice in c_mod.input_vars + c_mod.diff_vars + c_mod.alg_vars:
        for t in time:
            if t != time.first():
                assert not _slice[t].fixed

    # Check for correct dof
    assert (degrees_of_freedom(c_mod) == 
            c_mod.n_iv*(c_mod._samples_per_horizon + 1))
    # The +1 is to account for the inputs at the initial conditions,
    # which maybe should be fixed...

    for con in activated_equalities_generator(c_mod):
        assert abs(value(con.body) - value(con.upper)) < 1e-6


def test_initialization_from_initial_conditions(nmpc):

    dof_before = degrees_of_freedom(nmpc.c_mod)

    nmpc.initialize_control_problem(strategy='initial_conditions')

    dof_after = degrees_of_freedom(nmpc.c_mod)
    assert dof_after == dof_before

    c_mod = nmpc.c_mod
    locator = c_mod.var_locator
    time = c_mod.time
    t0 = time.first()

    # Check that expected value copying was performed
    for _slice in c_mod.diff_vars + c_mod.alg_vars + c_mod.deriv_vars:
        for t in time:
            assert _slice[t] == _slice[time.first()]

    # Expect only violated equalities to be accumulation equations
    # ^ This is false, as equalities involving inputs could be violated too
    for con in activated_equalities_generator(c_mod):
        # If the equality does not contain any inputs, it should
        # only be violated if it is an accumulation equation
        if abs(value(con.body) - value(con.upper)) > 1e-6:
            if not any([locator[id(v)].category == 'input' 
                        for v in identify_variables(con.expr)]):
                assert 'accumulation' in con.local_name


def test_solve_control_problem(nmpc):

    c_mod = nmpc.c_mod

    init_obj_value = value(c_mod.tracking_objective.expr)
    nmpc.solve_control_problem()
    final_obj_value = value(c_mod.tracking_objective.expr)

    # Not always true because initial model might not be feasible
    assert final_obj_value < init_obj_value

    for con in activated_equalities_generator(c_mod):
        assert abs(value(con.body) - value(con.upper)) < 1e-6

    for var in unfixed_variables_generator(c_mod):
        if var.lb is not None:
            assert var.lb - var.value < 1e-6
        if var.ub is not None:
            assert var.value - var.ub < 1e-6


def test_inject_inputs(nmpc):
    
    c_mod = nmpc.c_mod
    p_mod = nmpc.p_mod
    sample_time = nmpc.sample_time
    time = p_mod.time
    nmpc.inject_inputs_into(p_mod, c_mod, t_src=2, t_tgt=0)
    # Here I am copying the inputs at the incorrect time
    # (t=2 instead of t=0.5, one sampling time) just to explicitly
    # test both functions inject_inputs_into and inject_inputs_into_plant

    for i, _slice in enumerate(p_mod.input_vars):
        for t in time:
            if t > sample_time or t == 0:
                continue
            assert _slice[t].value == c_mod.input_vars[i][2].value

    nmpc.inject_inputs_into_plant(0)
    for i, _slice in enumerate(p_mod.input_vars):
        for t in time:
            if t > sample_time or t == 0:
                continue
            assert _slice[t].value == c_mod.input_vars[i][sample_time].value


def test_simulate_over_range(nmpc):

    p_mod = nmpc.p_mod
    c_mod = nmpc.c_mod
    time = p_mod.time
    sample_time = nmpc.sample_time

#    was_violated = {id(con):
#            abs(value(con.body)-value(con.upper))>=1e-6 
#            for con in activated_equalities_generator(p_mod)}
    # ^ Can't calculate value because many variables are not initialized

    assert degrees_of_freedom(p_mod) == 0
    nmpc.simulate_over_range(p_mod, 0, 3)
    assert degrees_of_freedom(p_mod) == 0

    for con in activated_equalities_generator(p_mod):
        if 'disc_eq' in con.local_name or 'balance' in con.local_name:
            # (know disc. and balance equations will be directly indexed
            # by time)
            index = con.index()
            if not type(index) is tuple:
                index = (index,)
            t_index = index[0]
            if t_index <= 3:
                # Equalities in simulated range should not be violated
                assert abs(value(con.body)-value(con.upper)) < 1e-6
#            else:
#                # Equalities outside simulated range should be violated
#                # if they were violated before
#                if was_violated[id(con)]:
#                    assert abs(value(con.body)-value(con.upper)) >= 1e-6

    nmpc.simulate_plant(0)

    # Check that plant simulation matches controller simulation.
    # Only valid because there is no noise or plant-model-mismatch
    # and plant/controller have the same time discretizations
    p_varlist = p_mod.diff_vars + p_mod.alg_vars + p_mod.deriv_vars
    c_varlist = c_mod.diff_vars + c_mod.alg_vars + c_mod.deriv_vars
    for i, pvar in enumerate(p_varlist):
        for t in time:
            if t > sample_time or t == 0:
                continue
            cvar = c_varlist[i]
            assert abs(pvar[t].value - cvar[t].value) < 1e-5


def test_initialize_from_previous(nmpc):
    c_mod = nmpc.c_mod
    time = c_mod.time
    sample_time = nmpc.sample_time

    assert nmpc.controller_solved

    c_varlist = c_mod.diff_vars + c_mod.alg_vars + c_mod.deriv_vars

    prev_values = [{t: _slice[t].value
                       for t in time}
                       for _slice in c_varlist]
    nmpc.initialize_from_previous_sample(c_mod)

    for i, cvar in enumerate(c_varlist):
        for t in time:
            if time.last() - t < sample_time or t == time.first():
                continue
            t_next = t + sample_time
            # Addition with ContinuousSet indices can result in rounding
            # errors. Here round to eighth decimal place.
            t_next = round(1e8*t_next)/1e8
            assert t_next in time
            assert cvar[t].value == prev_values[i][t_next]


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

    nmpc = test_constructor()

    m_steady = make_model(steady=True)

    # Validate attributes and lack of bounds in controller model.

    # Test helper functions for validation
    # - make sure they can catch errors

    # Next major test is to test the addition of set point
    # steady model/solve -> weights -> objective function
    
    test_validate_setpoint(nmpc, m_steady.fs)

    nmpc.s_mod = m_steady.fs

    test_construct_objective_weight_matrices(nmpc)

    test_add_objective_function(nmpc)

    test_add_pwc_constraints(nmpc)

    test_initialization_by_simulation(nmpc)

    test_initialization_from_initial_conditions(nmpc)

    test_solve_control_problem(nmpc)

    test_inject_inputs(nmpc)

    test_simulate_over_range(nmpc)

    test_initialize_from_previous(nmpc)
