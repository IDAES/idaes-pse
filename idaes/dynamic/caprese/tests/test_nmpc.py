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
from pytest import approx
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
from idaes.dynamic.caprese import nmpc
from idaes.dynamic.caprese.nmpc import *
from idaes.dynamic.caprese.util import *
import idaes.logger as idaeslog
from cstr_for_testing import make_model
import random

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
solver_available = SolverFactory('ipopt').available()
if solver_available:
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes'}
else:
    solver = None


def assert_categorization(model):
    init_input_set = ComponentSet([model.mixer.S_inlet.flow_rate[0],
                                   model.mixer.E_inlet.flow_rate[0]])

    init_deriv_list = [model.cstr.control_volume.energy_accumulation[0, 'aq']]
    init_diff_list = [model.cstr.control_volume.energy_holdup[0, 'aq']]
    init_fixed_list = [model.cstr.control_volume.volume[0],
                       model.mixer.E_inlet.temperature[0],
                       model.mixer.S_inlet.temperature[0]]

    init_ic_list = [model.cstr.control_volume.energy_holdup[0, 'aq']]

    init_alg_list = [
        model.cstr.outlet.flow_rate[0],
        model.cstr.outlet.temperature[0],
        model.cstr.inlet.flow_rate[0],
        model.cstr.inlet.temperature[0],
        model.mixer.outlet.flow_rate[0],
        model.mixer.outlet.temperature[0]
        ]

    for j in model.properties.component_list:
        init_deriv_list.append(
                model.cstr.control_volume.material_accumulation[0, 'aq', j])
        init_diff_list.append(
                model.cstr.control_volume.material_holdup[0, 'aq', j])
        
        init_fixed_list.append(model.mixer.E_inlet.conc_mol[0, j])
        init_fixed_list.append(model.mixer.S_inlet.conc_mol[0, j])

        init_ic_list.append(
                model.cstr.control_volume.material_holdup[0, 'aq', j])

        init_alg_list.extend([
            model.cstr.outlet.conc_mol[0, j],
            model.cstr.outlet.flow_mol_comp[0, j],
            model.cstr.inlet.conc_mol[0, j],
            model.cstr.inlet.flow_mol_comp[0, j],
            model.cstr.control_volume.rate_reaction_generation[0, 'aq', j],
            model.mixer.outlet.conc_mol[0, j],
            model.mixer.outlet.flow_mol_comp[0, j],
            model.mixer.E_inlet.flow_mol_comp[0, j],
            model.mixer.S_inlet.flow_mol_comp[0, j]
            ])

    for r in model.reactions.rate_reaction_idx:
        init_alg_list.extend([
            model.cstr.control_volume.reactions[0].reaction_coef[r],
            model.cstr.control_volume.reactions[0].reaction_rate[r],
            model.cstr.control_volume.rate_reaction_extent[0, r]
            ])

    init_deriv_set = ComponentSet(init_deriv_list)
    init_diff_set = ComponentSet(init_diff_list)
    init_fixed_set = ComponentSet(init_fixed_list)
    init_ic_set = ComponentSet(init_ic_list)
    init_alg_set = ComponentSet(init_alg_list)

    assert model._NMPC_NAMESPACE.input_vars.n_vars == len(init_input_set)
    for v in model._NMPC_NAMESPACE.input_vars:
        assert v[0] in init_input_set

    assert model._NMPC_NAMESPACE.deriv_vars.n_vars == len(init_deriv_set)
    for v in model._NMPC_NAMESPACE.deriv_vars:
        assert v[0] in init_deriv_set

    assert len(model._NMPC_NAMESPACE.diff_vars) == len(init_deriv_set)
    for v in model._NMPC_NAMESPACE.diff_vars:
        assert v[0] in init_diff_set

    assert len(model._NMPC_NAMESPACE.fixed_vars) == len(init_fixed_set)
    for v in model._NMPC_NAMESPACE.fixed_vars:
        assert v[0] in init_fixed_set

    assert len(model._NMPC_NAMESPACE.alg_vars) == len(init_alg_set)
    for v in model._NMPC_NAMESPACE.alg_vars:
        assert v[0] in init_alg_set

    assert len(model._NMPC_NAMESPACE.ic_vars) == len(init_ic_set)
    for v in model._NMPC_NAMESPACE.ic_vars:
        assert v[0] in init_ic_set

    assert len(model._NMPC_NAMESPACE.scalar_vars) == 0

    for var in model._NMPC_NAMESPACE.deriv_vars:
        assert len(var) == len(model._NMPC_NAMESPACE.get_time())
        assert var.index_set() is model._NMPC_NAMESPACE.get_time()
    for var in model._NMPC_NAMESPACE.alg_vars:
        assert len(var) == len(model._NMPC_NAMESPACE.get_time())
        assert var.index_set() is model._NMPC_NAMESPACE.get_time()


@pytest.fixture(scope='session')
def nmpc():
    # This tests the same model constructed in the test_nmpc_constructor_1 file
    m_plant = make_model(horizon=6, ntfe=60, ntcp=2)
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
    sample_time = 0.5
    # Six samples per horizon, five elements per sample

    initial_plant_inputs = [m_plant.fs.mixer.S_inlet.flow_rate[0],
                            m_plant.fs.mixer.E_inlet.flow_rate[0]]

    nmpc = NMPCSim(m_plant.fs, m_plant.fs.time,
            m_controller.fs, m_controller.fs.time, 
            inputs_at_t0=initial_plant_inputs,
            solver=solver, outlvl=idaeslog.DEBUG, 
            sample_time=sample_time)
    # IPOPT output looks a little weird solving for initial conditions here...
    # has non-zero dual infeasibility, iteration 1 has a non-zero
    # regularization coefficient. (Would love to debug this with a
    # transparent NLP solver...)

    assert hasattr(nmpc, 'p_mod')
    assert hasattr(nmpc, 'p_mod_time')
    assert hasattr(nmpc, 'c_mod')
    assert hasattr(nmpc, 'c_mod_time')
    assert hasattr(nmpc, 'sample_time')

    p_mod = nmpc.p_mod
    p_time = nmpc.p_mod_time
    c_mod = nmpc.c_mod
    c_time = nmpc.c_mod_time

    assert hasattr(p_mod, '_NMPC_NAMESPACE')
    assert hasattr(c_mod, '_NMPC_NAMESPACE')

    assert hasattr(p_mod._NMPC_NAMESPACE, 'diff_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'deriv_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'alg_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'input_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'fixed_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'scalar_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'ic_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'n_diff_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'n_deriv_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'n_input_vars')
    assert hasattr(p_mod._NMPC_NAMESPACE, 'n_alg_vars')

    assert hasattr(c_mod._NMPC_NAMESPACE, 'diff_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'deriv_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'alg_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'input_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'fixed_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'scalar_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'ic_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'n_diff_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'n_deriv_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'n_input_vars')
    assert hasattr(c_mod._NMPC_NAMESPACE, 'n_alg_vars')

    # Check that variables have been categorized properly
    assert_categorization(c_mod)
    assert_categorization(p_mod)

    return nmpc


def test_calculate_full_state_setpoint(nmpc):
    c_mod = nmpc.c_mod

    # Deactivate tracking objective from previous tests
    #c_mod._NMPC_NAMESPACE.tracking_objective.deactivate()
    
    set_point = [(c_mod.cstr.outlet.conc_mol[0, 'P'], 0.4),
                 (c_mod.cstr.outlet.conc_mol[0, 'S'], 0.0),
                 (c_mod.cstr.control_volume.energy_holdup[0, 'aq'], 300),
                 (c_mod.mixer.E_inlet.flow_rate[0], 0.1),
                 (c_mod.mixer.S_inlet.flow_rate[0], 2.0)]

    c_mod.mixer.E_inlet.flow_rate[0].fix(0.1)
    c_mod.mixer.S_inlet.flow_rate[0].fix(2.0)

    weight_tolerance = 5e-7
    weight_override = [
            (c_mod.mixer.E_inlet.flow_rate[0.], 20.),
            (c_mod.mixer.S_inlet.flow_rate[0.], 2.),
            (c_mod.cstr.control_volume.energy_holdup[0., 'aq'], 0.1),
            (c_mod.cstr.outlet.conc_mol[0., 'P'], 1.),
            (c_mod.cstr.outlet.conc_mol[0., 'S'], 1.),
            ]

    nmpc.calculate_full_state_setpoint(set_point,
            objective_weight_tolerance=weight_tolerance,
            objective_weight_override=weight_override)

    assert hasattr(c_mod._NMPC_NAMESPACE, 'user_setpoint')
    user_setpoint = c_mod._NMPC_NAMESPACE.user_setpoint
    assert hasattr(c_mod._NMPC_NAMESPACE, 'user_setpoint_weights')
    user_setpoint_weights = c_mod._NMPC_NAMESPACE.user_setpoint_weights
    assert hasattr(c_mod._NMPC_NAMESPACE, 'user_setpoint_vars')
    user_setpoint_vars = c_mod._NMPC_NAMESPACE.user_setpoint_vars

    for i, var in enumerate(user_setpoint_vars):
        if var.local_name.startswith('conc'):
            assert user_setpoint_weights[i] == 1.
        elif var.local_name.startswith('energy'):
            assert user_setpoint_weights[i] == 0.1
        elif var.local_name.startswith('E_'):
            assert user_setpoint_weights[i] == 20.
        elif var.local_name.startswith('S_'):
            assert user_setpoint_weights[i] == 2.
        
    alg_vars = c_mod._NMPC_NAMESPACE.alg_vars
    diff_vars = c_mod._NMPC_NAMESPACE.diff_vars
    input_vars = c_mod._NMPC_NAMESPACE.input_vars
    categories = [
            VariableCategory.DIFFERENTIAL,
            VariableCategory.ALGEBRAIC,
            VariableCategory.DERIVATIVE,
            VariableCategory.INPUT,
            ]
    category_dict = c_mod._NMPC_NAMESPACE.category_dict
    for categ in categories:
        group = category_dict[categ]
        # Assert that setpoint has been populated with non-None values
        assert not any([sp is None for sp in group.setpoint])
        # Assert that setpoint (target) and reference (initial) values are 
        # different in some way
        assert not all([sp == ref for sp, ref in 
            zip(group.setpoint, group.reference)])
        # Assert that initial and reference values are the same
        assert all([ref == var[0].value for ref, var in 
            zip(group.reference, group.varlist)])

def test_add_setpoint_to_controller(nmpc):
    c_mod = nmpc.c_mod
    weight_override = []
    for j in c_mod.properties.component_list:
        weight_override.append(
                (c_mod.cstr.control_volume.material_holdup[0, 'aq', j], 1))
    weight_override.append(
            (c_mod.cstr.control_volume.energy_holdup[0, 'aq'], 0.1))
    weight_override.append((c_mod.mixer.E_inlet.flow_rate[0], 2))
    weight_override.append((c_mod.mixer.S_inlet.flow_rate[0], 0.2))

    state_categories = [VariableCategory.DIFFERENTIAL]
    nmpc.add_setpoint_to_controller(
            objective_weight_override=weight_override,
            objective_state_categories=state_categories,
            time_resolution_option=TimeResolutionOption.FINITE_ELEMENTS,
            objective_name='test_objective')

    diff_vars = c_mod._NMPC_NAMESPACE.diff_vars
    input_vars = c_mod._NMPC_NAMESPACE.input_vars
    diff_weights = diff_vars.weights
    diff_sp = diff_vars.setpoint
    input_weights = input_vars.weights
    input_sp = input_vars.setpoint
    for i, var in enumerate(c_mod._NMPC_NAMESPACE.diff_vars):
        if var[0].local_name.startswith('material_holdup'):
            assert diff_weights[i] == 1
        elif var[0].local_name.startswith('energy_holdup'):
            assert diff_weights[i] == 0.1
    for i, var in enumerate(c_mod._NMPC_NAMESPACE.input_vars):
        if var[0].local_name.startswith('E'):
            assert input_weights[i] == 2.
        elif var[0].local_name.startswith('S'):
            assert input_weights[i] == 0.2

    # Then assert that objective has correct value
    time = list(c_mod._NMPC_NAMESPACE.get_time().get_finite_elements())
    time.pop(0)

    obj_state_term = sum(sum(diff_weights[i]*(var[t] - diff_sp[i])**2
                        for i, var in enumerate(diff_vars))
                        for t in time)

    obj_control_term = sum(sum(input_weights[i]*(var[time[k]] - 
                        var[time[k-1]])**2
                        for i, var in enumerate(input_vars))
                        for k in range(1, len(time)))

    obj_expr = obj_state_term + obj_control_term

    assert hasattr(c_mod._NMPC_NAMESPACE, 'test_objective')
    assert (value(obj_expr) == 
            approx(value(c_mod._NMPC_NAMESPACE.test_objective.expr), 1e-6))

    c_mod._NMPC_NAMESPACE.test_objective.deactivate()

def test_construct_objective_weights(nmpc):
    
    c_mod = nmpc.c_mod
    dynamic_weight_tol = 5e-7
    dynamic_weight_overwrite = \
            [(nmpc.c_mod.cstr.control_volume.energy_holdup[0, 'aq'], 0.1)]

    nmpc.construct_objective_weights(nmpc.c_mod,
            objective_weight_override=dynamic_weight_overwrite,
            objective_weight_tolerance=dynamic_weight_tol)

    diff_weights = c_mod._NMPC_NAMESPACE.diff_vars

    # Validate that weights are as expected
    category_dict = c_mod._NMPC_NAMESPACE.category_dict
    for categ, group in category_dict.items():
        if categ == VariableCategory.SCALAR or categ == VariableCategory.FIXED:
            continue
        for i, var in enumerate(group):
            if var[0].local_name.startswith('energy_holdup'):
                assert (group.weights[i] == dynamic_weight_overwrite[0][1])
            else:
                diff = max(abs(group.reference[i]-group.setpoint[i]),
                        dynamic_weight_tol)
                assert (group.weights[i] == 1/diff)


def test_add_objective_function(nmpc):

    c_mod = nmpc.c_mod
    nmpc.add_objective_function(c_mod,
            control_penalty_type=ControlPenaltyType.ACTION,
            time_resolution_option=TimeResolutionOption.SAMPLE_POINTS,
            name='tracking_objective')

    # Validate that something called 'tracking_objective' has been added
    assert hasattr(c_mod._NMPC_NAMESPACE, 'tracking_objective')
    assert value(c_mod._NMPC_NAMESPACE.tracking_objective.expr) > 0

    diff_vars = c_mod._NMPC_NAMESPACE.diff_vars
    input_vars = c_mod._NMPC_NAMESPACE.input_vars
    diff_weights = diff_vars.weights
    diff_sp = diff_vars.setpoint
    input_weights = input_vars.weights
    input_sp = input_vars.setpoint

    time = c_mod._NMPC_NAMESPACE.sample_points

    obj_state_term = sum(sum(diff_weights[i]*(var[t] - diff_sp[i])**2
                        for i, var in enumerate(diff_vars))
                        for t in time)

    obj_control_term = sum(sum(input_weights[i]*(var[time[k]] - 
                        var[time[k-1]])**2
                        for i, var in enumerate(input_vars))
                        for k in range(1, len(time)))

    obj_expr = obj_state_term + obj_control_term

    assert (value(obj_expr) == 
            approx(value(c_mod._NMPC_NAMESPACE.tracking_objective.expr), 1e-6))
    # Controller model has not been initialized yet, so value of
    # objective function may not be meaningful


def test_constrain_control_inputs_piecewise_constant(nmpc):
    sample_time = 0.5
    nmpc.constrain_control_inputs_piecewise_constant(sample_time=sample_time)

    c_mod = nmpc.c_mod

    assert nmpc.sample_time == sample_time
    assert nmpc.c_mod._NMPC_NAMESPACE.samples_per_horizon == 6

    # Test that components were added
    assert hasattr(c_mod._NMPC_NAMESPACE, 'pwc_constraint')

    # Test that constraints have the correct indexing set
    n_sample = int(c_mod.time.last()/sample_time)
    sample_points = [sample_time*i
            for i in range(1, n_sample+1)]
    # By convention, sample_points omits time.first()

    assert (sample_points == nmpc.c_mod._NMPC_NAMESPACE.sample_points)

    for t in sample_points:
        for i in range(c_mod._NMPC_NAMESPACE.n_input_vars):
            assert (t, i) not in c_mod._NMPC_NAMESPACE.pwc_constraint
            # ^ tuple because pwc_constraint is now indexed by time and the location
            # into the input list

    # Rough test that the constraints are correct - contain the correct 
    # variables.
    time = c_mod._NMPC_NAMESPACE.get_time()
    for i, t in enumerate(c_mod._NMPC_NAMESPACE.get_time()):
        if t not in sample_points:
            t_next = time[i+2]
            var_in_0 = [id(v) for v in 
                identify_variables(c_mod._NMPC_NAMESPACE.pwc_constraint[t, 0].expr)]
            var_in_1 = [id(v) for v in 
                identify_variables(c_mod._NMPC_NAMESPACE.pwc_constraint[t, 1].expr)]
            assert len(var_in_0) == 2
            assert len(var_in_1) == 2
            assert (id(c_mod._NMPC_NAMESPACE.input_vars.varlist[0][t]) 
                    in var_in_0)
            assert (id(c_mod._NMPC_NAMESPACE.input_vars.varlist[0][t_next]) 
                    in var_in_0)
            assert (id(c_mod._NMPC_NAMESPACE.input_vars.varlist[1][t]) 
                    in var_in_1)
            assert (id(c_mod._NMPC_NAMESPACE.input_vars.varlist[1][t_next]) 
                    in var_in_1)


@pytest.mark.skipif(not solver_available, reason='IPOPT is not available')
def test_initialization_by_time_element(nmpc):

    nmpc.initialize_control_problem(
            control_init_option=ControlInitOption.BY_TIME_ELEMENT)

    c_mod = nmpc.c_mod
    time = c_mod.time
    input_vars = c_mod._NMPC_NAMESPACE.input_vars
    diff_vars = c_mod._NMPC_NAMESPACE.diff_vars
    alg_vars = c_mod._NMPC_NAMESPACE.alg_vars

    # Validate that model has been correctly unfixed.
    # (At least at non-initial time)
    for _slice in input_vars.varlist + diff_vars.varlist + alg_vars.varlist:
        for t in time:
            if t != time.first():
                assert not _slice[t].fixed

    # Check for correct dof
    assert (degrees_of_freedom(c_mod) == 
            c_mod._NMPC_NAMESPACE.n_input_vars*
            c_mod._NMPC_NAMESPACE.samples_per_horizon)
    # The +1 is to account for the inputs at the initial conditions,
    # which maybe should be fixed...
    # ^These inputs will now remain fixed, as time.first() is not
    # a sample point.

    for con in activated_equalities_generator(c_mod):
        # Don't require pwc constraints to be satisfied,
        # as they were not active during initialization
        if not con.local_name.startswith('pwc'):
            assert value(con.body) == approx(value(con.upper), abs=1e-6)


def test_initialization_from_initial_conditions(nmpc):

    dof_before = degrees_of_freedom(nmpc.c_mod)

    nmpc.initialize_control_problem(
            control_init_option=ControlInitOption.FROM_INITIAL_CONDITIONS)

    dof_after = degrees_of_freedom(nmpc.c_mod)
    assert dof_after == dof_before

    c_mod = nmpc.c_mod
    locator = c_mod._NMPC_NAMESPACE.var_locator
    time = c_mod._NMPC_NAMESPACE.get_time()
    t0 = time.first()

    deriv_vars = c_mod._NMPC_NAMESPACE.deriv_vars
    diff_vars = c_mod._NMPC_NAMESPACE.diff_vars
    alg_vars = c_mod._NMPC_NAMESPACE.alg_vars

    # Check that expected value copying was performed
    for _slice in diff_vars.varlist + alg_vars.varlist + deriv_vars.varlist:
        for t in time:
            assert _slice[t] == _slice[t0]

    # Expect only violated equalities to be accumulation equations
    # ^ This is false, as equalities involving inputs could be violated too
    for con in activated_equalities_generator(c_mod):
        # If the equality does not contain any inputs, it should
        # only be violated if it is an accumulation equation
        if abs(value(con.body) - value(con.upper)) > 1e-6:
            if not any([locator[v].category == VariableCategory.INPUT 
                        for v in identify_variables(con.expr)]):
                assert 'accumulation' in con.local_name

@pytest.mark.skipif(not solver_available, reason='IPOPT unavailable')
def test_solve_control_problem(nmpc):
    c_mod = nmpc.c_mod

    init_obj_value = value(c_mod._NMPC_NAMESPACE.tracking_objective.expr)
    nmpc.solve_control_problem()
    final_obj_value = value(c_mod._NMPC_NAMESPACE.tracking_objective.expr)

    # Not always true because initial model might not be feasible
    assert final_obj_value < init_obj_value

    for con in activated_equalities_generator(c_mod):
        assert abs(value(con.body) - value(con.upper)) < 1e-6

    for var in unfixed_variables_generator(c_mod):
        if var.lb is not None:
            assert var.lb - var.value < 1e-6
        if var.ub is not None:
            assert var.value - var.ub < 1e-6


def test_inject_control_inputs(nmpc):

    c_mod = nmpc.c_mod
    p_mod = nmpc.p_mod
    sample_time = nmpc.sample_time
    time = p_mod.time
   # nmpc.inject_inputs_into(p_mod, c_mod, t_src=2, t_tgt=0)
   # # Here I am copying the inputs at the incorrect time
   # # (t=2 instead of t=0.5, one sampling time) just to explicitly
   # # test both functions inject_inputs_into and inject_control_inputs_into_plant

   # for i, _slice in enumerate(p_mod.input_vars):
   #     for t in time:
   #         if t > sample_time or t == 0:
   #             continue
   #         assert _slice[t].value == c_mod.input_vars[i][2].value

    nmpc.inject_control_inputs_into_plant(0, add_input_noise=False)
    for i, _slice in enumerate(p_mod._NMPC_NAMESPACE.input_vars):
        for t in time:
            if t > sample_time or t == 0:
                continue
            assert _slice[t].value == \
                c_mod._NMPC_NAMESPACE.input_vars.varlist[i][sample_time].value


def test_initialize_by_element_in_range(nmpc):

    p_mod = nmpc.p_mod
    c_mod = nmpc.c_mod
    time = p_mod.time
    sample_time = nmpc.sample_time

#    was_violated = {id(con):
#            abs(value(con.body)-value(con.upper))>=1e-6 
#            for con in activated_equalities_generator(p_mod)}
    # ^ Can't calculate value because many variables are not initialized

    assert degrees_of_freedom(p_mod) == 0
    initialize_by_element_in_range(p_mod, time, 0, 3,
            dae_vars=p_mod._NMPC_NAMESPACE.dae_vars,
            time_linking_vars=p_mod._NMPC_NAMESPACE.diff_vars.varlist)
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
    p_varlist = (p_mod._NMPC_NAMESPACE.diff_vars.varlist + 
                 p_mod._NMPC_NAMESPACE.alg_vars.varlist + 
                 p_mod._NMPC_NAMESPACE.deriv_vars.varlist)
    c_varlist = (c_mod._NMPC_NAMESPACE.diff_vars.varlist + 
                 c_mod._NMPC_NAMESPACE.alg_vars.varlist + 
                 c_mod._NMPC_NAMESPACE.deriv_vars.varlist)
    for i, pvar in enumerate(p_varlist):
        for t in time:
            if t > sample_time or t == 0:
                continue
            cvar = c_varlist[i]
            assert pvar[t].value == approx(cvar[t].value, abs=1e-5)


def test_calculate_error_between_states(nmpc):
    c_mod = nmpc.c_mod
    p_mod = nmpc.p_mod
    error1 = nmpc.calculate_error_between_states(c_mod, p_mod, 0.5, 0.5)
    c_varlist = c_mod._NMPC_NAMESPACE.diff_vars.varlist
    p_varlist = p_mod._NMPC_NAMESPACE.diff_vars.varlist
    for c_var, p_var in zip(c_varlist, p_varlist):
        assert c_var[0.5].value == approx(p_var[0.5].value, abs=1e-5)
    assert error1 == approx(0., abs=1e-5)


def test_initialize_from_previous(nmpc):
    c_mod = nmpc.c_mod
    time = c_mod.time
    sample_time = nmpc.sample_time

    assert nmpc.controller_solved

    c_varlist = (c_mod._NMPC_NAMESPACE.diff_vars.varlist + 
                 c_mod._NMPC_NAMESPACE.alg_vars.varlist + 
                 c_mod._NMPC_NAMESPACE.deriv_vars.varlist)

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


def test_transfer_current_plant_state_to_controller(nmpc):
    c_mod = nmpc.c_mod
    random.seed(12345)
    nmpc.transfer_current_plant_state_to_controller(3., add_plant_noise=True)

    p_ic_vars = nmpc.p_mod._NMPC_NAMESPACE.controller_ic_vars
    c_ic_vars = nmpc.c_mod._NMPC_NAMESPACE.ic_vars

    for p_var, c_var in zip(p_ic_vars, c_ic_vars):
        assert p_var[3.].value == approx(c_var[0.].value, rel=0.1)


#def test_calculate_full_state_setpoint(nmpc):
#    c_mod = nmpc.c_mod
#
#    # Deactivate tracking objective from previous tests
#    c_mod._NMPC_NAMESPACE.tracking_objective.deactivate()
#    
#    set_point = [(c_mod.cstr.outlet.conc_mol[0, 'P'], 0.4),
#                 (c_mod.cstr.outlet.conc_mol[0, 'S'], 0.0),
#                 (c_mod.cstr.control_volume.energy_holdup[0, 'aq'], 300),
#                 (c_mod.mixer.E_inlet.flow_rate[0], 0.1),
#                 (c_mod.mixer.S_inlet.flow_rate[0], 2.0)]
#
#    c_mod.mixer.E_inlet.flow_rate[0].fix(0.1)
#    c_mod.mixer.S_inlet.flow_rate[0].fix(2.0)
#
#    weight_tolerance = 5e-7
#    weight_override = [
#            (c_mod.mixer.E_inlet.flow_rate[0.], 20.),
#            (c_mod.mixer.S_inlet.flow_rate[0.], 2.),
#            (c_mod.cstr.control_volume.energy_holdup[0., 'aq'], 0.1),
#            (c_mod.cstr.outlet.conc_mol[0., 'P'], 1.),
#            (c_mod.cstr.outlet.conc_mol[0., 'S'], 1.),
#            ]
#
#    nmpc.calculate_full_state_setpoint(set_point,
#            objective_weight_tolerance=weight_tolerance,
#            objective_weight_override=weight_override)
#
#    assert hasattr(c_mod._NMPC_NAMESPACE, 'user_setpoint')
#    user_setpoint = c_mod._NMPC_NAMESPACE.user_setpoint
#    assert hasattr(c_mod._NMPC_NAMESPACE, 'user_setpoint_weights')
#    user_setpoint_weights = c_mod._NMPC_NAMESPACE.user_setpoint_weights
#    assert hasattr(c_mod._NMPC_NAMESPACE, 'user_setpoint_vars')
#    user_setpoint_vars = c_mod._NMPC_NAMESPACE.user_setpoint_vars
#
#    for i, var in enumerate(user_setpoint_vars):
#        if var.local_name.startswith('conc'):
#            assert user_setpoint_weights[i] == 1.
#        elif var.local_name.startswith('energy'):
#            assert user_setpoint_weights[i] == 0.1
#        elif var.local_name.startswith('E_'):
#            assert user_setpoint_weights[i] == 20.
#        elif var.local_name.startswith('S_'):
#            assert user_setpoint_weights[i] == 2.
#        
#    alg_vars = c_mod._NMPC_NAMESPACE.alg_vars
#    diff_vars = c_mod._NMPC_NAMESPACE.diff_vars
#    input_vars = c_mod._NMPC_NAMESPACE.input_vars
#    categories = [
#            VariableCategory.DIFFERENTIAL,
#            VariableCategory.ALGEBRAIC,
#            VariableCategory.DERIVATIVE,
#            VariableCategory.INPUT,
#            ]
#    category_dict = c_mod._NMPC_NAMESPACE.category_dict
#    for categ in categories:
#        group = category_dict[categ]
#        # Assert that setpoint has been populated with non-None values
#        assert not any([sp is None for sp in group.setpoint])
#        # Assert that setpoint (target) and reference (initial) values are 
#        # different in some way
#        assert not all([sp == ref for sp, ref in 
#            zip(group.setpoint, group.reference)])
#        # Assert that initial and reference values are the same
#        assert all([ref == var[0].value for ref, var in 
#            zip(group.reference, group.varlist)])


@pytest.fixture
def dynamic_cstr_model():
    return make_model().fs


@pytest.fixture
def steady_cstr_model():
    return make_model(steady=True).fs


