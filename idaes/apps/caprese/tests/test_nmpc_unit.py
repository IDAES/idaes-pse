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
Test for Caprese's module for NMPC.
"""

import pyomo.environ as aml
import pyomo.dae as dae
import pyomo.network as pyn
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom,
        activated_equalities_generator, unfixed_variables_generator)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.apps.caprese import nmpc
from idaes.apps.caprese.nmpc import *
from idaes.apps.caprese.util import *
from idaes.apps.caprese.common.config import (
        VariableCategory,
        ElementInitializationInputOption,
        ControlInitOption,
        )
from idaes.apps.caprese.tests.test_model import (
        make_model, 
        make_small_model,
        initialize_t0,
        copy_values_forward,
        )
import idaes.logger as idaeslog
import random
import pytest

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
solver_available = SolverFactory('ipopt').available()
if solver_available:
    solver = SolverFactory('ipopt')
#    solver.options = {'tol': 1e-8,
#                      'mu_init': 1e-8,
#                      'bound_push': 1e-8,
#                      'halt_on_ampl_error': 'yes'}
# Want to run with default options.
else:
    solver = None


class TestNMPCSim(object):
    def make_nmpc(self, sample_time=0.5, horizon=1, nfe=2):
        plant = make_model(horizon=horizon, nfe=nfe)
        plant_time = plant.time
        plant_inputs = [plant.flow_in[plant_time.first()]]
        controller = make_model(horizon=horizon, nfe=nfe)
        controller_time = controller.time
        nmpc = NMPCSim(
                plant_model=plant,
                plant_time_set=plant_time,
                controller_model=controller,
                controller_time_set=controller_time,
                inputs_at_t0=plant_inputs,
                sample_time=sample_time,
                )
        return nmpc
    
    @pytest.mark.unit
    def test_constructor(self):
        """
        Categorization is tested rigorously in the test_nmpc_constructor
        files. This test is very basic by comparison
        """
        nmpc = self.make_nmpc()
        assert not nmpc.controller_solved
        assert nmpc.current_plant_time == 0
        assert hasattr(nmpc.controller, nmpc.namespace_name)
        assert hasattr(nmpc.plant, nmpc.namespace_name)
        controller_namespace = getattr(nmpc.controller, nmpc.namespace_name)
        plant_namespace = getattr(nmpc.plant, nmpc.namespace_name)

        locator_name = 'var_locator'
        category_name = 'category_dict'
        assert hasattr(controller_namespace, locator_name)
        assert hasattr(plant_namespace, locator_name)
        assert hasattr(controller_namespace, category_name)
        assert hasattr(plant_namespace, category_name)

    @pytest.mark.unit
    def test_add_namespace_to(self):
        nmpc = self.make_nmpc()
        model = make_small_model()
        nmpc.add_namespace_to(model, model.time)
        assert hasattr(model, nmpc.namespace_name)

    @pytest.mark.unit
    def test_validate_sample_time(self):
        nmpc = self.make_nmpc()
        model = make_small_model()
        sample_time = 0.5
        nmpc.add_namespace_to(model, model.time)
        nmpc.validate_sample_time(sample_time, model)
        namespace = getattr(model, nmpc.namespace_name)
        assert hasattr(namespace, 'sample_points')
        assert hasattr(namespace, 'fe_per_sample')

        sample_point_set = set(namespace.sample_points)
        # NOTE: My implementation has wavered between whether or not 
        # time.first() is a sample point. I now assert, hopefully for
        # good, that it is.
        for p in [0.0, 0.5, 1]:
            assert p in sample_point_set

    @pytest.mark.unit
    def test_validate_slices(self):
        nmpc = self.make_nmpc()
        plant_namespace = getattr(nmpc.plant, nmpc.namespace_name)
        controller_namespace = getattr(nmpc.controller, nmpc.namespace_name)
        src_slices = [
                plant_namespace.diff_vars[0],
                plant_namespace.alg_vars[0],
                ]
        src_names = [_slice[0].name for _slice in src_slices]
        tgt_names = [
                nmpc.controller.conc[0,'A'].name,
                nmpc.controller.flow_out[0].name,
                ]
        # This will be violated if component construction order (and
        # therefore categorization) is changed.
        assert src_names == tgt_names
        tgt_slices = nmpc.validate_slices(
                nmpc.controller, 
                nmpc.plant, 
                nmpc.plant_time, 
                src_slices)
        tgt_slice_names = [_slice[0].name for _slice in tgt_slices]
        assert tgt_names == tgt_slice_names

    @pytest.mark.unit
    def test_validate_fixedness(self):
        nmpc = self.make_nmpc()
        nmpc.validate_fixedness(nmpc.plant, nmpc.controller)
        nmpc.plant.conc[0,'A'].unfix()
        with pytest.raises(AssertionError):
            nmpc.validate_fixedness(nmpc.plant)
        nmpc.plant.conc[0,'A'].fix()
        t1 = nmpc.plant_time.get_finite_elements()[1]
        nmpc.plant.conc_in[t1, 'A'].unfix()
        with pytest.raises(AssertionError):
            nmpc.validate_fixedness(nmpc.plant)

    @pytest.mark.unit
    def test_transfer_current_plant_state_to_controller(self):
        nmpc = self.make_nmpc()
        ts = nmpc.plant_time.first() + nmpc.sample_time
        nmpc.plant.conc[ts,'A'].set_value(4)
        nmpc.plant.conc[ts,'B'].set_value(5)
        nmpc.transfer_current_plant_state_to_controller(ts, 
                add_plant_noise=False)
        t0 = nmpc.controller_time.first()
        assert nmpc.controller.conc[t0,'A'].value == 4
        assert nmpc.controller.conc[t0,'B'].value == 5

        random.seed(123)
        nmpc.transfer_current_plant_state_to_controller(ts, 
                add_plant_noise=True,
                noise_weights=[1,1])
        assert nmpc.controller.conc[t0,'A'].value != 4
        assert nmpc.controller.conc[t0,'B'].value != 5

    @pytest.mark.unit
    def test_inject_control_inputs_into_plant(self):
        nmpc = self.make_nmpc()
        n_samples = 2
        sample_inputs = {1: 3, 2: 6}
        plant_sample_points = [nmpc.plant_time.first() + i*nmpc.sample_time
                for i in range(n_samples+1)]
        plant_sample = {i+1: [t for t in nmpc.plant_time 
            if t > plant_sample_points[i] and t <= plant_sample_points[i+1]]
            for i in range(n_samples)}
        tc = nmpc.controller_time.first() + nmpc.sample_time

        for i in range(n_samples):
            nmpc.controller.flow_in[tc].set_value(sample_inputs[i+1])
            nmpc.inject_control_inputs_into_plant(plant_sample_points[i])
            for t in plant_sample[i+1]:
                assert nmpc.plant.flow_in[t] == sample_inputs[i+1]
            nmpc.inject_control_inputs_into_plant(plant_sample_points[i],
                    add_input_noise=True,
                    noise_weights=[1])
            for t in plant_sample[i+1]:
                assert nmpc.plant.flow_in[t] != sample_inputs[i+1]

    @pytest.mark.unit
    def test_has_consistent_initial_conditions(self):
        nmpc = self.make_nmpc()
        with pytest.raises(ValueError):
            # Model has not been properly initialized
            nmpc.has_consistent_initial_conditions(nmpc.plant)

        t0 = nmpc.plant_time.first()
        nmpc.plant.rate[t0,:].set_value(0.0)
        nmpc.plant.flow_out[t0].set_value(nmpc.plant.flow_in[t0].value)
        for j in nmpc.plant.components:
            calculate_variable_from_constraint(
                    nmpc.plant.dcdt[t0,j],
                    nmpc.plant.material_balance[t0,j],
                    )
        assert nmpc.has_consistent_initial_conditions(nmpc.plant)

        nmpc.plant.flow_out[t0].set_value(0.0)
        assert not nmpc.has_consistent_initial_conditions(nmpc.plant)

    @pytest.mark.unit
    def test_solve_consistent_initial_conditions(self):
        # solve, assert has_consistent..., assert same dof
        nmpc = self.make_nmpc()
        dof_prior = degrees_of_freedom(nmpc.plant)

        # Just going to assume this model needs no additional
        # initialization before solving for initial conditions.
        nmpc.solve_consistent_initial_conditions(nmpc.plant)
        assert nmpc.has_consistent_initial_conditions(nmpc.plant)
        dof_post = degrees_of_freedom(nmpc.plant)
        assert dof_prior == dof_post

        nmpc.controller.flow_in.setlb(0.0)
        nmpc.controller.flow_out.setlb(0.0)
        nmpc.controller.flow_in.setub(10.0)
        nmpc.controller.flow_out.setub(10.0)
        nmpc.solve_consistent_initial_conditions(nmpc.controller,
                strip_bounds=False)
        assert nmpc.has_consistent_initial_conditions(nmpc.controller)

    @pytest.mark.unit
    def test_calculate_full_state_setpoint(self):
        # come up with setpoint, solve, check data structures,
        # make sure values of setpoint make sense.
        nmpc = self.make_nmpc()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        t0 = nmpc.controller_time.first()
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        categories = [VariableCategory.DIFFERENTIAL, VariableCategory.ALGEBRAIC]
        for categ in categories:
            group = namespace.category_dict[categ]
            for i, var in enumerate(group):
                assert group.setpoint[i] is not None

        assert namespace.diff_vars.setpoint[0] == pytest.approx(3.75)
        assert namespace.diff_vars.setpoint[1] == pytest.approx(1.25)
        assert namespace.alg_vars.setpoint[0] == pytest.approx(3.0)
        assert namespace.alg_vars.setpoint[1] == pytest.approx(-3.75)
        assert namespace.alg_vars.setpoint[2] == pytest.approx(3.75)
        assert namespace.input_vars.setpoint[0] == pytest.approx(3.0)

    @pytest.mark.unit
    def test_solve_setpoint(self):
        # create objective function on model
        # get dof
        # call solve_setpoint
        # compare dof
        # compare values
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        mod = nmpc.controller
        mod.objective = Objective(expr=
                (mod.flow_in[t0] - 3.0)**2 + (mod.conc[t0,'A'] - 5)**2)
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)

        dof_prior = degrees_of_freedom(mod)
        nmpc.solve_setpoint(require_steady=False, outlvl=idaeslog.DEBUG)
        dof_post = degrees_of_freedom(mod)

        assert dof_prior == dof_post
        assert namespace.diff_vars[0][t0].value == pytest.approx(5.0, abs=1e-3)
        assert namespace.diff_vars[1][t0].value == pytest.approx(7.11, abs=1e-3)
        assert namespace.alg_vars[0][t0].value == pytest.approx(3.0, abs=1e-3)
        assert namespace.alg_vars[1][t0].value == pytest.approx(-5.0, abs=1e-3)
        assert namespace.alg_vars[2][t0].value == pytest.approx(5.0, abs=1e-3)
        assert namespace.input_vars[0][t0].value == pytest.approx(3.0, abs=1e-3)

    @pytest.mark.unit
    def test_add_setpoint_to_controller(self):
        # This test just creates an objective then makes sure it has the 
        # correct variables. More detailed tests will be performed for
        # construct_objective_weights and add_objective_function.
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        categories = [VariableCategory.DIFFERENTIAL]
        nmpc.add_setpoint_to_controller(objective_state_categories=categories)
        assert hasattr(namespace, 'tracking_objective')
        obj_variables = ComponentSet(
                identify_variables(namespace.tracking_objective.expr))

        n_samples = len(namespace.sample_points) - 1
        assert len(obj_variables) == (len(namespace.diff_vars) + 
                len(namespace.input_vars))*n_samples
        for var in namespace.diff_vars.varlist + namespace.input_vars.varlist:
            for ts in namespace.sample_points:
                if ts == t0:
                    continue
                assert var[ts] in obj_variables

    @pytest.mark.unit
    def test_set_reference_values_from_initial(self):
        nmpc = self.make_nmpc()
        time = nmpc.plant_time
        t0 = time.first()
        namespace = getattr(nmpc.plant, nmpc.namespace_name)
        diff_vars = namespace.diff_vars
        nmpc.set_reference_values_from_initial(diff_vars)
        for i, var in enumerate(diff_vars):
            assert var[t0].value == diff_vars.reference[i]

    @pytest.mark.unit
    def test_construct_objective_weights(self):
        """
        Populates setpoint values,

        """
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        
        reference = namespace.diff_vars.reference
        setpoint = namespace.diff_vars.setpoint
        n_diff_vars = len(namespace.diff_vars)
        predicted_diff_weights = [
                1./abs(reference[i] - setpoint[i]) for i in range(n_diff_vars)
                ]
        tol = 1e-5
        nmpc.construct_objective_weights(
                nmpc.controller,
                categories=[
                    VariableCategory.DIFFERENTIAL, 
                    VariableCategory.INPUT,
                    ],
                objective_weight_tolerance=tol,
                )
        for pred, act in zip(predicted_diff_weights,
                namespace.diff_vars.weights):
            assert pred == act
        assert namespace.input_vars.weights[0] == 1/tol

        override = [
                (nmpc.controller.rate[t0,'A'], 2),
                (nmpc.controller.rate[t0,'B'], 2),
                (nmpc.controller.flow_out[t0], 2),
                ]
        nmpc.construct_objective_weights(
                nmpc.controller,
                categories=[
                    VariableCategory.ALGEBRAIC,
                    ],
                objective_weight_override=override,
                )
        for w in namespace.alg_vars.weights:
            assert w == 2

    @pytest.mark.unit
    def test_add_objective_function(self):
        # populate setpoints,
        # construct weights (with override everywhere)
        # add objective function
        # compare objective function to prediction
        # for diff+error, diff+action, diff+alg+action setpoints
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        copy_values_forward(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        override = [
                (nmpc.controller.conc[t0,'A'], 1),
                (nmpc.controller.conc[t0,'B'], 1),
                (nmpc.controller.rate[t0,'A'], 1),
                (nmpc.controller.rate[t0,'B'], 1),
                (nmpc.controller.flow_out[t0], 1),
                (nmpc.controller.flow_in[t0], 1),
                ]
        nmpc.construct_objective_weights(
                nmpc.controller,
                objective_weight_override=override,
                )
        pred_obj = {i: 0. for i in range(1, 4)}
        sample_points = namespace.sample_points[1:]
        setpoint = namespace.diff_vars.setpoint
        for i, v in enumerate(namespace.diff_vars):
            for t in sample_points:
                pred_obj[1] += (v[t] - setpoint[i])**2
                pred_obj[2] += (v[t] - setpoint[i])**2
                pred_obj[3] += (v[t] - setpoint[i])**2
        setpoint = namespace.alg_vars.setpoint
        for i, v in enumerate(namespace.alg_vars):
            for t in sample_points:
                pred_obj[3] += (v[t] - setpoint[i])**2
        setpoint = namespace.input_vars.setpoint
        for i, v in enumerate(namespace.input_vars):
            for t in sample_points:
                t_prev = t - nmpc.sample_time
                pred_obj[1] += (v[t] - setpoint[i])**2
                pred_obj[2] += (v[t] - v[t_prev])**2
                pred_obj[3] += (v[t] - v[t_prev])**2

        nmpc.add_objective_function(
                nmpc.controller,
                control_penalty_type=ControlPenaltyType.ERROR,
                objective_state_categories=[
                    VariableCategory.DIFFERENTIAL,
                    ]
                )
        assert aml.value(namespace.objective.expr == pred_obj[1])
        namespace.del_component(namespace.objective)

        nmpc.add_objective_function(
                nmpc.controller,
                control_penalty_type=ControlPenaltyType.ACTION,
                objective_state_categories=[
                    VariableCategory.DIFFERENTIAL,
                    ]
                )
        assert aml.value(namespace.objective.expr == pred_obj[2])
        namespace.del_component(namespace.objective)

        nmpc.add_objective_function(
                nmpc.controller,
                control_penalty_type=ControlPenaltyType.ACTION,
                objective_state_categories=[
                    VariableCategory.DIFFERENTIAL,
                    VariableCategory.ALGEBRAIC,
                    ]
                )
        assert aml.value(namespace.objective.expr == pred_obj[3])
        namespace.del_component(namespace.objective)

    @pytest.mark.unit
    def test_set_bounds_from_initial(self):
        # change bounds on some VarGroup
        # copy forward in time
        # check at all t > t0
        nmpc = self.make_nmpc()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        time = nmpc.controller_time
        t0 = time.first()
        for var in namespace.input_vars:
            var[t0].domain = aml.NonNegativeReals
            var[t0].setub(10.)
        for var in namespace.diff_vars:
            var[t0].domain = aml.NonNegativeReals
            var[t0].setub(100.)
        namespace.diff_vars.varlist[0].setlb(1.)
        namespace.diff_vars.varlist[1].setlb(-1.)

        for var in namespace.scalar_vars:
            var.setlb(0.1)
            var.setub(11.)
        nmpc.set_bounds_from_initial(namespace.input_vars)
        nmpc.set_bounds_from_initial(namespace.diff_vars)
        nmpc.set_bounds_from_initial(namespace.scalar_vars)
        # ^ Won't actually do anything for scalar vars...
        for t in time:
            if t == t0:
                continue
            for var in namespace.input_vars:
                assert var[t].ub == var[t0].ub
                assert var[t].lb == 0.0
            for var in namespace.diff_vars:
                assert var[t].ub == var[t0].ub
            assert namespace.diff_vars.varlist[0][t].lb == 1.0
            assert namespace.diff_vars.varlist[1][t].lb == 0.0

    @pytest.mark.unit
    def test_constrain_control_inputs_piecewise_constant(self):
        # construct nmpc
        # constrain_control_inputs
        # assert hasattr pwc_constraint_list, pwc_constraint
        # for each, assert has correct variables
        nmpc = self.make_nmpc()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        time = nmpc.controller_time
        t0 = time.first()
        sample_points = namespace.sample_points
        sample_point_set = set(sample_points)
        nmpc.constrain_control_inputs_piecewise_constant()
        assert hasattr(namespace, 'pwc_constraint')
        assert hasattr(namespace, 'pwc_constraint_list')
        for var, con in zip(namespace.input_vars, 
                namespace.pwc_constraint_list):
            for t in time:
                if t in sample_point_set:
                    continue
                t_next = time.next(t)
                incident_vars = ComponentSet(
                        identify_variables(con[t].expr))
                assert var[t] in incident_vars
                assert var[t_next] in incident_vars

    @pytest.mark.unit
    def test_initialize_control_problem(self):
        # Each initialization strategy should be tested in more detail
        # in other functions.
        # This function should test that each strategy can be passed to
        # the user-facing function, and that the values left in the
        # controller model are somewhat reasonable
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        namespace = getattr(nmpc.controller, nmpc.namespace_name)
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        copy_values_forward(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        override = [
                (nmpc.controller.conc[t0,'A'], 1),
                (nmpc.controller.conc[t0,'B'], 1),
                (nmpc.controller.rate[t0,'A'], 1),
                (nmpc.controller.rate[t0,'B'], 1),
                (nmpc.controller.flow_out[t0], 1),
                (nmpc.controller.flow_in[t0], 1),
                ]
        nmpc.add_setpoint_to_controller(objective_weight_override=override)
        nmpc.constrain_control_inputs_piecewise_constant()

        nmpc.initialize_control_problem(
                control_init_option=ControlInitOption.FROM_INITIAL_CONDITIONS)
        nmpc.initialize_control_problem(
                control_init_option=ControlInitOption.BY_TIME_ELEMENT)
        with pytest.raises(RuntimeError):
            nmpc.initialize_control_problem(
                    control_init_option=ControlInitOption.FROM_PREVIOUS)

    @pytest.mark.unit
    def test_initialize_by_solving_elements(self):
        # Make nmpc
        # add setpoint to controller
        # add pwc_constraints
        # call initialize_by_solving_elements
        # assert that values are as expected.
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        tl = time.last()
        controller = nmpc.controller
        namespace = getattr(controller, nmpc.namespace_name)
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        copy_values_forward(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        override = [
                (nmpc.controller.conc[t0,'A'], 1),
                (nmpc.controller.conc[t0,'B'], 1),
                (nmpc.controller.rate[t0,'A'], 1),
                (nmpc.controller.rate[t0,'B'], 1),
                (nmpc.controller.flow_out[t0], 1),
                (nmpc.controller.flow_in[t0], 1),
                ]
        nmpc.add_setpoint_to_controller(objective_weight_override=override)
        nmpc.constrain_control_inputs_piecewise_constant()
        nmpc.controller.flow_in[:].set_value(2.0)
        nmpc.initialize_by_solving_elements(
                controller,
                time,
                input_type=ElementInitializationInputOption.INITIAL,
                )
        input_vars = namespace.input_vars
        diff_vars = namespace.diff_vars
        deriv_vars = namespace.deriv_vars
        assert input_vars[0][tl].value == 2.0
        assert diff_vars[0][tl].value == pytest.approx(3.185595567867036)
        assert diff_vars[1][tl].value == pytest.approx(1.1532474073395755)
        assert deriv_vars[0][tl].value == pytest.approx(0.44321329639889284)
        assert deriv_vars[1][tl].value == pytest.approx(0.8791007531878847)

        nmpc.initialize_by_solving_elements(
                controller,
                time,
                input_type=ElementInitializationInputOption.SET_POINT,
                )
        for t in time:
            if t != t0:
                assert input_vars[0][t] == 3.0
        assert diff_vars[0][tl].value == pytest.approx(3.7037037037037037)
        assert diff_vars[1][tl].value == pytest.approx(1.0746896480968502)
        assert deriv_vars[0][tl].value == pytest.approx(0.1851851851851849)
        assert deriv_vars[1][tl].value == pytest.approx(0.47963475941315314)

    @pytest.mark.unit
    def test_initialize_from_previous_sample(self):
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        tl = time.last()
        controller = nmpc.controller
        namespace = getattr(controller, nmpc.namespace_name)
        sample_points = namespace.sample_points
        t1 = sample_points[1]
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        # One way to deal with this pitfall is to add bounds to inlet flow.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        copy_values_forward(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        override = [
                (nmpc.controller.conc[t0,'A'], 1),
                (nmpc.controller.conc[t0,'B'], 1),
                (nmpc.controller.rate[t0,'A'], 1),
                (nmpc.controller.rate[t0,'B'], 1),
                (nmpc.controller.flow_out[t0], 1),
                (nmpc.controller.flow_in[t0], 1),
                ]
        nmpc.add_setpoint_to_controller(objective_weight_override=override)
        nmpc.constrain_control_inputs_piecewise_constant()
        nmpc.solve_control_problem()
        
        diff_vars = namespace.diff_vars
        input_vars = namespace.input_vars
        sample_time = nmpc.sample_time
        n_samples = namespace.samples_per_horizon
        cat_dict = namespace.category_dict
        categories = [
                VariableCategory.DIFFERENTIAL,
                VariableCategory.INPUT,
                ]
        expected = {categ: {s: [{t-sample_time: 
            cat_dict[categ][i][t].value for t in time
            if sample_points[s] < t and t <= sample_points[s+1]}
            for i in range(len(cat_dict[categ]))]
            for s in range(1, n_samples)}
            for categ in categories}
        for categ in categories:
            expected[categ][n_samples] = [{t: cat_dict[categ].setpoint[i] 
                for t in time if sample_points[n_samples-1] < t 
                and t <= sample_points[n_samples]}
                for i in range(len(cat_dict[categ]))]

        nmpc.initialize_from_previous_sample(nmpc.controller)

        for s in range(1, n_samples):
            interval = [t for t in time 
                    if sample_points[s-1] < t and t <= sample_points[s]]
            for categ in categories:
                for i, var in enumerate(cat_dict[categ]):
                    for t in interval:
                        assert var[t].value == expected[categ][s][i][t]

    @pytest.mark.unit
    def test_initialize_from_initial_conditions(self):
        # make nmpc
        # add setpoint to controller
        # call initialize
        # assert that values are same as at t0
        nmpc = self.make_nmpc()
        controller = nmpc.controller
        time = nmpc.controller_time
        t0 = time.first()
        namespace = getattr(controller, nmpc.namespace_name)
        sample_points = namespace.sample_points
        n_samples = namespace.samples_per_horizon
        samples = {i: [t for t in time 
            if sample_points[i-1] < t and t <= sample_points[i]]
            for i in range(1, n_samples)}
        sample_val_dict = {
                1: 2.5,
                2: 1.5,
                }
        for i in samples:
            for t in samples[i]:
                nmpc.controller.flow_in[t].set_value(sample_val_dict[i])
        initialize_t0(nmpc.controller)
        nmpc.initialize_from_initial_conditions(nmpc.controller)
        categories = [
                VariableCategory.DIFFERENTIAL,
                VariableCategory.ALGEBRAIC,
                VariableCategory.DERIVATIVE,
                ]
        cat_dict = namespace.category_dict
        for categ in categories:
            for var in cat_dict[categ]:
                for t in time:
                    assert var[t].value == var[t0].value
        # Assert inputs were not changed
        for i in samples:
            for t in samples[i]:
                var = nmpc.controller.flow_in[t]
                assert var.value == sample_val_dict[i]

    @pytest.mark.unit
    def test_solve_control_problem(self):
        # make nmpc
        # add setpoint
        # initialize, from ICs probably
        # solve control problem
        # assert values are as expected.
        # Could probably even do this without initialization
        nmpc = self.make_nmpc()
        time = nmpc.controller_time
        t0 = time.first()
        tl = time.last()
        controller = nmpc.controller
        namespace = getattr(controller, nmpc.namespace_name)
        t1 = namespace.sample_points[1]
        setpoint = [
                (nmpc.controller.flow_in[t0], 3.0),
                ]
        # Note: without this initialization, the setpoint arrives at a
        # local optimum with inlet flow ~= 0.
        # One way to deal with this pitfall is to add bounds to inlet flow.
        nmpc.controller.flow_in[:].set_value(3.0)
        initialize_t0(nmpc.controller)
        copy_values_forward(nmpc.controller)
        nmpc.calculate_full_state_setpoint(setpoint, outlvl=idaeslog.DEBUG)
        override = [
                (nmpc.controller.conc[t0,'A'], 1),
                (nmpc.controller.conc[t0,'B'], 1),
                (nmpc.controller.rate[t0,'A'], 1),
                (nmpc.controller.rate[t0,'B'], 1),
                (nmpc.controller.flow_out[t0], 1),
                (nmpc.controller.flow_in[t0], 1),
                ]
        nmpc.add_setpoint_to_controller(objective_weight_override=override)
        nmpc.constrain_control_inputs_piecewise_constant()
        nmpc.solve_control_problem()
        input_vars = namespace.input_vars
        diff_vars = namespace.diff_vars
        assert input_vars[0][t1].value == pytest.approx(3.192261151432352)
        assert input_vars[0][tl].value == pytest.approx(2.9818775607191648)
        assert diff_vars[0][tl].value == pytest.approx(3.7101450012850137)
        assert diff_vars[1][tl].value == pytest.approx(1.0898406680173942)
        assert nmpc.controller_solved

    @pytest.mark.unit
    def test_simulate_plant(self):
        # populate plant with some random inputs
        # simulate
        # assert values are as expected
        # assert 0 dof
        nmpc = self.make_nmpc()
        time = nmpc.plant_time
        namespace = getattr(nmpc.plant, nmpc.namespace_name)
        sample_points = namespace.sample_points
        t1 = sample_points[1]
        n_samples = namespace.samples_per_horizon
        samples = {i: [t for t in time 
            if sample_points[i-1] < t and t <= sample_points[i]]
            for i in range(1, n_samples)}
        for t in samples[1]:
            nmpc.plant.flow_in[t].set_value(2.5)
        nmpc.simulate_plant(t_start=nmpc.plant_time.first())
        assert degrees_of_freedom(nmpc.plant) == 0
        input_vars = namespace.input_vars
        diff_vars = namespace.diff_vars
        assert input_vars[0][t1].value == 2.5
        assert diff_vars[0][t1].value == pytest.approx(3.0155642023346307)
        assert diff_vars[1][t1].value == pytest.approx(0.5914009717947231)

    @pytest.mark.unit
    def test_calculate_error_between_states(self):
        # make nmpc
        # artificially change values
        # calculate error.
        nmpc = self.make_nmpc()
        initialize_t0(nmpc.controller)
        copy_values_forward(nmpc.controller)
        initialize_t0(nmpc.controller)
        copy_values_forward(nmpc.controller)

        nmpc.controller.conc[1, 'A'].set_value(5)
        nmpc.controller.conc[1, 'B'].set_value(5)
        nmpc.plant.conc[1, 'A'].set_value(3)
        nmpc.plant.conc[1, 'B'].set_value(3)
        Q_matrix = [1,1]
        error = nmpc.calculate_error_between_states(
                nmpc.controller,
                nmpc.plant,
                1,
                1,
                Q_matrix=Q_matrix,
                )
        assert error == 1*2**2 + 1*2**2
