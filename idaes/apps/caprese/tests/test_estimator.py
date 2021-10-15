#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
""" Tests for the estimator model subclass of block
"""

import pytest
import pyomo.environ as pyo

from idaes.apps.caprese.categorize import (
        categorize_dae_variables_and_constraints,
        )

from idaes.apps.caprese.tests.test_simple_model import (
        make_model,
        initialize_t0,
        )

from idaes.apps.caprese.estimator import (
        _EstimatorBlockData,
        EstimatorBlock,
        SimpleEstimatorBlock,
        IndexedEstimatorBlock,
        )
from idaes.apps.caprese.common.config import (
        VariableCategory,
        ConstraintCategory,
        )
from idaes.apps.caprese.dynamic_var import (
        DiffVar,
        DerivVar,
        AlgVar,
        InputVar,
        ModelDisturbanceVar,
        MeasurementErrorVar,
        ActualMeasurementVar,
        )
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.common.collections import ComponentSet
from pyomo.dae.flatten import (
        flatten_components_along_sets,
        flatten_dae_components,
        )

__author__ = "Kuan-Han Lin"

solver_available = pyo.SolverFactory('ipopt').available()
if solver_available:
    solver = pyo.SolverFactory('ipopt')
else:
    solver = None
    
    
class TestEstimatorBlock(object):
    
    @pytest.mark.unit
    def test_construct(self, sample_time=0.5, horizon=1., nfe=2):
        model = make_model(horizon=horizon, nfe=nfe)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[0,'A']]
        estimator = EstimatorBlock(
                model=model,
                time=time,
                inputs=inputs,
                measurements=measurements,
                sample_time=sample_time,
                )
        
        assert type(estimator) is SimpleEstimatorBlock
        assert isinstance(estimator, EstimatorBlock)
        assert isinstance(estimator, _EstimatorBlockData)
        
        estimator.construct()
        assert estimator[None] is estimator
        assert estimator.mod is model
        assert estimator.time is time
        assert estimator._sample_time
        assert all(i1 is i2 for i1, i2 in zip(estimator._inputs, inputs))
        assert all(i1 is i2 for i1, i2 in zip(estimator._measurements, measurements))
        assert hasattr(estimator, "category_dict")
        assert hasattr(estimator, "con_category_dict")
        assert hasattr(estimator, "MHE_VARS_CONS_BLOCK")
        MHEBlock = estimator.MHE_VARS_CONS_BLOCK
        assert hasattr(MHEBlock, "MEASUREMENT_SET")
        assert len(MHEBlock.MEASUREMENT_SET.ordered_data()) == len(measurements)
        assert hasattr(MHEBlock, "DIFFERENTIAL_SET")
        assert len(MHEBlock.DIFFERENTIAL_SET.ordered_data()) == \
            len(estimator.category_dict[VariableCategory.DIFFERENTIAL])
        
    @pytest.mark.unit
    def test_construct_indexed(self, sample_time=0.5):#, horizon=1., nfe=2):
        estimator_set = pyo.Set(initialize=[0, 1, 2])
        estimator_set.construct()
        horizon_map = {0: 1., 1: 3., 2: 5.}
        nfe_map = {0: 2, 1: 6, 2: 10}
        model_map = {i: make_model(horizon_map[i], nfe_map[i])
                for i in estimator_set}
        time_map = {i: model_map[i].time for i in estimator_set}
        inputs_map = {i: [model_map[i].flow_in[0]] for i in estimator_set}
        measurements_map = {i: [model_map[i].conc[0, 'A']] for i in estimator_set}
        sample_time_map = {i: sample_time for i in estimator_set}
        
        estimator = EstimatorBlock(
                estimator_set,
                model=model_map,
                time=time_map,
                inputs=inputs_map,
                measurements=measurements_map,
                sample_time=sample_time_map,
                )
        
        assert type(estimator) is IndexedEstimatorBlock
        assert isinstance(estimator, IndexedEstimatorBlock)
        estimator.construct()
        assert all(b.parent_component() is estimator for b in estimator.values())

        for i in estimator_set:
            assert i in estimator
        for i, c in estimator.items():
            assert c.mod is model_map[i]
            assert c.time is time_map[i]
            assert c._sample_time is sample_time_map[i]
            assert all(i1 is i2 for i1, i2 in zip(c._inputs, inputs_map[i]))
            assert all(i1 is i2 for i1, i2 in zip(c._measurements, measurements_map[i]))
            assert hasattr(c, "category_dict")
            assert hasattr(c, "con_category_dict")
            assert hasattr(c, "MHE_VARS_CONS_BLOCK")
            MHEBlock = c.MHE_VARS_CONS_BLOCK
            assert hasattr(MHEBlock, "MEASUREMENT_SET")
            assert len(MHEBlock.MEASUREMENT_SET.ordered_data()) == len(measurements_map[i])
            assert hasattr(MHEBlock, "DIFFERENTIAL_SET")
            assert len(MHEBlock.DIFFERENTIAL_SET.ordered_data()) == \
                len(c.category_dict[VariableCategory.DIFFERENTIAL])

    @pytest.mark.unit
    def make_estimator(self, sample_time=0.5, horizon=1., nfe=2):
        model = make_model(horizon=horizon, nfe=nfe)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[0,'A']]
        estimator = EstimatorBlock(
                model=model,
                time=time,
                inputs=inputs,
                measurements=measurements,
                sample_time=sample_time,
                )
        estimator.construct()
        estimator.set_sample_time(sample_time)
        return estimator

    @pytest.mark.unit
    def test_use_user_given_var_con_categ_dict(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        sample_time = 0.5
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[t0,'A']]

        scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=pyo.Var)
        scalar_cons, dae_cons = flatten_dae_components(model, time, ctype = pyo.Constraint)
        input_ids = [id(vart0) for vart0 in inputs]
        input_vars = [var for var in dae_vars if id(var[t0]) in input_ids]
        measurement_ids = [id(vart0) for vart0 in measurements]
        measurement_vars = [var for var in dae_vars if id(var[t0]) in measurement_ids]
        category_dict, con_category_dict = categorize_dae_variables_and_constraints(
                                                                                    model,
                                                                                    dae_vars,
                                                                                    dae_cons,
                                                                                    time,
                                                                                    input_vars = input_vars,
                                                                                    )

        # Provide the category_dict when constructing EstimatorBlock
        estimator = EstimatorBlock(
                model=model,
                time=time,
                category_dict={None: category_dict},
                con_category_dict={None: con_category_dict},
                inputs=input_vars,
                measurements=measurement_vars,
                sample_time=sample_time,
                )
        estimator.construct()
        
        # Make sure _category_dict is not None
        assert estimator._category_dict
        
        # Make sure dae variables are saved in estimator.dae_vars
        vart0_id_dae_vars = [id(var[t0]) for var in estimator.dae_vars]
        for categ, varlist in estimator.category_dict.items():
            if categ not in [VariableCategory.MEASUREMENT,
                             VariableCategory.ACTUALMEASUREMENT,
                             VariableCategory.MEASUREMENTERROR,
                             VariableCategory.MODELDISTURBANCE]:
                for var in varlist:
                    assert id(var[t0]) in vart0_id_dae_vars
                    vart0_id_dae_vars.remove(id(var[t0]))
        assert vart0_id_dae_vars == []
        
        # Make sure con_category_dict is not None
        assert estimator._con_category_dict
        assert estimator.con_category_dict == estimator._con_category_dict
        
    @pytest.mark.unit
    def test_categorize_var_con_for_MHE(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        mod = estimator.mod
        
        # line253 - 287: make sure variables are categoried in correct cateogries.
        pred_diff_vars = [mod.conc[t0, "A"], mod.conc[t0, "B"]]
        diff_vars = estimator.category_dict[VariableCategory.DIFFERENTIAL]
        diff_vart0_ids = [id(var[t0]) for var in diff_vars]
        assert len(diff_vars) == len(pred_diff_vars)
        for var in pred_diff_vars:
            assert id(var) in diff_vart0_ids
            diff_vart0_ids.remove(id(var))
        assert diff_vart0_ids == []
        
        pred_alg_vars = [mod.flow_out[t0], mod.rate[t0, "A"], mod.rate[t0, "B"]]
        alg_vars = estimator.category_dict[VariableCategory.ALGEBRAIC]
        alg_vart0_ids = [id(var[t0]) for var in alg_vars]
        assert len(alg_vars) == len(pred_alg_vars)
        for var in pred_alg_vars:
            assert id(var) in alg_vart0_ids
            alg_vart0_ids.remove(id(var))
        assert alg_vart0_ids == []
        
        pred_deri_vars = [mod.dcdt[t0, "A"], mod.dcdt[t0, "B"]]
        deri_vars = estimator.category_dict[VariableCategory.DERIVATIVE]
        deri_vart0_ids = [id(var[t0]) for var in deri_vars]
        assert len(deri_vars) == len(pred_deri_vars)
        for var in pred_deri_vars:
            assert id(var) in deri_vart0_ids
            deri_vart0_ids.remove(id(var))
        assert deri_vart0_ids == []
        
        pred_input_vars = [mod.flow_in[t0]]
        input_vars = estimator.category_dict[VariableCategory.INPUT]
        input_vart0_ids = [id(var[t0]) for var in input_vars]
        assert len(input_vars) == len(pred_input_vars)
        for var in pred_input_vars:
            assert id(var) in input_vart0_ids
            input_vart0_ids.remove(id(var))
        assert input_vart0_ids == []
        
        
        # line293 - 322: make sure constraints are categoried in correct cateogries.
        tlast = time.last()
        
        pred_diff_cons = [mod.material_balance[tlast, "A"],
                          mod.material_balance[tlast, "B"]]
        diff_cons = estimator.con_category_dict[ConstraintCategory.DIFFERENTIAL]
        diff_contlast_ids = [id(con[tlast]) for con in diff_cons]
        assert len(diff_cons) == len(pred_diff_cons)
        for con in pred_diff_cons:
            assert id(con) in diff_contlast_ids
            diff_contlast_ids.remove(id(con))
        assert diff_contlast_ids == []
        
        pred_alg_cons = [mod.flow_eqn[tlast],
                         mod.rate_eqn[tlast, "A"],
                         mod.rate_eqn[tlast, "B"]]
        alg_cons = estimator.con_category_dict[ConstraintCategory.ALGEBRAIC]
        alg_contlast_ids = [id(con[tlast]) for con in alg_cons]
        assert len(alg_cons) == len(pred_alg_cons)
        for con in pred_alg_cons:
            assert id(con) in alg_contlast_ids
            alg_contlast_ids.remove(id(con))
        assert alg_contlast_ids == []
        
        pred_disc_cons = [mod.dcdt_disc_eq[tlast, "A"], 
                          mod.dcdt_disc_eq[tlast, "B"]]
        disc_cons = estimator.con_category_dict[ConstraintCategory.DISCRETIZATION]
        disc_contlast_ids = [id(con[tlast]) for con in disc_cons]
        assert len(disc_cons) == len(pred_disc_cons)
        for con in pred_disc_cons:
            assert id(con) in disc_contlast_ids
            disc_contlast_ids.remove(id(con))
        assert disc_contlast_ids == []
        
    @pytest.mark.unit
    def test_reference_var_category_dict(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        assert VariableCategory.DIFFERENTIAL in estimator.category_dict
        assert VariableCategory.ALGEBRAIC in estimator.category_dict
        assert VariableCategory.DERIVATIVE in estimator.category_dict
        assert VariableCategory.INPUT in estimator.category_dict
        
        CATEGORY_MAP = {VariableCategory.DIFFERENTIAL: DiffVar,
                        VariableCategory.ALGEBRAIC: AlgVar,
                        VariableCategory.DERIVATIVE: DerivVar,
                        VariableCategory.INPUT: InputVar}
        
        # Make sure variables are references and their ctype are correct.
        for categ, ctype in CATEGORY_MAP.items():
            varlist = estimator.category_dict[categ]
            for var in varlist:
                assert var.is_reference()
                assert var.ctype == ctype
        
    @pytest.mark.unit
    def test_add_sample_point_set(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        assert hasattr(estimator, 'SAMPLEPOINT_SET')
        assert all(i1 is i2 for i1, i2 
                   in zip(estimator.SAMPLEPOINT_SET.ordered_data(), estimator.sample_points))
        
    @pytest.mark.unit
    def test_add_actual_measurement_param(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        MHEBlock = estimator.MHE_VARS_CONS_BLOCK
        assert hasattr(MHEBlock, "ACTUAL_MEASUREMENT_BLOCK")
        acemeablock = MHEBlock.ACTUAL_MEASUREMENT_BLOCK
        assert acemeablock.dim() == 1
        assert acemeablock.index_set() == estimator.MEASUREMENT_SET
        
        for bind in estimator.MEASUREMENT_SET:
            curr_block = acemeablock[bind]
            assert hasattr(curr_block, "actual_measurement")
            assert curr_block.actual_measurement.index_set() == estimator.SAMPLEPOINT_SET
            var = curr_block.actual_measurement
            for sp in estimator.SAMPLEPOINT_SET:
                assert var[sp].value == 0.0
                assert var[sp].fixed
                
    @pytest.mark.unit
    def test_add_measurement_error(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        MHEBlock = estimator.MHE_VARS_CONS_BLOCK
        assert hasattr(MHEBlock, "MEASUREMENT_ERROR_BLOCK")
        meaerrblock = MHEBlock.MEASUREMENT_ERROR_BLOCK
        assert meaerrblock.dim() == 1
        assert meaerrblock.index_set() == estimator.MEASUREMENT_SET
        
        for bind in estimator.MEASUREMENT_SET:
            curr_block = meaerrblock[bind]
            assert hasattr(curr_block, "measurement_error")
            assert curr_block.measurement_error.index_set() == estimator.SAMPLEPOINT_SET
            var = curr_block.measurement_error
            for sp in estimator.SAMPLEPOINT_SET:
                assert var[sp].value == 0.0
                     
    @pytest.mark.unit
    def test_add_model_disturbance(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        MHEBlock = estimator.MHE_VARS_CONS_BLOCK
        assert hasattr(MHEBlock, "MODEL_DISTURBANCE_BLOCK")
        moddisblock = MHEBlock.MODEL_DISTURBANCE_BLOCK
        assert moddisblock.dim() == 1
        assert moddisblock.index_set() == estimator.DIFFERENTIAL_SET
        
        for bind in estimator.MEASUREMENT_SET:
            curr_block = moddisblock[bind]
            assert hasattr(curr_block, "model_disturbance")
            assert curr_block.model_disturbance.index_set() == estimator.SAMPLEPOINT_SET
            var = curr_block.model_disturbance
            for sp in estimator.SAMPLEPOINT_SET:
                if sp == t0:
                    assert var[sp].fixed
                assert var[sp].value == 0.0
                
    @pytest.mark.unit
    def test_add_MHE_vars_to_category_dict(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        # Make sure category dict contains these categories
        assert VariableCategory.DIFFERENTIAL in estimator.category_dict
        assert VariableCategory.ALGEBRAIC in estimator.category_dict
        assert VariableCategory.DERIVATIVE in estimator.category_dict
        assert VariableCategory.INPUT in estimator.category_dict

        assert VariableCategory.ACTUALMEASUREMENT in estimator.category_dict
        assert VariableCategory.MEASUREMENTERROR in estimator.category_dict
        assert VariableCategory.MODELDISTURBANCE in estimator.category_dict
        
        
        # Make sure MHE categories contains variables we expect
        MHEBlock = estimator.MHE_VARS_CONS_BLOCK
        
        actmeablock = MHEBlock.ACTUAL_MEASUREMENT_BLOCK
        pred_actmea_vars = ComponentSet((actmeablock[ind].actual_measurement[t0]
                                         for ind in estimator.MEASUREMENT_SET))
        actmea_vars = estimator.category_dict[VariableCategory.ACTUALMEASUREMENT]
        for var in actmea_vars:
            assert var.is_reference()
            assert var.ctype == ActualMeasurementVar
            assert var._attr == "actualmeasurement"
            assert var[t0] in pred_actmea_vars
            pred_actmea_vars.remove(var[t0])
        assert pred_actmea_vars._data == {}
        
        meaerrblock = MHEBlock.MEASUREMENT_ERROR_BLOCK
        pred_meaerr_vars = ComponentSet((meaerrblock[ind].measurement_error[t0]
                                         for ind in estimator.MEASUREMENT_SET))
        meaerr_vars = estimator.category_dict[VariableCategory.MEASUREMENTERROR]
        for var in meaerr_vars:
            assert var.is_reference()
            assert var.ctype == MeasurementErrorVar
            assert var._attr == "measurementerror"
            assert var[t0] in pred_meaerr_vars
            pred_meaerr_vars.remove(var[t0])
        assert pred_meaerr_vars._data == {}
        
        moddisblock = MHEBlock.MODEL_DISTURBANCE_BLOCK
        pred_moddis_vars = ComponentSet((moddisblock[ind].model_disturbance[t0]
                                         for ind in estimator.DIFFERENTIAL_SET))
        moddis_vars = estimator.category_dict[VariableCategory.MODELDISTURBANCE]
        for var in moddis_vars:
            assert var.is_reference()
            assert var.ctype == ModelDisturbanceVar
            assert var._attr == "modeldisturbance"
            assert var[t0] in pred_moddis_vars
            pred_moddis_vars.remove(var[t0])
        assert pred_moddis_vars._data == {}
        
    @pytest.mark.unit
    def test_add_mea_moddis_componentmap(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        moddis_block = estimator.MODELDISTURBANCE_BLOCK
        diffvar_map_moddis = estimator.diffvar_map_moddis
        for ind, var in enumerate(estimator.differential_vars):
            assert diffvar_map_moddis[var[t0]] == moddis_block[ind].var   
        derivar_map_moddis = estimator.derivar_map_moddis
        for ind, var in enumerate(estimator.derivative_vars):
            assert derivar_map_moddis[var[t0]] == moddis_block[ind].var
        
            
        meaerr_block = estimator.MEASUREMENTERROR_BLOCK
        meavar_map_meaerr = estimator.meavar_map_meaerr
        for ind, var in enumerate(estimator.measurement_vars):
            assert meavar_map_meaerr[var[t0]] == meaerr_block[ind].var
            
        actmea_block = estimator.ACTUALMEASUREMENT_BLOCK
        meavar_map_actmea = estimator.meavar_map_actmea
        for ind, var in enumerate(estimator.measurement_vars):
            assert meavar_map_actmea[var[t0]] == actmea_block[ind].var
        
    @pytest.mark.unit
    def test_add_measurement_constraint(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        MHEBlock = estimator.MHE_VARS_CONS_BLOCK
        
        assert hasattr(MHEBlock, "MEASUREMENT_CONSTRAINT_BLOCK") 
        meacon_block = MHEBlock.MEASUREMENT_CONSTRAINT_BLOCK
        assert meacon_block.dim() == 1
        assert meacon_block.index_set() == estimator.MEASUREMENT_SET
        
        actmea_block = estimator.MHE_VARS_CONS_BLOCK.ACTUAL_MEASUREMENT_BLOCK
        meaerr_block = estimator.MHE_VARS_CONS_BLOCK.MEASUREMENT_ERROR_BLOCK
        mea_block = estimator.MEASUREMENT_BLOCK
        for bind in estimator.MEASUREMENT_SET:
            assert len(list(meacon_block[bind].component_objects(pyo.Constraint))) == 1
            assert hasattr(meacon_block[bind], "con_mea_err")
            curr_con = meacon_block[bind].con_mea_err
            assert curr_con.index_set() == time
            assert all(i1 is i2 for i1, i2 in zip(curr_con.keys(), estimator.SAMPLEPOINT_SET))
            for t in time: 
                if t in estimator.SAMPLEPOINT_SET:
                    pred_expr = actmea_block[bind].actual_measurement[t] ==\
                                (mea_block[bind].var[t] + meaerr_block[bind].measurement_error[t])
                    assert curr_con[t].expr.to_string() == pred_expr.to_string()
            
    @pytest.mark.unit
    def test_add_disturbance_to_differential_cons(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        MHEBlock = estimator.MHE_VARS_CONS_BLOCK
        
        assert hasattr(MHEBlock, "DISTURBED_DIFFERENTIAL_CONSTRAINT_BLOCK") 
        ddc_block = MHEBlock.DISTURBED_DIFFERENTIAL_CONSTRAINT_BLOCK
        assert ddc_block.dim() == 1
        assert len(list(ddc_block.keys())) == len(estimator.con_category_dict[ConstraintCategory.DIFFERENTIAL])
        
        def map_t_to_sp(tind):
            for item in estimator.sample_point_indices:
                if item >= tind:
                    return estimator.time[item]
        
        diff_cons = estimator.con_category_dict[ConstraintCategory.DIFFERENTIAL]
        moddis_block = estimator.MODELDISTURBANCE_BLOCK
        for bind in estimator.DIFFERENTIAL_SET:
            curr_block = ddc_block[bind]
            assert len(list(curr_block.component_objects(pyo.Constraint))) == 1
            assert hasattr(curr_block, "disturbed_diff_con")
            curr_con = curr_block.disturbed_diff_con
            assert curr_con.index_set() == time
            for t in time:
                tind = time.find_nearest_index(t)
                sp = map_t_to_sp(tind)
                pred_expr = diff_cons[bind][t].body - moddis_block[bind].var[sp] == 0.0
                assert curr_con[t].expr.to_string() == pred_expr.to_string()
                assert not diff_cons[bind][t].active
        
                
    @pytest.mark.unit
    def test_add_steady_state_objective_for_MHE_initialization(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        desiredss = [
                (estimator.mod.flow_in[t0], 3.0),
                ]
        weights = [
                (estimator.mod.flow_in[t0], 2.0),
                ]
        estimator.mod.flow_in[:].set_value(3.0)
        initialize_t0(estimator.mod)
        
        estimator.add_single_time_optimization_objective(desiredss, weights)

        assert hasattr(estimator, 'single_time_optimization_objective')

        pred_obj_expr = (2.0*(estimator.mod.flow_in[t0] - 3.0)**2)
        obj_expr = estimator.single_time_optimization_objective.expr
        assert pyo.value(pred_obj_expr) == pyo.value(obj_expr)
        assert pred_obj_expr.to_string() == obj_expr.to_string()

        estimator.del_component(estimator.single_time_optimization_objective)

        desiredss = [
                (estimator.mod.flow_in[t0], 3.0),
                (estimator.mod.conc[t0,'A'], 1.5),
                ]
        weights = [
                (estimator.mod.flow_in[t0], 1.0),
                (estimator.mod.conc[t0,'A'], 5.0),
                ]
        estimator.add_single_time_optimization_objective(desiredss, weights)
        pred_obj_expr = (
                1.0*(estimator.mod.flow_in[t0] - 3.0)**2 + 
                5.0*(estimator.mod.conc[t0,'A'] - 1.5)**2
                )
        obj_expr = estimator.single_time_optimization_objective.expr
        assert pyo.value(pred_obj_expr) == pyo.value(obj_expr)
        assert pred_obj_expr.to_string() == obj_expr.to_string()
        
    @pytest.mark.component
    @pytest.mark.skipif(not solver_available, reason='IPOPT is not available')
    def test_solve_steady_state_for_MHE_initialization(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        desiredss = [
                (estimator.mod.flow_in[t0], 3.0),
                ]
        weights = [
                (estimator.mod.flow_in[t0], 1.0),
                ]
        estimator.add_single_time_optimization_objective(desiredss, weights)
        estimator.mod.flow_in[:].set_value(3.0)
        initialize_t0(estimator.mod)

        dof_prior = degrees_of_freedom(estimator)
        estimator.solve_single_time_optimization(solver, 
                                            ic_type = "differential_var",
                                            require_steady = True,
                                            load_setpoints = False,
                                            restore_ic_input_after_solve = False,
                                            isMHE_block = True,)
        dof_post = degrees_of_freedom(estimator)

        assert dof_prior == dof_post

        assert estimator.differential_vars[0][t0].value == \
                pytest.approx(3.75, abs=1e-3)
        assert estimator.differential_vars[1][t0].value == \
                pytest.approx(1.25, abs=1e-3)
        assert estimator.algebraic_vars[0][t0].value == \
                pytest.approx(3.0, abs=1e-3)
        assert estimator.algebraic_vars[1][t0].value == \
                pytest.approx(-3.75, abs=1e-3)
        assert estimator.algebraic_vars[2][t0].value == \
                pytest.approx(3.75, abs=1e-3)
        assert estimator.input_vars[0][t0].value == \
                pytest.approx(3.0, abs=1e-3)
                
        # Make sure all undisturbed differential equations are deactivated
        diff_cons = estimator.con_category_dict[ConstraintCategory.DIFFERENTIAL]
        assert not all(con[t].active for con in diff_cons for t in time)
        # Make sure MHE block is active
        assert estimator.MHE_VARS_CONS_BLOCK.active
        # Make sure derivatives are not fixed
        deri_vars = estimator.category_dict[VariableCategory.DERIVATIVE]
        assert not all(var[t].fixed for var in deri_vars for t in time)
        # Make sure all inputs are fixed
        input_vars = estimator.category_dict[VariableCategory.INPUT]
        assert all(var[t].fixed for var in input_vars for t in time)
        # Make sure differential variables at t0 are fixed
        diff_vars = estimator.category_dict[VariableCategory.DIFFERENTIAL]
        assert all(var[t0].fixed for var in diff_vars)
                
    @pytest.mark.unit
    def test_initialize_actualmeasurements_at_t0(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        # Re-initialize measurements and actualmeasurements to make sure they are different before assertion
        estimator.vectors.measurement[...].set_value(20.)
        estimator.vectors.actualmeasurement[...].set_value(100.)
        
        estimator.initialize_actualmeasurements_at_t0()
        
        for ind in estimator.ACTUALMEASUREMENT_SET:
            var_ind = (ind, t0)
            actmea_var = estimator.vectors.actualmeasurement
            assert actmea_var[var_ind].value == 20.
            
    @pytest.mark.component
    @pytest.mark.skipif(not solver_available, reason='IPOPT is not available')
    def test_initialize_past_info_with_steady_state(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        desiredss = [
                (estimator.mod.flow_in[t0], 3.0),
                ]
        weights = [
                (estimator.mod.flow_in[t0], 1.0),
                ]
        
        estimator.mod.flow_in[:].set_value(3.0)
        initialize_t0(estimator.mod)
        
        estimator.initialize_past_info_with_steady_state(desiredss,
                                                         weights, 
                                                         solver)
        
        for tp in time:
            assert estimator.differential_vars[0][tp].value == \
                    pytest.approx(3.75, abs=1e-3)
            assert estimator.differential_vars[1][tp].value == \
                    pytest.approx(1.25, abs=1e-3)
            assert estimator.algebraic_vars[0][tp].value == \
                    pytest.approx(3.0, abs=1e-3)
            assert estimator.algebraic_vars[1][tp].value == \
                    pytest.approx(-3.75, abs=1e-3)
            assert estimator.algebraic_vars[2][tp].value == \
                    pytest.approx(3.75, abs=1e-3)
            assert estimator.input_vars[0][tp].value == \
                    pytest.approx(3.0, abs=1e-3)
                    
        for sp in estimator.sample_points:
            curr_var = estimator.vectors.actualmeasurement[0, sp]
            assert curr_var.value == 3.75
            assert curr_var.fixed
                
    @pytest.mark.unit
    def test_add_noise_minimize_objective(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        # Re-initialize noises so they are not zero for objective evalutaion
        estimator.vectors.modeldisturbance[...].set_value(0.5)
        estimator.vectors.measurementerror[...].set_value(2.5)
        
        
        model_disturbance_weights = [
                                    (estimator.mod.conc[t0,'A'], 0.1),
                                    (estimator.mod.conc[t0,'B'], 0.2),
                                    ]
        measurement_error_weights = [(estimator.mod.conc[0,'A'], 10.),]
        estimator.add_noise_minimize_objective(model_disturbance_weights,
                                               measurement_error_weights,
                                               givenform = "weight")
        
        assert hasattr(estimator, 'noise_minimize_objective')
        
        MD_block = estimator.MHE_VARS_CONS_BLOCK.MODEL_DISTURBANCE_BLOCK
        ME_block = estimator.MHE_VARS_CONS_BLOCK.MEASUREMENT_ERROR_BLOCK
        pred_obj_expr = sum(0.1*(MD_block[0].model_disturbance[t])**2 for t in estimator.sample_points if t!=t0) + \
                        sum(0.2*(MD_block[1].model_disturbance[t])**2 for t in estimator.sample_points if t!=t0) + \
                        sum(10.*(ME_block[0].measurement_error[t])**2 for t in estimator.sample_points)
        obj_expr = estimator.noise_minimize_objective.expr
        assert pyo.value(pred_obj_expr) == pyo.value(obj_expr)
        assert pyo.value(obj_expr) > 0
        estimator.del_component("noise_minimize_objective")
        
        
        model_disturbance_variances = [
                                        (estimator.mod.conc[t0,'A'], 10.),
                                        (estimator.mod.conc[t0,'B'], 20.),
                                        ]
        measurement_error_variances = [(estimator.mod.conc[0,'A'], 0.1),]
        estimator.add_noise_minimize_objective(model_disturbance_variances,
                                               measurement_error_variances,
                                               givenform = "variance")
        
        assert hasattr(estimator, 'noise_minimize_objective')
        
        pred_obj_expr = sum(1/10.*(MD_block[0].model_disturbance[t])**2 for t in estimator.sample_points if t!=t0) + \
                        sum(1/20.*(MD_block[1].model_disturbance[t])**2 for t in estimator.sample_points if t!=t0) + \
                        sum(1/0.1*(ME_block[0].measurement_error[t])**2 for t in estimator.sample_points)
        obj_expr = estimator.noise_minimize_objective.expr
        assert pyo.value(pred_obj_expr) == pyo.value(obj_expr)
        assert pyo.value(obj_expr) > 0
        
    @pytest.mark.unit
    def test_check_var_con_dof(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        estimator.check_var_con_dof(skip_dof_check = False)
        
        input_vars = estimator.vectors.input
        for ind in input_vars.keys():
            assert input_vars[ind].fixed
            
        diffvars = estimator.vectors.differential
        for ind in diffvars.keys():
            assert not diffvars[ind].fixed
            
        for ind in estimator.MODELDISTURBANCE_SET:
            assert estimator.vectors.modeldisturbance[ind, 0].fixed
            assert pyo.value(estimator.vectors.modeldisturbance[ind, 0]) == 0.0
            
    @pytest.mark.unit
    def test_load_inputs_for_MHE(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        
        # Re-initialize inputs so they are not 100.
        estimator.vectors.input[...].set_value(0.5)
        
        input = [100.]
        estimator.load_inputs_into_last_sample(input)
        
        spi = estimator.sample_point_indices
        check_t_list = [tp for tp in time if tp > time[spi[-2]] and tp <= time[spi[-1]]]
        for tp in check_t_list:
            assert pyo.value(estimator.vectors.input[0, tp]) == 100.
    
    @pytest.mark.unit
    def test_generate_estimates_at_time(self):
        estimator = self.make_estimator()
        time = estimator.time
        t0 = time.first()
        tlast = time.last()
        
        # Re-set values for differential variables
        estimator.vectors.differential[0,tlast].set_value(105.)
        estimator.vectors.differential[1,tlast].set_value(205.)

        estimates = estimator.generate_estimates_at_time(tlast)
        assert estimates == [105., 205.]

    @pytest.mark.unit
    def test_load_measurements(self):
        blk = self.make_estimator()
        time = blk.time
        t0 = time.first()
        vals = [0.25]
        blk.load_measurements(vals, target = "measurement", timepoint = t0)        
        for b, val in zip(blk.MEASUREMENT_BLOCK.values(), vals):
            assert b.var[t0].value == val
            
        vals2 = [0.75]
        t_last = time.last()
        blk.load_measurements(vals2, target = "actualmeasurement", timepoint = t_last)        
        for b, val in zip(blk.ACTUALMEASUREMENT_BLOCK.values(), vals2):
            assert b.var[t_last].value == val

abc = TestEstimatorBlock()
abc.test_solve_steady_state()