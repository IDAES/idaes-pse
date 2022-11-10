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
""" Tests for the dynamic model subclass of block
"""

import pyomo.environ as pyo
import pyomo.dae as dae
import pyomo.network as pyn
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.core.base.block import _BlockData, SubclassOf
from pyomo.dae.flatten import flatten_dae_components

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    activated_equalities_generator,
    unfixed_variables_generator,
)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.apps.caprese.tests.test_simple_model import (
    make_model,
    make_small_model,
    initialize_t0,
    copy_values_forward,
)
from idaes.apps.caprese.dynamic_block import (
    DynamicBlock,
    SimpleDynamicBlock,
    IndexedDynamicBlock,
    _DynamicBlockData,
)
from idaes.apps.caprese.categorize import categorize_dae_variables
from idaes.apps.caprese.common.config import (
    VariableCategory,
    InputOption,
)

VC = VariableCategory
from idaes.apps.caprese.nmpc_var import (
    NmpcVar,
    DiffVar,
    DerivVar,
    AlgVar,
    InputVar,
    FixedVar,
    MeasuredVar,
)
import idaes.logger as idaeslog
import random
import pytest

__author__ = "Robert Parker"

solver_available = pyo.SolverFactory("ipopt").available()
if solver_available:
    solver = pyo.SolverFactory("ipopt")
else:
    solver = None


class TestDynamicBlock(object):
    @pytest.mark.unit
    def test_init_simple(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[0, "A"], model.conc[0, "B"]]
        block = DynamicBlock(
            model=model,
            time=time,
            inputs=inputs,
            measurements=measurements,
        )
        # Assert that we have the correct type
        assert type(block) is SimpleDynamicBlock
        assert isinstance(block, DynamicBlock)
        assert isinstance(block, _DynamicBlockData)

        block.construct()
        # Assert that we behave like a simple block
        assert block[None] is block
        assert all(b is block for b in block[:])

        # Assert that input attributes have been processed correctly
        assert block.mod is model
        assert block.time is time
        assert all(i1 is i2 for i1, i2 in zip(block._inputs, inputs))
        assert all(i1 is i2 for i1, i2 in zip(block._measurements, measurements))

        # Assert that utility attributes have been added
        assert hasattr(block, "category_dict")
        assert hasattr(block, "vardata_map")
        assert hasattr(block, "measurement_vars")
        assert hasattr(block, "differential_vars")
        assert hasattr(block, "algebraic_vars")
        assert hasattr(block, "derivative_vars")
        assert hasattr(block, "input_vars")
        assert hasattr(block, "fixed_vars")

        subblocks = [
            block.mod,
            block.vectors,
            block.DIFFERENTIAL_BLOCK,
            block.ALGEBRAIC_BLOCK,
            block.INPUT_BLOCK,
            block.FIXED_BLOCK,
            block.DERIVATIVE_BLOCK,
            block.MEASUREMENT_BLOCK,
        ]

        block_objects = ComponentSet(
            block.component_objects(pyo.Block, descend_into=False)
        )
        # Assert that subblocks have been added
        assert len(subblocks) == len(block_objects)
        for b in subblocks:
            assert b in block_objects

        # Assert that we can add variables and constraints to the block
        block.v = pyo.Var(initialize=3)
        block.c = pyo.Constraint(expr=block.v == 5)
        assert block.v.value == 3
        assert block.v in ComponentSet(identify_variables(block.c.expr))

    @pytest.mark.unit
    def test_init_indexed(self):
        block_set = pyo.Set(initialize=[0, 1, 2])
        block_set.construct()
        horizon_map = {0: 1.0, 1: 3.0, 2: 5.0}
        nfe_map = {0: 2, 1: 6, 2: 10}
        model_map = {i: make_model(horizon_map[i], nfe_map[i]) for i in block_set}
        time_map = {i: model_map[i].time for i in block_set}
        inputs_map = {i: [model_map[i].flow_in[0]] for i in block_set}
        measurements_map = {
            i: [model_map[i].conc[0, "A"], model_map[i].conc[0, "B"]] for i in block_set
        }
        # Construct block with a dict for each of its arguments
        block = DynamicBlock(
            block_set,
            model=model_map,
            time=time_map,
            inputs=inputs_map,
            measurements=measurements_map,
        )
        # Make sure we have the right type
        assert type(block) is IndexedDynamicBlock
        assert isinstance(block, DynamicBlock)

        block.construct()
        assert all(b.parent_component() is block for b in block.values())

        # Check __contains__
        for i in block_set:
            assert i in block

        # Check attributes and subblocks of each data object
        for i, b in block.items():
            assert b.mod is model_map[i]
            assert b.time is time_map[i]
            assert all(i1 is i2 for i1, i2 in zip(b._inputs, inputs_map[i]))
            assert all(i1 is i2 for i1, i2 in zip(b._measurements, measurements_map[i]))

            assert hasattr(b, "category_dict")
            assert hasattr(b, "vardata_map")
            assert hasattr(b, "measurement_vars")
            assert hasattr(b, "differential_vars")
            assert hasattr(b, "algebraic_vars")
            assert hasattr(b, "derivative_vars")
            assert hasattr(b, "input_vars")
            assert hasattr(b, "fixed_vars")

            subblocks = [
                b.mod,
                b.vectors,
                b.DIFFERENTIAL_BLOCK,
                b.ALGEBRAIC_BLOCK,
                b.INPUT_BLOCK,
                b.FIXED_BLOCK,
                b.DERIVATIVE_BLOCK,
                b.MEASUREMENT_BLOCK,
            ]

            block_objects = ComponentSet(
                b.component_objects(pyo.Block, descend_into=False)
            )
            assert len(subblocks) == len(block_objects)
            for sb in subblocks:
                assert sb in block_objects

            b.v = pyo.Var(initialize=3)
            b.c = pyo.Constraint(expr=b.v == 5)
            assert b.v.value == 3
            assert b.v in ComponentSet(identify_variables(b.c.expr))

    @pytest.mark.unit
    def test_init_rule(self):
        block_set = pyo.Set(initialize=range(3))
        block_set.construct()
        # Create same maps as before
        horizon_map = {0: 1, 1: 3, 2: 5}
        nfe_map = {0: 2, 1: 6, 2: 10}
        model_map = {
            i: make_model(horizon=horizon_map[i], nfe=nfe_map[i]) for i in block_set
        }

        # Create rule to construct DynamicBlock with
        def dynamic_block_rule(b, i):
            model = model_map[i]
            time = model.time
            t0 = time.first()
            inputs = [model.flow_in[t0]]
            measurements = [model.conc[0, "A"], model.conc[0, "B"]]

            # Won't be obvious that these attrs need to be set if
            # constructing from a rule
            b.mod = model
            super(_BlockData, b).__setattr__("time", time)
            b._inputs = inputs
            b._measurements = measurements

        # Create DynamicBlock from a rule
        block = DynamicBlock(block_set, rule=dynamic_block_rule)
        assert type(block) is IndexedDynamicBlock
        assert isinstance(block, DynamicBlock)

        block.construct()

        # Make sure iterating over block.values works as expected
        assert all(b.parent_component() is block for b in block.values())

        # Make sure __contains__ works
        for i in block_set:
            assert i in block

        # Assert correct attributes and subblocks
        for i, b in block.items():
            assert b.mod is model_map[i]
            assert b.time is model_map[i].time
            t0 = b.time.first()
            assert all(
                i1 is i2 for i1, i2 in zip(b._inputs, [model_map[i].flow_in[t0]])
            )
            assert all(
                i1 is i2
                for i1, i2 in zip(
                    b._measurements,
                    [model_map[i].conc[t0, "A"], model_map[i].conc[t0, "B"]],
                )
            )

            assert hasattr(b, "category_dict")
            assert hasattr(b, "vardata_map")
            assert hasattr(b, "measurement_vars")
            assert hasattr(b, "differential_vars")
            assert hasattr(b, "algebraic_vars")
            assert hasattr(b, "derivative_vars")
            assert hasattr(b, "input_vars")
            assert hasattr(b, "fixed_vars")

            subblocks = [
                b.mod,
                b.vectors,
                b.DIFFERENTIAL_BLOCK,
                b.ALGEBRAIC_BLOCK,
                b.INPUT_BLOCK,
                b.FIXED_BLOCK,
                b.DERIVATIVE_BLOCK,
                b.MEASUREMENT_BLOCK,
            ]

            block_objects = ComponentSet(
                b.component_objects(pyo.Block, descend_into=False)
            )
            assert len(subblocks) == len(block_objects)
            for sb in subblocks:
                assert sb in block_objects

            b.v = pyo.Var(initialize=3)
            b.c = pyo.Constraint(expr=b.v == 5)
            assert b.v.value == 3
            assert b.v in ComponentSet(identify_variables(b.c.expr))

    @pytest.mark.unit
    def test_construct(self):
        m = make_small_model()
        time = m.time
        t0 = time.first()
        helper = DynamicBlock(
            model=m,
            time=time,
            inputs=[m.flow_in[0]],
            measurements=[m.conc[0, "A"], m.conc[0, "B"]],
        )
        helper.construct()
        assert hasattr(helper, "category_dict")
        assert hasattr(helper, "vardata_map")
        assert hasattr(helper, "measurement_vars")
        assert hasattr(helper, "differential_vars")
        assert hasattr(helper, "algebraic_vars")
        assert hasattr(helper, "derivative_vars")
        assert hasattr(helper, "input_vars")
        assert hasattr(helper, "fixed_vars")

        # Make sure category dict contains variables we expect
        assert VariableCategory.DIFFERENTIAL in helper.category_dict
        assert VariableCategory.ALGEBRAIC in helper.category_dict
        assert VariableCategory.DERIVATIVE in helper.category_dict
        assert VariableCategory.INPUT in helper.category_dict
        assert VariableCategory.FIXED in helper.category_dict

        pred_diff_vars = ComponentSet(
            (
                m.conc[t0, "A"],
                m.conc[t0, "B"],
            )
        )
        diff_vars = list(helper.component_objects(DiffVar))
        assert len(pred_diff_vars) == len(diff_vars)
        for var in diff_vars:
            assert var[t0] in pred_diff_vars

        pred_alg_vars = ComponentSet(
            (
                m.flow_out[t0],
                m.rate[t0, "A"],
                m.rate[t0, "B"],
            )
        )
        alg_vars = list(helper.component_objects(AlgVar))
        assert len(pred_alg_vars) == len(alg_vars)
        for var in alg_vars:
            assert var[t0] in pred_alg_vars

        pred_input_vars = ComponentSet((m.flow_in[t0],))
        input_vars = list(helper.component_objects(InputVar))
        assert len(pred_input_vars) == len(input_vars)
        for var in input_vars:
            assert var[t0] in pred_input_vars

        pred_fixed_vars = ComponentSet(
            (
                m.conc_in[t0, "A"],
                m.conc_in[t0, "B"],
            )
        )
        fixed_vars = list(helper.component_objects(FixedVar))
        assert len(pred_fixed_vars) == len(fixed_vars)
        for var in fixed_vars:
            assert var[t0] in pred_fixed_vars

        pred_deriv_vars = ComponentSet(
            (
                m.dcdt[t0, "A"],
                m.dcdt[t0, "B"],
            )
        )
        deriv_vars = list(helper.component_objects(DerivVar))
        assert len(pred_deriv_vars) == len(deriv_vars)
        for var in deriv_vars:
            assert var[t0] in pred_deriv_vars

    @pytest.mark.unit
    def test_category_blocks(self):
        m = make_small_model()
        time = m.time
        t0 = time.first()
        helper = DynamicBlock(
            model=m,
            time=time,
            inputs=[m.flow_in[0]],
            measurements=[m.conc[0, "A"], m.conc[0, "B"]],
        )
        helper.construct()

        # Test that NmpcVectors and category blocks behave
        # as we expect
        diff_vars = helper.vectors.differential
        diff_var_set = helper.DIFFERENTIAL_SET * m.time
        assert diff_vars.index_set() == diff_var_set
        # And that they contain the same variables as the corresponding
        # lists on our helper block.
        for b, v in zip(helper.DIFFERENTIAL_BLOCK.values(), helper.differential_vars):
            assert b.var is v

        helper.vectors.derivative.set_setpoint(0.0)
        for var in helper.DERIVATIVE_BLOCK[:].var:
            assert var.setpoint == 0.0
        for b, v in zip(helper.DERIVATIVE_BLOCK.values(), helper.derivative_vars):
            assert b.var is v

        sp = [1, 2, 3]
        helper.vectors.algebraic.set_setpoint(sp)
        for i, j in zip(helper.ALGEBRAIC_SET, sp):
            assert helper.ALGEBRAIC_BLOCK[i].var.setpoint == j

        alg_vars = list(helper.component_objects(AlgVar))
        pred_alg_vars = ComponentSet(
            (
                m.flow_out[t0],
                m.rate[t0, "A"],
                m.rate[t0, "B"],
            )
        )
        assert len(pred_alg_vars) == len(alg_vars)
        for var in alg_vars:
            assert var[t0] in pred_alg_vars
        for b, v in zip(helper.ALGEBRAIC_BLOCK.values(), helper.algebraic_vars):
            assert b.var is v

        assert list(helper.vectors.algebraic.get_setpoint()) == sp

        vals = (10, 11)
        diff_vars.values = vals
        for var, val in zip(helper.DIFFERENTIAL_BLOCK[:].var, vals):
            for v in var[:]:
                assert v.value == val

        _slice = helper.DIFFERENTIAL_BLOCK[:].var
        diff_var_vals = diff_vars.values
        for var, vals in zip(_slice, diff_var_vals):
            for v, vl in zip(var[:], vals):
                assert v.value == vl

    @pytest.mark.unit
    def test_add_references(self):
        m = make_small_model()
        time = m.time
        t0 = time.first()
        blk = DynamicBlock(
            model=m,
            time=time,
            inputs=[m.flow_in[0]],
            measurements=[m.conc[0, "A"], m.conc[0, "B"]],
        )
        blk.construct()

        attrs = [
            "differential",
            "algebraic",
            "derivative",
            "input",
            "fixed",
            "measurement",
        ]

        for attr in attrs:
            # Make sure we've added the expected attribute
            assert hasattr(blk.vectors, attr)
            vec = getattr(blk.vectors, attr)

            # These "vectors" should be two-dimensional;
            # indexed by a coordinate (index into a list) and by time.
            assert vec.dim() == 2
            assert list(vec.index_set().subsets())[1] is time

            # Make sure we can use the set/get setpoint methods.
            # This should be tested more extensively elsewhere.
            setpoints = list(range(len(list(vec[:, t0]))))
            vec.set_setpoint(setpoints)
            assert list(vec.get_setpoint()) == setpoints

            # Make sure the underlying var has the ctype we expect
            # (check using the attribute/component name)
            for var in vec._generate_referenced_vars():
                assert var.ctype._attr == attr

    @pytest.mark.unit
    def test_validate_sample_time(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[t0, "A"], model.conc[t0, "B"]]
        blk = DynamicBlock(
            model=model,
            time=time,
            inputs=inputs,
            measurements=measurements,
        )
        blk.construct()
        blk.validate_sample_time(0.5)
        assert hasattr(blk, "sample_points")
        assert hasattr(blk, "fe_per_sample")
        assert hasattr(blk, "sample_point_indices")

        sample_point_set = set(blk.sample_points)
        sample_point_indices = set(blk.sample_point_indices)
        for p in [0.0, 0.5, 1.0]:
            assert p in sample_point_set
        for i in [1, 3, 5]:
            assert i in sample_point_indices
        assert len(sample_point_set) == 3
        assert len(sample_point_indices) == 3

        with pytest.raises(ValueError, match=r".*integer divider.*"):
            blk.validate_sample_time(0.6)

        with pytest.raises(ValueError, match=r"Could not find a time point.*"):
            blk.validate_sample_time(1 / 3.0)

        with pytest.raises(
            ValueError, match=r".*tolerance is larger than.*not.*unique.*"
        ):
            blk.validate_sample_time(0.5, tolerance=0.09)
            # min spacing in continuous set: 0.166667

        blk.validate_sample_time(0.5, tolerance=0.08)
        sample_point_set = set(blk.sample_points)
        for p in [0.0, 0.5, 1.0]:
            assert p in sample_point_set
        for i in [1, 3, 5]:
            assert i in sample_point_indices
        assert len(sample_point_set) == 3
        assert len(sample_point_indices) == 3

    @pytest.mark.unit
    def test_set_sample_time(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[t0, "A"], model.conc[t0, "B"]]
        blk = DynamicBlock(
            model=model,
            time=time,
            inputs=inputs,
            measurements=measurements,
        )
        blk.construct()

        blk.set_sample_time(1.0)
        assert blk.sample_points == [0.0, 1.0]

        blk.set_sample_time(0.5)
        assert blk.sample_points == [0.0, 0.5, 1.0]

    @pytest.mark.unit
    def make_block(self, sample_time=0.5, horizon=1, nfe=2):
        model = make_model(horizon=horizon, nfe=nfe)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[t0, "A"], model.conc[t0, "B"]]
        scalar_vars, dae_vars = flatten_dae_components(
            model,
            time,
            pyo.Var,
        )
        category_dict = categorize_dae_variables(
            dae_vars, time, inputs, measurements=measurements
        )
        dyn_block = DynamicBlock(
            model=model,
            time=time,
            category_dict={None: category_dict},
            # inputs=inputs,
            # measurements=measurements,
        )
        dyn_block.construct()
        dyn_block.set_sample_time(sample_time)
        return dyn_block

    @pytest.mark.unit
    def test_init_sample_to_setpoint(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()

        # Initialize all data objects with some consistent value
        # so I can make sure certain data objects are skipped
        # by the subsequent initialization
        init_val = -1.0
        blk.vectors.differential[...].set_value(init_val)
        blk.vectors.algebraic[...].set_value(init_val)
        blk.vectors.input[...].set_value(init_val)
        blk.vectors.derivative[...].set_value(init_val)

        # Note that I have no pointer from vectors.differential to
        # the list diff_vars... I can construct an identical list with
        # the _generate_referenced_vars method, but that is O(n)...
        # Not sure yet if this will be a problem/inconvenience...
        blk.vectors.differential.set_setpoint(range(len(blk.differential_vars)))
        blk.vectors.algebraic.set_setpoint(range(len(blk.algebraic_vars)))
        blk.vectors.input.set_setpoint(range(len(blk.input_vars)))
        blk.vectors.derivative.set_setpoint(range(len(blk.derivative_vars)))

        # Possible sample indices: {1, 2}
        sample_points = blk.sample_points
        blk.initialize_sample_to_setpoint(2)
        vectors = [
            blk.vectors.differential,
            blk.vectors.algebraic,
            blk.vectors.input,
            blk.vectors.derivative,
        ]
        for vec in vectors:
            for i, t in vec:
                if t <= sample_points[1]:
                    assert vec[i, t].value == init_val
                else:
                    assert vec[i, t].value == i

        blk.vectors.differential.set_setpoint(
            2 * i for i in range(len(blk.differential_vars))
        )
        blk.vectors.algebraic.set_setpoint(
            2 * i for i in range(len(blk.algebraic_vars))
        )
        blk.vectors.input.set_setpoint(2 * i for i in range(len(blk.input_vars)))
        blk.vectors.derivative.set_setpoint(
            2 * i for i in range(len(blk.derivative_vars))
        )

        blk.initialize_sample_to_setpoint(1)
        for vec in vectors:
            for i, t in vec:
                if t == t0:
                    assert vec[i, t].value == init_val
                elif t <= sample_points[1]:
                    assert vec[i, t].value == 2 * i
                else:
                    assert vec[i, t].value == i

    @pytest.mark.unit
    def test_init_to_setpoint(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()

        # Initialize all data objects with some consistent value
        # so I can make sure certain data objects are skipped
        # by the subsequent initialization
        init_val = -1.0
        blk.vectors.differential[...].set_value(init_val)
        blk.vectors.algebraic[...].set_value(init_val)
        blk.vectors.input[...].set_value(init_val)
        blk.vectors.derivative[...].set_value(init_val)

        # Note that I have no pointer from vectors.differential to
        # the list diff_vars... I can construct an identical list with
        # the _generate_referenced_vars method, but that is O(n)...
        # Not sure yet if this will be a problem/inconvenience...
        blk.vectors.differential.set_setpoint(range(len(blk.differential_vars)))
        blk.vectors.algebraic.set_setpoint(range(len(blk.algebraic_vars)))
        blk.vectors.input.set_setpoint(range(len(blk.input_vars)))
        blk.vectors.derivative.set_setpoint(range(len(blk.derivative_vars)))

        blk.initialize_to_setpoint()
        vectors = [
            blk.vectors.differential,
            blk.vectors.algebraic,
            blk.vectors.input,
            blk.vectors.derivative,
        ]
        for vec in vectors:
            for i, t in vec:
                if t == t0:
                    assert vec[i, t].value == init_val
                else:
                    assert vec[i, t].value == i

        new_sp = 7.0
        blk.vectors.algebraic.set_setpoint(new_sp)
        blk.initialize_to_setpoint(ctype=AlgVar)
        vectors = [
            blk.vectors.differential,
            blk.vectors.input,
            blk.vectors.derivative,
        ]
        for vec in vectors:
            for i, t in vec:
                if t == t0:
                    assert vec[i, t].value == init_val
                else:
                    assert vec[i, t].value == i

        vec = blk.vectors.algebraic
        for i, t in vec:
            if t == t0:
                assert vec[i, t].value == init_val
            else:
                assert vec[i, t].value == new_sp

    @pytest.mark.unit
    def test_init_sample_to_initial(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()

        # So I can tell when something wasn't initialized:
        init_val = -1.0  # `init_val` name might be confusing...
        blk.vectors.differential[...].set_value(init_val)
        blk.vectors.algebraic[...].set_value(init_val)
        blk.vectors.input[...].set_value(init_val)
        blk.vectors.derivative[...].set_value(init_val)

        vectors = [
            blk.vectors.differential,
            blk.vectors.algebraic,
            blk.vectors.derivative,
        ]
        # Set initial conditions to "index"
        for vec in vectors:
            for i, var in enumerate(vec[:, t0]):
                var.set_value(i)

        # Possible sample indices: {1, 2}
        sample_points = blk.sample_points
        blk.initialize_sample_to_initial(2)
        for vec in vectors:
            # We initialized sample 2 to the initial point in that
            # sample. This was just `init_val`
            for i, t in vec:
                if t != t0:
                    assert vec[i, t].value == init_val

        t1 = sample_points[1]
        for vec in vectors:
            for i, var in enumerate(vec[:, t1]):
                var.set_value(i)
        blk.initialize_sample_to_initial(2)
        for vec in vectors:
            for i, t in vec:
                if t == t0:
                    assert vec[i, t].value == i
                elif t < t1:
                    assert vec[i, t].value == init_val
                else:
                    assert vec[i, t].value == i

        blk.initialize_sample_to_initial(1)
        for vec in vectors:
            for i, t in vec:
                assert vec[i, t].value == i

    @pytest.mark.unit
    def test_init_to_initial_conditions(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()

        # So I can tell when something wasn't initialized:
        init_val = -1.0  # `init_val` name might be confusing...
        blk.vectors.differential[...].set_value(init_val)
        blk.vectors.algebraic[...].set_value(init_val)
        blk.vectors.input[...].set_value(init_val)
        blk.vectors.derivative[...].set_value(init_val)

        # Don't have a concise way to set value at t0
        # from a vector.
        # E.g. vectors.differential[:,t0].set_value(range(len(diff_vars)))
        # or vectors.differential[:,t0] = range(len(diff_vars))
        # The latter requires handling iterables in a slice setitem call
        vectors = [
            blk.vectors.differential,
            blk.vectors.algebraic,
            blk.vectors.input,
            blk.vectors.derivative,
        ]
        # Set initial conditions to "index"
        for vec in vectors:
            for i, var in enumerate(vec[:, t0]):
                var.set_value(i)

        blk.initialize_to_initial_conditions(ctype=(DiffVar, AlgVar))
        for vec in (blk.vectors.differential, blk.vectors.algebraic):
            for i, t in vec:
                assert vec[i, t].value == i
        for vec in (blk.vectors.input, blk.vectors.derivative):
            for i, t in vec:
                if t == t0:
                    assert vec[i, t].value == i
                else:
                    assert vec[i, t].value == init_val
                    # value that was used to initialize, not the initial
                    # condition value...

    @pytest.mark.component
    @pytest.mark.skipif(not solver_available, reason="IPOPT is not available")
    def test_init_by_solving_elements(self):
        blk = self.make_block()
        time = blk.time
        blk.mod.flow_in[:].set_value(3.0)
        initialize_t0(blk.mod)
        copy_values_forward(blk.mod)
        blk.mod.flow_in[:].set_value(2.0)
        blk.initialize_by_solving_elements(solver)

        t0 = time.first()
        tl = time.last()
        vectors = blk.vectors
        assert vectors.input[0, tl].value == 2.0
        assert vectors.differential[0, tl].value == pytest.approx(3.185595567867036)
        assert vectors.differential[1, tl].value == pytest.approx(1.1532474073395755)
        assert vectors.derivative[0, tl].value == pytest.approx(0.44321329639889284)
        assert vectors.derivative[1, tl].value == pytest.approx(0.8791007531878847)

        vectors.input.set_setpoint(3.0)
        blk.initialize_by_solving_elements(solver, input_option=InputOption.SETPOINT)
        for t in time:
            if t == t0:
                assert vectors.input[0, t].value == 2.0
            else:
                assert vectors.input[0, t].value == 3.0
        assert vectors.differential[0, tl].value == pytest.approx(3.7037037037037037)
        assert vectors.differential[1, tl].value == pytest.approx(1.0746896480968502)
        assert vectors.derivative[0, tl].value == pytest.approx(0.1851851851851849)
        assert vectors.derivative[1, tl].value == pytest.approx(0.47963475941315314)

    @pytest.mark.component
    @pytest.mark.skipif(not solver_available, reason="IPOPT is not available")
    def test_init_samples_by_element(self):
        blk = self.make_block()
        time = blk.time
        blk.mod.flow_in[:].set_value(3.0)
        initialize_t0(blk.mod)
        copy_values_forward(blk.mod)
        blk.mod.flow_in[:].set_value(2.0)
        blk.set_sample_time(0.5)
        blk.initialize_samples_by_element((1, 2), solver)

        t0 = time.first()
        tl = time.last()
        vectors = blk.vectors
        assert vectors.input[0, tl].value == 2.0
        assert vectors.differential[0, tl].value == pytest.approx(3.185595567867036)
        assert vectors.differential[1, tl].value == pytest.approx(1.1532474073395755)
        assert vectors.derivative[0, tl].value == pytest.approx(0.44321329639889284)
        assert vectors.derivative[1, tl].value == pytest.approx(0.8791007531878847)

        blk.vectors.input.set_setpoint(3.0)
        blk.initialize_samples_by_element(1, solver, input_option=InputOption.SETPOINT)
        blk.initialize_samples_by_element(2, solver, input_option=InputOption.SETPOINT)
        for t in time:
            if t == t0:
                assert vectors.input[0, t].value == 2.0
            else:
                assert vectors.input[0, t].value == 3.0
        assert vectors.differential[0, tl].value == pytest.approx(3.7037037037037037)
        assert vectors.differential[1, tl].value == pytest.approx(1.0746896480968502)
        assert vectors.derivative[0, tl].value == pytest.approx(0.1851851851851849)
        assert vectors.derivative[1, tl].value == pytest.approx(0.47963475941315314)

    @pytest.mark.unit
    def test_advance_by(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()
        tl = time.last()

        for t in time:
            blk.vectors.differential[:, t].set_value(t)
            blk.vectors.algebraic[:, t].set_value(t)
            blk.vectors.derivative[:, t].set_value(t)
            blk.vectors.input[:, t].set_value(t)
            blk.vectors.fixed[:, t].set_value(t)

        ctypes = (DiffVar, DerivVar, AlgVar, InputVar, FixedVar)

        shift = tl - t0  # 1.0, two samples
        blk.advance_by_time(shift)
        for t in time:
            if t == t0:
                for v in blk.component_objects(ctypes):
                    assert v[t].value == tl
            else:
                for v in blk.component_objects(ctypes):
                    assert v[t].value == t

        ctypes_to_shift = (DiffVar, DerivVar)
        ctypes_to_not_shift = (AlgVar, InputVar, FixedVar)
        shift = (tl - t0) / 2  # 0.5, one sample
        idx_of_last_sample = time.find_nearest_index(tl - shift)
        time_of_last_sample = time.at(idx_of_last_sample)
        blk.advance_by_time(shift, ctype=ctypes_to_shift)
        for t in time:
            if t <= shift:
                for v in blk.component_objects(ctypes_to_shift):
                    assert v[t].value == t + shift
            else:
                for v in blk.component_objects(ctypes_to_shift):
                    assert v[t].value == t
            if t == t0:
                for v in blk.component_objects(ctypes_to_not_shift):
                    assert v[t].value == tl
            else:
                for v in blk.component_objects(ctypes_to_not_shift):
                    assert v[t].value == t

    @pytest.mark.unit
    def test_advance_one_sample(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()
        tl = time.last()

        for t in time:
            blk.vectors.differential[:, t].set_value(t)
            blk.vectors.algebraic[:, t].set_value(t)
            blk.vectors.derivative[:, t].set_value(t)
            blk.vectors.input[:, t].set_value(t)
            blk.vectors.fixed[:, t].set_value(t)

        ctypes = (DiffVar, DerivVar, AlgVar, InputVar, FixedVar)

        shift = tl - t0  # 1.0, two samples
        blk.set_sample_time(shift)
        blk.advance_one_sample()
        for t in time:
            if t == t0:
                for v in blk.component_objects(ctypes):
                    assert v[t].value == tl
            else:
                for v in blk.component_objects(ctypes):
                    assert v[t].value == t

        ctypes_to_shift = (DiffVar, DerivVar)
        ctypes_to_not_shift = (AlgVar, InputVar, FixedVar)
        shift = (tl - t0) / 2  # 0.5, one sample
        idx_of_last_sample = time.find_nearest_index(tl - shift)
        time_of_last_sample = time.at(idx_of_last_sample)
        blk.set_sample_time(shift)
        blk.advance_one_sample(ctype=ctypes_to_shift)
        for t in time:
            if t <= shift:
                for v in blk.component_objects(ctypes_to_shift):
                    assert v[t].value == t + shift
            else:
                for v in blk.component_objects(ctypes_to_shift):
                    assert v[t].value == t
            if t == t0:
                for v in blk.component_objects(ctypes_to_not_shift):
                    assert v[t].value == tl
            else:
                for v in blk.component_objects(ctypes_to_not_shift):
                    assert v[t].value == t

    @pytest.mark.unit
    def test_generate_time_in_sample(self):
        blk = self.make_block()
        time = blk.time
        ts = blk.sample_points[1]
        i_s = time.find_nearest_index(ts)

        time_in_sample = list(blk.generate_time_in_sample(ts))
        i_prev = time.find_nearest_index(ts - blk.sample_time)
        pred_time_in_sample = [time.at(i) for i in range(i_prev + 1, i_s + 1)]
        assert time_in_sample == pred_time_in_sample

        time_in_sample = list(blk.generate_time_in_sample(ts, include_t0=True))
        pred_time_in_sample = [time.at(i) for i in range(i_prev, i_s + 1)]
        assert time_in_sample == pred_time_in_sample

        i_0 = 1
        t0 = time.at(i_0)
        ts = time.last()
        time_in_sample = list(blk.generate_time_in_sample(ts, t0=t0))
        pred_time_in_sample = list(t for t in time if t != t0)
        assert time_in_sample == pred_time_in_sample

    @pytest.mark.unit
    def test_get_data_from_sample(self):
        blk = self.make_block()
        time = blk.time
        ts = blk.sample_points[2]

        n_diff = len(blk.DIFFERENTIAL_SET)
        n_input = len(blk.INPUT_SET)
        n_alg = len(blk.ALGEBRAIC_SET)

        for (i, t), var in blk.vectors.differential.items():
            var.set_value(i * t)
        for (i, t), var in blk.vectors.input.items():
            var.set_value((n_diff + i) * t)
        for (i, t), var in blk.vectors.algebraic.items():
            var.set_value((n_diff + n_input + i) * t)

        # Test default variables with out including t0
        data = blk.get_data_from_sample(ts)
        # By default extract values from differential and input variables
        data_time = list(blk.generate_time_in_sample(ts))

        assert len(data) == n_diff + n_input
        for var in blk.component_objects((DiffVar, InputVar)):
            _slice = var.referent
            cuid = pyo.ComponentUID(_slice)
            assert cuid in data

        for i, b in blk.DIFFERENTIAL_BLOCK.items():
            cuid = pyo.ComponentUID(b.var.referent)
            values = list(b.var[t].value for t in data_time)
            assert values == data[cuid]
        for i, b in blk.INPUT_BLOCK.items():
            cuid = pyo.ComponentUID(b.var.referent)
            values = list(b.var[t].value for t in data_time)
            assert values == data[cuid]

        # Test default variables, including the data point at t0
        data = blk.get_data_from_sample(ts, include_t0=True)
        data_time = list(blk.generate_time_in_sample(ts, include_t0=True))
        assert len(data) == n_diff + n_input
        for var in blk.component_objects((DiffVar, InputVar)):
            _slice = var.referent
            cuid = pyo.ComponentUID(_slice)
            assert cuid in data
        for i, b in blk.DIFFERENTIAL_BLOCK.items():
            cuid = pyo.ComponentUID(b.var.referent)
            values = list(b.var[t].value for t in data_time)
            assert values == data[cuid]
        for i, b in blk.INPUT_BLOCK.items():
            cuid = pyo.ComponentUID(b.var.referent)
            values = list(b.var[t].value for t in data_time)
            assert values == data[cuid]

        data = blk.get_data_from_sample(ts, variables=(VariableCategory.ALGEBRAIC,))
        data_time = list(blk.generate_time_in_sample(ts))

        assert len(data) == n_alg
        for var in blk.component_objects(AlgVar):
            cuid = pyo.ComponentUID(var.referent)
            assert cuid in data

        for i, b in blk.ALGEBRAIC_BLOCK.items():
            cuid = pyo.ComponentUID(b.var.referent)
            values = list(b.var[t].value for t in data_time)
            assert values == data[cuid]

    @pytest.mark.unit
    def test_set_variance(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()

        variance_list = [
            (var[t0], 0.05) for var in blk.component_objects(SubclassOf(NmpcVar))
        ]
        blk.set_variance(variance_list)

        for var in blk.DIFFERENTIAL_BLOCK[:].var:
            assert var.variance == 0.05
        for var in blk.INPUT_BLOCK[:].var:
            assert var.variance == 0.05
        for var in blk.ALGEBRAIC_BLOCK[:].var:
            assert var.variance == 0.05
        for var in blk.MEASUREMENT_BLOCK[:].var:
            assert var.variance == 0.05
        for var in blk.FIXED_BLOCK[:].var:
            assert var.variance == 0.05
        for var in blk.DERIVATIVE_BLOCK[:].var:
            assert var.variance == 0.05

    @pytest.mark.unit
    def test_generate_inputs_at_time(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()
        t1 = time.at(2)
        vals = list(0.1 * i for i in blk.INPUT_SET)
        for var, val in zip(blk.vectors.input[:, t1], vals):
            var.set_value(val)
        inputs = list(blk.generate_inputs_at_time(t1))
        assert inputs == vals

    @pytest.mark.unit
    def test_generate_measurements_at_time(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()
        ts = blk.sample_points[1]
        vals = list(0.1 * i for i in blk.MEASUREMENT_SET)
        for var, val in zip(blk.vectors.measurement[:, ts], vals):
            var.set_value(val)
        meas = list(blk.generate_measurements_at_time(ts))
        assert meas == vals

    @pytest.mark.unit
    def test_inject_inputs(self):
        blk = self.make_block()
        time = blk.time
        vals = list(0.25 * i for i in blk.INPUT_SET)
        blk.inject_inputs(vals)
        for b, val in zip(blk.INPUT_BLOCK.values(), vals):
            for t in time:
                assert b.var[t].value == val

    @pytest.mark.unit
    def test_load_measurements(self):
        blk = self.make_block()
        time = blk.time
        t0 = time.first()
        vals = list(0.25 * i for i in blk.MEASUREMENT_SET)
        blk.load_measurements(vals)
        for b, val in zip(blk.MEASUREMENT_BLOCK.values(), vals):
            assert b.var[t0].value == val

    @pytest.mark.unit
    def test_categories_only_measurement_input(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in]
        measurements = [
            pyo.Reference(model.conc[:, "A"]),
            pyo.Reference(model.conc[:, "B"]),
        ]
        category_dict = {
            VC.INPUT: inputs,
            VC.MEASUREMENT: measurements,
        }
        db = DynamicBlock(
            model=model,
            time=time,
            category_dict={None: category_dict},
        )
        db.construct()
        db.set_sample_time(0.5)

        db.vectors.input[:, t0].set_value(1.1)
        db.vectors.measurement[:, t0].set_value(2.2)
        assert model.flow_in[t0].value == 1.1
        assert model.conc[t0, "A"].value == 2.2
        assert model.conc[t0, "B"].value == 2.2

    @pytest.mark.component
    def test_initialize_only_measurement_input(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in]
        measurements = [
            pyo.Reference(model.conc[:, "A"]),
            pyo.Reference(model.conc[:, "B"]),
        ]
        category_dict = {
            VC.INPUT: inputs,
            VC.MEASUREMENT: measurements,
        }
        db = DynamicBlock(
            model=model,
            time=time,
            category_dict={None: category_dict},
        )
        db.construct()
        db.set_sample_time(0.5)

        db.mod.flow_in[:].set_value(3.0)
        initialize_t0(db.mod)
        copy_values_forward(db.mod)
        db.mod.flow_in[:].set_value(2.0)

        # Don't need to know any of the special categories to initialize
        # by element. This is only because we have an implicit discretization.
        db.initialize_by_solving_elements(solver)

        t0 = time.first()
        tl = time.last()
        vectors = db.vectors
        assert vectors.input[0, tl].value == 2.0
        assert vectors.measurement[0, tl].value == pytest.approx(3.185595567867036)
        assert vectors.measurement[1, tl].value == pytest.approx(1.1532474073395755)
        assert model.dcdt[tl, "A"].value == pytest.approx(0.44321329639889284)
        assert model.dcdt[tl, "B"].value == pytest.approx(0.8791007531878847)

    @pytest.mark.unit
    def test_extra_category(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in]
        measurements = [
            pyo.Reference(model.conc[:, "A"]),
            pyo.Reference(model.conc[:, "B"]),
        ]
        disturbances = [
            pyo.Reference(model.conc_in[:, "A"]),
            pyo.Reference(model.conc_in[:, "B"]),
        ]
        category_dict = {
            VC.INPUT: inputs,
            VC.MEASUREMENT: measurements,
            VC.DISTURBANCE: disturbances,
        }
        db = DynamicBlock(
            model=model,
            time=time,
            category_dict={None: category_dict},
        )
        db.construct()
        db.set_sample_time(0.5)

        db.vectors.input[:, t0].set_value(1.1)
        db.vectors.measurement[:, t0].set_value(2.2)
        db.vectors.disturbance[:, t0].set_value(3.3)
        assert model.flow_in[t0].value == 1.1
        assert model.conc[t0, "A"].value == 2.2
        assert model.conc[t0, "B"].value == 2.2
        assert model.conc_in[t0, "A"].value == 3.3
        assert model.conc_in[t0, "B"].value == 3.3

    @pytest.mark.unit
    def test_empty_category(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in]
        measurements = [
            pyo.Reference(model.conc[:, "A"]),
            pyo.Reference(model.conc[:, "B"]),
        ]
        category_dict = {
            VC.INPUT: inputs,
            VC.MEASUREMENT: measurements,
            VC.ALGEBRAIC: [],
        }
        db = DynamicBlock(
            model=model,
            time=time,
            category_dict={None: category_dict},
        )
        db.construct()
        db.set_sample_time(0.5)

        # Categories with no variables are removed.
        # If they were retained, trying to iterate
        # over, e.g., vectors.algebraic[:, :] would
        # fail due to inconsistent dimension.
        assert VC.ALGEBRAIC not in db.category_dict
