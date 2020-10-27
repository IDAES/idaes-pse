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
"""

import pyomo.environ as aml
import pyomo.dae as dae
import pyomo.network as pyn
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.core.base.block import _BlockData

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom,
        activated_equalities_generator, unfixed_variables_generator)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.apps.caprese.tests.test_simple_model import (
        make_model, 
        make_small_model,
        initialize_t0,
        copy_values_forward,
        )
from idaes.apps.caprese.model import (
        DynamicBlock,
        SimpleDynamicBlock,
        IndexedDynamicBlock,
        _DynamicBlockData,
        )
from idaes.apps.caprese.common.config import VariableCategory
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

class TestDynamicBlock(object):

    def test_construct_simple(self):
        model = make_model(horizon=1, nfe=2)
        time = model.time
        t0 = time.first()
        inputs = [model.flow_in[t0]]
        measurements = [model.conc[0,'A'], model.conc[0,'B']]
        block = DynamicBlock(
                model=model,
                time=time,
                inputs=inputs,
                measurements=measurements,
                )
        assert type(block) is SimpleDynamicBlock
        assert isinstance(block, DynamicBlock)
        assert isinstance(block, _DynamicBlockData)

        block.construct()
        assert block[None] is block
        assert all(b is block for b in block[:])

        assert block.mod is model
        assert block.time is time
        assert all(i1 is i2 for i1, i2 in zip(block._inputs, inputs))
        assert all(i1 is i2 for i1, i2 in zip(block._measurements, measurements))

        assert hasattr(block, 'category_dict')
        assert hasattr(block, 'measured_vars')
        assert hasattr(block, 'vardata_map')

        subblocks = [
                block.mod,
                block.vectors,
                block.DIFFERENTIAL_BLOCK,
                block.ALGEBRAIC_BLOCK,
                block.INPUT_BLOCK,
                block.FIXED_BLOCK,
                block.DERIVATIVE_BLOCK,
                ]

        block_objects = ComponentSet(
                block.component_objects(aml.Block, descend_into=False))
        assert len(subblocks) == len(block_objects)
        for b in subblocks:
            assert b in block_objects

        block.v = aml.Var(initialize=3)
        block.c = aml.Constraint(expr=block.v==5)
        assert block.v.value == 3
        assert block.v in ComponentSet(identify_variables(block.c.expr))

    def test_construct_indexed(self):
        block_set = aml.Set(initialize=[0,1,2])
        block_set.construct()
        horizon_map = {0: 1., 1: 3., 2: 5.}
        nfe_map = {0: 2, 1: 6, 2: 10}
        model_map = {i: make_model(horizon_map[i], nfe_map[i])
                for i in block_set}
        time_map = {i: model_map[i].time for i in block_set}
        inputs_map = {i: [model_map[i].flow_in[0]] for i in block_set}
        measurements_map = {
                i: [model_map[i].conc[0,'A'], model_map[i].conc[0,'B']]
                for i in block_set
                }
        block = DynamicBlock(
                block_set,
                model=model_map,
                time=time_map,
                inputs=inputs_map,
                measurements=measurements_map,
                )
        assert type(block) is IndexedDynamicBlock
        assert isinstance(block, DynamicBlock)

        block.construct()
        # Why is block._data empty here?
        # Are blocks constructed sparse by default?
        assert all(b.parent_component() is block for b in block.values())

        for i in block_set:
            assert i in block

        for i, b in block.items():
            assert b.mod is model_map[i]
            assert b.time is time_map[i]
            assert all(i1 is i2 for i1, i2 in zip(b._inputs, inputs_map[i]))
            assert all(i1 is i2 for i1, i2 in 
                    zip(b._measurements, measurements_map[i]))

            assert hasattr(b, 'category_dict')
            assert hasattr(b, 'measured_vars')
            assert hasattr(b, 'vardata_map')

            subblocks = [
                    b.mod,
                    b.vectors,
                    b.DIFFERENTIAL_BLOCK,
                    b.ALGEBRAIC_BLOCK,
                    b.INPUT_BLOCK,
                    b.FIXED_BLOCK,
                    b.DERIVATIVE_BLOCK,
                    ]

            block_objects = ComponentSet(
                    b.component_objects(aml.Block, descend_into=False))
            assert len(subblocks) == len(block_objects)
            for sb in subblocks:
                assert sb in block_objects

            b.v = aml.Var(initialize=3)
            b.c = aml.Constraint(expr=b.v==5)
            assert b.v.value == 3
            assert b.v in ComponentSet(identify_variables(b.c.expr))

    def test_construct_rule(self):
        block_set = aml.Set(initialize=range(3))
        block_set.construct()
        horizon_map = {0: 1, 1: 3, 2: 5}
        nfe_map = {0: 2, 1: 6, 2: 10}
        model_map = {
                i: make_model(horizon=horizon_map[i], nfe=nfe_map[i])
                for i in block_set
                }

        def dynamic_block_rule(b, i):
            model = model_map[i]
            time = model.time
            t0 = time.first()
            inputs = [model.flow_in[t0]]
            measurements = [model.conc[0,'A'], model.conc[0,'B']]

            # Won't be obvious that these attrs need to be set if
            # constructing from a rule
            b.mod = model
            super(_BlockData, b).__setattr__('time', time)
            b._inputs = inputs
            b._measurements = measurements

        block = DynamicBlock(block_set, rule=dynamic_block_rule)
        assert type(block) is IndexedDynamicBlock
        assert isinstance(block, DynamicBlock)

        block.construct()
        # Why is block._data empty here?
        # Are blocks constructed sparse by default?
        assert all(b.parent_component() is block for b in block.values())

        for i in block_set:
            assert i in block

        for i, b in block.items():
            assert b.mod is model_map[i]
            assert b.time is model_map[i].time
            t0 = b.time.first()
            assert all(i1 is i2 for i1, i2 in zip(b._inputs, 
                [model_map[i].flow_in[t0]]))
            assert all(i1 is i2 for i1, i2 in zip(b._measurements, 
                [model_map[i].conc[t0,'A'], model_map[i].conc[t0,'B']]))

            assert hasattr(b, 'category_dict')
            assert hasattr(b, 'measured_vars')
            assert hasattr(b, 'vardata_map')

            subblocks = [
                    b.mod,
                    b.vectors,
                    b.DIFFERENTIAL_BLOCK,
                    b.ALGEBRAIC_BLOCK,
                    b.INPUT_BLOCK,
                    b.FIXED_BLOCK,
                    b.DERIVATIVE_BLOCK,
                    ]

            block_objects = ComponentSet(
                    b.component_objects(aml.Block, descend_into=False))
            assert len(subblocks) == len(block_objects)
            for sb in subblocks:
                assert sb in block_objects

            b.v = aml.Var(initialize=3)
            b.c = aml.Constraint(expr=b.v==5)
            assert b.v.value == 3
            assert b.v in ComponentSet(identify_variables(b.c.expr))

    def test_init(self):
        m = make_small_model()
        time = m.time
        t0 = time.first()
        helper = DynamicBlock(
                model=m,
                time=time,
                inputs=[m.flow_in[0]],
                measurements=[m.conc[0,'A'], m.conc[0,'B']])
        helper.construct()
        assert hasattr(helper, 'category_dict')
        assert hasattr(helper, 'measured_vars')
        assert hasattr(helper, 'vardata_map')

        # Make sure category dict contains variables we expect
        assert VariableCategory.DIFFERENTIAL in helper.category_dict
        assert VariableCategory.ALGEBRAIC in helper.category_dict
        assert VariableCategory.DERIVATIVE in helper.category_dict
        assert VariableCategory.INPUT in helper.category_dict
        assert VariableCategory.FIXED in helper.category_dict
    
        pred_diff_vars = ComponentSet((
                m.conc[t0,'A'],
                m.conc[t0,'B'],
                ))
        diff_vars = list(helper.component_objects(DiffVar))
        assert len(pred_diff_vars) == len(diff_vars)
        for var in diff_vars:
            assert var[t0] in pred_diff_vars

        pred_alg_vars = ComponentSet((
                m.flow_out[t0],
                m.rate[t0,'A'],
                m.rate[t0,'B'],
                ))
        alg_vars = list(helper.component_objects(AlgVar))
        assert len(pred_alg_vars) == len(alg_vars)
        for var in alg_vars:
            assert var[t0] in pred_alg_vars

        pred_input_vars = ComponentSet((
                m.flow_in[t0],
                ))
        input_vars = list(helper.component_objects(InputVar))
        assert len(pred_input_vars) == len(input_vars)
        for var in input_vars:
            assert var[t0] in pred_input_vars

        pred_fixed_vars = ComponentSet((
                m.conc_in[t0,'A'],
                m.conc_in[t0,'B'],
                ))
        fixed_vars = list(helper.component_objects(FixedVar))
        assert len(pred_fixed_vars) == len(fixed_vars)
        for var in fixed_vars:
            assert var[t0] in pred_fixed_vars

        pred_deriv_vars = ComponentSet((
                m.dcdt[t0,'A'],
                m.dcdt[t0,'B'],
                ))
        deriv_vars = list(helper.component_objects(DerivVar))
        assert len(pred_deriv_vars) == len(deriv_vars)
        for var in deriv_vars:
            assert var[t0] in pred_deriv_vars

