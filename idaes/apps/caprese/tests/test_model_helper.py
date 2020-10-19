from pyomo.common.collections import ComponentSet
from idaes.apps.caprese.model import *
from idaes.apps.caprese.common.config import VariableCategory
from idaes.apps.caprese.nmpc_var import *
from idaes.apps.caprese.tests.test_simple_model import (
        make_model,
        make_small_model,
        )

class TestDynamicHelper(object):

    def test_init(self):
        m = make_small_model()
        time = m.time
        t0 = time.first()
        helper = DynamicModelHelper(m, time, inputs=[m.flow_in[0]])
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
        diff_vars = helper.category_dict[VariableCategory.DIFFERENTIAL]
        assert len(pred_diff_vars) == len(diff_vars)
        for var in diff_vars:
            assert var[t0] in pred_diff_vars

        pred_alg_vars = ComponentSet((
                m.flow_out[t0],
                m.rate[t0,'A'],
                m.rate[t0,'B'],
                ))
        alg_vars = helper.category_dict[VariableCategory.ALGEBRAIC]
        assert len(pred_alg_vars) == len(alg_vars)
        for var in alg_vars:
            assert var[t0] in pred_alg_vars

        pred_input_vars = ComponentSet((
                m.flow_in[t0],
                ))
        input_vars = helper.category_dict[VariableCategory.INPUT]
        assert len(pred_input_vars) == len(input_vars)
        for var in input_vars:
            assert var[t0] in pred_input_vars

        pred_fixed_vars = ComponentSet((
                m.conc_in[t0,'A'],
                m.conc_in[t0,'B'],
                ))
        fixed_vars = helper.category_dict[VariableCategory.FIXED]
        assert len(pred_fixed_vars) == len(fixed_vars)
        for var in fixed_vars:
            assert var[t0] in pred_fixed_vars

        pred_deriv_vars = ComponentSet((
                m.dcdt[t0,'A'],
                m.dcdt[t0,'B'],
                ))
        deriv_vars = helper.category_dict[VariableCategory.DERIVATIVE]
        assert len(pred_deriv_vars) == len(deriv_vars)
        for var in deriv_vars:
            assert var[t0] in pred_deriv_vars

    def test_category_blocks(self):
        m = make_small_model()
        time = m.time
        t0 = time.first()
        helper = DynamicModelHelper(m, time, inputs=[m.flow_in[0]])

        diff_vars = helper.namespace.vectors.diff
        diff_var_set = helper.namespace.DIFFERENTIAL_SET*m.time
        assert diff_vars.index_set() == diff_var_set

        helper.namespace.vectors.deriv.set_setpoint(0.0)
        for var in helper.namespace.DERIVATIVE_BLOCK[:].var:
            assert var.setpoint == 0.

        sp = [1,2,3]
        helper.namespace.vectors.alg.set_setpoint(sp)
        for i,j in zip(helper.namespace.ALGEBRAIC_SET, sp):
            assert helper.namespace.ALGEBRAIC_BLOCK[i].var.setpoint == j

        alg_vars = list(m.component_objects(AlgVar))
        pred_alg_vars = ComponentSet((
                m.flow_out[t0],
                m.rate[t0,'A'],
                m.rate[t0,'B'],
                ))
        assert len(pred_alg_vars) == len(alg_vars)
        for var in alg_vars:
            assert var[t0] in pred_alg_vars

        assert list(helper.namespace.vectors.alg.get_setpoint()) == sp
