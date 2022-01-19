import pyomo.common.unittest as unittest
import pytest

import pyomo.environ as pyo
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from idaes.apps.nmpc.cost_expressions import (
    get_tracking_cost_from_constant_setpoint,
)

@pytest.mark.unit
class TestTrackingCost(unittest.TestCase):

    def test_tracking_cost_no_weights(self):
        m = pyo.ConcreteModel()
        m.time = pyo.Set(initialize=[1, 2, 3])
        m.v1 = pyo.Var(m.time, initialize={i: 1*i for i in m.time})
        m.v2 = pyo.Var(m.time, initialize={i: 2*i for i in m.time})

        setpoint_data = {
            str(pyo.ComponentUID(m.v1)): 3.0,
            str(pyo.ComponentUID(m.v2)): 4.0,
        }

        m.tracking_expr = get_tracking_cost_from_constant_setpoint(
            [m.v1, m.v2],
            m.time,
            setpoint_data,
        )

        var_sets = {
            i: ComponentSet(identify_variables(m.tracking_expr[i]))
            for i in m.time
        }
        for i in m.time:
            self.assertIn(m.v1[i], var_sets[i])
            self.assertIn(m.v2[i], var_sets[i])
            pred_value = (1*i - 3)**2 + (2*i - 4)**2
            self.assertEqual(pred_value, pyo.value(m.tracking_expr[i]))

    def test_tracking_cost_with_weights(self):
        m = pyo.ConcreteModel()
        m.time = pyo.Set(initialize=[1, 2, 3])
        m.v1 = pyo.Var(m.time, initialize={i: 1*i for i in m.time})
        m.v2 = pyo.Var(m.time, initialize={i: 2*i for i in m.time})

        setpoint_data = {
            str(pyo.ComponentUID(m.v1)): 3.0,
            str(pyo.ComponentUID(m.v2)): 4.0,
        }
        weight_data = {
            str(pyo.ComponentUID(m.v1)): 0.1,
            str(pyo.ComponentUID(m.v2)): 0.5,
        }

        m.tracking_expr = get_tracking_cost_from_constant_setpoint(
            [m.v1, m.v2],
            m.time,
            setpoint_data,
            weight_data=weight_data,
        )

        var_sets = {
            i: ComponentSet(identify_variables(m.tracking_expr[i]))
            for i in m.time
        }
        for i in m.time:
            self.assertIn(m.v1[i], var_sets[i])
            self.assertIn(m.v2[i], var_sets[i])
            pred_value = 0.1*(1*i - 3)**2 + 0.5*(2*i - 4)**2
            self.assertAlmostEqual(pred_value, pyo.value(m.tracking_expr[i]))

    def test_exceptions(self):
        m = pyo.ConcreteModel()
        m.time = pyo.Set(initialize=[1, 2, 3])
        m.v1 = pyo.Var(m.time, initialize={i: 1*i for i in m.time})
        m.v2 = pyo.Var(m.time, initialize={i: 2*i for i in m.time})

        setpoint_data = {
            str(pyo.ComponentUID(m.v1)): 3.0,
        }
        weight_data = {
            str(pyo.ComponentUID(m.v1)): 0.1,
        }

        with self.assertRaisesRegex(KeyError, "Setpoint data"):
            m.tracking_expr = get_tracking_cost_from_constant_setpoint(
                [m.v1, m.v2],
                m.time,
                setpoint_data,
            )

        setpoint_data = {
            str(pyo.ComponentUID(m.v1)): 3.0,
            str(pyo.ComponentUID(m.v2)): 4.0,
        }

        with self.assertRaisesRegex(KeyError, "Tracking weight"):
            m.tracking_expr = get_tracking_cost_from_constant_setpoint(
                [m.v1, m.v2],
                m.time,
                setpoint_data,
                weight_data=weight_data,
            )
