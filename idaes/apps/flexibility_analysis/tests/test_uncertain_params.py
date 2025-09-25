#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
from idaes.apps.flexibility_analysis import _check_dependencies
import pyomo.environ as pe
import unittest
from idaes.apps.flexibility_analysis.uncertain_params import _replace_uncertain_params
from idaes.apps.flexibility_analysis.indices import _ConIndex, _VarIndex
from pyomo.core.expr.compare import compare_expressions
import pytest


@pytest.mark.unit
class TestReplaceUncertainParams(unittest.TestCase):
    def test_replace_mutable_parameter(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.y = pe.Var()
        m.p = pe.Param(mutable=True)
        m.c1 = pe.Constraint(expr=(m.x + m.p, 0))
        m.c2 = pe.Constraint(expr=(None, m.y - m.x, 0))
        m.c3 = pe.Constraint(expr=(-1, m.y + m.p, 1))

        nominal_values = pe.ComponentMap()
        nominal_values[m.p] = 2.3

        param_bounds = pe.ComponentMap()
        param_bounds[m.p] = (1, 4)

        _replace_uncertain_params(
            m=m,
            uncertain_params=[m.p],
            param_nominal_values=nominal_values,
            param_bounds=param_bounds,
        )

        self.assertEqual(len(m.unc_cons), 4)

        self.assertEqual(len(m.unc_param_vars), 1)
        v_ndx = _VarIndex(m.p, None)
        self.assertEqual(m.unc_param_vars[v_ndx].lb, 1)
        self.assertEqual(m.unc_param_vars[v_ndx].ub, 4)
        self.assertEqual(m.unc_param_vars[v_ndx].value, 2.3)

        self.assertFalse(m.c1.active)
        c_ndx = _ConIndex(m.c1, None)
        self.assertEqual(m.unc_cons[c_ndx].lower, 0)
        self.assertEqual(m.unc_cons[c_ndx].upper, 0)
        self.assertTrue(
            compare_expressions(m.unc_cons[c_ndx].body, m.x + m.unc_param_vars[v_ndx])
        )

        self.assertFalse(m.c2.active)
        c_ndx = _ConIndex(m.c2, "ub")
        self.assertEqual(m.unc_cons[c_ndx].lower, None)
        self.assertEqual(m.unc_cons[c_ndx].upper, 0)
        self.assertTrue(compare_expressions(m.unc_cons[c_ndx].body, m.y - m.x))

        self.assertFalse(m.c3.active)
        c_ndx = _ConIndex(m.c3, "lb")
        self.assertEqual(m.unc_cons[c_ndx].lower, -1)
        self.assertEqual(m.unc_cons[c_ndx].upper, None)
        self.assertTrue(
            compare_expressions(m.unc_cons[c_ndx].body, m.y + m.unc_param_vars[v_ndx])
        )

        c_ndx = _ConIndex(m.c3, "ub")
        self.assertEqual(m.unc_cons[c_ndx].lower, None)
        self.assertEqual(m.unc_cons[c_ndx].upper, 1)
        self.assertTrue(
            compare_expressions(m.unc_cons[c_ndx].body, m.y + m.unc_param_vars[v_ndx])
        )

    def test_replace_var(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.p = pe.Var()
        m.c1 = pe.Constraint(expr=(m.x + m.p, 0))

        nominal_values = pe.ComponentMap()
        nominal_values[m.p] = 2.3

        param_bounds = pe.ComponentMap()
        param_bounds[m.p] = (1, 4)

        _replace_uncertain_params(
            m=m,
            uncertain_params=[m.p],
            param_nominal_values=nominal_values,
            param_bounds=param_bounds,
        )

        self.assertFalse(m.c1.active)
        c_ndx = _ConIndex(m.c1, None)
        v_ndx = _VarIndex(m.p, None)
        self.assertEqual(m.unc_cons[c_ndx].lower, 0)
        self.assertEqual(m.unc_cons[c_ndx].upper, 0)
        self.assertTrue(
            compare_expressions(m.unc_cons[c_ndx].body, m.x + m.unc_param_vars[v_ndx])
        )
        self.assertEqual(m.unc_param_vars[v_ndx].lb, 1)
        self.assertEqual(m.unc_param_vars[v_ndx].ub, 4)
        self.assertEqual(m.unc_param_vars[v_ndx].value, 2.3)

    def test_non_constant_var(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.p = pe.Param(initialize=1, mutable=True)
        m.x.setlb(m.p)

        nominal_values = pe.ComponentMap()
        nominal_values[m.p] = 2.3

        param_bounds = pe.ComponentMap()
        param_bounds[m.p] = (1, 4)

        with self.assertRaises(ValueError):
            _replace_uncertain_params(
                m=m,
                uncertain_params=[m.p],
                param_nominal_values=nominal_values,
                param_bounds=param_bounds,
            )

    def test_non_constant_var2(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.p = pe.Param(initialize=1, mutable=True)
        m.x.setub(m.p)

        nominal_values = pe.ComponentMap()
        nominal_values[m.p] = 2.3

        param_bounds = pe.ComponentMap()
        param_bounds[m.p] = (1, 4)

        with self.assertRaises(ValueError):
            _replace_uncertain_params(
                m=m,
                uncertain_params=[m.p],
                param_nominal_values=nominal_values,
                param_bounds=param_bounds,
            )
