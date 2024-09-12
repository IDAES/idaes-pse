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
from idaes.apps.flexibility_analysis.var_utils import (
    get_all_unfixed_variables,
    get_used_unfixed_variables,
    BoundsManager,
    _remove_var_bounds,
    _apply_var_bounds,
)
import pytest


@pytest.mark.unit
class TestGetVariables(unittest.TestCase):
    def test_get_all_unfixed_variables(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.y = pe.Var()
        m.x_ref = pe.Reference(m.x)
        m.y.fix(3)

        unfixed_vars = get_all_unfixed_variables(m)
        self.assertEqual(len(unfixed_vars), 1)
        self.assertIn(m.x, unfixed_vars)
        self.assertNotIn(m.y, unfixed_vars)

    def test_get_used_unfixed_variables(self):
        m = pe.ConcreteModel()
        m.x = pe.Var([1, 2, 3, 4, 5, 6])
        m.x[3].fix(1)
        m.x[4].fix(1)
        m.x[6].fix(1)
        m.c1 = pe.Constraint(expr=m.x[1] == m.x[3])
        m.obj = pe.Objective(expr=m.x[5] + m.x[6])
        uuf_vars = get_used_unfixed_variables(m)
        self.assertEqual(len(uuf_vars), 2)
        self.assertIn(m.x[1], uuf_vars)
        self.assertIn(m.x[5], uuf_vars)


@pytest.mark.unit
class TestBounds(unittest.TestCase):
    def test_bounds_manager1(self):
        m = pe.ConcreteModel()
        m.x = pe.Var(bounds=(-1, 1))

        bm = BoundsManager(m)
        bm.save_bounds()
        m.x.setlb(-2)
        self.assertEqual(m.x.lb, -2)
        bm.pop_bounds()
        self.assertEqual(m.x.lb, -1)

    def test_remove_var_bounds(self):
        m = pe.ConcreteModel()
        m.x = pe.Var(bounds=(-1, 1))
        m.y = pe.Var(domain=pe.NonNegativeReals)

        bm = BoundsManager(m)
        bm.save_bounds()

        _remove_var_bounds(m)
        self.assertIsNone(m.x.lb)
        self.assertIsNone(m.x.ub)
        self.assertIsNone(m.y.lb)
        self.assertIsNone(m.y.ub)

        bm.pop_bounds()
        self.assertEqual(m.x.lb, -1)
        self.assertEqual(m.x.ub, 1)
        self.assertEqual(m.y.lb, 0)
        self.assertIsNone(m.y.ub)

    def test_remove_var_bounds_exception(self):
        m = pe.ConcreteModel()
        m.x = pe.Var(domain=pe.Binary)
        with self.assertRaises(ValueError):
            _remove_var_bounds(m)

    def test_apply_var_bounds(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.y = pe.Var(bounds=(-1, 1))

        new_bounds = pe.ComponentMap()
        new_bounds[m.x] = (-5, 5)
        new_bounds[m.y] = (0, 2)

        _apply_var_bounds(new_bounds)
        self.assertEqual(m.x.lb, -5)
        self.assertEqual(m.x.ub, 5)
        self.assertEqual(m.y.lb, 0)
        self.assertEqual(m.y.ub, 1)
