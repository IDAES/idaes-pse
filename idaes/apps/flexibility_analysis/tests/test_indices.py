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
import unittest
import pyomo.environ as pe
from idaes.apps.flexibility_analysis.indices import _VarIndex, _ConIndex
import pytest


@pytest.mark.unit
class TestIndices(unittest.TestCase):
    def test_var_index(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.y = pe.Var()

        vi1 = _VarIndex(m.x, "lb")

        self.assertIs(vi1.var, m.x)
        self.assertEqual(vi1.bound, "lb")

        vi2 = _VarIndex(m.x, "lb")
        self.assertEqual(vi1, vi2)
        self.assertEqual(hash(vi1), hash(vi2))

        vi2 = _VarIndex(m.x, "ub")
        self.assertNotEqual(vi1, vi2)

        vi2 = _VarIndex(m.y, "lb")
        self.assertNotEqual(vi1, vi2)

        self.assertEqual(str(vi1), "('x', 'lb')")
        self.assertEqual(repr(vi1), "('x', 'lb')")

    def test_con_index(self):
        m = pe.ConcreteModel()
        m.x = pe.Var()
        m.y = pe.Var()
        m.c1 = pe.Constraint(expr=m.x == m.y)
        m.c2 = pe.Constraint(expr=m.x == 2 * m.y)

        ci1 = _ConIndex(m.c1, "lb")

        self.assertIs(ci1.con, m.c1)
        self.assertEqual(ci1.bound, "lb")

        ci2 = _ConIndex(m.c1, "lb")
        self.assertEqual(ci1, ci2)
        self.assertEqual(hash(ci1), hash(ci2))

        ci2 = _ConIndex(m.c1, "ub")
        self.assertNotEqual(ci1, ci2)

        ci2 = _ConIndex(m.c2, "lb")
        self.assertNotEqual(ci1, ci2)

        self.assertEqual(str(ci1), "('c1', 'lb')")
        self.assertEqual(repr(ci1), "('c1', 'lb')")

        vi1 = _VarIndex(m.x, "lb")
        self.assertNotEqual(ci1, vi1)
        self.assertNotEqual(vi1, ci1)
