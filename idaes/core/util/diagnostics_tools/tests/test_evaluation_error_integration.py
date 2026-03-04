#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains tests for integrating evaluation error detection into
the DiagnosticsToolbox.
"""

from io import StringIO
import math
from unittest import TestCase

import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Objective,
    Var,
    Param,
    Integers,
    log,
    tan,
    asin,
    acos,
    sqrt,
)

from idaes.core.util.diagnostics_tools.diagnostics_toolbox import (
    DiagnosticsToolbox,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


class TestEvalErrorDetection(TestCase):
    @pytest.mark.unit
    def test_div(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == 1 / m.x)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(0)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w, "c: Potential division by 0 in 1/x; Denominator bounds are (0, inf)"
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w, "c: Potential division by 0 in 1/x; Denominator bounds are (-1, inf)"
        )

    @pytest.mark.unit
    def test_pow1(self):
        m = ConcreteModel()
        m.x = Var(bounds=(None, None))
        m.y = Var()
        m.p = Param(initialize=2, mutable=True)
        m.c = Constraint(expr=m.y == m.x**m.p)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.p.value = 2.5
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**p; base bounds are (-inf, inf); exponent bounds are (2.5, 2.5)",
        )

        m.x.setlb(1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

    @pytest.mark.unit
    def test_pow2(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.p = Var(domain=Integers)
        m.c = Constraint(expr=m.y == m.x**m.p)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**p; base bounds are (-inf, inf); exponent bounds are (-inf, inf)",
        )

    @pytest.mark.unit
    def test_pow3(self):
        m = ConcreteModel()
        m.x = Var(bounds=(0, None))
        m.y = Var()
        m.p = Var(bounds=(0, None))
        m.c = Constraint(expr=m.y == m.x**m.p)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

    @pytest.mark.unit
    def test_pow4(self):
        m = ConcreteModel()
        m.x = Var(bounds=(0, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == m.x ** (-2))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**(-2); base bounds are (0, inf); exponent bounds are (-2, -2)",
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**(-2); base bounds are (-1, inf); exponent bounds are (-2, -2)",
        )

    @pytest.mark.unit
    def test_pow5(self):
        m = ConcreteModel()
        m.x = Var(bounds=(0, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == m.x ** (-2.5))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**(-2.5); base bounds are (0, inf); exponent bounds are (-2.5, -2.5)",
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**(-2.5); base bounds are (-1, inf); exponent bounds are (-2.5, -2.5)",
        )

    @pytest.mark.unit
    def test_log(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == log(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(0)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential log of a non-positive number in log(x); Argument bounds are (0, inf)",
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential log of a non-positive number in log(x); Argument bounds are (-1, inf)",
        )

    @pytest.mark.unit
    def test_tan(self):
        m = ConcreteModel()
        m.x = Var(bounds=(-math.pi / 4, math.pi / 4))
        m.y = Var()
        m.c = Constraint(expr=m.y == tan(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-math.pi)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: tan(x) may evaluate to -inf or inf; Argument bounds are (-3.141592653589793, 0.7853981633974483)",
        )

    @pytest.mark.unit
    def test_asin(self):
        m = ConcreteModel()
        m.x = Var(bounds=(-0.5, 0.5))
        m.y = Var()
        m.c = Constraint(expr=m.y == asin(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of asin outside [-1, 1] in asin(x); Argument bounds are (-inf, 0.5)",
        )

        m.x.setlb(-0.5)
        m.x.setub(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of asin outside [-1, 1] in asin(x); Argument bounds are (-0.5, inf)",
        )

    @pytest.mark.unit
    def test_acos(self):
        m = ConcreteModel()
        m.x = Var(bounds=(-0.5, 0.5))
        m.y = Var()
        m.c = Constraint(expr=m.y == acos(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of acos outside [-1, 1] in acos(x); Argument bounds are (-inf, 0.5)",
        )

        m.x.setlb(-0.5)
        m.x.setub(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of acos outside [-1, 1] in acos(x); Argument bounds are (-0.5, inf)",
        )

    @pytest.mark.unit
    def test_sqrt(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == sqrt(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential square root of a negative number in sqrt(x); Argument bounds are (-1, inf)",
        )

    @pytest.mark.unit
    def test_display(self):
        stream = StringIO()
        m = ConcreteModel()
        m.x = Var()
        m.y = Var()
        m.obj = Objective(expr=m.x**2 + m.y**2.5)
        m.c1 = Constraint(expr=m.y >= log(m.x))
        m.c2 = Constraint(expr=m.y >= (m.x - 1) ** 2.5)
        m.c3 = Constraint(expr=m.x - 1 >= 0)
        dtb = DiagnosticsToolbox(m)
        dtb.display_potential_evaluation_errors(stream=stream)
        expected = "====================================================================================\n3 WARNINGS\n\n    c1: Potential log of a non-positive number in log(x); Argument bounds are (-inf, inf)\n    c2: Potential evaluation error in (x - 1)**2.5; base bounds are (-inf, inf); exponent bounds are (2.5, 2.5)\n    obj: Potential evaluation error in y**2.5; base bounds are (-inf, inf); exponent bounds are (2.5, 2.5)\n\n====================================================================================\n"
        got = stream.getvalue()
        exp_list = expected.split("\n")
        got_list = got.split("\n")
        self.assertEqual(len(exp_list), len(got_list))
        for _exp, _got in zip(exp_list, got_list):
            self.assertEqual(_exp, _got)
