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
This module contains tests for the constraint term analysis tools.
"""

import re
import pytest

from pyomo.environ import (
    ConcreteModel,
    value,
    Var,
    Param,
    exp,
    Expr_if,
    log,
    log10,
    sin,
    cos,
    tan,
    asin,
    acos,
    atan,
    sqrt,
    sinh,
    cosh,
    tanh,
    asinh,
    acosh,
    atanh,
)
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.core import expr as EXPR

from idaes.models.properties.modular_properties.eos.ceos_common import (
    cubic_roots_available,
    CubicThermoExpressions,
    CubicType as CubicEoS,
)
from idaes.core.util.diagnostics_tools.constraint_term_analysis import (
    ConstraintTermAnalysisVisitor,
)
from idaes.models.properties import iapws95

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


class TestConstraintTermAnalysisVisitor:
    @pytest.mark.unit
    def test_sum_combinations(self):
        # Check method to generate sums of all combinations of terms
        # excludes single term sums
        terms = [1, 2, 3, 4, 5]
        visitor = ConstraintTermAnalysisVisitor(max_canceling_terms=None)
        sums = [i for i in visitor._generate_combinations(terms)]

        print(sums)

        expected = [
            ((0, 1), (1, 2)),
            ((0, 1), (2, 3)),
            ((0, 1), (3, 4)),
            ((0, 1), (4, 5)),
            ((1, 2), (2, 3)),
            ((1, 2), (3, 4)),
            ((1, 2), (4, 5)),
            ((2, 3), (3, 4)),
            ((2, 3), (4, 5)),
            ((3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3)),
            ((0, 1), (1, 2), (3, 4)),
            ((0, 1), (1, 2), (4, 5)),
            ((0, 1), (2, 3), (3, 4)),
            ((0, 1), (2, 3), (4, 5)),
            ((0, 1), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4)),
            ((1, 2), (2, 3), (4, 5)),
            ((1, 2), (3, 4), (4, 5)),
            ((2, 3), (3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3), (3, 4)),
            ((0, 1), (1, 2), (2, 3), (4, 5)),
            ((0, 1), (1, 2), (3, 4), (4, 5)),
            ((0, 1), (2, 3), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5)),
        ]

        assert sums == expected

    @pytest.mark.unit
    def test_sum_combinations_limited_default(self):
        # Check method to generate sums of all combinations of terms
        # excludes single term sums
        terms = [1, 2, 3, 4, 5]
        visitor = ConstraintTermAnalysisVisitor()
        sums = [i for i in visitor._generate_combinations(terms)]

        expected = expected = [
            ((0, 1), (1, 2)),
            ((0, 1), (2, 3)),
            ((0, 1), (3, 4)),
            ((0, 1), (4, 5)),
            ((1, 2), (2, 3)),
            ((1, 2), (3, 4)),
            ((1, 2), (4, 5)),
            ((2, 3), (3, 4)),
            ((2, 3), (4, 5)),
            ((3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3)),
            ((0, 1), (1, 2), (3, 4)),
            ((0, 1), (1, 2), (4, 5)),
            ((0, 1), (2, 3), (3, 4)),
            ((0, 1), (2, 3), (4, 5)),
            ((0, 1), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4)),
            ((1, 2), (2, 3), (4, 5)),
            ((1, 2), (3, 4), (4, 5)),
            ((2, 3), (3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3), (3, 4)),
            ((0, 1), (1, 2), (2, 3), (4, 5)),
            ((0, 1), (1, 2), (3, 4), (4, 5)),
            ((0, 1), (2, 3), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4), (4, 5)),
            # 5 term combination should be excluded
        ]

        assert sums == expected

    @pytest.mark.unit
    def test_check_sum_cancellations(self):
        terms = [1, -2, 3, -4, 5]
        visitor = ConstraintTermAnalysisVisitor()
        cancellations = visitor._check_sum_cancellations(terms)

        # We expect to canceling combinations
        # Results should be a list with each entry being a tuple containing a
        # canceling combination
        # In turn, each element of the tuple should be a 2-tuple with the
        # position of the term in the input list and its value
        expected = [((0, 1), (2, 3), (3, -4)), ((0, 1), (1, -2), (3, -4), (4, 5))]

        assert cancellations == expected

    @pytest.mark.unit
    def test_int(self):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=7)

        assert vv == [7]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_float(self):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=7.7)

        assert vv == [7.7]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_scalar_param(self):
        m = ConcreteModel()
        m.scalar_param = Param(initialize=1, mutable=True)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.scalar_param
        )

        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_indexed_param(self):
        m = ConcreteModel()
        m.indexed_param = Param(["a", "b"], initialize=1, mutable=True)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.indexed_param["a"]
        )
        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_scalar_var(self):
        m = ConcreteModel()
        m.scalar_var = Var(initialize=1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.scalar_var
        )

        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_indexed_var(self):
        m = ConcreteModel()
        m.indexed_var = Var(["a", "b"], initialize=1)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.indexed_var["a"]
        )
        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_scalar_var_fixed(self):
        m = ConcreteModel()
        m.scalar_var = Var(initialize=1)
        m.scalar_var.fix()
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.scalar_var
        )

        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_indexed_var_fixed(self):
        m = ConcreteModel()
        m.indexed_var = Var(["a", "b"], initialize=1)
        m.indexed_var.fix()

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.indexed_var["a"]
        )
        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_equality_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e-7)
        m.v2 = Var(initialize=1e7)

        expr = m.v1 == m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_equality_expr_constant(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e-7)
        m.v2 = Var(initialize=1e7)

        # Fix v1, not constant yet as v2 still free
        m.v1.fix()

        expr = m.v1 == m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Fix v2, now constant
        m.v2.fix()

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v1["a"].set_value(1e-7)

        expr = sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [1e-7, 1e7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_sum_expr_constant(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v1["a"].set_value(1e-7)
        m.v1.fix()

        expr = sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [1e-7, 1e7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=3)
        m.v3 = Var(initialize=5.0000001)
        m.v4 = Var(initialize=5)

        return m

    @pytest.mark.unit
    def test_product_expr(self, model):
        m = ConcreteModel()
        expr = model.v1 * model.v2 * model.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(30.0000006, rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_product_sum_expr(self, model):
        expr = (model.v1 + model.v2) * (model.v3 + model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx((2 + 3) * (5.0000001 + 5), rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_product_sum_expr_w_negation(self, model):
        expr = (model.v1 + model.v2) * (model.v3 - model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(0.0000005, rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_expr(self, model):
        expr = model.v1 / model.v2 / model.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(2 / 3 / 5.0000001, rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_sum_expr(self, model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=(model.v1 + model.v2) / (model.v3 + model.v4)
        )

        assert vv == [
            pytest.approx(2 / (5.0000001 + 5), rel=1e-8),
            pytest.approx(3 / (5.0000001 + 5), rel=1e-8),
        ]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_sum_expr_w_negation(self, model):
        expr = (model.v1 - model.v2) / (model.v3 - model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [
            pytest.approx(2 / 0.0000001, rel=1e-8),
            pytest.approx(-3 / 0.0000001, rel=1e-8),
        ]
        assert len(mm) == 0
        # Check for cancellation should be deferred
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_sum_expr_w_zero_denominator(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=3)
        m.v3 = Var(initialize=5)
        m.v4 = Var(initialize=5)

        expr = (m.v1 - m.v2) / (m.v3 - m.v4)
        with pytest.raises(
            ZeroDivisionError,
            match=re.escape(
                "Error in ConstraintTermAnalysisVisitor: found division with denominator of 0 "
                "((v1 - v2)/(v3 - v4))."
            ),
        ):
            ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

    @pytest.mark.unit
    def test_division_sum_expr_w_negation_deferred(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=5.0000001)
        m.v4 = Var(initialize=5)

        expr = ((m.v1 - m.v2) / (m.v3 - m.v4)) ** 2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(0, abs=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_expr_error(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=0)

        with pytest.raises(
            ZeroDivisionError,
            match=re.escape(
                "Error in ConstraintTermAnalysisVisitor: found division with "
                "denominator of 0 (v1/v2)."
            ),
        ):
            ConstraintTermAnalysisVisitor().walk_expression(expr=m.v1 / m.v2)

    @pytest.mark.unit
    def test_pow_expr(self, model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=model.v1**model.v2
        )

        assert vv == [8]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_pow_sum_expr(self, model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=(model.v1 + model.v2) ** (model.v3 + model.v4)
        )

        assert vv == [pytest.approx(5**10.0000001, rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_pow_sum_expr_w_negation(self, model):
        expr = (model.v1 + model.v2) ** (model.v3 - model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx((2 + 3) ** (0.0000001), rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.fixture(scope="class")
    def func_model(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=0.99999)
        m.v3 = Var(initialize=-100)

        m.v2.fix()

        return m

    @pytest.mark.unit
    def test_negation_expr(self, func_model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=-func_model.v1
        )

        assert vv == [-1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_negation_sum_expr(self, func_model):
        expr = -(func_model.v1 - func_model.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        # Checking the cancellation should be deferred
        # Expect to get two values back
        assert vv == [pytest.approx(-1, rel=1e-8), pytest.approx(0.99999, rel=1e-8)]
        assert len(mm) == 0
        # Check is deferred, so no cancellations identified
        assert len(cc) == 0
        assert not k
        assert not tr

    # acosh has bounds that don't fit with other functions - we will assume we caught enough
    func_list = [
        exp,
        log,
        log10,
        sqrt,
        sin,
        cos,
        tan,
        asin,
        acos,
        atan,
        sinh,
        cosh,
        tanh,
        asinh,
        atanh,
    ]
    func_error_list = [log, log10, sqrt, asin, acos, acosh, atanh]

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_list)
    def test_func_expr(self, func_model, func):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=func(func_model.v2)
        )

        assert vv == [pytest.approx(value(func(0.99999)), rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_list)
    def test_func_sum_expr(self, func_model, func):
        expr = func(func_model.v1 - func_model.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(value(func(0.00001)), rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_list)
    def test_func_sum_expr_w_negation(self, func_model, func):
        expr = func(func_model.v1 - func_model.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(value(func(0.00001)), rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_error_list)
    def test_func_expr_error(self, func_model, func):
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Error in ConstraintTermAnalysisVisitor: error evaluating "
            ),
        ):
            ConstraintTermAnalysisVisitor().walk_expression(expr=func(func_model.v3))

    @pytest.mark.unit
    def test_expr_if(self, model):
        model.exprif = Expr_if(
            IF=model.v1,
            THEN=model.v1 + model.v2,
            ELSE=model.v3 - model.v4,
        )

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=model.exprif
        )

        assert vv == [pytest.approx(5, rel=1e-8)]
        assert len(mm) == 0
        assert model.exprif in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.skipif(
        not AmplInterface.available(), reason="pynumero_ASL is not available"
    )
    @pytest.mark.skipif(not cubic_roots_available, reason="Cubic roots not available")
    def test_ext_func(self):
        # Use the cubic root external function to test
        m = ConcreteModel()
        m.a = Var(initialize=1)
        m.b = Var(initialize=1)

        m.expr_write = CubicThermoExpressions(m)
        Z = m.expr_write.z_liq(eos=CubicEoS.PR, A=m.a, B=m.b)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=Z)

        assert vv == [pytest.approx(-2.1149075414767577, rel=1e-8)]
        assert len(mm) == 0
        assert Z in cc
        assert len(cc) == 1
        assert not k
        assert not tr

        # Check that model state did not change
        assert value(m.a) == 1
        assert value(m.b) == 1
        assert value(Z) == pytest.approx(-2.1149075414767577, rel=1e-8)

    @pytest.mark.unit
    def test_equality_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)

        expr = m.v2 == sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 3e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_inequality_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)

        expr = m.v2 <= sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 3e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_compound_equality_expr_1(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)

        expr = 6 * m.v2 == 8 * sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-6e-7, 2.4e8]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_ranged_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e7)
        m.v2 = Var(initialize=1e-7)
        m.v3 = Var(initialize=1e7)

        m.expr = EXPR.RangedExpression(args=(m.v1, m.v2, m.v3), strict=True)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=m.expr)

        assert vv == [-1e7, 1e-7, -1e7]
        assert m.expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Fix v1 and v2 to make first two terms constant
        m.v1.fix()
        m.v2.fix()
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=m.expr)

        # Should not be flagged as constant due to v3
        assert vv == [-1e7, 1e-7, -1e7]
        assert m.expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Fix v3 to make all terms constant
        m.v3.fix()
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=m.expr)

        # Should now be constant
        assert vv == [-1e7, 1e-7, -1e7]
        assert m.expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_compound_equality_expr_2(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)
        m.v3 = Var(initialize=1e3)

        # Set this small so we get two mismatched warnings
        m.v1["a"].set_value(1e-7)

        expr1 = sum(m.v1[i] for i in m.v1)
        expr = 6 * m.v2 == 8 * expr1 + m.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-6e-7, pytest.approx(8 * (2e7 + 1e-7) + 1000, rel=1e-8)]
        assert expr in mm
        assert expr1 in mm
        assert len(mm) == 2
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.v1 - m.v2
        )

        assert vv == [2, -2]
        assert len(mm) == 0
        # We do not check cancellation at the sum, so this should be empty
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_expr_w_canceling_sum(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=3)

        expr = m.v3 * (m.v1 - m.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [6, -6]
        assert len(mm) == 0
        # Check for cancellation should be deferred
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_expr_w_deferred_canceling_sum(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=3)

        expr = (m.v3 * (m.v1 - m.v2)) ** 2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0]
        assert len(mm) == 0
        # We should get a warning about canceling sums here
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

        # Check for tolerance of sum cancellation
        m.v2.set_value(2.00022)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor(
            term_cancellation_tolerance=1e-4
        ).walk_expression(expr=expr)

        assert vv == [pytest.approx((3 * -0.00022) ** 2, rel=1e-8)]
        assert len(mm) == 0
        # This should pass as the difference is greater than tol
        assert len(cc) == 0
        assert not k
        assert not tr

        # Change tolerance so it should identify cancellation
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor(
            term_cancellation_tolerance=1e-3
        ).walk_expression(expr=expr)

        assert vv == [pytest.approx((3 * -0.00022) ** 2, rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_equality_expr_safe(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e7)
        m.v2 = Var(initialize=1e7)

        # This is a standard constraint form, so should have no issues despite cancellation
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=0 == m.v1 - m.v2
        )

        assert vv == [0, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_equality_expr_zero_term(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e7)
        m.v2 = Var(initialize=1e7)
        m.v3 = Var(initialize=0)

        expr = m.v3 == m.v1 - m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        # No canceling terms as v3=0 is ignored thus we have a=b form
        assert len(cc) == 0
        assert not k
        assert not tr

        # Set v3 above zero tolerance
        m.v3.set_value(1e-4)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)
        assert vv == [-1e-4, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_equality_expr_compound(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b"], initialize=5e6)
        m.v2 = Var(initialize=1e7)
        m.v3 = Var(initialize=0)

        expr = m.v3 == sum(v for v in m.v1.values()) - m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        # No canceling terms as v3=0 is ignored thus we have a=b form
        assert len(cc) == 0
        assert not k
        assert not tr

        # Set v3 above zero tolerance
        m.v3.set_value(1e-4)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)
        assert vv == [-1e-4, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    # Double check for a+eps=c form gets flagged in some way
    @pytest.mark.unit
    def test_canceling_equality_expr_canceling_sides(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=1e-8)

        expr1 = m.v1 + m.v3
        expr = m.v2 == expr1
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, pytest.approx(1 + 1e-8, abs=1e-8)]
        assert expr1 in mm
        assert len(mm) == 1
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    # Check to make sure simple linking constraints are not flagged as canceling
    @pytest.mark.unit
    def test_linking_equality_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)

        expr = m.v1 == m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, 1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_linking_equality_expr_compound(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=1)

        expr = m.v1 == m.v2 * m.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, 1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.component
    def test_external_function_w_string_argument(self):
        m = ConcreteModel()
        m.properties = iapws95.Iapws95ParameterBlock()
        m.state = m.properties.build_state_block([0])

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.state[0].enth_mol
        )

        assert vv == [pytest.approx(1.1021387e-2, rel=1e-6)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

        # Test nested external functions
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.state[0].temperature
        )

        assert vv == [pytest.approx(270.4877, rel=1e-6)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_zero_tolerance(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e-12)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=1e-12)

        expr1 = m.v1 - m.v3
        expr2 = expr1 + 1
        expr = m.v2 == expr2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, pytest.approx(1, abs=1e-8)]
        # We expect no mismatches, as smallest terms are below zero tolerance
        assert len(mm) == 0
        # We expect no canceling terms, as v1 and v3 are ignored (zero tolerance), leaving
        # 1 == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Set v1 and v3 above zero tolerance
        m.v1.set_value(1e-8)
        m.v3.set_value(1e-8)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, pytest.approx(1, abs=1e-8)]
        assert expr2 in mm
        assert len(mm) == 1
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    # Test to make sure scaled constraints are not flagged as issues
    @pytest.mark.unit
    def test_scaled_equality_expr_product(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=2)

        expr = 0 == m.v3 * (m.v1 - m.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, 0]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_scaled_equality_expr_division(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=2)

        expr = 0 == (m.v1 - m.v2) / m.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, 0]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_mole_fraction_constraint(self):
        m = ConcreteModel()
        m.mole_frac_a = Var(initialize=0.5)
        m.flow_a = Var(initialize=100)
        m.flow_b = Var(initialize=100)

        expr = m.mole_frac_a * (m.flow_a + m.flow_b) == m.flow_a
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(-100, rel=1e-5), pytest.approx(100, rel=1e-5)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_mole_fraction_constraint_trace(self):
        m = ConcreteModel()
        m.mole_frac_a = Var(initialize=0.999810)
        m.flow_a = Var(initialize=122.746)
        m.flow_b = Var(initialize=0.0233239)

        expr = m.mole_frac_a * (m.flow_a + m.flow_b) == m.flow_a
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [
            pytest.approx(-122.746, rel=1e-5),
            pytest.approx(122.746, rel=1e-5),
        ]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr
