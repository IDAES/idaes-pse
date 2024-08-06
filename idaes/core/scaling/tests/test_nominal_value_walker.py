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
"""
This module contains tests for the NominalValueExtractionVisitor.
"""
import math
import logging
import pytest

import pyomo.environ as pyo
from pyomo.contrib.pynumero.asl import AmplInterface

from idaes.core.scaling.util import NominalValueExtractionVisitor, set_scaling_factor
from idaes.models.properties.modular_properties.eos.ceos_common import (
    cubic_roots_available,
    CubicThermoExpressions,
    CubicType as CubicEoS,
)

__author__ = "Andrew Lee"


class TestNominalValueExtractionVisitor:
    @pytest.fixture(scope="class")
    def m(self):
        m = pyo.ConcreteModel()
        m.set = pyo.Set(initialize=["a", "b", "c"])

        return m

    @pytest.mark.unit
    def test_int(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=7) == [7]

    @pytest.mark.unit
    def test_float(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=7.7) == [7.7]

    @pytest.mark.unit
    def test_negative_float(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=-7.7) == [-7.7]

    @pytest.mark.unit
    def test_zero(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=0) == [0]

    @pytest.mark.unit
    def test_true(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=True) == [1]

    @pytest.mark.unit
    def test_false(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=False) == [0]

    @pytest.mark.unit
    def test_scalar_param_no_scale(self, m):
        m.scalar_param = pyo.Param(initialize=1, mutable=True)
        assert NominalValueExtractionVisitor().walk_expression(expr=m.scalar_param) == [
            1
        ]

    @pytest.mark.unit
    def test_scalar_param_w_scale(self, m):
        m.scalar_param = pyo.Param(default=12, mutable=True)
        set_scaling_factor(m.scalar_param, 1 / 10)
        assert NominalValueExtractionVisitor().walk_expression(expr=m.scalar_param) == [
            12
        ]

    @pytest.mark.unit
    def test_indexed_param_w_scale(self, m):
        m.indexed_param = pyo.Param(m.set, initialize=1, mutable=True)
        set_scaling_factor(m.indexed_param["a"], 1 / 13)
        set_scaling_factor(m.indexed_param["b"], 1 / 14)
        set_scaling_factor(m.indexed_param["c"], 1 / 15)

        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_param["a"]
        ) == [1]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_param["b"]
        ) == [1]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_param["c"]
        ) == [1]

    @pytest.mark.unit
    def test_scalar_var_no_scale(self, m):
        m.scalar_var = pyo.Var(initialize=10)
        # Should use current value
        assert NominalValueExtractionVisitor().walk_expression(expr=m.scalar_var) == [
            10
        ]

    @pytest.mark.unit
    def test_scalar_var_w_scale(self, m):
        set_scaling_factor(m.scalar_var, 1 / 21)
        assert NominalValueExtractionVisitor().walk_expression(expr=m.scalar_var) == [
            21
        ]

    @pytest.mark.unit
    def test_var_neg_bounds(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(bounds=(-1000, 0))
        set_scaling_factor(m.var, 1 / 4)

        # Expect nominal value to be negative
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [-4]

    @pytest.mark.unit
    def test_var_neg_upper_bound(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(bounds=(None, -2000))
        set_scaling_factor(m.var, 1 / 4)

        # Expect nominal value to be negative
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [-4]

    @pytest.mark.unit
    def test_var_neg_domain(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(domain=pyo.NegativeReals)

        set_scaling_factor(m.var, 1 / 4)
        # Expect nominal value to be negative
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [-4]

    @pytest.mark.unit
    def test_var_neg_value(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(initialize=-1)
        set_scaling_factor(m.var, 1 / 4)

        # Expect nominal value to be negative
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [-4]

    @pytest.mark.unit
    def test_var_fixed_value(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(initialize=-1)
        m.var.fix()
        set_scaling_factor(m.var, 1 / 4)

        # Nominal value should be value
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [-1]

    @pytest.mark.unit
    def test_var_pos_bounds(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(bounds=(0, 1000))
        set_scaling_factor(m.var, 1 / 4)

        # Expect nominal value to be positive
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [4]

    @pytest.mark.unit
    def test_var_pos_lower_bound(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(bounds=(1000, None))
        set_scaling_factor(m.var, 1 / 4)

        # Expect nominal value to be positive
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [4]

    @pytest.mark.unit
    def test_var_pos_domain(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(domain=pyo.PositiveReals)
        set_scaling_factor(m.var, 1 / 4)

        # Expect nominal value to be positive
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [4]

    @pytest.mark.unit
    def test_var_pos_value(self):
        m = pyo.ConcreteModel()
        m.var = pyo.Var(initialize=1)
        set_scaling_factor(m.var, 1 / 4)

        # Expect nominal value to be positive
        assert NominalValueExtractionVisitor().walk_expression(expr=m.var) == [4]

    @pytest.mark.unit
    def test_indexed_var_no_scale(self, m):
        m.indexed_var = pyo.Var(m.set, initialize=1)
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["a"]
        ) == [1]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["b"]
        ) == [1]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["c"]
        ) == [1]

    @pytest.mark.unit
    def test_indexed_var_w_scale(self, m):
        set_scaling_factor(m.indexed_var["a"], 1 / 22)
        set_scaling_factor(m.indexed_var["b"], 1 / 23)
        set_scaling_factor(m.indexed_var["c"], 1 / 24)

        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["a"]
        ) == [22]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["b"]
        ) == [23]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["c"]
        ) == [24]

    @pytest.mark.unit
    def test_indexed_var_w_scale_partial_fixed(self, m):
        m.indexed_var["a"].fix(20)
        set_scaling_factor(m.indexed_var["a"], 1 / 22)
        set_scaling_factor(m.indexed_var["b"], 1 / 23)
        set_scaling_factor(m.indexed_var["c"], 1 / 24)

        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["a"]
        ) == [20]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["b"]
        ) == [23]
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.indexed_var["c"]
        ) == [24]

        # Clean up for future tests
        m.indexed_var["a"].unfix()

    @pytest.mark.unit
    def test_equality_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.scalar_var == m.indexed_var["a"]
        ) == [21, 22]

    @pytest.mark.unit
    def test_inequality_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.scalar_var <= m.indexed_var["a"]
        ) == [21, 22]

    @pytest.mark.unit
    def test_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=sum(m.indexed_var[i] for i in m.set)
        ) == [22, 23, 24]

    @pytest.mark.unit
    def test_additive_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.scalar_var + m.indexed_var["a"] + m.scalar_param
        ) == [21, 22, 12]

    @pytest.mark.unit
    def test_additive_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.scalar_var + m.indexed_var["a"] - m.scalar_param
        ) == [21, 22, -12]

    @pytest.mark.unit
    def test_product_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.scalar_var * m.indexed_var["a"] * m.scalar_param
        ) == [21 * 22 * 12]

    @pytest.mark.unit
    def test_product_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=(m.scalar_var + m.indexed_var["a"])
            * (m.scalar_param + m.indexed_var["b"])
        ) == [21 * 12, 21 * 23, 22 * 12, 22 * 23]

    @pytest.mark.unit
    def test_product_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=(m.scalar_var - m.indexed_var["a"])
            * (m.scalar_param - m.indexed_var["b"])
        ) == [21 * 12, -21 * 23, -22 * 12, 22 * 23]

    @pytest.mark.unit
    def test_division_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.scalar_var / m.indexed_var["a"] / m.scalar_param
        ) == [21 / 22 / 12]

    @pytest.mark.unit
    def test_division_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=(m.scalar_var + m.indexed_var["a"])
            / (m.scalar_param + m.indexed_var["b"])
        ) == [(21 + 22) / (12 + 23)]

    @pytest.mark.unit
    def test_division_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=(m.scalar_var - m.indexed_var["a"])
            / (m.scalar_param - m.indexed_var["b"])
        ) == [(21 - 22) / (12 - 23)]

    @pytest.mark.unit
    def test_division_expr_error(self, m, caplog):
        caplog.set_level(logging.DEBUG, logger="idaes.core.scaling")
        assert NominalValueExtractionVisitor().walk_expression(
            expr=1 / (m.scalar_var - 21)
        ) == [1]

        expected = "Nominal value of 0 found in denominator of division expression. "
        "Assigning a value of 1. You should check you scaling factors and models to "
        "ensure there are no values of 0 that can appear in these functions."

        assert expected in caplog.text

    @pytest.mark.unit
    def test_pow_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.scalar_var ** m.indexed_var["a"]
        ) == pytest.approx([21**22], rel=1e-12)

    @pytest.mark.unit
    def test_pow_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=(m.scalar_var + m.indexed_var["a"])
            ** (m.scalar_param + m.indexed_var["b"])
        ) == [
            pytest.approx((21 + 22) ** (12 + 23), rel=1e-12),
        ]

    @pytest.mark.unit
    def test_pow_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=(m.scalar_var - m.indexed_var["a"])
            ** (m.scalar_param - m.indexed_var["b"])
        ) == [
            pytest.approx(abs(21 - 22) ** (12 - 23), rel=1e-12),
        ]

    @pytest.mark.unit
    def test_negation_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=-m.scalar_var) == [
            -21
        ]

    @pytest.mark.unit
    def test_negation_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=-(m.scalar_var + m.indexed_var["a"])
        ) == [-21, -22]

    @pytest.mark.unit
    def test_log_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.log(m.scalar_var)
        ) == [pytest.approx(math.log(21), rel=1e-12)]

    @pytest.mark.unit
    def test_log_expr_error(self, m):
        with pytest.raises(
            ValueError,
            match="Evaluation error occurred when getting nominal value in log expression "
            "with input 0.0. You should check you scaling factors and model to "
            "address any numerical issues or scale this constraint manually.",
        ):
            assert NominalValueExtractionVisitor().walk_expression(
                expr=pyo.log(m.scalar_var - 21)
            )

    @pytest.mark.unit
    def test_log_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.log(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.log(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_log_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.log(-m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.log(-21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_log10_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.log10(m.scalar_var)
        ) == [pytest.approx(math.log10(21), rel=1e-12)]

    @pytest.mark.unit
    def test_log10_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.log10(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.log10(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_log10_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.log10(-m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.log10(-21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_log10_expr_error(self, m):
        with pytest.raises(
            ValueError,
            match="Evaluation error occurred when getting nominal value in log10 expression "
            "with input 0.0. You should check you scaling factors and model to "
            "address any numerical issues or scale this constraint manually.",
        ):
            assert NominalValueExtractionVisitor().walk_expression(
                expr=pyo.log10(m.scalar_var - 21)
            )

    @pytest.mark.unit
    def test_sqrt_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sqrt(m.scalar_var)
        ) == [pytest.approx(21**0.5, rel=1e-12)]

    @pytest.mark.unit
    def test_sqrt_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sqrt(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx((21 + 22) ** 0.5, rel=1e-12)]

    @pytest.mark.unit
    def test_sqrt_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sqrt(-m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx((-21 + 22) ** 0.5, rel=1e-12)]

    @pytest.mark.unit
    def test_sqrt_expr_error(self, m):
        with pytest.raises(
            ValueError,
            match="Evaluation error occurred when getting nominal value in sqrt expression "
            "with input -21.0. You should check you scaling factors and model to "
            "address any numerical issues or scale this constraint manually.",
        ):
            assert NominalValueExtractionVisitor().walk_expression(
                expr=pyo.sqrt(-m.scalar_var)
            )

    @pytest.mark.unit
    def test_sin_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sin(m.scalar_var)
        ) == [pytest.approx(math.sin(21), rel=1e-12)]

    @pytest.mark.unit
    def test_sin_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sin(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.sin(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_sin_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sin(m.scalar_var - m.indexed_var["a"])
        ) == [pytest.approx(math.sin(21 - 22), rel=1e-12)]

    @pytest.mark.unit
    def test_cos_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.cos(m.scalar_var)
        ) == [pytest.approx(math.cos(21), rel=1e-12)]

    @pytest.mark.unit
    def test_cos_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.cos(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.cos(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_cos_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.cos(m.scalar_var - m.indexed_var["a"])
        ) == [pytest.approx(math.cos(21 - 22), rel=1e-12)]

    @pytest.mark.unit
    def test_tan_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.tan(m.scalar_var)
        ) == [pytest.approx(math.tan(21), rel=1e-12)]

    @pytest.mark.unit
    def test_tan_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.tan(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.tan(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_tan_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.tan(m.scalar_var - m.indexed_var["a"])
        ) == [pytest.approx(math.tan(21 - 22), rel=1e-12)]

    @pytest.mark.unit
    def test_sinh_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sinh(m.scalar_var)
        ) == [pytest.approx(math.sinh(21), rel=1e-12)]

    @pytest.mark.unit
    def test_sinh_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sinh(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.sinh(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_sinh_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.sinh(m.scalar_var - m.indexed_var["a"])
        ) == [pytest.approx(math.sinh(21 - 22), rel=1e-12)]

    @pytest.mark.unit
    def test_cosh_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.cosh(m.scalar_var)
        ) == [pytest.approx(math.cosh(21), rel=1e-12)]

    @pytest.mark.unit
    def test_cosh_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.cosh(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.cosh(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_cosh_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.cosh(m.scalar_var - m.indexed_var["a"])
        ) == [pytest.approx(math.cosh(21 - 22), rel=1e-12)]

    @pytest.mark.unit
    def test_tanh_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.tanh(m.scalar_var)
        ) == [pytest.approx(math.tanh(21), rel=1e-12)]

    @pytest.mark.unit
    def test_tanh_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.tanh(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.tanh(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_tanh_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.tanh(m.scalar_var - m.indexed_var["a"])
        ) == [pytest.approx(math.tanh(21 - 22), rel=1e-12)]

    @pytest.mark.unit
    def test_asin_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(expr=pyo.asin(1)) == [
            pytest.approx(math.asin(1), rel=1e-12)
        ]

    @pytest.mark.unit
    def test_asin_sum_expr(self, m):
        m.scalar_param.set_value(0.5)
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.asin(0.5 + m.scalar_param)
        ) == [pytest.approx(math.asin(1), rel=1e-12)]

    @pytest.mark.unit
    def test_asin_sum_expr_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.asin(0.5 - m.scalar_param)
        ) == [pytest.approx(math.asin(0), rel=1e-12)]

    @pytest.mark.unit
    def test_acos_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.acos(m.scalar_param)
        ) == [pytest.approx(math.acos(0.5), rel=1e-12)]

    @pytest.mark.unit
    def test_acos_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.acos(0.5 + m.scalar_param)
        ) == [pytest.approx(math.acos(1), rel=1e-12)]

    @pytest.mark.unit
    def test_acos_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.acos(0.5 - m.scalar_param)
        ) == [pytest.approx(math.acos(0), rel=1e-12)]

    @pytest.mark.unit
    def test_asinh_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.asinh(m.scalar_var)
        ) == [pytest.approx(math.asinh(21), rel=1e-12)]

    @pytest.mark.unit
    def test_asinh_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.asinh(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.asinh(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_asinh_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.asinh(m.scalar_var - m.indexed_var["a"])
        ) == [pytest.approx(math.asinh(21 - 22), rel=1e-12)]

    @pytest.mark.unit
    def test_acosh_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.acosh(m.scalar_var)
        ) == [pytest.approx(math.acosh(21), rel=1e-12)]

    @pytest.mark.unit
    def test_acosh_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.acosh(m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.acosh(21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_acosh_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.acosh(-m.scalar_var + m.indexed_var["a"])
        ) == [pytest.approx(math.acosh(-21 + 22), rel=1e-12)]

    @pytest.mark.unit
    def test_atanh_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.atanh(m.scalar_param)
        ) == [pytest.approx(math.atanh(0.5), rel=1e-12)]

    @pytest.mark.unit
    def test_atanh_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.atanh(0.4 + m.scalar_param)
        ) == [pytest.approx(math.atanh(0.9), rel=1e-12)]

    @pytest.mark.unit
    def test_atanh_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.atanh(0.4 - m.scalar_param)
        ) == [pytest.approx(math.atanh(-0.1), rel=1e-12)]

    @pytest.mark.unit
    def test_exp_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.exp(m.scalar_param)
        ) == [pytest.approx(math.exp(0.5), rel=1e-12)]

    @pytest.mark.unit
    def test_exp_sum_expr(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.exp(0.4 + m.scalar_param)
        ) == [pytest.approx(math.exp(0.9), rel=1e-12)]

    @pytest.mark.unit
    def test_exp_sum_expr_w_negation(self, m):
        assert NominalValueExtractionVisitor().walk_expression(
            expr=pyo.exp(-0.4 + m.scalar_param)
        ) == [pytest.approx(math.exp(0.1), rel=1e-12)]

    @pytest.mark.unit
    def test_expr_if(self, m):
        m.exprif = pyo.Expr_if(
            IF=m.scalar_param,
            THEN=m.indexed_var["a"],
            ELSE=m.indexed_var["b"] + m.indexed_var["c"],
        )

        assert NominalValueExtractionVisitor().walk_expression(expr=m.exprif) == [
            22,
            23,
            24,
        ]

    @pytest.mark.unit
    def test_expr_if_w_negation(self, m):
        m.exprif = pyo.Expr_if(
            IF=m.scalar_param,
            THEN=m.indexed_var["a"],
            ELSE=m.indexed_var["b"] - m.indexed_var["c"],
        )

        assert NominalValueExtractionVisitor().walk_expression(expr=m.exprif) == [
            22,
            23,
            -24,
        ]

    @pytest.mark.unit
    @pytest.mark.skipif(
        not AmplInterface.available(), reason="pynumero_ASL is not available"
    )
    @pytest.mark.skipif(not cubic_roots_available, reason="Cubic roots not available")
    def test_ext_func(self):
        # Use the cubic root external function to test
        m = pyo.ConcreteModel()
        m.a = pyo.Var(initialize=1)
        m.b = pyo.Var(initialize=1)

        set_scaling_factor(m.a, 1 / 2)
        set_scaling_factor(m.b, 1 / 4)

        m.expr_write = CubicThermoExpressions(m)
        Z = m.expr_write.z_liq(eos=CubicEoS.PR, A=m.a, B=m.b)

        expected_mag = -9.489811292072448
        assert NominalValueExtractionVisitor().walk_expression(expr=Z) == [
            pytest.approx(expected_mag, rel=1e-8)
        ]

        # Check that model state did not change
        assert pyo.value(m.a) == 1
        assert pyo.value(m.b) == 1
        assert pyo.value(Z) == pytest.approx(-2.1149075414767577, rel=1e-8)

        # Now, change the actual state to the expected magnitudes and confirm result
        m.a.set_value(2)
        m.b.set_value(4)
        assert pyo.value(Z) == pytest.approx(expected_mag, rel=1e-8)

    @pytest.mark.unit
    def test_Expression(self, m):
        m.expression = pyo.Expression(
            expr=m.scalar_param ** (sum(m.indexed_var[i] for i in m.set))
        )

        assert NominalValueExtractionVisitor().walk_expression(expr=m.expression) == [
            0.5 ** (22 + 23 + 24)
        ]

    @pytest.mark.unit
    def test_constraint(self, m):
        m.constraint = pyo.Constraint(expr=m.scalar_var == m.expression)

        assert NominalValueExtractionVisitor().walk_expression(
            expr=m.constraint.expr
        ) == [21, 0.5 ** (22 + 23 + 24)]
