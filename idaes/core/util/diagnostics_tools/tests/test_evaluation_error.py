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
This module contains tests for the evaluation error checker.
"""

import pytest
from math import inf

from pyomo.environ import (
    ConcreteModel,
    Var,
    Binary,
    Integers,
    log,
    log10,
    tan,
    sin,
    cos,
    asin,
    acos,
    sqrt,
    exp,
)
from pyomo.common.config import ConfigDict, ConfigValue

from idaes.core.util.diagnostics_tools.evaluation_error import (
    _get_bounds_with_inf,
    _check_eval_error_division,
    _check_eval_error_pow,
    _check_eval_error_log,
    _check_eval_error_tan,
    _check_eval_error_asin,
    _check_eval_error_acos,
    _check_eval_error_sqrt,
    _check_eval_error_unary,
    EvalErrorWalker,
)

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee"


@pytest.fixture
def config():
    """Create a basic config for testing"""
    cfg = ConfigDict()
    cfg.declare(
        "warn_for_evaluation_error_at_bounds",
        ConfigValue(default=False, domain=bool),
    )
    return cfg


@pytest.fixture
def config_with_bounds_warning():
    """Create a config that warns for errors at bounds"""
    cfg = ConfigDict()
    cfg.declare(
        "warn_for_evaluation_error_at_bounds",
        ConfigValue(default=True, domain=bool),
    )
    return cfg


class TestGetBoundsWithInf:
    """Test the _get_bounds_with_inf function"""

    @pytest.mark.unit
    def test_bounded_variable(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        lb, ub = _get_bounds_with_inf(m.x)
        assert lb == 1
        assert ub == 10

    @pytest.mark.unit
    def test_unbounded_variable(self):
        m = ConcreteModel()
        m.x = Var()
        lb, ub = _get_bounds_with_inf(m.x)
        assert lb == -inf
        assert ub == inf

    @pytest.mark.unit
    def test_lower_bounded_only(self):
        m = ConcreteModel()
        m.x = Var(bounds=(0, None))
        lb, ub = _get_bounds_with_inf(m.x)
        assert lb == 0
        assert ub == inf

    @pytest.mark.unit
    def test_upper_bounded_only(self):
        m = ConcreteModel()
        m.x = Var(bounds=(None, 100))
        lb, ub = _get_bounds_with_inf(m.x)
        assert lb == -inf
        assert ub == 100

    @pytest.mark.unit
    def test_expression_with_bounds(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 5))
        m.y = Var(bounds=(2, 3))
        expr = m.x + m.y
        lb, ub = _get_bounds_with_inf(expr)
        assert lb == 3
        assert ub == 8


class TestCheckEvalErrorDivision:
    """Test the _check_eval_error_division function"""

    @pytest.mark.unit
    def test_safe_division(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(2, 5))
        expr = m.x / m.y
        warn_list = []
        _check_eval_error_division(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_division_by_zero_possible(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(-1, 1))
        expr = m.x / m.y
        warn_list = []
        _check_eval_error_division(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential division by 0" in warn_list[0]

    @pytest.mark.unit
    def test_division_by_zero_at_bounds_no_warning(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(0, 5))
        expr = m.x / m.y
        warn_list = []
        _check_eval_error_division(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_division_by_zero_at_bounds_with_warning(self, config_with_bounds_warning):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(0, 5))
        expr = m.x / m.y
        warn_list = []
        _check_eval_error_division(expr, warn_list, config_with_bounds_warning)
        assert len(warn_list) == 1
        assert "Potential division by 0" in warn_list[0]

    @pytest.mark.unit
    def test_division_negative_denominator_range(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(-5, -1))
        expr = m.x / m.y
        warn_list = []
        _check_eval_error_division(expr, warn_list, config)
        assert len(warn_list) == 0


class TestCheckEvalErrorPow:
    """Test the _check_eval_error_pow function"""

    @pytest.mark.unit
    def test_safe_power_positive_base(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(0.5, 2))
        expr = m.x**m.y
        warn_list = []
        _check_eval_error_pow(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_negative_base_fractional_exponent(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-5, -1))
        m.y = Var(bounds=(0.5, 2))
        expr = m.x**m.y
        warn_list = []
        _check_eval_error_pow(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential evaluation error" in warn_list[0]

    @pytest.mark.unit
    def test_integer_exponent_zero_base(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 5))
        m.n = Var(domain=Integers, bounds=(2, 4))
        expr = m.x**m.n
        warn_list = []
        _check_eval_error_pow(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_binary_exponent(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 5))
        m.b = Var(domain=Binary)
        expr = m.x**m.b
        warn_list = []
        _check_eval_error_pow(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_fixed_integer_exponent(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 5))
        m.x.fix(2.0)
        m.y = Var(bounds=(0, 10))
        expr = m.y**2
        warn_list = []
        _check_eval_error_pow(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_negative_base_integer_exponent(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-5, -1))
        m.n = Var(domain=Integers, bounds=(2, 4))
        expr = m.x**m.n
        warn_list = []
        _check_eval_error_pow(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_nonnegative_base_nonnegative_exponent(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 5))
        m.y = Var(bounds=(0, 2))
        expr = m.x**m.y
        warn_list = []
        _check_eval_error_pow(expr, warn_list, config)
        assert len(warn_list) == 0


class TestCheckEvalErrorLog:
    """Test the _check_eval_error_log function"""

    @pytest.mark.unit
    def test_safe_log(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        expr = log(m.x)
        warn_list = []
        _check_eval_error_log(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_log_of_negative(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-5, -1))
        expr = log(m.x)
        warn_list = []
        _check_eval_error_log(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential log of a non-positive number" in warn_list[0]

    @pytest.mark.unit
    def test_log_of_zero_at_bound_no_warning(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 10))
        expr = log(m.x)
        warn_list = []
        _check_eval_error_log(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_log_of_zero_at_bound_with_warning(self, config_with_bounds_warning):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 10))
        expr = log(m.x)
        warn_list = []
        _check_eval_error_log(expr, warn_list, config_with_bounds_warning)
        assert len(warn_list) == 1
        assert "Potential log of a non-positive number" in warn_list[0]

    @pytest.mark.unit
    def test_log_crossing_zero(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-1, 1))
        expr = log(m.x)
        warn_list = []
        _check_eval_error_log(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential log of a non-positive number" in warn_list[0]


class TestCheckEvalErrorTan:
    """Test the _check_eval_error_tan function"""

    @pytest.mark.unit
    def test_safe_tan(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 1))
        expr = tan(m.x)
        warn_list = []
        _check_eval_error_tan(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_tan_unbounded(self, config):
        m = ConcreteModel()
        m.x = Var()
        expr = tan(m.x)
        warn_list = []
        _check_eval_error_tan(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "may evaluate to -inf or inf" in warn_list[0]


class TestCheckEvalErrorAsin:
    """Test the _check_eval_error_asin function"""

    @pytest.mark.unit
    def test_safe_asin(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-1, 1))
        expr = asin(m.x)
        warn_list = []
        _check_eval_error_asin(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_asin_outside_domain(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-2, 2))
        expr = asin(m.x)
        warn_list = []
        _check_eval_error_asin(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential evaluation of asin outside [-1, 1]" in warn_list[0]

    @pytest.mark.unit
    def test_asin_above_domain(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 2))
        expr = asin(m.x)
        warn_list = []
        _check_eval_error_asin(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential evaluation of asin outside [-1, 1]" in warn_list[0]

    @pytest.mark.unit
    def test_asin_below_domain(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-2, 0))
        expr = asin(m.x)
        warn_list = []
        _check_eval_error_asin(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential evaluation of asin outside [-1, 1]" in warn_list[0]


class TestCheckEvalErrorAcos:
    """Test the _check_eval_error_acos function"""

    @pytest.mark.unit
    def test_safe_acos(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-1, 1))
        expr = acos(m.x)
        warn_list = []
        _check_eval_error_acos(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_acos_outside_domain(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-2, 2))
        expr = acos(m.x)
        warn_list = []
        _check_eval_error_acos(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential evaluation of acos outside [-1, 1]" in warn_list[0]

    @pytest.mark.unit
    def test_acos_above_domain(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 1.5))
        expr = acos(m.x)
        warn_list = []
        _check_eval_error_acos(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential evaluation of acos outside [-1, 1]" in warn_list[0]


class TestCheckEvalErrorSqrt:
    """Test the _check_eval_error_sqrt function"""

    @pytest.mark.unit
    def test_safe_sqrt(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 10))
        expr = sqrt(m.x)
        warn_list = []
        _check_eval_error_sqrt(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_sqrt_of_negative(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-5, -1))
        expr = sqrt(m.x)
        warn_list = []
        _check_eval_error_sqrt(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential square root of a negative number" in warn_list[0]

    @pytest.mark.unit
    def test_sqrt_crossing_zero(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-1, 1))
        expr = sqrt(m.x)
        warn_list = []
        _check_eval_error_sqrt(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential square root of a negative number" in warn_list[0]

    @pytest.mark.unit
    def test_sqrt_strictly_positive(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        expr = sqrt(m.x)
        warn_list = []
        _check_eval_error_sqrt(expr, warn_list, config)
        assert len(warn_list) == 0


class TestCheckEvalErrorUnary:
    """Test the _check_eval_error_unary function"""

    @pytest.mark.unit
    def test_log_function(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-1, 1))
        expr = log(m.x)
        warn_list = []
        _check_eval_error_unary(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential log of a non-positive number" in warn_list[0]

    @pytest.mark.unit
    def test_log10_function(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-1, 1))
        expr = log10(m.x)
        warn_list = []
        _check_eval_error_unary(expr, warn_list, config)
        assert len(warn_list) == 1
        assert "Potential log of a non-positive number" in warn_list[0]

    @pytest.mark.unit
    def test_safe_unary_function(self, config):
        """Test a unary function that doesn't have special handling"""
        m = ConcreteModel()
        m.x = Var(bounds=(0, 1))
        expr = sin(m.x)
        warn_list = []
        _check_eval_error_unary(expr, warn_list, config)
        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_exp_function(self, config):
        """Test exp which should not generate warnings"""
        m = ConcreteModel()
        m.x = Var(bounds=(-10, 10))
        expr = exp(m.x)
        warn_list = []
        _check_eval_error_unary(expr, warn_list, config)
        assert len(warn_list) == 0


class TestEvalErrorWalker:
    """Test the EvalErrorWalker class"""

    @pytest.mark.unit
    def test_walker_division_error(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(-1, 1))
        expr = m.x / m.y

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 1
        assert "Potential division by 0" in warn_list[0]

    @pytest.mark.unit
    def test_walker_power_error(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-5, -1))
        m.y = Var(bounds=(0.5, 2))
        expr = m.x**m.y

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 1
        assert "Potential evaluation error" in warn_list[0]

    @pytest.mark.unit
    def test_walker_log_error(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-5, -1))
        expr = log(m.x)

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 1
        assert "Potential log of a non-positive number" in warn_list[0]

    @pytest.mark.unit
    def test_walker_sqrt_error(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-5, -1))
        expr = sqrt(m.x)

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 1
        assert "Potential square root of a negative number" in warn_list[0]

    @pytest.mark.unit
    def test_walker_asin_error(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-2, 2))
        expr = asin(m.x)

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 1
        assert "Potential evaluation of asin outside [-1, 1]" in warn_list[0]

    @pytest.mark.unit
    def test_walker_acos_error(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-2, 2))
        expr = acos(m.x)

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 1
        assert "Potential evaluation of acos outside [-1, 1]" in warn_list[0]

    @pytest.mark.unit
    def test_walker_multiple_errors(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(-1, 1))
        m.y = Var(bounds=(-1, 1))
        m.z = Var(bounds=(1, 10))
        # Expression with multiple potential errors
        expr = log(m.x) + m.z / m.y

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 2
        assert any("Potential log" in w for w in warn_list)
        assert any("Potential division by 0" in w for w in warn_list)

    @pytest.mark.unit
    def test_walker_no_errors(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(2, 5))
        expr = m.x + m.y * 2 + exp(m.x)

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 0

    @pytest.mark.unit
    def test_walker_complex_expression(self, config):
        m = ConcreteModel()
        m.x = Var(bounds=(0, 10))
        m.y = Var(bounds=(-1, 1))
        m.z = Var(bounds=(1, 5))
        # Complex expression with nested operations
        expr = sqrt(m.x) / m.y + log(m.z) + m.x**2

        walker = EvalErrorWalker(config)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        # Should find division by zero error
        assert len(warn_list) == 1
        assert "Potential division by 0" in warn_list[0]

    @pytest.mark.unit
    def test_walker_with_bounds_warning_config(self, config_with_bounds_warning):
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(0, 5))
        expr = m.x / m.y

        walker = EvalErrorWalker(config_with_bounds_warning)
        walker.walk_expression(expr)
        warn_list = walker._warn_list

        assert len(warn_list) == 1
        assert "Potential division by 0" in warn_list[0]

    @pytest.mark.unit
    def test_walker_exit_node_return(self, config):
        """Test that exitNode returns the warning list"""
        m = ConcreteModel()
        m.x = Var(bounds=(1, 10))
        m.y = Var(bounds=(-1, 1))
        expr = m.x / m.y

        walker = EvalErrorWalker(config)
        result = walker.walk_expression(expr)

        # The walk_expression should return the final result from exitNode
        assert isinstance(result, list)
        assert len(result) == 1
        assert "Potential division by 0" in result[0]
