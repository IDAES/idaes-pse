#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for math util methods.
"""

import pytest
from pyomo.environ import ConcreteModel, Param, Var, value
from idaes.core.util.math import (
    smooth_abs,
    smooth_minmax,
    smooth_min,
    smooth_max,
    safe_sqrt,
    safe_log,
)

__author__ = "Andrew Lee"


@pytest.fixture(scope="module")
def simple_model():
    """Build a simple model for testing."""
    m = ConcreteModel()
    m.a = Var(initialize=4.0)
    m.b = Var(initialize=-4.0)
    m.e = Param(default=1e-4)

    return m


@pytest.mark.unit
def test_smooth_abs_maths():
    # Test basic smooth_abs functionalliy
    assert smooth_abs(4, 0) == 4.0
    assert smooth_abs(-4, 0) == 4.0
    assert smooth_abs(10.0, 0.0) == 10.0
    assert smooth_abs(-10.0, 0.0) == 10.0

    assert smooth_abs(2, 1e-4) == pytest.approx(2.0, abs=1e-4)
    assert smooth_abs(-2, 1e-4) == pytest.approx(2.0, abs=1e-4)
    assert smooth_abs(10) == pytest.approx(10.0, abs=1e-4)
    assert smooth_abs(-10) == pytest.approx(10.0, abs=1e-4)


@pytest.mark.unit
def test_smooth_abs_expr(simple_model):
    # Test that smooth_abs works with Pyomo components
    assert value(smooth_abs(simple_model.a, 0)) == 4.0
    assert value(smooth_abs(simple_model.b, 0)) == 4.0

    assert value(smooth_abs(simple_model.a, simple_model.e)) == pytest.approx(
        4.0, abs=1e-4
    )
    assert value(smooth_abs(simple_model.b, simple_model.e)) == pytest.approx(
        4.0, abs=1e-4
    )


@pytest.mark.unit
def test_smooth_abs_a_errors():
    # Test that smooth_abs returns meaningful errors when given invalid arg
    with pytest.raises(TypeError):
        smooth_abs("foo")
    with pytest.raises(TypeError):
        smooth_abs([1, 2, 3])


@pytest.mark.unit
def test_smooth_abs_eps_errors():
    # Test that smooth_abs returns meaningful errors when given invalid eps
    with pytest.raises(TypeError):
        smooth_abs(1.0, "a")
    with pytest.raises(TypeError):
        smooth_abs(1.0, [1, 2, 3])


@pytest.mark.unit
def test_smooth_minmax_maths():
    # Test basic smooth_minmax functionality
    assert smooth_minmax(1, 2, 0, sense="max") == 2
    assert smooth_minmax(1, 2, 0, sense="min") == 1
    assert smooth_minmax(5.0, 3, 0.0, sense="max") == 5
    assert smooth_minmax(5.0, 3, 0.0, sense="min") == 3

    assert smooth_minmax(2.0, 12.0, 1e-4, "max") == pytest.approx(12.0, abs=1e-4)
    assert smooth_minmax(2.0, 12.0, 1e-4, "min") == pytest.approx(2.0, abs=1e-4)
    assert smooth_minmax(32.0, 12.0, sense="max") == pytest.approx(32.0, abs=1e-4)
    assert smooth_minmax(32.0, 12.0, sense="min") == pytest.approx(12.0, abs=1e-4)


@pytest.mark.unit
def test_smooth_minmax_default_sense():
    # Test that smooth_minmax defaults to maximise
    assert (
        smooth_minmax(
            1,
            2,
            0,
        )
        == 2
    )


@pytest.mark.unit
def test_smooth_minmax_expr(simple_model):
    # Test that smooth_minmax works with Pyomo components
    assert value(smooth_minmax(simple_model.a, simple_model.b, 0, sense="max")) == 4.0
    assert value(smooth_minmax(simple_model.a, simple_model.b, 0, sense="min")) == -4.0

    assert value(
        smooth_minmax(simple_model.a, simple_model.b, sense="max")
    ) == pytest.approx(4.0, abs=1e-4)
    assert value(
        smooth_minmax(simple_model.a, simple_model.b, sense="min")
    ) == pytest.approx(-4.0, abs=1e-4)

    assert value(
        smooth_minmax(simple_model.a, simple_model.b, simple_model.e, sense="max")
    ) == pytest.approx(4.0, abs=1e-4)
    assert value(
        smooth_minmax(simple_model.a, simple_model.b, simple_model.e, sense="min")
    ) == pytest.approx(-4.0, abs=1e-4)


@pytest.mark.unit
def test_smooth_abs_ab_errors():
    # Test that smooth_abs returns meaningful errors when given invalid args
    with pytest.raises(TypeError):
        smooth_abs("foo", 1)
    with pytest.raises(TypeError):
        smooth_abs(3, [1, 2, 3])


@pytest.mark.unit
def test_smooth_minmax_eps_errors():
    # Test that smooth_abs returns meaningful errors when given invalid eps
    with pytest.raises(TypeError):
        smooth_minmax(1.0, 1.0, "foo")
    with pytest.raises(TypeError):
        smooth_minmax(1.0, 1.0, [1, 2, 3])


@pytest.mark.unit
def test_smooth_minmax_sense_errors():
    # Test that smooth_abs returns meaningful errors when given invalid sense
    with pytest.raises(ValueError):
        smooth_minmax(1.0, 1.0, sense="foo")
    with pytest.raises(ValueError):
        smooth_minmax(1.0, 1.0, sense=1.0)
    with pytest.raises(ValueError):
        smooth_minmax(1.0, 1.0, sense=[1.0])


@pytest.mark.unit
def test_smooth_max(simple_model):
    # Test that smooth_max gives correct values
    assert smooth_max(3.0, 12.0) == pytest.approx(12.0, abs=1e-4)
    assert value(
        smooth_max(simple_model.a, simple_model.b, simple_model.e)
    ) == pytest.approx(4.0, abs=1e-4)


@pytest.mark.unit
def test_smooth_min(simple_model):
    # Test that smooth_min gives correct values
    assert smooth_min(3.0, 12.0) == pytest.approx(3.0, abs=1e-4)
    assert value(
        smooth_min(simple_model.a, simple_model.b, simple_model.e)
    ) == pytest.approx(-4.0, abs=1e-4)


@pytest.mark.unit
def test_smooth_min(simple_model):
    # Test that smooth_min gives correct values
    assert safe_sqrt(4) == pytest.approx(2.0, abs=1e-4)
    assert safe_sqrt(0, eps=1e-6) == pytest.approx(0.0, abs=1e-3)
    assert safe_sqrt(-4) == pytest.approx(0.0, abs=1e-4)
    assert value(safe_sqrt(simple_model.a, simple_model.e)) == pytest.approx(
        2.0, abs=1e-4
    )


@pytest.mark.unit
def test_smooth_min(simple_model):
    # Test that smooth_min gives correct values
    assert safe_log(4) == pytest.approx(1.386294, abs=1e-4)
    assert safe_log(0) < -5
    assert safe_log(-4) < -5
    assert value(safe_log(simple_model.a, simple_model.e)) == pytest.approx(
        1.386294, abs=1e-4
    )
