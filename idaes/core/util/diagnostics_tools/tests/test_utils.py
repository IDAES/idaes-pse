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
This module contains test for the diagnostics utility functions.
"""

import re

import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Set,
    Var,
)
from pyomo.common.collections import ComponentSet

from idaes.core.util.diagnostics_tools.utils import (
    var_in_block,
    vars_fixed_to_zero,
    vars_near_zero,
    vars_violating_bounds,
    vars_with_none_value,
    vars_with_extreme_values,
    check_parallel_jacobian,
    extreme_jacobian_rows,
    extreme_jacobian_columns,
    extreme_jacobian_entries,
)
from idaes.core.scaling import set_scaling_factor
from idaes.core.scaling.util import get_jacobian


@pytest.fixture
def model():
    m = ConcreteModel()
    m.b = Block()

    m.v1 = Var()
    m.v2 = Var()
    m.v3 = Var()
    m.v4 = Var()

    m.v1.fix(0)
    m.v2.fix(3)
    m.v3.set_value(0)

    return m


@pytest.mark.unit
def test_var_in_block(model):
    assert var_in_block(model.v1, model)
    assert not var_in_block(model.v1, model.b)


@pytest.mark.unit
def test_vars_fixed_to_zero(model):
    zero_vars = vars_fixed_to_zero(model)
    assert isinstance(zero_vars, ComponentSet)
    assert len(zero_vars) == 1
    for i in zero_vars:
        assert i is model.v1


@pytest.mark.unit
def test_vars_near_zero(model):
    model.v3.set_value(1e-5)

    near_zero_vars = vars_near_zero(model, variable_zero_value_tolerance=1e-5)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 2
    for i in near_zero_vars:
        assert i.local_name in ["v1", "v3"]

    near_zero_vars = vars_near_zero(model, variable_zero_value_tolerance=1e-6)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 1
    for i in near_zero_vars:
        assert i is model.v1

    set_scaling_factor(model.v3, 1e5)
    near_zero_vars = vars_near_zero(model, variable_zero_value_tolerance=1e-5)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 1
    for i in near_zero_vars:
        assert i is model.v1

    near_zero_vars = vars_near_zero(model, variable_zero_value_tolerance=1)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 2
    for i in near_zero_vars:
        assert i.local_name in ["v1", "v3"]


@pytest.mark.unit
def test_vars_with_none_value(model):
    none_value = vars_with_none_value(model)

    assert isinstance(none_value, ComponentSet)
    assert len(none_value) == 1
    for i in none_value:
        assert i is model.v4


@pytest.mark.unit
def test_vars_with_bounds_issues(model):
    model.v1.setlb(2)
    model.v1.setub(6)
    model.v2.setlb(0)
    model.v2.setub(10)
    model.v4.set_value(10)
    model.v4.setlb(0)
    model.v4.setub(1)

    bounds_issue = vars_violating_bounds(model, tolerance=0)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 2
    for i in bounds_issue:
        assert i.local_name in ["v1", "v4"]

    m = ConcreteModel()
    m.v = Var(initialize=-1e-4, bounds=(0, 1))

    # Just outside lower bound
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1

    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1

    # Just inside lower bound
    m.v.set_value(1e-4)

    set_scaling_factor(m.v, 1, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    # Just inside upper bound
    m.v.set_value(1 - 1e-4)

    set_scaling_factor(m.v, 1, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    # Just outside upper bound
    m.v.set_value(1 + 1e-4)

    set_scaling_factor(m.v, 1, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1

    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1


@pytest.mark.unit
def test_vars_with_extreme_values():
    m = ConcreteModel()
    m.v1 = Var(initialize=1e-12)  # below zero
    m.v2 = Var(initialize=1e-8)  # small
    m.v3 = Var(initialize=1e-4)
    m.v4 = Var(initialize=1e0)
    m.v5 = Var(initialize=1e4)
    m.v6 = Var(initialize=1e8)
    m.v7 = Var(initialize=1e12)  # large

    xvars = vars_with_extreme_values(m, large=1e9, small=1e-7, zero=1e-10)

    assert len(xvars) == 2
    for i in xvars:
        assert i.name in ["v2", "v7"]

    set_scaling_factor(m.v1, 1e3)
    set_scaling_factor(m.v2, 1e6)
    set_scaling_factor(m.v4, 1e-12)
    set_scaling_factor(m.v6, 1e3)
    set_scaling_factor(m.v7, 1e-12)

    xvars = vars_with_extreme_values(m, large=1e9, small=1e-7, zero=1e-10)

    assert len(xvars) == 2
    for i in xvars:
        assert i.name in ["v1", "v6"]


class TestExtremeJacobianMethods:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.I = Set(initialize=[i for i in range(5)])

        m.x = Var(m.I, initialize=1.0)

        diag = [1e7, 1, 10, 0.1, 1e-7]
        out = [1, 1, 1, 1, 1]

        @m.Constraint(m.I)
        def dummy_eqn(b, i):
            if i == 0:
                # Off-diagonal element so that extreme
                # rows and extreme columns are different
                rhs = b.x[4]
            else:
                rhs = 0
            return out[i] == diag[i] * b.x[i] + rhs

        return m

    @pytest.mark.unit
    def test_extreme_jacobian_rows(self, model):
        m = model

        def assert_unscaled_jacobian_correct(m, scaled=False):
            jac, nlp = get_jacobian(m, include_scaling_factors=scaled)
            out = extreme_jacobian_rows(jac, nlp)
            assert type(out) == list
            assert len(out) == 2
            assert out[0][0] == pytest.approx(1e7)
            assert out[0][1] is m.dummy_eqn[0]
            assert out[1][0] == pytest.approx(1e-7)
            assert out[1][1] is m.dummy_eqn[4]

        # No scaling
        assert_unscaled_jacobian_correct(m, scaled=False)
        assert_unscaled_jacobian_correct(m, scaled=True)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_rows(jac, nlp, large=1e8)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e-7)
        assert out[0][1] is m.dummy_eqn[4]

        out = extreme_jacobian_rows(jac, nlp, small=1e-8)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e7)
        assert out[0][1] is m.dummy_eqn[0]

        out = extreme_jacobian_rows(jac, nlp, large=1e8, small=1e-8)
        assert len(out) == 0

        # Constraint scaling
        set_scaling_factor(m.dummy_eqn[0], 1e-7)
        set_scaling_factor(m.dummy_eqn[4], 1e7)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_rows(jac, nlp)
        assert len(out) == 0

        # Variable scaling
        set_scaling_factor(m.x[1], 1e7)
        set_scaling_factor(m.x[2], 1e-7)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_rows(jac, nlp)
        assert len(out) == 2
        assert out[0][0] == pytest.approx(1e-7)
        assert out[0][1] is m.dummy_eqn[1]
        assert out[1][0] == pytest.approx(1e8)
        assert out[1][1] is m.dummy_eqn[2]

        # More constraint scaling
        set_scaling_factor(m.dummy_eqn[1], 1e7)
        set_scaling_factor(m.dummy_eqn[2], 1e-8)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_rows(jac, nlp)
        assert len(out) == 0

    @pytest.mark.unit
    def test_extreme_jacobian_columns(self, model):
        m = model

        def assert_unscaled_jacobian_correct(m, scaled=False):
            jac, nlp = get_jacobian(m, include_scaling_factors=scaled)
            out = extreme_jacobian_columns(jac, nlp)
            assert type(out) == list
            assert len(out) == 1
            assert out[0][0] == pytest.approx(1e7)
            assert out[0][1] is m.x[0]

        # No scaling
        assert_unscaled_jacobian_correct(m, scaled=True)
        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_columns(jac, nlp, large=1e8)
        assert len(out) == 0

        out = extreme_jacobian_columns(jac, nlp, small=1e-8)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e7)
        assert out[0][1] is m.x[0]

        out = extreme_jacobian_columns(jac, nlp, large=1e8, small=1e-8)
        assert len(out) == 0

        # Constraint scaling
        set_scaling_factor(m.dummy_eqn[0], 1e-7)
        set_scaling_factor(m.dummy_eqn[4], 1e7)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_columns(jac, nlp)
        assert len(out) == 0

        # Variable scaling
        set_scaling_factor(m.x[1], 1e7)
        set_scaling_factor(m.x[2], 1e-7)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_columns(jac, nlp)
        assert len(out) == 2
        assert out[0][0] == pytest.approx(1e-7)
        assert out[0][1] is m.x[1]
        assert out[1][0] == pytest.approx(1e8)
        assert out[1][1] is m.x[2]

        # More constraint scaling
        set_scaling_factor(m.dummy_eqn[1], 1e7)
        set_scaling_factor(m.dummy_eqn[2], 1e-8)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_columns(jac, nlp)
        assert len(out) == 0

    @pytest.mark.unit
    def test_extreme_jacobian_entries(self, model):
        m = model

        def assert_unscaled_jacobian_correct(m, scaled=False):
            jac, nlp = get_jacobian(m, include_scaling_factors=scaled)
            out = extreme_jacobian_entries(jac, nlp)
            assert type(out) == list
            assert len(out) == 2
            assert out[0][0] == pytest.approx(1e7)
            assert out[0][1] is m.dummy_eqn[0]
            assert out[0][2] is m.x[0]
            assert out[1][0] == pytest.approx(1e-7)
            assert out[1][1] is m.dummy_eqn[4]
            assert out[1][2] is m.x[4]

        # No scaling
        assert_unscaled_jacobian_correct(m, scaled=False)
        assert_unscaled_jacobian_correct(m, scaled=True)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_entries(jac, nlp, large=1e8)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e-7)
        assert out[0][1] is m.dummy_eqn[4]
        assert out[0][2] is m.x[4]

        out = extreme_jacobian_entries(jac, nlp, small=1e-8)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e7)
        assert out[0][1] is m.dummy_eqn[0]
        assert out[0][2] is m.x[0]

        out = extreme_jacobian_entries(jac, nlp, zero=1e-6)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e7)
        assert out[0][1] is m.dummy_eqn[0]
        assert out[0][2] is m.x[0]

        # Constraint scaling
        set_scaling_factor(m.dummy_eqn[0], 1e-7)
        set_scaling_factor(m.dummy_eqn[4], 1e7)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_entries(jac, nlp)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e-7)
        assert out[0][1] is m.dummy_eqn[0]
        assert out[0][2] is m.x[4]

        out = extreme_jacobian_entries(jac, nlp, zero=1e-6)
        assert len(out) == 0

        out = extreme_jacobian_entries(jac, nlp, small=1e-8)
        assert len(out) == 0

        # Variable scaling
        set_scaling_factor(m.x[1], 1e7)
        set_scaling_factor(m.x[2], 1e-7)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_entries(jac, nlp)
        assert len(out) == 3
        assert out[0][0] == pytest.approx(1e-7)
        assert out[0][1] is m.dummy_eqn[0]
        assert out[0][2] is m.x[4]
        assert out[1][0] == pytest.approx(1e-7)
        assert out[1][1] is m.dummy_eqn[1]
        assert out[1][2] is m.x[1]
        assert out[2][0] == pytest.approx(1e8)
        assert out[2][1] is m.dummy_eqn[2]
        assert out[2][2] is m.x[2]

        # More constraint scaling
        set_scaling_factor(m.dummy_eqn[1], 1e7)
        set_scaling_factor(m.dummy_eqn[2], 1e-8)

        assert_unscaled_jacobian_correct(m, scaled=False)

        jac, nlp = get_jacobian(m)
        out = extreme_jacobian_entries(jac, nlp)
        assert len(out) == 1
        assert out[0][0] == pytest.approx(1e-7)
        assert out[0][1] is m.dummy_eqn[0]
        assert out[0][2] is m.x[4]


class TestCheckParallelJacobian:
    @pytest.mark.unit
    def test_invalid_direction(self):
        m = ConcreteModel()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Unrecognised value for direction (foo). " "Must be 'row' or 'column'."
            ),
        ):
            check_parallel_jacobian(m, direction="foo")

    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.v1 = Var(initialize=1e-8)
        m.v2 = Var()
        m.v3 = Var()
        m.v4 = Var()

        m.c1 = Constraint(expr=m.v1 == m.v2 - 0.99999 * m.v4)
        m.c2 = Constraint(expr=m.v1 + 1.00001 * m.v4 == 1e-8 * m.v3)
        m.c3 = Constraint(expr=1e8 * (m.v1 + m.v4) + 1e10 * m.v2 == 1e-6 * m.v3)
        m.c4 = Constraint(expr=-m.v1 == -0.99999 * (m.v2 - m.v4))

        return m

    @pytest.mark.unit
    def test_rows(self, model):
        assert check_parallel_jacobian(model, direction="row") == [(model.c1, model.c4)]
        assert check_parallel_jacobian(model, direction="row", tolerance=0) == []

    @pytest.mark.unit
    def test_columns(self, model):
        pcol = check_parallel_jacobian(model, direction="column")

        expected = [
            ("v1", "v2"),
            ("v1", "v4"),
            ("v2", "v4"),
        ]

        for i in pcol:
            assert tuple(sorted([i[0].name, i[1].name])) in expected
