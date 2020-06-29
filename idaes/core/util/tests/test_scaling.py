##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
This module contains tests for scaling.
"""

import pytest
import pyomo.environ as pyo
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.scaling import (
    ScalingBasis,
    calculate_scaling_factors,
    badly_scaled_var_generator,
    grad_fd,
    constraint_fd_autoscale,
    scale_single_constraint,
    scale_constraints
)

__author__ = "John Eslick, Tim Bartholomew"


class TestCalculateScalingFactors():
    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=2, bounds=(1, 7))
        m.y = pyo.Var(initialize=3)
        m.z = pyo.Var([1, 2, 3], initialize=0)
        m.b1 = pyo.Block()
        m.b1.c1 = pyo.Constraint(expr=m.z[1] == m.x + m.y)
        m.b1.c2 = pyo.Constraint(expr=m.z[2] == m.x * m.y)
        m.e1 = pyo.Expression(expr=m.y * m.x)
        m.c3 = pyo.Constraint(expr=m.z[3] == m.e1)

        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_expression = pyo.Suffix()
        m.nominal_value = pyo.Suffix()
        m.b1.scaling_expression = pyo.Suffix()

        m.nominal_value[m.x] = 10
        m.scaling_factor[m.y] = 1 / 11
        m.scaling_factor[m.e1] = 1 / 200
        m.scaling_expression[m.z[1]] = 1 / m.y
        m.scaling_expression[m.z[2]] = 1 / (2 * m.x * m.y)
        m.b1.scaling_expression[m.b1.c1] = 1 / m.x
        m.b1.scaling_expression[m.b1.c2] = 1 / (m.x * m.y)
        m.scaling_expression[m.c3] = 1 / m.e1
        return m

    @pytest.mark.unit
    def test_value_basis(self, model):
        # Test value based scaling factor calculations
        calculate_scaling_factors(model, basis=ScalingBasis.Value)
        # check that scaling factors are correctly calculated based on value
        assert model.scaling_factor[model.z[1]] == pytest.approx(1 / 3)
        assert model.scaling_factor[model.z[2]] == pytest.approx(1 / 12)
        assert model.b1.scaling_factor[model.b1.c1] == pytest.approx(1 / 2)
        assert model.b1.scaling_factor[model.b1.c2] == pytest.approx(1 / 6)
        assert model.scaling_factor[model.c3] == pytest.approx(1 / 6)

    @pytest.mark.unit
    def test_inverse_basis(self, model):
        # Test inverse based scaling factor calculations
        calculate_scaling_factors(model, basis=ScalingBasis.InverseVarScale)
        assert model.scaling_factor[model.z[1]] == pytest.approx(1 / 11)
        assert model.scaling_factor[model.z[2]] == pytest.approx(1 / 220)
        assert model.b1.scaling_factor[model.b1.c1] == pytest.approx(1 / 10)
        assert model.b1.scaling_factor[model.b1.c2] == pytest.approx(1 / 110)
        assert model.scaling_factor[model.c3] == pytest.approx(1 / 200)

    @pytest.mark.unit
    def test_var_basis(self, model):
        # Test scaling factor based calculation
        calculate_scaling_factors(model, basis=ScalingBasis.VarScale)
        assert model.scaling_factor[model.z[1]] == pytest.approx(11)
        assert model.scaling_factor[model.z[2]] == pytest.approx(55)
        assert model.b1.scaling_factor[model.b1.c1] == pytest.approx(10)
        assert model.b1.scaling_factor[model.b1.c2] == pytest.approx(110)
        assert model.scaling_factor[model.c3] == pytest.approx(200)

    @pytest.mark.unit
    def test_lower_basis(self, model):
        # Test scaling factor based on lower bound, falling back on 1 where no
        # lower bound exists
        calculate_scaling_factors(model, basis=ScalingBasis.Lower)
        assert model.scaling_factor[model.z[1]] == pytest.approx(1 / 1)
        assert model.scaling_factor[model.z[2]] == pytest.approx(1 / 2)
        assert model.b1.scaling_factor[model.b1.c1] == pytest.approx(1 / 1)
        assert model.b1.scaling_factor[model.b1.c2] == pytest.approx(1 / 1)
        assert model.scaling_factor[model.c3] == pytest.approx(1)

    @pytest.mark.unit
    def test_upper_basis(self, model):
        # Test scaling factor based on upper bound, falling back on 1 where no
        # upper bound exists
        calculate_scaling_factors(model, basis=ScalingBasis.Upper)
        assert model.scaling_factor[model.z[1]] == pytest.approx(1 / 1)
        assert model.scaling_factor[model.z[2]] == pytest.approx(1 / 14)
        assert model.b1.scaling_factor[model.b1.c1] == pytest.approx(1 / 7)
        assert model.b1.scaling_factor[model.b1.c2] == pytest.approx(1 / 7)
        assert model.scaling_factor[model.c3] == pytest.approx(1)

    @pytest.mark.unit
    def test_mixed_basis(self, model):
        # Test scaling factor based on midpoint of bounds and falling back on
        # scale factor
        calculate_scaling_factors(model, basis=[ScalingBasis.Mid,
                                                ScalingBasis.InverseVarScale])
        assert model.scaling_factor[model.z[1]] == pytest.approx(1 / 11)
        assert model.scaling_factor[model.z[2]] == pytest.approx(1 / 88)
        assert model.b1.scaling_factor[model.b1.c1] == pytest.approx(1 / 4)
        assert model.b1.scaling_factor[model.b1.c2] == pytest.approx(1 / 44)
        assert model.scaling_factor[model.c3] == pytest.approx(1 / 200)

    @pytest.mark.unit
    def test_expressions(self, model):
        # Check that the scaling factor calculations didn't change the scaling
        # factor expressions
        assert pyo.value(model.scaling_expression[model.z[1]]) == \
               pytest.approx(1 / 3)
        assert pyo.value(model.scaling_expression[model.z[2]]) == \
               pytest.approx(1 / 12)
        assert pyo.value(model.b1.scaling_expression[model.b1.c1]) == \
               pytest.approx(1 / 2)
        assert pyo.value(model.b1.scaling_expression[model.b1.c2]) == \
               pytest.approx(1 / 6)
        assert model.scaling_factor[model.c3] == pytest.approx(1 / 200)


@pytest.mark.unit
def test_find_badly_scaled_vars():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=1e6)
    m.y = pyo.Var(initialize=1e-8)
    m.z = pyo.Var(initialize=1e-20)
    m.b = pyo.Block()
    m.b.w = pyo.Var(initialize=1e10)

    a = [id(v) for v, sv in badly_scaled_var_generator(m)]
    assert id(m.x) in a
    assert id(m.y) in a
    assert id(m.b.w) in a
    assert id(m.z) not in a

    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    m.b.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    m.scaling_factor[m.x] = 1e-6
    m.scaling_factor[m.y] = 1e6
    m.scaling_factor[m.z] = 1
    m.b.scaling_factor[m.b.w] = 1e-5

    a = [id(v) for v, sv in badly_scaled_var_generator(m)]
    assert id(m.x) not in a
    assert id(m.y) not in a
    assert id(m.b.w) in a
    assert id(m.z) not in a


@pytest.mark.unit
def test_grad_fd():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=1e6)
    m.y = pyo.Var(initialize=1e8)
    m.z = pyo.Var(initialize=1e4)
    m.e = pyo.Expression(expr=100 * (m.x / m.z)**2)

    sfx = 1e-6
    sfy = 1e-7
    sfz = 1e-4

    m.c1 = pyo.Constraint(expr=0 == 10 * m.x - 4 * m.y + 5)
    m.c2 = pyo.Constraint(expr=0 == m.e + 5)

    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    m.scaling_expression = pyo.Suffix()
    m.scaling_factor[m.x] = sfx
    m.scaling_factor[m.y] = sfy
    m.scaling_factor[m.z] = sfz
    m.scaling_expression[m.c2] = 1/(100 * (m.x / m.z)**2)
    calculate_scaling_factors(m)

    g, v = grad_fd(m.c1, scaled=False)
    if id(v[0]) == id(m.x):
        ix = 0
        iz = 1
    else:
        iz = 0
        ix = 1
    assert g[ix] == pytest.approx(10, rel=1e-2)
    assert g[iz] == pytest.approx(-4, rel=1e-2)

    g, v = grad_fd(m.c1, scaled=True)
    if id(v[0]) == id(m.x):
        ix = 0
        iz = 1
    else:
        iz = 0
        ix = 1
    assert g[ix] == pytest.approx(10/sfx, rel=1e-2)
    assert g[iz] == pytest.approx(-4/sfy, rel=1e-2)

    g, v = grad_fd(m.c2, scaled=False)
    if id(v[0]) == id(m.x):
        ix = 0
        iz = 1
    else:
        iz = 0
        ix = 1
    assert g[ix] == pytest.approx(pyo.value(200*m.x/m.z**2), rel=1e-2)
    assert g[iz] == pytest.approx(pyo.value(-200*m.x**2/m.z**3), rel=1e-2)
    assert g[ix] == pytest.approx(2, rel=1e-2)
    assert g[iz] == pytest.approx(-200, rel=1e-3)

    g, v = grad_fd(c=m.c2, scaled=True)
    if id(v[0]) == id(m.x):
        ix = 0
        iz = 1
    else:
        iz = 0
        ix = 1
    # makesure the replacements and subs in calculating the scaled gradient
    # didn't have the side affect of changing anything.  Especailly the named
    # expression
    assert pyo.value(m.x) == 1e6
    assert pyo.value(m.z) == 1e4
    assert pyo.value(m.e) == pytest.approx(pyo.value(100 * (m.x / m.z)**2))
    assert g[ix] == pytest.approx(
        pyo.value(200*(sfz/sfx)**2*sfx*m.x/((sfz*m.z)**2)) *
        m.scaling_factor[m.c2], rel=1e-2, )
    assert g[iz] == pytest.approx(
        pyo.value(-200*(sfz/sfx)**2*(sfx*m.x)**2/(sfz*m.z)**3) *
        m.scaling_factor[m.c2], rel=1e-2, )
    assert g[ix] == pytest.approx(2, rel=1e-2)
    assert g[iz] == pytest.approx(-2, rel=1e-3)

    # change the constraint scale factor to get the gradient back over 100 and
    # make sure the auto scale correctly accounts for existing scale factor
    m.scaling_factor[m.c2] = 1e-1
    constraint_fd_autoscale(m.c2)
    assert m.scaling_factor[m.c2] == pytest.approx(5e-5, rel=1e-2)


class TestScaleSingleConstraint():
    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=500)
        m.c1 = pyo.Constraint(expr=m.x <= 1e3)
        m.c2 = pyo.Constraint(expr=m.x == 1e3)
        m.c3 = pyo.Constraint(expr=m.x >= 1e3)
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.c1] = 1 / 1e3
        m.scaling_factor[m.c2] = 1 / 1e3
        m.scaling_factor[m.c3] = 1 / 1e3
        return m

    @pytest.mark.unit
    def test_unscaled_constraints(self, model):
        assert model.c1.lower is None
        assert model.c1.body is model.x
        assert model.c1.upper.value == pytest.approx(1e3)
        assert model.c2.lower.value == pytest.approx(1e3)
        assert model.c2.body is model.x
        assert model.c2.upper.value == pytest.approx(1e3)
        assert model.c3.lower.value == pytest.approx(1e3)
        assert model.c3.body is model.x
        assert model.c3.upper is None

    @pytest.mark.unit
    def test_not_constraint(self, model):
        with pytest.raises(TypeError):
            scale_single_constraint(model.x)

    @pytest.mark.unit
    def test_less_than_constraint(self, model):
        scale_single_constraint(model.c1)
        assert model.c1.lower is None
        assert model.c1.body() == pytest.approx(model.x.value / 1e3)
        assert model.c1.upper.value == pytest.approx(1)

    @pytest.mark.unit
    def test_equality_constraint(self, model):
        scale_single_constraint(model.c2)
        assert model.c2.lower.value == pytest.approx(1)
        assert model.c2.body() == pytest.approx(model.x.value / 1e3)
        assert model.c2.upper.value == pytest.approx(1)

    @pytest.mark.unit
    def test_greater_than_constraint(self, model):
        scale_single_constraint(model.c3)
        assert model.c3.lower.value == pytest.approx(1)
        assert model.c3.body() == pytest.approx(model.x.value / 1e3)
        assert model.c3.upper is None

    @pytest.mark.unit
    def test_scaling_factor_and_expression_replacement(self, model):
        model.c4 = pyo.Constraint(expr=model.x <= 1e6)
        model.scaling_factor[model.c4] = 1e-6
        model.scaling_expression = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        model.scaling_expression[model.c4] = 1 / model.x
        scale_single_constraint(model.c4)
        assert model.c4.upper.value == pytest.approx(1)
        assert model.scaling_factor[model.c4] == pytest.approx(1)
        assert model.scaling_expression[model.c4] == pytest.approx(1)

    @pytest.fixture(scope="class")
    def model2(self):
        m = pyo.ConcreteModel()
        m.y = pyo.Var()
        m.c = pyo.Constraint(expr=m.y <= 1e3)
        return m

    @pytest.mark.unit
    def test_no_scaling_factor_suffix(self, model2):
        with pytest.raises(ConfigurationError):
            scale_single_constraint(model2.c)

    @pytest.mark.unit
    def test_no_scaling_factor(self, model2):
        model2.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        scale_single_constraint(model2.c)
        assert model2.c.upper.value == pytest.approx(1e3)


class TestScaleConstraints():
    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=1e3)
        m.y = pyo.Var(initialize=1e6)
        m.c1 = pyo.Constraint(expr=m.x == 1e3)
        m.c2 = pyo.Constraint(expr=m.y == 1e6)
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.c1] = 1e-3
        m.scaling_factor[m.c2] = 1e-6

        m.b1 = pyo.Block()
        m.b1.c1 = pyo.Constraint(expr=m.x <= 1e9)
        m.b1.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.b1.scaling_factor[m.b1.c1] = 1e-9

        m.b1.b2 = pyo.Block()
        m.b1.b2.c1 = pyo.Constraint(expr=m.x <= 1e12)
        m.b1.b2.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.b1.b2.scaling_factor[m.b1.b2.c1] = 1e-12

        return m

    @pytest.mark.unit
    def test_scale_one_block(self, model):
        scale_constraints(model, descend_into=False)
        # scaled
        assert model.c1.lower.value == pytest.approx(1)
        assert model.c1.body() == pytest.approx(model.x.value / 1e3)
        assert model.c1.upper.value == pytest.approx(1)
        assert model.c2.lower.value == pytest.approx(1)
        assert model.c2.body() == pytest.approx(model.y.value / 1e6)
        assert model.c2.upper.value == pytest.approx(1)
        # unscaled
        assert model.b1.c1.upper.value == pytest.approx(1e9)
        assert model.b1.b2.c1.upper.value == pytest.approx(1e12)

    @pytest.mark.unit
    def test_scale_model(self, model):
        scale_constraints(model)
        assert model.c1.upper.value == pytest.approx(1)
        assert model.b1.c1.upper.value == pytest.approx(1)
        assert model.b1.b2.c1.upper.value == pytest.approx(1)
