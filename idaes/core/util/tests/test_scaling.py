##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
This module contains tests for scaling expressions.
"""

import pytest
import pyomo.environ as pyo
from idaes.core.util.scaling import (
    ScalingBasis,
    calculate_scaling_factors,
)

__author__ = "John Eslick"


import pytest
import pyomo.environ as pyo
from idaes.core.util.scaling import (
    ScalingBasis,
    calculate_scaling_factors,
    badly_scaled_var_generator,
    grad_fd,
    constraint_fd_autoscale,
)

__author__ = "John Eslick"

def test_1():
    # Test scalar variables
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=2, bounds=(1,7))
    m.y = pyo.Var(initialize=3)
    m.z = pyo.Var([1,2,3], initialize=0)
    m.b1 = pyo.Block()
    m.b1.c1 = pyo.Constraint(expr=m.z[1]==m.x + m.y)
    m.b1.c2 = pyo.Constraint(expr=m.z[2]==m.x*m.y)
    m.e1 = pyo.Expression(expr=m.y*m.x)
    m.c3 = pyo.Constraint(expr=m.z[3]==m.e1)

    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    m.scaling_expression = pyo.Suffix()
    m.nominal_value = pyo.Suffix()
    m.b1.scaling_expression = pyo.Suffix()

    m.nominal_value[m.x] = 10
    m.scaling_factor[m.y] = 1/11
    m.scaling_factor[m.e1] = 1/200
    m.scaling_expression[m.z[1]] = 1/m.y
    m.scaling_expression[m.z[2]] = 1/(2*m.x*m.y)
    m.b1.scaling_expression[m.b1.c1] = 1/m.x
    m.b1.scaling_expression[m.b1.c2] = 1/(m.x*m.y)
    m.scaling_expression[m.c3] = 1/m.e1

    # Test value based scaling factor calculations
    calculate_scaling_factors(m, basis=ScalingBasis.Value)
    # check that scaling factors are correctly calculated based on value
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/3)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/12)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/2)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/6)
    assert m.scaling_factor[m.c3] == pytest.approx(1/6)

    # Test scaling factor based calculation
    calculate_scaling_factors(m, basis=ScalingBasis.InverseVarScale)
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/11)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/220)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/10)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/110)
    assert m.scaling_factor[m.c3] == pytest.approx(1/200)

    # Test scaling factor based calculation
    calculate_scaling_factors(m, basis=ScalingBasis.VarScale)
    assert m.scaling_factor[m.z[1]] == pytest.approx(11)
    assert m.scaling_factor[m.z[2]] == pytest.approx(55)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(10)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(110)
    assert m.scaling_factor[m.c3] == pytest.approx(200)

    # Test scaling factor based on midpoint of bounds and falling back on scale
    # factor
    calculate_scaling_factors(m, basis=[ScalingBasis.Mid, ScalingBasis.InverseVarScale])
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/11)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/88)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/4)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/44)
    assert m.scaling_factor[m.c3] == pytest.approx(1/200)

    # Test scaling factor based on lower bound, falling back on 1 where no lower
    # bound exists
    calculate_scaling_factors(m, basis=ScalingBasis.Lower)
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/1)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/2)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/1)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/1)
    assert m.scaling_factor[m.c3] == pytest.approx(1)

    # Test scaling factor based on upper bound, falling back on 1 where no lower
    # bound exists
    calculate_scaling_factors(m, basis=ScalingBasis.Upper)
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/1)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/14)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/7)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/7)
    assert m.scaling_factor[m.c3] == pytest.approx(1)

    # Check that the scaling factor calculations didn't change the scaling
    # factor expressions
    assert pyo.value(m.scaling_expression[m.z[1]]) == pytest.approx(1/3)
    assert pyo.value(m.scaling_expression[m.z[2]]) == pytest.approx(1/12)
    assert pyo.value(m.b1.scaling_expression[m.b1.c1]) == pytest.approx(1/2)
    assert pyo.value(m.b1.scaling_expression[m.b1.c2]) == pytest.approx(1/6)
    assert m.scaling_factor[m.c3] == pytest.approx(1)

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
        pyo.value(200*(sfz/sfx)**2*sfx*m.x/((sfz*m.z)**2))*m.scaling_factor[m.c2],
        rel=1e-2,
    )
    assert g[iz] == pytest.approx(
        pyo.value(-200*(sfz/sfx)**2*(sfx*m.x)**2/(sfz*m.z)**3)*m.scaling_factor[m.c2],
        rel=1e-2,
    )
    assert g[ix] == pytest.approx(2, rel=1e-2)
    assert g[iz] == pytest.approx(-2, rel=1e-3)

    # change the constraint scale factor to get the gradient back over 100 and
    # make sure the auto scale correctly accounts for existing scale factor
    m.scaling_factor[m.c2] = 1e-1
    constraint_fd_autoscale(m.c2)
    assert m.scaling_factor[m.c2] == pytest.approx(5e-5, rel=1e-2)
