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
    apply_scaling,
)

__author__ = "John Eslick"

def test_1():
    # Test scalar variables
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=2, bounds=(1,7))
    m.y = pyo.Var(initialize=3)
    m.z = pyo.Var([1,2], initialize=0)
    m.b1 = pyo.Block()
    m.b1.c1 = pyo.Constraint(expr=m.z[1]==m.x + m.y)
    m.b1.c2 = pyo.Constraint(expr=m.z[2]==m.x*m.y)

    m.scaling_factor = pyo.Suffix()
    m.scaling_expression = pyo.Suffix()
    m.b1.scaling_expression = pyo.Suffix()

    m.scaling_factor[m.x] = 1/10
    m.scaling_factor[m.y] = 1/11
    m.scaling_expression[m.z[1]] = 1/m.y
    m.scaling_expression[m.z[2]] = 1/(2*m.x*m.y)
    m.b1.scaling_expression[m.b1.c1] = 1/m.x
    m.b1.scaling_expression[m.b1.c2] = 1/(m.x*m.y)

    # Test value based scaling factor calculations
    apply_scaling(m, basis=ScalingBasis.Value)
    # check that scaling factors are correctly calculated based on value
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/3)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/12)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/2)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/6)

    # Test scaling factor based calculation
    apply_scaling(m, basis=ScalingBasis.InverseVarScale)
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/11)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/220)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/10)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/110)

    # Test scaling factor based calculation
    apply_scaling(m, basis=ScalingBasis.VarScale)
    assert m.scaling_factor[m.z[1]] == pytest.approx(11)
    assert m.scaling_factor[m.z[2]] == pytest.approx(55)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(10)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(110)

    # Test scaling factor based on midpoint of bounds and falling back on scale
    # factor
    apply_scaling(m, basis=[ScalingBasis.Mid, ScalingBasis.InverseVarScale])
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/11)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/88)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/4)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/44)

    # Test scaling factor based on lower bound, falling back on 1 where no lower
    # bound exists
    apply_scaling(m, basis=ScalingBasis.Lower)
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/1)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/2)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/1)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/1)

    # Test scaling factor based on upper bound, falling back on 1 where no lower
    # bound exists
    apply_scaling(m, basis=ScalingBasis.Upper)
    assert m.scaling_factor[m.z[1]] == pytest.approx(1/1)
    assert m.scaling_factor[m.z[2]] == pytest.approx(1/14)
    assert m.b1.scaling_factor[m.b1.c1] == pytest.approx(1/7)
    assert m.b1.scaling_factor[m.b1.c2] == pytest.approx(1/7)

    # Check that the scaling factor calculations didn't change the scaling
    # factor expressions
    assert pyo.value(m.scaling_expression[m.z[1]]) == pytest.approx(1/3)
    assert pyo.value(m.scaling_expression[m.z[2]]) == pytest.approx(1/12)
    assert pyo.value(m.b1.scaling_expression[m.b1.c1]) == pytest.approx(1/2)
    assert pyo.value(m.b1.scaling_expression[m.b1.c2]) == pytest.approx(1/6)
