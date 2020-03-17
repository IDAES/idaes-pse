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
This module contains tests for the variable replace transformation.
"""

import pytest
import pyomo.environ as pyo
import idaes.core.plugins
from idaes.core.util.model_statistics import degrees_of_freedom


__author__ = "John Eslick"


@pytest.fixture
def model():
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3, 4], initialize=2)
    m.y = pyo.Var(initialize=3)
    m.z = pyo.Var([1, 2, 3, 4], initialize=0)
    m.e = pyo.Expression(expr=m.x[1] + m.x[2])
    m.z.fix()
    m.z[4] = 4.0

    # In the first round this becomes x[1].fix(m.z[1])
    m.c1 = pyo.Constraint(expr=m.z[1] == m.x[1])
    # In the first round this beconmes y.fix(m.z[2])
    m.c2 = pyo.Constraint(expr=m.z[2] == m.y)
    # In the second round this becomes x[2].fix(m.z[3]-m.x[1])
    m.c3 = pyo.Constraint(expr=m.z[3] == m.e)
    # This constraint is nonlinear but m.x[1] will be fixed
    m.c4 = pyo.Constraint(expr=m.z[4] == m.x[3] ** 2 + m.x[1])
    # This constraint is eliminated by a substitution, x[3] and x[4] are not
    # fixed, so it's a linear equation with two variables.  This should be in the
    # second round after x[1] becomes fixed in the first
    m.c5 = pyo.Constraint(expr=m.z[4] == m.x[4] + m.x[1] + m.x[3])

    # Solution
    m.x[1] = 0
    m.x[2] = 0
    m.x[3] = 2
    m.x[4] = 2
    m.y = 0

    # make a dict to check variable values against (vars and named expressions)
    m.sol = {}
    for v in m.component_data_objects((pyo.Var, pyo.Expression)):
        m.sol[id(v)] = pyo.value(v)
    return m


def test_model(model):
    m = model
    assert degrees_of_freedom(m) == 0
    assert len([c for c in m.component_data_objects(pyo.Constraint, active=True)]) == 5


def test_transform(model):
    m = model
    elim = pyo.TransformationFactory("simple_equality_eliminator")
    elim.apply_to(m)
    assert degrees_of_freedom(m) == 0

    assert len([c for c in m.component_data_objects(pyo.Constraint, active=True)]) == 1
    assert m.x[1].fixed
    assert m.x[2].fixed
    assert m.y.fixed


def test_reverse(model):
    m = model
    elim = pyo.TransformationFactory("simple_equality_eliminator")
    elim.apply_to(m)

    # mess up the values for all the substituted vars
    for i, v in elim._subs_map.items():
        v = 1e12
    elim.revert()

    for v in m.component_data_objects((pyo.Var, pyo.Expression)):
        assert pyo.value(v) == m.sol[id(v)]

    assert len([c for c in m.component_data_objects(pyo.Constraint, active=True)]) == 5
    assert not m.x[1].fixed
    assert not m.x[2].fixed
    assert not m.y.fixed

def test_bounds(model):
    m = model

    elim = pyo.TransformationFactory("simple_equality_eliminator")

    m.x[3].setlb(-1)
    m.x[3].setub(10)
    m.x[4].setlb(1)
    m.x[4].setub(12)

    elim.apply_to(m)

    # The bounds here come from constraint c5.  z4 is 4, x1 is 0, and x4 and x5
    # are not fixed, so you get 4 = x3 + x4.  I'm not sure if x3 or x4 will be
    # replaced, so I calculated bounds on both.

    # For x3 the bounds that come from x4 are 4 - x4.lb and 4 - x4.ub or (-8 and 3)
    # so the tightest set of bounds from x3 and x4 are (-1, 3)

    # If instead x3 is replaced, the bounds on x4 that come from x3 are
    # 4 - 10 and 4 + 1 or (-6, 5).  The tightest set of bounds on x4 is (1, 5)
    if id(m.x[3]) not in elim._subs_map:
        assert m.x[3].lb == -1
        assert m.x[3].ub == 3
    else:
        assert m.x[4].lb == 1
        assert m.x[4].ub == 5
