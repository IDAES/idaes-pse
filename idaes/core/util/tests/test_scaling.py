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
This module contains tests for scaling.
"""

import pytest
import pyomo.environ as pyo
import pyomo.dae as dae
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.relational_expr import (
    EqualityExpression,
    InequalityExpression,
    RangedExpression,
)
from pyomo.network import Port, Arc
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import number_activated_objectives
import idaes.core.util.scaling as sc
import logging
from idaes.core.solvers.tests.test_petsc import dae_with_non_time_indexed_constraint

__author__ = "John Eslick, Tim Bartholomew"


@pytest.mark.unit
def test_none_left_mult():
    with pytest.raises(
        TypeError, match="unsupported operand type\(s\) for \*: 'int' and 'NoneType'"
    ):
        assert sc.__none_left_mult(4, None) is None
    with pytest.raises(
        TypeError, match="unsupported operand type\(s\) for \*: 'float' and 'NoneType'"
    ):
        assert sc.__none_left_mult(4.0, None) is None
    assert sc.__none_left_mult(None, 4) is None
    assert sc.__none_left_mult(3, 4) == 12


@pytest.mark.unit
def test_scale_constraint():
    m = pyo.ConcreteModel()
    m.x = pyo.Var()
    m.y = pyo.Var()

    m.c_eq = pyo.Constraint(expr=m.x == m.y)
    m.c_ineq = pyo.Constraint(expr=m.x <= m.y)
    m.c_range = pyo.Constraint(expr=(0, m.x + m.y, 1))

    sc.__scale_constraint(m.c_eq, 2)
    assert isinstance(m.c_eq.expr, EqualityExpression)
    sc.__scale_constraint(m.c_ineq, 0.5)
    assert isinstance(m.c_ineq.expr, InequalityExpression)
    sc.__scale_constraint(m.c_range, 10)
    assert isinstance(m.c_range.expr, RangedExpression)


@pytest.mark.unit
def test_scale_arcs(caplog):
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3, 4])
    m.y = pyo.Var([1, 2, 3, 4])
    m.z = pyo.Var([1, 2, 3, 4])
    m.w = pyo.Var([1, 2, 3, 4])

    def rule_schmequality(port, name, index_set):
        # Copied from Pyomo's implementation of Port.Equality
        #  ___________________________________________________________________________
        #
        #  Pyomo: Python Optimization Modeling Objects
        #  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
        #  Under the terms of Contract DE-NA0003525 with National Technology and
        #  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
        #  rights in this software.
        #  This software is distributed under the 3-clause BSD License.
        #  ___________________________________________________________________________
        #
        #  This module was originally developed as part of the PyUtilib project
        #  Copyright (c) 2008 Sandia Corporation.
        #  This software is distributed under the BSD License.
        #  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
        #  the U.S. Government retains certain rights in this software.
        #  ___________________________________________________________________________
        for arc in port.arcs(active=True):
            Port._add_equality_constraint(arc, name, index_set)

    m.p1 = Port()
    m.p1.add(m.x[1], name="x")
    m.p1.add(m.y[1], name="y")
    m.p1.add(m.z[1], name="z", rule=Port.Extensive)
    m.p1.add(m.w[1], name="w", rule=rule_schmequality)

    m.p = Port([2, 3, 4])
    for i in m.p.keys():
        m.p[i].add(m.x[i], name="x")
        m.p[i].add(m.y[i], name="y")
        m.p[i].add(m.z[i], name="z", rule=Port.Extensive)
        m.p[i].add(m.w[i], name="w", rule=rule_schmequality)

    def arc_rule(b, i):
        if i == 1:
            return (m.p1, m.p[2])
        elif i == 2:
            return (m.p[3], m.p[4])

    m.arcs = Arc([1, 2], rule=arc_rule)

    sc.set_scaling_factor(m.x, 10)
    sc.set_scaling_factor(m.y, 20)
    sc.set_scaling_factor(m.z, 0.1)
    sc.set_scaling_factor(m.w, 0.05)
    sc.set_scaling_factor(m.x[1], 5)

    # make sure there is no error if the scaling is done with unexpanded arcs
    sc.scale_arc_constraints(m)

    # expand and make sure it works
    caplog.set_level(logging.WARNING)
    caplog.clear()
    pyo.TransformationFactory("network.expand_arcs").apply_to(m)
    sc.scale_arc_constraints(m)

    logrec = caplog.records[0]
    assert logrec.levelno == logging.WARNING
    assert "Variable z on Port p1" in logrec.message
    assert "on arcs[1] will" in logrec.message
    logrec = caplog.records[1]
    assert logrec.levelno == logging.WARNING
    assert "Variable w on Port p1" in logrec.message
    assert "on arcs[1] will" in logrec.message
    logrec = caplog.records[2]
    assert logrec.levelno == logging.WARNING
    assert "Variable z on Port p[3]" in logrec.message
    assert "on arcs[2] will" in logrec.message
    logrec = caplog.records[3]
    assert logrec.levelno == logging.WARNING
    assert "Variable w on Port p[3]" in logrec.message
    assert "on arcs[2] will" in logrec.message

    m.x[1] = 1
    m.x[2] = 2
    m.x[3] = 3
    m.x[4] = 4
    m.y[1] = 11
    m.y[2] = 12
    m.y[3] = 13
    m.y[4] = 14

    # for all the arc constraints the differnce is 1 the scale factor is the
    # smallest scale factor for variables in a constraint.  Make sure the
    # constraints are scaled as expected.
    assert abs(m.arcs_expanded[1].x_equality.body()) == 5
    assert abs(m.arcs_expanded[2].x_equality.body()) == 10
    assert abs(m.arcs_expanded[1].y_equality.body()) == 20
    assert abs(m.arcs_expanded[2].y_equality.body()) == 20


@pytest.mark.unit
def test_map_scaling_factor(caplog):
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3, 4])
    sc.set_scaling_factor(m.x[1], 11)
    sc.set_scaling_factor(m.x[2], 12)
    sc.set_scaling_factor(m.x[3], 13)
    caplog.set_level(logging.WARNING)
    caplog.clear()
    assert sc.map_scaling_factor(m.x.values(), warning=True) == 1
    logrec = caplog.records[0]
    assert logrec.levelno == logging.WARNING
    assert "scaling factor" in logrec.message

    assert sc.map_scaling_factor(m.x.values(), func=max) == 13
    assert sc.map_scaling_factor(m.x.values(), default=20) == 11
    with pytest.raises(TypeError):
        sc.map_scaling_factor(m.x.values(), default=None)

    # test min_scaling_factor with just calls map_scaling_factor
    assert sc.min_scaling_factor(m.x.values()) == 1
    assert sc.min_scaling_factor(m.x.values(), default=14) == 11


@pytest.mark.unit
def test_propagate_indexed_scaling():
    m = pyo.ConcreteModel()
    m.b = pyo.Block()
    m.a = pyo.Var()
    m.x = pyo.Var([1, 2, 3], initialize=1e6)
    m.y = pyo.Var([1, 2, 3], initialize=1e-8)
    m.z = pyo.Var([1, 2, 3], initialize=1e-20)

    @m.Constraint([1, 2, 3])
    def c1(b, i):
        return m.x[i] == 0

    @m.Constraint([1, 2, 3])
    def c2(b, i):
        return m.y[i] == 0

    m.b.w = pyo.Var([1, 2, 3], initialize=1e10)
    m.b.c1 = pyo.Constraint(expr=m.b.w[1] == 0)
    m.b.c2 = pyo.Constraint(expr=m.b.w[2] == 0)

    sc.set_scaling_factor(m.a, 104)
    sc.set_scaling_factor(m.b.c1, 14)

    # Set suffix directly since set_scaling_factor also sets data objects
    m.scaling_factor[m.x] = 11
    m.scaling_factor[m.y] = 13
    m.b.scaling_factor[m.b.w] = 16
    m.scaling_factor[m.c1] = 14

    for i in [1, 2, 3]:
        assert sc.get_scaling_factor(m.x[i]) is None
        assert sc.get_scaling_factor(m.y[i]) is None
        assert sc.get_scaling_factor(m.z[i]) is None
        assert sc.get_scaling_factor(m.b.w[i]) is None
        assert sc.get_scaling_factor(m.c1[i]) is None
        assert sc.get_scaling_factor(m.c2[i]) is None
    assert sc.get_scaling_factor(m.x) == 11
    assert sc.get_scaling_factor(m.y) == 13
    assert sc.get_scaling_factor(m.z) is None
    assert sc.get_scaling_factor(m.b.w) == 16
    assert sc.get_scaling_factor(m.c1) == 14
    assert sc.get_scaling_factor(m.c2) is None

    sc.propagate_indexed_component_scaling_factors(m)
    for i in [1, 2, 3]:
        assert sc.get_scaling_factor(m.x[i]) == 11
        assert sc.get_scaling_factor(m.y[i]) == 13
        assert sc.get_scaling_factor(m.z[i]) is None
        assert sc.get_scaling_factor(m.b.w[i]) == 16
        assert sc.get_scaling_factor(m.c1[i]) == 14
        assert sc.get_scaling_factor(m.c2[i]) is None


@pytest.mark.unit
def test_calculate_scaling_factors():
    r"""This tests the method to find and execute calculate_scaling_factors
    methods, here we make sure they are found and run in a right order.  The
    contents of specific calculate_scaling_factors methods will have to be
    tested in the models they belong to.

        m          To be correct:
       / \           f and g are before e
      a   b          c and d are before a
     / \   \         e is before b
    c   d   e        a and b are before m
           / \
          f   g
    """
    o = []  # list of compoent names in the order their calculate_scaling_factors
    # method is called
    def rule(blk):
        # This rule for building a block just adds a calculate scaling factor
        # function to the block, which adds the block name to the list o.  Then
        # by looking at the order of entries in o, we can see the order in which
        # the calculate_scaling_factors methods where called.
        def ca():
            o.append(blk.name)

        blk.calculate_scaling_factors = ca

    # Create blocks with the tree structure shown in the docstring
    m = pyo.ConcreteModel(name="m", rule=rule)
    m.a = pyo.Block(rule=rule)
    m.b = pyo.Block(rule=rule)
    m.a.c = pyo.Block(rule=rule)
    m.a.d = pyo.Block(rule=rule)
    m.b.e = pyo.Block(rule=rule)
    m.b.e.f = pyo.Block(rule=rule)
    m.b.e.g = pyo.Block(rule=rule)
    # Execute the calculate scaling factors methods.
    sc.calculate_scaling_factors(m)
    # We should iterate through the blocks in construction order so we can
    # depend on the order that the calculate_scaling_factors() are called
    # being deterministic, so we just need to check the specific expected order,
    # even though there are additional "correct" orders.
    assert tuple(o) == ("a.c", "a.d", "a", "b.e.f", "b.e.g", "b.e", "b", "m")


@pytest.mark.unit
def test_set_get_unset(caplog):
    """Make sure the Jacobian from Pynumero matches expectation.  This is
    mostly to ensure we understand the interface and catch if things change.
    """
    m = pyo.ConcreteModel()
    m.x = pyo.Var()
    m.z = pyo.Var([1, 2, 3, 4])
    m.c1 = pyo.Constraint(expr=0 == m.x)

    @m.Constraint([1, 2, 3, 4])
    def c2(b, i):
        return b.z[i] == 0

    m.ex = pyo.Expression(expr=m.x)

    sc.set_scaling_factor(m.z, 10)
    assert sc.get_scaling_factor(m.z) == 10
    assert sc.get_scaling_factor(m.z[1]) == 10
    assert sc.get_scaling_factor(m.z[2]) == 10
    assert sc.get_scaling_factor(m.z[3]) == 10
    assert sc.get_scaling_factor(m.z[4]) == 10
    sc.unset_scaling_factor(m.z)
    assert sc.get_scaling_factor(m.z) is None
    assert sc.get_scaling_factor(m.z[1]) is None
    assert sc.get_scaling_factor(m.z[2]) is None
    assert sc.get_scaling_factor(m.z[3]) is None
    assert sc.get_scaling_factor(m.z[4]) is None
    sc.set_scaling_factor(m.z, 10, data_objects=False)
    assert sc.get_scaling_factor(m.z) == 10
    assert sc.get_scaling_factor(m.z[1]) is None
    assert sc.get_scaling_factor(m.z[2]) is None
    assert sc.get_scaling_factor(m.z[3]) is None
    assert sc.get_scaling_factor(m.z[4]) is None
    sc.set_scaling_factor(m.z, 10)
    sc.unset_scaling_factor(m.z, data_objects=False)
    assert sc.get_scaling_factor(m.z) is None
    assert sc.get_scaling_factor(m.z[1]) == 10
    assert sc.get_scaling_factor(m.z[2]) == 10
    assert sc.get_scaling_factor(m.z[3]) == 10
    assert sc.get_scaling_factor(m.z[4]) == 10
    sc.unset_scaling_factor(m.z)

    caplog.set_level(logging.WARNING)
    caplog.clear()
    assert sc.get_scaling_factor(m.z[1], warning=True) is None
    assert sc.get_scaling_factor(m.z[1], warning=True, default=1) == 1
    for i in [0, 1]:  # two calls should be two log records
        logrec = caplog.records[i]
        assert logrec.levelno == logging.WARNING
        assert "scaling factor" in logrec.message

    # This one is a bit of a mystery, what do you really expect if you provide
    # a default and ask for an exception.  I'll guess if you provide a default,
    # your don't really want an exception.
    assert sc.get_scaling_factor(m.z[1], exception=True, default=1) == 1

    caplog.clear()
    with pytest.raises(KeyError):
        sc.get_scaling_factor(m.z[1], exception=True)
    logrec = caplog.records[0]
    assert logrec.levelno == logging.ERROR
    assert "scaling factor" in logrec.message

    # Okay it's pretty well tested, but make sure it works for constraints and
    # expressions

    sc.set_scaling_factor(m.x, 11)
    sc.set_scaling_factor(m.ex, 2)
    sc.set_scaling_factor(m.c1, 3)
    sc.set_scaling_factor(m.c2, 4)

    assert sc.get_scaling_factor(m.x) == 11
    assert sc.get_scaling_factor(m.ex) == 2
    assert sc.get_scaling_factor(m.c1) == 3
    assert sc.get_scaling_factor(m.c2) == 4
    assert sc.get_scaling_factor(m.c2[1]) == 4
    assert sc.get_scaling_factor(m.c2[2]) == 4
    assert sc.get_scaling_factor(m.c2[3]) == 4
    assert sc.get_scaling_factor(m.c2[4]) == 4

    # Check the underlying suffix
    assert m.scaling_factor[m.x] == 11
    assert m.scaling_factor[m.ex] == 2
    assert m.scaling_factor[m.c1] == 3
    assert m.scaling_factor[m.c2] == 4
    assert m.scaling_factor[m.c2[1]] == 4
    assert m.scaling_factor[m.c2[2]] == 4
    assert m.scaling_factor[m.c2[3]] == 4
    assert m.scaling_factor[m.c2[4]] == 4

    # Test overwriting scaling factors
    sc.set_scaling_factor(m.x, 13, overwrite=False)
    sc.set_scaling_factor(m.ex, 17, overwrite=False)
    # z doesn't have a scaling factor, make sure that 19 gets set as planned
    # to everything besides z[1], whose factor is not overwritten
    sc.set_scaling_factor(m.z[1], 29, overwrite=False)
    sc.set_scaling_factor(m.z, 19, overwrite=False)
    sc.set_scaling_factor(m.c1, 7, overwrite=False)
    sc.set_scaling_factor(m.c2, 23, overwrite=False)

    assert sc.get_scaling_factor(m.x) == 11
    assert sc.get_scaling_factor(m.ex) == 2
    assert sc.get_scaling_factor(m.z) == 19  # z was unset at the beginning of this
    assert sc.get_scaling_factor(m.c1) == 3
    assert sc.get_scaling_factor(m.c2) == 4

    assert m.scaling_factor[m.x] == 11
    assert m.scaling_factor[m.c1] == 3
    assert m.scaling_factor[m.z] == 19
    assert m.scaling_factor[m.c1] == 3
    assert m.scaling_factor[m.c2] == 4

    for i in range(1, 5):
        if i == 1:
            assert sc.get_scaling_factor(m.z[i]) == 29
            assert m.scaling_factor[m.z[i]] == 29
        else:
            assert sc.get_scaling_factor(m.z[i]) == 19
            assert m.scaling_factor[m.z[i]] == 19
        assert sc.get_scaling_factor(m.c2[i]) == 4
        assert m.scaling_factor[m.c2[i]] == 4

    sc.set_scaling_factor(m.x, 13, overwrite=True)
    sc.set_scaling_factor(m.ex, 17, overwrite=True)
    sc.set_scaling_factor(m.c1, 7, overwrite=True)
    # Make sure that overwrite properly overwrites subcomponents too
    sc.set_scaling_factor(m.z, 31, overwrite=True)
    sc.set_scaling_factor(m.c2, 23, overwrite=True)

    assert sc.get_scaling_factor(m.x) == 13
    assert sc.get_scaling_factor(m.ex) == 17
    assert sc.get_scaling_factor(m.z) == 31
    assert sc.get_scaling_factor(m.c1) == 7
    assert sc.get_scaling_factor(m.c2) == 23

    assert m.scaling_factor[m.x] == 13
    assert m.scaling_factor[m.ex] == 17
    assert m.scaling_factor[m.z] == 31
    assert m.scaling_factor[m.c1] == 7
    assert m.scaling_factor[m.c2] == 23

    for i in range(1, 5):
        assert sc.get_scaling_factor(m.z[i]) == 31
        assert sc.get_scaling_factor(m.c2[i]) == 23

        assert m.scaling_factor[m.z[i]] == 31
        assert m.scaling_factor[m.c2[i]] == 23

    # Test the ambiguous case where a parent component has a scaling factor, the subcomponents do not,
    # and set_scaling_factor is called with data_objects=True and overwrite=False. In this case,
    # avoid setting subcomponent scaling factors if the parent component has a scaling factor.
    sc.unset_scaling_factor(m.z, data_objects=True)
    sc.set_scaling_factor(m.z, 31, data_objects=False)
    sc.set_scaling_factor(m.z, 37, data_objects=True, overwrite=False)

    assert sc.get_scaling_factor(m.z) == 31
    assert m.scaling_factor[m.z] == 31
    for i in range(1, 5):
        assert sc.get_scaling_factor(m.z[i]) is None


@pytest.mark.unit
def test_set_and_get_scaling_factor():
    m = pyo.ConcreteModel()
    m.x = pyo.Var()
    m.z = pyo.Var([1, 2, 3, 4])
    m.c1 = pyo.Constraint(expr=0 == m.x)

    @m.Constraint([1, 2, 3, 4])
    def c2(b, i):
        return b.z[i] == 0

    m.ex = pyo.Expression(expr=m.x)

    # Make sure exception is raised for indexed components regardless of whether they have scaling factors set
    with pytest.raises(AttributeError) as exception:
        sc.set_and_get_scaling_factor(m.z, 10)
    assert (
        exception.value.args[0]
        == "Ambiguous which scaling factor to return for indexed component z."
    )

    sc.set_scaling_factor(m.z, 11)

    with pytest.raises(AttributeError) as exception:
        sc.set_and_get_scaling_factor(m.z, 10)
    assert (
        exception.value.args[0]
        == "Ambiguous which scaling factor to return for indexed component z."
    )

    for j in range(1, 5):
        assert sc.set_and_get_scaling_factor(m.z[j], 10) == 11
        assert sc.get_scaling_factor(m.z[j]) == 11

    sc.unset_scaling_factor(m.z[1])
    assert sc.set_and_get_scaling_factor(m.z[1], 10) == 10
    assert sc.get_scaling_factor(m.z[1]) == 10

    assert sc.set_and_get_scaling_factor(m.x, 7) == 7
    assert sc.get_scaling_factor(m.x) == 7
    assert sc.set_and_get_scaling_factor(m.x, 13) == 7
    assert sc.get_scaling_factor(m.x) == 7

    assert sc.set_and_get_scaling_factor(m.c1, 17) == 17
    assert sc.get_scaling_factor(m.c1) == 17
    assert sc.set_and_get_scaling_factor(m.c1, 19) == 17
    assert sc.get_scaling_factor(m.c1) == 17

    sc.set_scaling_factor(m.ex, 23)
    assert sc.set_and_get_scaling_factor(m.ex, 29) == 23
    assert sc.get_scaling_factor(m.ex) == 23
    sc.unset_scaling_factor(m.ex)
    assert sc.set_and_get_scaling_factor(m.ex, 29) == 29
    assert sc.get_scaling_factor(m.ex) == 29


@pytest.mark.unit
def test_find_badly_scaled_vars():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=1e6)
    m.y = pyo.Var(initialize=1e-8)
    m.z = pyo.Var(initialize=1e-20)
    m.b = pyo.Block()
    m.b.w = pyo.Var(initialize=1e10)

    a = [id(v) for v, sv in sc.badly_scaled_var_generator(m)]
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

    a = [id(v) for v, sv in sc.badly_scaled_var_generator(m)]
    assert id(m.x) not in a
    assert id(m.y) not in a
    assert id(m.b.w) in a
    assert id(m.z) not in a


@pytest.mark.unit
def test_find_unscaled_vars_and_constraints():
    m = pyo.ConcreteModel()
    m.b = pyo.Block()
    m.x = pyo.Var(initialize=1e6)
    m.y = pyo.Var(initialize=1e-8)
    m.z = pyo.Var(initialize=1e-20)
    m.c1 = pyo.Constraint(expr=m.x == 0)
    m.c2 = pyo.Constraint(expr=m.y == 0)
    m.b.w = pyo.Var([1, 2, 3], initialize=1e10)
    m.b.c1 = pyo.Constraint(expr=m.b.w[1] == 0)
    m.b.c2 = pyo.Constraint(expr=m.b.w[2] == 0)
    m.c3 = pyo.Constraint(expr=m.z == 0)

    sc.set_scaling_factor(m.x, 1)
    sc.set_scaling_factor(m.b.w[1], 2)
    sc.set_scaling_factor(m.c1, 1)
    sc.set_scaling_factor(m.b.c1, 1)
    sc.constraint_scaling_transform(m.c3, 1)

    a = [id(v) for v in sc.unscaled_variables_generator(m)]
    # Make sure we pick up the right variales
    assert id(m.x) not in a
    assert id(m.y) in a
    assert id(m.z) in a
    assert id(m.b.w[1]) not in a
    assert id(m.b.w[2]) in a
    assert id(m.b.w[3]) in a
    assert len(a) == 4  # make sure we didn't pick up any other random stuff

    a = [id(v) for v in sc.unscaled_constraints_generator(m)]
    assert id(m.c1) not in a
    assert id(m.b.c1) not in a
    assert id(m.c2) in a
    assert id(m.b.c2) in a
    assert id(m.c3) not in a
    assert len(a) == 2  # make sure we didn't pick up any other random stuff


class TestSingleConstraintScalingTransform:
    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.x = pyo.Var(initialize=500)
        m.c1 = pyo.Constraint(expr=m.x <= 1e3)
        m.c2 = pyo.Constraint(expr=m.x == 1e3)
        m.c3 = pyo.Constraint(expr=m.x >= 1e3)
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
            sc.constraint_scaling_transform(model.x, 1)

    @pytest.mark.unit
    def test_less_than_constraint(self, model):
        sc.constraint_scaling_transform(model.c1, 1e-3)
        assert model.c1.lower is None
        assert model.c1.body() == pytest.approx(model.x.value / 1e3)
        assert model.c1.upper.value == pytest.approx(1)
        assert sc.get_scaling_factor(model.c1) is None
        sc.constraint_scaling_transform_undo(model.c1)
        assert model.c1.lower is None
        assert model.c1.body() == pytest.approx(model.x.value)
        assert model.c1.upper.value == pytest.approx(1e3)

    @pytest.mark.unit
    def test_equality_constraint(self, model):
        sc.constraint_scaling_transform(model.c2, 1e-3)
        # Do transformation again to be sure that, despite being called twice in
        # a row, the transformation only happens once.
        sc.constraint_scaling_transform(model.c2, 1e-3)
        assert model.c2.lower.value == pytest.approx(1)
        assert model.c2.body() == pytest.approx(model.x.value / 1e3)
        assert model.c2.upper.value == pytest.approx(1)
        assert sc.get_constraint_transform_applied_scaling_factor(model.c2) == 1e-3

        # Check overwrite protection
        sc.constraint_scaling_transform(model.c2, 5, overwrite=False)
        assert model.c2.lower.value == pytest.approx(1)
        assert model.c2.body() == pytest.approx(model.x.value / 1e3)
        assert model.c2.upper.value == pytest.approx(1)
        assert sc.get_constraint_transform_applied_scaling_factor(model.c2) == 1e-3

        sc.constraint_scaling_transform_undo(model.c2)
        assert sc.get_constraint_transform_applied_scaling_factor(model.c2) is None
        assert model.c2.lower.value == pytest.approx(1e3)
        assert model.c2.body() == pytest.approx(model.x.value)
        assert model.c2.upper.value == pytest.approx(1e3)

    @pytest.mark.unit
    def test_greater_than_constraint(self, model):
        sc.constraint_scaling_transform(model.c3, 1e-3)
        assert sc.get_scaling_factor(model.c3) == None
        assert model.c3.lower.value == pytest.approx(1)
        assert model.c3.body() == pytest.approx(model.x.value / 1e3)
        assert model.c3.upper is None
        sc.constraint_scaling_transform_undo(model.c3)
        assert model.c3.lower.value == pytest.approx(1e3)
        assert model.c3.body() == pytest.approx(model.x.value)


class TestScaleSingleConstraint:
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


class TestScaleConstraintsPynumero:
    def model(self):
        m = pyo.ConcreteModel()
        x = m.x = pyo.Var(initialize=1e3)
        y = m.y = pyo.Var(initialize=1e6)
        z = m.z = pyo.Var(initialize=1e4)
        m.c1 = pyo.Constraint(expr=0 == -x * y + z)
        m.c2 = pyo.Constraint(expr=0 == 3 * x + 4 * y + 2 * z)
        m.c3 = pyo.Constraint(expr=0 <= z**3)
        return m

    @pytest.mark.unit
    def test_jacobian(self):
        """Make sure the Jacobian from Pynumero matches expectation.  This is
        mostly to ensure we understand the interface and catch if things change.
        """
        m = self.model()
        assert number_activated_objectives(m) == 0
        jac, jac_scaled, nlp = sc.constraint_autoscale_large_jac(m, no_scale=True)
        assert number_activated_objectives(m) == 0

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

        assert jac[c1_row, x_col] == pytest.approx(-1e6)
        assert jac[c1_row, y_col] == pytest.approx(-1e3)
        assert jac[c1_row, z_col] == pytest.approx(1)

        assert jac[c2_row, x_col] == pytest.approx(3)
        assert jac[c2_row, y_col] == pytest.approx(4)
        assert jac[c2_row, z_col] == pytest.approx(2)

        assert jac[c3_row, z_col] == pytest.approx(3e8)

        # Make sure scaling factors don't affect the result
        sc.set_scaling_factor(m.c1, 1e-6)
        sc.set_scaling_factor(m.x, 1e-3)
        sc.set_scaling_factor(m.y, 1e-6)
        sc.set_scaling_factor(m.z, 1e-4)
        jac, jac_scaled, nlp = sc.constraint_autoscale_large_jac(m, no_scale=True)
        assert jac[c1_row, x_col] == pytest.approx(-1e6)
        # Check the scaled jacobian calculation
        assert jac_scaled[c1_row, x_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, z_col] == pytest.approx(0.01)

    @pytest.mark.unit
    def test_scale_no_var_scale(self):
        """Make sure the Jacobian from Pynumero matches expectation.  This is
        mostly to ensure we understand the interface and catch if things change.
        """
        m = self.model()
        jac, jac_scaled, nlp = sc.constraint_autoscale_large_jac(m)

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

        assert jac_scaled[c1_row, x_col] == pytest.approx(-100)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-0.1)
        assert jac_scaled[c1_row, z_col] == pytest.approx(1e-4)
        assert m.scaling_factor[m.c1] == pytest.approx(1e-4)

        assert jac_scaled[c2_row, x_col] == pytest.approx(3)
        assert jac_scaled[c2_row, y_col] == pytest.approx(4)
        assert jac_scaled[c2_row, z_col] == pytest.approx(2)

        assert jac_scaled[c3_row, z_col] == pytest.approx(3e2)

    @pytest.mark.unit
    def test_scale_with_var_scale(self):
        """Make sure the Jacobian from Pynumero matches expectation.  This is
        mostly to ensure we understand the interface and catch if things change.
        """
        m = self.model()
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.x] = 1e-3
        m.scaling_factor[m.y] = 1e-6
        m.scaling_factor[m.z] = 1e-4

        jac, jac_scaled, nlp = sc.constraint_autoscale_large_jac(m)

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

        assert jac_scaled[c1_row, x_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-1000)
        assert jac_scaled[c1_row, z_col] == pytest.approx(1e-2)
        assert m.scaling_factor[m.c1] == pytest.approx(1e-6)

        assert jac_scaled[c2_row, x_col] == pytest.approx(0.075)
        assert jac_scaled[c2_row, y_col] == pytest.approx(100)
        assert jac_scaled[c2_row, z_col] == pytest.approx(0.5)
        assert m.scaling_factor[m.c2] == pytest.approx(2.5e-5)

        assert jac_scaled[c3_row, z_col] == pytest.approx(3e6)
        assert m.scaling_factor[m.c3] == pytest.approx(1e-6)

    @pytest.mark.unit
    def test_scale_with_ignore_var_scale(self):
        """Make sure the Jacobian from Pynumero matches expectation.  This is
        mostly to ensure we understand the interface and catch if things change.
        """
        m = self.model()
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.x] = 1e-3
        m.scaling_factor[m.y] = 1e-6
        m.scaling_factor[m.z] = 1e-4

        jac, jac_scaled, nlp = sc.constraint_autoscale_large_jac(
            m, ignore_variable_scaling=True
        )

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

        assert jac_scaled[c1_row, x_col] == pytest.approx(-100)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-0.1)
        assert jac_scaled[c1_row, z_col] == pytest.approx(1e-4)
        assert m.scaling_factor[m.c1] == pytest.approx(1e-4)

        assert jac_scaled[c2_row, x_col] == pytest.approx(3)
        assert jac_scaled[c2_row, y_col] == pytest.approx(4)
        assert jac_scaled[c2_row, z_col] == pytest.approx(2)
        assert m.scaling_factor[m.c2] == pytest.approx(1)

        assert jac_scaled[c3_row, z_col] == pytest.approx(3e2)
        assert m.scaling_factor[m.c3] == pytest.approx(1e-6)

    @pytest.mark.unit
    def test_scale_with_ignore_constraint_scale(self):
        """Make sure ignore_constraint_scaling ignores given scaling factors."""
        m = pyo.ConcreteModel()
        m.a = pyo.Var([1, 2], initialize=1)
        m.c = pyo.Constraint(expr=(0, m.a[1] + m.a[2], 1))
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.c] = 1e6

        jac, jac_scaled, nlp = sc.constraint_autoscale_large_jac(
            m, ignore_constraint_scaling=True
        )
        assert m.scaling_factor[m.c] == pytest.approx(1)

    @pytest.mark.unit
    def test_condition_number(self):
        """Calculate the condition number of the Jacobian"""
        m = self.model()
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.x] = 1e-3
        m.scaling_factor[m.y] = 1e-6
        m.scaling_factor[m.z] = 1e-4
        m.scaling_factor[m.c1] = 1e-6
        m.scaling_factor[m.c2] = 1e-6
        m.scaling_factor[m.c3] = 1e-12

        n = sc.jacobian_cond(m, scaled=True)
        assert n == pytest.approx(500, abs=200)
        n = sc.jacobian_cond(m, scaled=False)
        assert n == pytest.approx(7.5e7, abs=5e6)

    @pytest.mark.unit
    def test_scale_with_ignore_var_scale_constraint_scale(self):
        """Make sure the Jacobian from Pynumero matches expectation.  This is
        mostly to ensure we understand the interface and catch if things change.
        """
        m = self.model()
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.scaling_factor[m.c1] = 1e-6
        m.scaling_factor[m.x] = 1e-3
        m.scaling_factor[m.y] = 1e-6
        m.scaling_factor[m.z] = 1e-4

        jac, jac_scaled, nlp = sc.constraint_autoscale_large_jac(
            m, ignore_variable_scaling=True
        )

        c1_row = nlp._condata_to_idx[m.c1]
        c2_row = nlp._condata_to_idx[m.c2]
        c3_row = nlp._condata_to_idx[m.c3]
        x_col = nlp._vardata_to_idx[m.x]
        y_col = nlp._vardata_to_idx[m.y]
        z_col = nlp._vardata_to_idx[m.z]

        assert jac_scaled[c1_row, x_col] == pytest.approx(-1)
        assert jac_scaled[c1_row, y_col] == pytest.approx(-1e-3)
        assert jac_scaled[c1_row, z_col] == pytest.approx(1e-6)
        assert m.scaling_factor[m.c1] == pytest.approx(1e-6)

        assert jac_scaled[c2_row, x_col] == pytest.approx(3)
        assert jac_scaled[c2_row, y_col] == pytest.approx(4)
        assert jac_scaled[c2_row, z_col] == pytest.approx(2)
        assert m.scaling_factor[m.c2] == pytest.approx(1)

        assert jac_scaled[c3_row, z_col] == pytest.approx(3e2)
        assert m.scaling_factor[m.c1] == pytest.approx(1e-6)


class TestScaleConstraints:
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


class TestCacheVars:
    @pytest.mark.unit
    def test_cache_vars(self):
        m = pyo.ConcreteModel()
        val1 = 1
        val2 = 2
        m.v1 = pyo.Var(initialize=val1)
        m.v2 = pyo.Var(initialize=val2)

        varlist = [m.v1, m.v2]
        varset = ComponentSet(varlist)

        with sc.CacheVars(varlist) as cache:
            assert cache.cache == [1, 2]
            for var in cache.vars:
                assert var in varset
            m.v1.set_value(11)
            m.v2.set_value(12)

        assert m.v1.value == val1
        assert m.v2.value == val2


class TestFlattenedScalingAssignment:
    def set_initial_scaling_factors(self, m):
        scaling_factor = m.scaling_factor
        for var in m.z.values():
            scaling_factor[var] = 0.1
        for var in m.u.values():
            scaling_factor[var] = 0.5

    @pytest.fixture(scope="class")
    def model(self):
        m = pyo.ConcreteModel()
        m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        m.time = dae.ContinuousSet(bounds=(0, 1))
        m.space = dae.ContinuousSet(bounds=(0, 1))
        m.z = pyo.Var(m.time, m.space)
        m.dz = dae.DerivativeVar(m.z, wrt=m.time)
        m.y = pyo.Var(m.time, m.space)
        m.u = pyo.Var(m.time)
        m.s = pyo.Var()

        def de_rule(m, t, x):
            return m.dz[t, x] == 5 * m.y[t, x] - 10 * m.z[t, x]

        m.de = pyo.Constraint(m.time, m.space, rule=de_rule)

        def ae_rule(m, t, x):
            return m.y[t, x] == 4 + m.z[t, x] ** 3

        m.ae = pyo.Constraint(m.time, m.space, rule=ae_rule)

        x0 = m.space.first()

        def ue_rule(m, t):
            return m.z[t, x0] == 2 * m.u[t]

        m.ue = pyo.Constraint(m.time, rule=ue_rule)

        tf, xf = m.time.last(), m.space.last()

        def se_rule(m):
            return m.z[tf, xf] == m.s

        m.se = pyo.Constraint(rule=se_rule)

        return m

    @pytest.mark.unit
    def test_scale_2d(self, model):
        m = model
        scaling_factor = m.scaling_factor
        self.set_initial_scaling_factors(m)

        assignment = [
            (m.y, m.ae),
            (m.dz, m.de),
        ]
        scaler = sc.FlattenedScalingAssignment(scaling_factor, assignment, (0, 0))

        y = scaler.get_representative_data_object(m.y)
        assert y is m.y[0, 0]

        for var in scaler.varlist:
            scaler.calculate_variable_scaling_factor(var)

        for var in m.y.values():
            assert scaling_factor[var] == pytest.approx(1 / (4 + 10**3))
        nominal_y = 1 / scaling_factor[y]

        for var in m.dz.values():
            assert scaling_factor[var] == pytest.approx(1 / (5 * nominal_y - 100))

        for con in scaler.conlist:
            scaler.set_constraint_scaling_factor(con)

        for index, con in m.ae.items():
            var = m.y[index]
            assert scaling_factor[con] == scaling_factor[var]

        for index, con in m.de.items():
            var = m.dz[index]
            assert scaling_factor[con] == scaling_factor[var]

        scaler.set_derivative_factor_from_state(m.dz)
        for index, dvar in m.dz.items():
            z = m.z[index]
            assert scaling_factor[z] == scaling_factor[dvar]

    @pytest.mark.unit
    def test_scale_1d(self, model):
        m = model
        scaling_factor = m.scaling_factor
        self.set_initial_scaling_factors(m)

        assignment = [
            (m.u, m.ue),
        ]
        scaler = sc.FlattenedScalingAssignment(scaling_factor, assignment, 0)

        u = scaler.get_representative_data_object(m.u)
        assert u is m.u[0]

        for var in scaler.varlist:
            scaler.calculate_variable_scaling_factor(var)
        for index, var in m.u.items():
            z = m.z[index, 0]
            assert scaling_factor[var] == pytest.approx(2 * scaling_factor[z])

        for con in scaler.conlist:
            scaler.set_constraint_scaling_factor(con)
        for index, con in m.ue.items():
            u = m.u[index]
            assert scaling_factor[con] == scaling_factor[u]

    @pytest.mark.unit
    def test_scale_0d(self, model):
        m = model
        scaling_factor = m.scaling_factor
        self.set_initial_scaling_factors(m)

        assignment = [
            (m.s, m.se),
            (m.y[0, 0], m.ae[0, 0]),
        ]
        scaler = sc.FlattenedScalingAssignment(scaling_factor, assignment, None)

        s = scaler.get_representative_data_object(m.s)
        y = scaler.get_representative_data_object(m.y[0, 0])
        assert s is m.s
        assert y is m.y[0, 0]

        for var in scaler.varlist:
            scaler.calculate_variable_scaling_factor(var)
        tf, xf = m.time.last(), m.space.last()
        assert scaling_factor[s] == scaling_factor[m.z[tf, xf]]

        assert scaling_factor[y] == pytest.approx(1 / (4 + 10**3))


@pytest.mark.unit
def test_extreme_jacobian_rows_and_columns():
    m = pyo.ConcreteModel()

    m.I = pyo.Set(initialize=[i for i in range(5)])

    m.x = pyo.Var(m.I, initialize=1.0)

    diag = [1e7, 1, 10, 0.1, 1e-7]
    out = [1, 1, 1, 1, 1]

    @m.Constraint(m.I)
    def dummy_eqn(b, i):
        return out[i] == diag[i] * m.x[i]

    out = sc.extreme_jacobian_rows(m)
    assert type(out) == list
    assert len(out) == 2
    assert out[0][0] == pytest.approx(1e7)
    assert out[0][1] is m.dummy_eqn[0]
    assert out[1][0] == pytest.approx(1e-7)
    assert out[1][1] is m.dummy_eqn[4]

    out = sc.extreme_jacobian_rows(m, large=1e8)
    assert len(out) == 1
    assert out[0][0] == pytest.approx(1e-7)
    assert out[0][1] is m.dummy_eqn[4]

    out = sc.extreme_jacobian_rows(m, small=1e-8)
    assert len(out) == 1
    assert out[0][0] == pytest.approx(1e7)
    assert out[0][1] is m.dummy_eqn[0]

    out = sc.extreme_jacobian_rows(m, large=1e8, small=1e-8)
    assert len(out) == 0

    out = sc.extreme_jacobian_columns(m)
    assert type(out) == list
    assert len(out) == 2
    assert out[0][0] == pytest.approx(1e7)
    assert out[0][1] is m.x[0]
    assert out[1][0] == pytest.approx(1e-7)
    assert out[1][1] is m.x[4]

    out = sc.extreme_jacobian_columns(m, large=1e8)
    assert len(out) == 1
    assert out[0][0] == pytest.approx(1e-7)
    assert out[0][1] is m.x[4]

    out = sc.extreme_jacobian_columns(m, small=1e-8)
    assert len(out) == 1
    assert out[0][0] == pytest.approx(1e7)
    assert out[0][1] is m.x[0]

    out = sc.extreme_jacobian_columns(m, large=1e8, small=1e-8)
    assert len(out) == 0

    sc.constraint_scaling_transform(m.dummy_eqn[0], 1e-7)
    sc.constraint_scaling_transform(m.dummy_eqn[4], 1e7)

    out = sc.extreme_jacobian_rows(m)
    assert len(out) == 0

    out = sc.extreme_jacobian_columns(m)
    assert len(out) == 0

    # Even if scaling factors are ignored, transformed constraints
    # remain transformed
    out = sc.extreme_jacobian_columns(m, scaled=False)
    assert len(out) == 0

    out = sc.extreme_jacobian_rows(m, scaled=False)
    assert len(out) == 0

    sc.set_scaling_factor(m.x[1], 1e7)
    sc.set_scaling_factor(m.x[2], 1e-7)

    out = sc.extreme_jacobian_columns(m, scaled=False)
    assert len(out) == 0

    out = sc.extreme_jacobian_rows(m, scaled=False)
    assert len(out) == 0

    out = sc.extreme_jacobian_columns(m)
    assert len(out) == 2
    assert out[0][0] == pytest.approx(1e-7)
    assert out[0][1] is m.x[1]
    assert out[1][0] == pytest.approx(1e8)
    assert out[1][1] is m.x[2]

    out = sc.extreme_jacobian_rows(m)
    assert len(out) == 2
    assert out[0][0] == pytest.approx(1e-7)
    assert out[0][1] is m.dummy_eqn[1]
    assert out[1][0] == pytest.approx(1e8)
    assert out[1][1] is m.dummy_eqn[2]

    sc.constraint_scaling_transform(m.dummy_eqn[1], 1e7)
    sc.constraint_scaling_transform(m.dummy_eqn[2], 1e-8)

    out = sc.extreme_jacobian_columns(m, scaled=False)
    assert len(out) == 2
    assert out[0][0] == pytest.approx(1e7)
    assert out[0][1] is m.x[1]
    assert out[1][0] == pytest.approx(1e-7)
    assert out[1][1] is m.x[2]

    out = sc.extreme_jacobian_rows(m, scaled=False)
    assert len(out) == 2
    assert out[0][0] == pytest.approx(1e7)
    assert out[0][1] is m.dummy_eqn[1]
    assert out[1][0] == pytest.approx(1e-7)
    assert out[1][1] is m.dummy_eqn[2]

    out = sc.extreme_jacobian_columns(m)
    assert len(out) == 0

    out = sc.extreme_jacobian_rows(m)
    assert len(out) == 0


def discretization_tester(transformation_method, scheme, t_skip, continuity_eqns=False):
    """Function to avoid repeated code in testing scaling different discretization methods.

    Args:
        transformation_method: Discretization method to use (presently "dae.finite_difference" or "dae.collocation")
        scheme: Discretization scheme to use
        t_skip: Times to skip in testing scaling of discretization and/or continuity equations
        continuity_eqns: Whether to look for continuity equations while testing. Right now implementation is based on
            the "LAGRANGE-LEGENDRE" scheme. If additional methods with continuity equations are added to Pyomo that
            behave differently, this conditional might have to be updated

    Returns:
        None
    """

    def approx(expected):
        return pytest.approx(expected, rel=1e-10)

    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint(
        nfe=3, transformation_method=transformation_method, scheme=scheme
    )

    for t in m.t:
        sc.set_scaling_factor(m.y[t, 1], 1)
        sc.set_scaling_factor(m.y[t, 2], 10)
        sc.set_scaling_factor(m.y[t, 3], 0.1)
        sc.set_scaling_factor(m.y[t, 4], 1000)
        sc.set_scaling_factor(m.y[t, 5], 13)
        # Default scaling factor for y[6] is 1

    sc.scale_time_discretization_equations(m, m.t, 1 / 7)

    for t in m.t:
        assert sc.get_scaling_factor(m.ydot[t, 1]) == approx(7)
        assert sc.get_scaling_factor(m.ydot[t, 2]) == approx(70)
        assert sc.get_scaling_factor(m.ydot[t, 3]) == approx(0.7)
        assert sc.get_scaling_factor(m.ydot[t, 4]) == approx(7000)
        assert sc.get_scaling_factor(m.ydot[t, 5]) == approx(7 * 13)
        assert sc.get_scaling_factor(m.ydot[t, 6]) == approx(7)

        # No discretization equations in t_skip
        if t in t_skip:
            continue
        if not continuity_eqns:
            assert sc.get_constraint_transform_applied_scaling_factor(
                m.ydot_disc_eq[t, 1]
            ) == approx(7)
            assert sc.get_constraint_transform_applied_scaling_factor(
                m.ydot_disc_eq[t, 2]
            ) == approx(70)
            assert sc.get_constraint_transform_applied_scaling_factor(
                m.ydot_disc_eq[t, 3]
            ) == approx(0.7)
            assert sc.get_constraint_transform_applied_scaling_factor(
                m.ydot_disc_eq[t, 4]
            ) == approx(7000)
            assert sc.get_constraint_transform_applied_scaling_factor(
                m.ydot_disc_eq[t, 5]
            ) == approx(7 * 13)
            assert sc.get_constraint_transform_applied_scaling_factor(
                m.ydot_disc_eq[t, 6]
            ) == approx(7)
        else:
            # For Lagrange-Legendre, continuity equations exist only on boundaries of finite elements (except t=0)
            # while discretization equations exist only in the interior of finite elements.
            try:
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.ydot_disc_eq[t, 1]
                ) == approx(7)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.ydot_disc_eq[t, 2]
                ) == approx(70)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.ydot_disc_eq[t, 3]
                ) == approx(0.7)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.ydot_disc_eq[t, 4]
                ) == approx(7000)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.ydot_disc_eq[t, 5]
                ) == approx(7 * 13)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.ydot_disc_eq[t, 6]
                ) == approx(7)
            except KeyError:
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.y_t_cont_eq[t, 1]
                ) == approx(1)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.y_t_cont_eq[t, 2]
                ) == approx(10)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.y_t_cont_eq[t, 3]
                ) == approx(0.1)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.y_t_cont_eq[t, 4]
                ) == approx(1000)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.y_t_cont_eq[t, 5]
                ) == approx(13)
                assert sc.get_constraint_transform_applied_scaling_factor(
                    m.y_t_cont_eq[t, 6]
                ) == approx(1)


@pytest.mark.unit
def test_scaling_discretization_equations_backward():
    discretization_tester("dae.finite_difference", "BACKWARD", [0], False)


@pytest.mark.unit
def test_scaling_discretization_equations_forward():
    discretization_tester("dae.finite_difference", "FORWARD", [180], False)


@pytest.mark.unit
def test_scaling_discretization_equations_lagrange_radau():
    discretization_tester("dae.collocation", "LAGRANGE-RADAU", [0], False)


@pytest.mark.unit
def test_scaling_discretization_equations_lagrange_legendre():
    discretization_tester("dae.collocation", "LAGRANGE-LEGENDRE", [0], True)


@pytest.mark.unit
def test_correct_set_identification():
    # Suggested by Robby Parker. The original implementation of scale_time_discretization_equations
    # used a Python function that tested for set equality instead of identity. As a result, trouble
    # could happen if time and some other indexing set had the same members and time appeared later
    # than that other indexing set.
    def approx(x):
        return pytest.approx(x, 1e-12)

    m = pyo.ConcreteModel()
    m.time = dae.ContinuousSet(initialize=[0, 1, 2])
    m.space = dae.ContinuousSet(initialize=[0, 1, 2])

    m.x = pyo.Var(m.space, m.time, initialize=0)
    m.xdot = dae.DerivativeVar(m.x, wrt=m.time)

    @m.Constraint(m.space, m.time)
    def diff_eqn(b, z, t):
        return b.xdot[z, t] == -b.x[z, t]

    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m, nfe=2, wrt=m.time, scheme="BACKWARD"
    )
    for i in range(3):
        sc.set_scaling_factor(m.x[0, i], 2)
        sc.set_scaling_factor(m.x[1, i], 3)
        sc.set_scaling_factor(m.x[2, i], 5)

    sc.scale_time_discretization_equations(m, m.time, 10)

    for i in range(3):
        assert sc.get_scaling_factor(m.xdot[0, i]) == approx(0.2)
        assert sc.get_scaling_factor(m.xdot[1, i]) == approx(0.3)
        assert sc.get_scaling_factor(m.xdot[2, i]) == approx(0.5)

    for i in range(1, 3):
        assert sc.get_constraint_transform_applied_scaling_factor(
            m.xdot_disc_eq[0, i]
        ) == approx(0.2)
        assert sc.get_constraint_transform_applied_scaling_factor(
            m.xdot_disc_eq[1, i]
        ) == approx(0.3)
        assert sc.get_constraint_transform_applied_scaling_factor(
            m.xdot_disc_eq[2, i]
        ) == approx(0.5)


@pytest.mark.unit
def test_indexed_blocks():
    m = pyo.ConcreteModel()
    m.time = dae.ContinuousSet(initialize=[0, 1, 2])
