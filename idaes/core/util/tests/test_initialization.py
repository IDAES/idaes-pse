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
Tests for math util methods.
"""

import pytest
from pyomo.environ import Block, ConcreteModel,  Constraint, Expression, \
                            Set, SolverFactory, Var, value
from pyomo.network import Arc, Port
from idaes.core.util.initialization import (propagate_state,
                                            solve_indexed_blocks)

__author__ = "Andrew Lee"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


def test_propagate_state():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    # Set values on first block
    m.b1.v1.value = 10
    m.b1.v2[1].value = 20
    m.b1.v2[2].value = 30

    # Make sure vars in block 2 haven't been changed accidentally
    assert m.b2.v1.value is None
    assert m.b2.v2[1].value is None
    assert m.b2.v2[2].value is None

    propagate_state(m.s1)

    # Check that values were propagated correctly
    assert m.b2.v1.value == m.b1.v1.value
    assert m.b2.v2[1].value == m.b1.v2[1].value
    assert m.b2.v2[2].value == m.b1.v2[2].value

    assert m.b1.v1.fixed is False
    assert m.b1.v2[1].fixed is False
    assert m.b1.v2[2].fixed is False
    assert m.b2.v1.fixed is False
    assert m.b2.v2[1].fixed is False
    assert m.b2.v2[2].fixed is False


def test_propagate_state_reverse():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    # Test reverse propogation - set values on second block
    m.b2.v1.value = 100
    m.b2.v2[1].value = 200
    m.b2.v2[2].value = 300

    # Make sure vars in block 1 haven't been changed accidentally
    assert m.b1.v1.value is None
    assert m.b1.v2[1].value is None
    assert m.b1.v2[2].value is None

    propagate_state(m.s1, direction="backward")

    # Check that values were propagated correctly
    assert m.b2.v1.value == m.b1.v1.value
    assert m.b2.v2[1].value == m.b1.v2[1].value
    assert m.b2.v2[2].value == m.b1.v2[2].value

    assert m.b1.v1.fixed is False
    assert m.b1.v2[1].fixed is False
    assert m.b1.v2[2].fixed is False
    assert m.b2.v1.fixed is False
    assert m.b2.v2[1].fixed is False
    assert m.b2.v2[2].fixed is False


def test_propagate_state_fixed():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    # Set values on first block
    m.b1.v1.value = 10
    m.b1.v2[1].value = 20
    m.b1.v2[2].value = 30

    # Make sure vars in block 2 haven't been changed accidentally
    assert m.b2.v1.value is None
    assert m.b2.v2[1].value is None
    assert m.b2.v2[2].value is None

    # Fix v1 in block 2
    m.b2.v1.fix(500)

    propagate_state(m.s1)

    # Check that values were propagated correctly
    assert m.b2.v1.value == 500
    assert m.b2.v2[1].value == m.b1.v2[1].value
    assert m.b2.v2[2].value == m.b1.v2[2].value

    assert m.b1.v1.fixed is False
    assert m.b1.v2[1].fixed is False
    assert m.b1.v2[2].fixed is False
    assert m.b2.v1.fixed is True
    assert m.b2.v2[1].fixed is False
    assert m.b2.v2[2].fixed is False


def test_propagate_state_Expression():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.e = Expression(expr=b.v1)

        b.p = Port()
        b.p.add(b.e, "E")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    with pytest.raises(TypeError):
        propagate_state(m.s1)


def test_propagate_state_invalid_stream():
    m = ConcreteModel()

    with pytest.raises(TypeError):
        propagate_state(m)


def test_propagate_state_invalid_direction():
    m = ConcreteModel()

    def block_rule(b):
        b.s = Set(initialize=[1, 2])
        b.v1 = Var()
        b.v2 = Var(b.s)

        b.p = Port()
        b.p.add(b.v1, "V1")
        b.p.add(b.v2, "V2")
        return

    m.b1 = Block(rule=block_rule)
    m.b2 = Block(rule=block_rule)

    m.s1 = Arc(source=m.b1.p, destination=m.b2.p)

    with pytest.raises(ValueError):
        propagate_state(m.s1, direction="foo")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve_indexed_block_list():
    # Create an indexed block and try to solve it
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3])

    def block_rule(b, x):
        b.v = Var(initialize=1.0)
        b.c = Constraint(expr=b.v == 2.0)
    m.b = Block(m.s, rule=block_rule)

    solve_indexed_blocks(solver=solver, blocks=[m.b])

    for i in m.s:
        assert value(m.b[i].v == 2.0)


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve_indexed_block_IndexedBlock():
    # Create an indexed block and try to solve it
    m = ConcreteModel()
    m.s = Set(initialize=[1, 2, 3])

    def block_rule(b, x):
        b.v = Var(initialize=1.0)
        b.c = Constraint(expr=b.v == 2.0)
    m.b = Block(m.s, rule=block_rule)

    solve_indexed_blocks(solver=solver, blocks=m.b)

    for i in m.s:
        assert value(m.b[i].v == 2.0)


def test_solve_indexed_block_error():
    # Try solve_indexed_block on non-block object
    with pytest.raises(TypeError):
        solve_indexed_blocks(solver=None, blocks=[1, 2, 3])
