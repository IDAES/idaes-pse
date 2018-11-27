##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for math util methods.
"""

import pytest
from pyomo.environ import Block, ConcreteModel,  Constraint, \
                            Set, SolverFactory, Var, value
from pyomo.network import Port
from idaes.core.util.initialization import (evaluate_variable_from_constraint,
                                            solve_indexed_blocks,
                                            fix_port, unfix_port)

__author__ = "Andrew Lee"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None


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


def test_fix_port_no_args():
    m = ConcreteModel()
    m.v = Var(initialize=1.0)
    m.c = Port(noruleinit=True)
    m.c.add(m.v, "v")

    fix_port(m.c, "v")

    assert m.v.fixed is True
    assert m.v.value == 1.0


def test_fix_port_value():
    m = ConcreteModel()
    m.v = Var(initialize=1.0)
    m.c = Port(noruleinit=True)
    m.c.add(m.v, "v")

    fix_port(m.c, "v", value=2.0)

    assert m.v.fixed is True
    assert m.v.value == 2.0


def test_fix_port_indexed():
    m = ConcreteModel()
    m.v1 = Var(initialize=1.0)
    m.v2 = Var(initialize=10.0)
    m.c = Port([1, 2], noruleinit=True)
    m.c[1].add(m.v1, "v1")
    m.c[2].add(m.v2, "v2")

    fix_port(m.c, "v1", port_idx=[1])

    assert m.v1.fixed is True
    assert m.v2.fixed is False
    assert m.v1.value == 1.0
    assert m.v2.value == 10.0


def test_fix_port_indexed_value():
    m = ConcreteModel()
    m.v1 = Var(initialize=1.0)
    m.v2 = Var(initialize=10.0)
    m.c = Port([1, 2], noruleinit=True)
    m.c[1].add(m.v1, "v1")
    m.c[2].add(m.v2, "v2")

    fix_port(m.c, "v2", port_idx=[2], value=30.0)

    assert m.v1.fixed is False
    assert m.v2.fixed is True
    assert m.v1.value == 1.0
    assert m.v2.value == 30.0


def test_fix_port_indexed_var():
    m = ConcreteModel()
    m.v = Var([1, 2], initialize=1.0)
    m.c = Port(noruleinit=True)
    m.c.add(m.v, "v")

    fix_port(m.c, "v", comp=1)

    assert m.v[1].fixed is True
    assert m.v[2].fixed is False
    assert m.v[1].value == 1.0
    assert m.v[2].value == 1.0


def test_fix_port_indexed_var_value():
    m = ConcreteModel()
    m.v = Var([1, 2], initialize=1.0)
    m.c = Port(noruleinit=True)
    m.c.add(m.v, "v")

    fix_port(m.c, "v", comp=1, value=2.0)

    assert m.v[1].fixed is True
    assert m.v[2].fixed is False
    assert m.v[1].value == 2.0
    assert m.v[2].value == 1.0


def test_fix_port_both_indexed():
    m = ConcreteModel()
    m.v1 = Var([1, 2], initialize=1.0)
    m.v2 = Var(initialize=1.0)
    m.c = Port([1, 2], noruleinit=True)
    m.c[1].add(m.v1, "v1")
    m.c[2].add(m.v2, "v2")

    fix_port(m.c, "v1", port_idx=[1], comp=1)

    assert m.v1[1].fixed is True
    assert m.v1[2].fixed is False
    assert m.v2.fixed is False
    assert m.v1[1].value == 1.0
    assert m.v1[2].value == 1.0
    assert m.v2.value == 1.0


def test_fix_port_both_indexed_value():
    m = ConcreteModel()
    m.v1 = Var([1, 2], initialize=1.0)
    m.v2 = Var(initialize=1.0)
    m.c = Port([1, 2], noruleinit=True)
    m.c[1].add(m.v1, "v1")
    m.c[2].add(m.v2, "v2")

    fix_port(m.c, "v1", port_idx=[1], comp=1, value=10.0)

    assert m.v1[1].fixed is True
    assert m.v1[2].fixed is False
    assert m.v2.fixed is False
    assert m.v1[1].value == 10.0
    assert m.v1[2].value == 1.0
    assert m.v2.value == 1.0


def test_unfix_port_no_args():
    m = ConcreteModel()
    m.v = Var(initialize=1.0)
    m.v.fix()
    m.c = Port(noruleinit=True)
    m.c.add(m.v, "v")

    unfix_port(m.c, "v")

    assert m.v.fixed is False
    assert m.v.value == 1.0


def test_unfix_port_indexed():
    m = ConcreteModel()
    m.v1 = Var(initialize=1.0)
    m.v2 = Var(initialize=10.0)
    m.v1.fix()
    m.v2.fix()
    m.c = Port([1, 2], noruleinit=True)
    m.c[1].add(m.v1, "v1")
    m.c[2].add(m.v2, "v2")

    unfix_port(m.c, "v1", port_idx=[1])

    assert m.v1.fixed is False
    assert m.v2.fixed is True
    assert m.v1.value == 1.0
    assert m.v2.value == 10.0


def test_unfix_port_indexed_var():
    m = ConcreteModel()
    m.v = Var([1, 2], initialize=1.0)
    m.v.fix()
    m.c = Port(noruleinit=True)
    m.c.add(m.v, "v")

    unfix_port(m.c, "v", comp=1)

    assert m.v[1].fixed is False
    assert m.v[2].fixed is True
    assert m.v[1].value == 1.0
    assert m.v[2].value == 1.0


def test_unfix_port_both_indexed():
    m = ConcreteModel()
    m.v1 = Var([1, 2], initialize=1.0)
    m.v2 = Var(initialize=1.0)
    m.v1.fix()
    m.v2.fix()
    m.c = Port([1, 2], noruleinit=True)
    m.c[1].add(m.v1, "v1")
    m.c[2].add(m.v2, "v2")

    unfix_port(m.c, "v1", port_idx=[1], comp=1)

    assert m.v1[1].fixed is False
    assert m.v1[2].fixed is True
    assert m.v2.fixed is True
    assert m.v1[1].value == 1.0
    assert m.v1[2].value == 1.0
    assert m.v2.value == 1.0


def test_evaluate_variable_from_constraint():
    m = ConcreteModel()
    m.a = Var(initialize=0.0)

    m.c1 = Constraint(expr=0 == m.a - 10.0)

    assert evaluate_variable_from_constraint(m.a, m.c1) == 10.0
    assert value(m.a) == 10.0

    m.c2 = Constraint(expr=20.0 == m.a)

    assert evaluate_variable_from_constraint(m.a, m.c2) == 20.0
    assert value(m.a) == 20.0

    m.c3 = Constraint(expr=0 == m.a + 30.0)

    assert evaluate_variable_from_constraint(m.a, m.c3) == -30.0
    assert value(m.a) == -30.0

    m.c4 = Constraint(expr=-40.0 == m.a)

    assert evaluate_variable_from_constraint(m.a, m.c4) == -40.0
    assert value(m.a) == -40.0

    m.c5 = Constraint(expr=-60.0 == m.a + 10.0)

    assert evaluate_variable_from_constraint(m.a, m.c5) == -70.0
    assert value(m.a) == -70.0

    m.c6 = Constraint(expr=-70.0 >= m.a + 10.0)

    assert evaluate_variable_from_constraint(m.a, m.c6) == -80.0
    assert value(m.a) == -80.0
