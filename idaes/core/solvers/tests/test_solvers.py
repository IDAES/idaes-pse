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
from pyomo.environ import SolverFactory
import pyomo.environ as pyo
import pytest
import idaes.core.plugins


def lp():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=3)
    m.y = pyo.Var(initialize=3)
    m.c1 = pyo.Constraint(expr=m.x >= 1)
    m.c2 = pyo.Constraint(expr=m.y >= 2)
    m.c3 = pyo.Constraint(expr=m.x <= 5)
    m.c4 = pyo.Constraint(expr=m.y <= 5)
    m.obj = pyo.Objective(expr=m.x + m.y)
    return m, 1

def milp():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(domain=pyo.Integers, initialize=3)
    m.y = pyo.Var(domain=pyo.Integers, initialize=3)
    m.c1 = pyo.Constraint(expr=m.x >= 0.5)
    m.c2 = pyo.Constraint(expr=m.y >= 1.5)
    m.c3 = pyo.Constraint(expr=m.x <= 5)
    m.c4 = pyo.Constraint(expr=m.y <= 5)
    m.obj = pyo.Objective(expr=m.x + m.y)
    return m, 1

def nlp():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=-0.1)
    m.y = pyo.Var(initialize=1)
    m.c = pyo.Constraint(expr=m.x >= 1)
    m.obj = pyo.Objective(expr=m.x**2 + m.y**2)
    return m, 1

def minlp():
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=-0.1)
    m.y = pyo.Var(initialize=1)
    m.i = pyo.Var(domain=pyo.Binary, initialize=1)
    m.c = pyo.Constraint(expr=m.x >= 1)
    m.obj = pyo.Objective(
        expr=m.i * (m.x**2 + m.y**2) + (1 - m.i) * 4 *(m.x**2 + m.y**2))
    return m, 1

@pytest.mark.unit
def test_couenne_available():
    if not SolverFactory('couenne').available():
        raise Exception("Could not find couenne.")

@pytest.mark.unit
def test_couenne_available():
    if not SolverFactory('bonmin').available():
        raise Exception("Could not find bonmin.")

@pytest.mark.unit
def test_sipopt_available():
    if not SolverFactory('ipopt_sens').available():
        raise Exception("Could not find ipopt_sens.")

@pytest.mark.unit
def test_cbc_available():
    if not SolverFactory('cbc').available():
        raise Exception("Could not find cbc.")

@pytest.mark.unit
def test_sipopt_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = nlp()
    solver = SolverFactory('ipopt_sens')
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)

@pytest.mark.unit
def test_bonmin_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = minlp()
    solver = SolverFactory('bonmin')
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)

@pytest.mark.unit
def test_couenne_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = minlp()
    solver = SolverFactory('couenne')
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)

@pytest.mark.unit
def test_cbc_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = milp()
    solver = SolverFactory('cbc')
    solver.solve(m)
    m.display()
    assert pytest.approx(x) == pyo.value(m.x)
