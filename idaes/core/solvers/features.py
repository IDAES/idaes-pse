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

from functools import lru_cache
import pyomo.environ as pyo
from pyomo.common.errors import ApplicationError

def lp():
    """This provides a simple LP model for solver testing.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for m.x
    """
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
    """This provides a simple MILP model for solver testing.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for m.x
    """
    m = pyo.ConcreteModel()
    m.x = pyo.Var(domain=pyo.Integers, initialize=3)
    m.y = pyo.Var(domain=pyo.Integers, initialize=3)
    m.c1 = pyo.Constraint(expr=m.x >= 0.5)
    m.c2 = pyo.Constraint(expr=m.y >= 1.5)
    m.c3 = pyo.Constraint(expr=m.x <= 5)
    m.c4 = pyo.Constraint(expr=m.y <= 5)
    m.obj = pyo.Objective(expr=m.x + m.y)
    return m, 1

def nle():
    """This provides a simple system of nonlinear equations model for solver
    testing.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for m.x
    """
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=-0.1)
    m.eq1 = pyo.Constraint(expr=m.x**2 == 1)
    return m, 1

def nlp():
    """This provides a simple NLP model for solver testing.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for m.x
    """
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=-0.1)
    m.y = pyo.Var(initialize=1)
    m.c = pyo.Constraint(expr=m.x >= 1)
    m.obj = pyo.Objective(expr=m.x**2 + m.y**2)
    return m, 1

def minlp():
    """This provides a simple MINLP model for solver testing.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for m.x and m.i
    """
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=-0.1)
    m.y = pyo.Var(initialize=1)
    m.i = pyo.Var(domain=pyo.Binary, initialize=0)
    m.c = pyo.Constraint(expr=m.x >= 1)
    m.obj = pyo.Objective(
        expr=m.i * (m.x**2 + m.y**2) + (1 - m.i) * 4 *(m.x**2 + m.y**2))
    return m, 1, 1

@lru_cache(maxsize=10)
def ipopt_has_linear_solver(linear_solver):
    """Check if IPOPT can use the specified linear solver.

    Args:
        linear_solver (str): linear solver in {"ma27", "ma57", "ma77", "ma86",
            "ma97", "pardiso", "pardisomkl", "spral", "wsmp", "mumps"} or other
            custom solver.

    Returns:
        (bool): True if Ipopt is available with the specified linear solver or False
            if either Ipopt or the linear solver is not available.
    """
    m, x = nlp()
    solver = pyo.SolverFactory('ipopt', options={"linear_solver": linear_solver})
    try:
        solver.solve(m)
    except ApplicationError:
        return False
    try:
        assert abs(x - pyo.value(m.x)) < 1e-8
    except AssertionError:
        return False # solver mysteriously doesn't work right
    return True
