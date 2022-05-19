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
import pyomo.dae as pyodae
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
    m.eq1 = pyo.Constraint(expr=m.x**3 == 1)
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
        expr=m.i * (m.x**2 + m.y**2) + (1 - m.i) * 4 * (m.x**2 + m.y**2)
    )
    return m, 1, 1


def dae(nfe=1):
    """This provides a DAE model for solver testing.

    The problem and expected result are from the problem given here:
    https://archimede.dm.uniba.it/~testset/report/chemakzo.pdf.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for y[1] to y[5] and y6
    """
    model = pyo.ConcreteModel(name="chemakzo")

    # Set problem parameter values
    model.k = pyo.Param([1, 2, 3, 4], initialize={1: 18.7, 2: 0.58, 3: 0.09, 4: 0.42})
    model.Ke = pyo.Param(initialize=34.4)
    model.klA = pyo.Param(initialize=3.3)
    model.Ks = pyo.Param(initialize=115.83)
    model.pCO2 = pyo.Param(initialize=0.9)
    model.H = pyo.Param(initialize=737)

    # Problem variables ydot = dy/dt,
    #    (dy6/dt is not explicitly in the equations, so only 5 ydots i.e.
    #    y6 is an algebraic variable and y1 to y5 are differential variables)
    model.t = pyodae.ContinuousSet(bounds=(0, 180))
    model.y = pyo.Var(model.t, [1, 2, 3, 4, 5], initialize=1.0)  #
    model.y6 = pyo.Var(model.t, initialize=1.0)  #
    model.ydot = pyodae.DerivativeVar(model.y, wrt=model.t)  # dy/dt
    model.r = pyo.Var(model.t, [1, 2, 3, 4, 5], initialize=1.0)
    model.Fin = pyo.Var(model.t, initialize=1.0)

    # Equations
    @model.Constraint(model.t)
    def eq_ydot1(b, t):
        return b.ydot[t, 1] == -2.0 * b.r[t, 1] + b.r[t, 2] - b.r[t, 3] - b.r[t, 4]

    @model.Constraint(model.t)
    def eq_ydot2(b, t):
        return b.ydot[t, 2] == -0.5 * b.r[t, 1] - b.r[t, 4] - 0.5 * b.r[t, 5] + b.Fin[t]

    @model.Constraint(model.t)
    def eq_ydot3(b, t):
        return b.ydot[t, 3] == b.r[t, 1] - b.r[t, 2] + b.r[t, 3]

    @model.Constraint(model.t)
    def eq_ydot4(b, t):
        return b.ydot[t, 4] == -b.r[t, 2] + b.r[t, 3] - 2.0 * b.r[t, 4]

    @model.Constraint(model.t)
    def eq_ydot5(b, t):
        return b.ydot[t, 5] == b.r[t, 2] - b.r[t, 3] + b.r[t, 5]

    @model.Constraint(model.t)
    def eq_y6(b, t):
        return 0 == b.Ks * b.y[t, 1] * b.y[t, 4] - b.y6[t]

    @model.Constraint(model.t)
    def eq_r1(b, t):
        return b.r[t, 1] == b.k[1] * b.y[t, 1] ** 4 * b.y[t, 2] ** 0.5

    @model.Constraint(model.t)
    def eq_r2(b, t):
        return b.r[t, 2] == b.k[2] * b.y[t, 3] * b.y[t, 4]

    @model.Constraint(model.t)
    def eq_r3(b, t):
        return b.r[t, 3] == b.k[2] / b.Ke * b.y[t, 1] * b.y[t, 5]

    @model.Constraint(model.t)
    def eq_r4(b, t):
        return b.r[t, 4] == b.k[3] * b.y[t, 1] * b.y[t, 4] ** 2

    @model.Constraint(model.t)
    def eq_r5(b, t):
        return b.r[t, 5] == b.k[4] * b.y6[t] ** 2 * b.y[t, 2] ** 0.5

    @model.Constraint(model.t)
    def eq_Fin(b, t):
        return b.Fin[t] == b.klA * (b.pCO2 / b.H - b.y[t, 2])

    # Set initial conditions and solve initial from the values of differential
    # variables.
    y0 = {1: 0.444, 2: 0.00123, 3: 0.0, 4: 0.007, 5: 0.0}  # initial differential vars
    for i in y0:
        model.y[0, i].fix(y0[i])

    discretizer = pyo.TransformationFactory("dae.finite_difference")
    discretizer.apply_to(model, nfe=nfe, scheme="BACKWARD")

    return (
        model,
        0.1150794920661702,
        0.1203831471567715e-2,
        0.1611562887407974,
        0.3656156421249283e-3,
        0.1708010885264404e-1,
        0.4873531310307455e-2,
    )


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
    solver = pyo.SolverFactory("ipopt", options={"linear_solver": linear_solver})
    try:
        solver.solve(m)
    except ApplicationError:
        return False
    try:
        assert abs(x - pyo.value(m.x)) < 1e-8
    except AssertionError:
        return False  # solver mysteriously doesn't work right
    return True
