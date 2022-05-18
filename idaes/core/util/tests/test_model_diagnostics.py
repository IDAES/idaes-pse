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
This module contains model diagnostic utility functions for use in IDAES (Pyomo) models.
"""

import pytest

# Need to update
import pyomo.environ as pyo
import numpy as np
import idaes.core.util.scaling as iscale

# TODO: Add pyomo.dae test case
"""
from pyomo.environ import TransformationFactory
from pyomo.dae import ContinuousSet, DerivativeVar
"""

# Need to update
from idaes.core.util.model_diagnostics import DegeneracyHunter

__author__ = "Alex Dowling, Douglas Allan"


@pytest.fixture()
def dummy_problem():
    m = pyo.ConcreteModel()

    m.I = pyo.Set(initialize=[i for i in range(5)])

    m.x = pyo.Var(m.I, initialize=1.0)

    diag = [100, 1, 10, 0.1, 5]
    out = [1, 1, 1, 1, 1]

    @m.Constraint(m.I)
    def dummy_eqn(b, i):
        return out[i] == diag[i] * m.x[i]

    m.obj = pyo.Objective(expr=0)
    return m


@pytest.fixture()
def u_exp():
    return np.array(
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]]
    )


@pytest.fixture()
def s_exp():
    return np.array([0.1, 1, 5, 10])


@pytest.fixture()
def v_exp():
    return np.array(
        [[0, 0, 0, 1, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]]
    ).T


@pytest.mark.unit
def test_dense_svd(dummy_problem, u_exp, s_exp, v_exp):
    m = dummy_problem
    dh = DegeneracyHunter(m)
    dh.svd_analysis(dense=True)
    assert dh.s == pytest.approx(s_exp, 1e-6, abs=1e-10)
    assert u_exp == pytest.approx(np.abs(dh.u), 1e-6, abs=1e-10)
    assert v_exp == pytest.approx(np.abs(dh.v), 1e-6, abs=1e-10)


@pytest.mark.unit
def test_sparse_svd(dummy_problem, u_exp, s_exp, v_exp):
    m = dummy_problem
    dh = DegeneracyHunter(m)
    dh.svd_analysis()
    assert dh.s == pytest.approx(s_exp, 1e-6, abs=1e-10)
    assert u_exp == pytest.approx(np.abs(dh.u), 1e-6, abs=1e-10)
    assert v_exp == pytest.approx(np.abs(dh.v), 1e-6, abs=1e-10)


@pytest.mark.unit
def test_scaling(dummy_problem, u_exp, s_exp, v_exp):
    ssf = iscale.set_scaling_factor
    cst = iscale.constraint_scaling_transform
    m = dummy_problem
    cst(m.dummy_eqn[0], 1e-2)
    ssf(m.x[1], 1e3)
    cst(m.dummy_eqn[1], 1e3)
    cst(m.dummy_eqn[2], 1e-1)
    ssf(m.x[3], 1e-3)
    cst(m.dummy_eqn[3], 1e2)
    cst(m.dummy_eqn[4], 0.2)

    dh = DegeneracyHunter(m)
    dh.svd_analysis(dense=False)
    assert dh.s == pytest.approx(np.ones((4,)), 1e-6)
    dh.svd_analysis(dense=True)
    assert dh.s == pytest.approx(np.ones((4,)), 1e-6)


@pytest.mark.unit
def test_underdetermined_variables_and_constraints(dummy_problem, capsys):
    m = dummy_problem
    dh = DegeneracyHunter(m)

    dh.svd_analysis(n_sv=3)
    captured = capsys.readouterr()
    assert captured.out == "Computing the 3 smallest singular value(s)\n"

    dh.underdetermined_variables_and_constraints()
    captured = capsys.readouterr()
    assert captured.out == (
        "Column:    Variable\n3: x[3]\n\nRow:    " "Constraint\n3: dummy_eqn[3]\n"
    )

    dh.underdetermined_variables_and_constraints(n_calc=3)
    captured = capsys.readouterr()
    assert captured.out == (
        "Column:    Variable\n4: x[4]\n\nRow:    " "Constraint\n4: dummy_eqn[4]\n"
    )
    with pytest.raises(
        ValueError,
        match="User wanted constraints and variables associated "
        "with the 4-th smallest singular value, "
        "but only 3 small singular values have been "
        "calculated. Run svd_analysis again and specify "
        "n_sv>=4.",
    ):
        dh.underdetermined_variables_and_constraints(n_calc=4)

    dh.underdetermined_variables_and_constraints(tol=2)
    captured = capsys.readouterr()
    assert captured.out == ("Column:    Variable\n\nRow:    Constraint\n")


@pytest.mark.unit
def test_underdetermined_calls_svd_analysis(dummy_problem, capsys):
    m = dummy_problem
    dh = DegeneracyHunter(m)

    dh.underdetermined_variables_and_constraints(n_calc=1)
    captured = capsys.readouterr()
    assert captured.out == (
        "Computing the 4 smallest singular value(s)\n"
        "Column:    Variable\n3: x[3]\n\nRow:    "
        "Constraint\n3: dummy_eqn[3]\n"
    )


@pytest.mark.unit
def test_sv_value_error(dummy_problem):
    m = dummy_problem
    dh = DegeneracyHunter(m)
    with pytest.raises(
        ValueError,
        match="For a 5 by 5 system, svd_analysis can compute at most 4 "
        "singular values and vectors, but 5 were called for.",
    ):
        dh.svd_analysis(n_sv=5)
    with pytest.raises(ValueError, match="Nonsense value for n_sv=-1 received."):
        dh.svd_analysis(n_sv=-1)


@pytest.mark.unit
def test_single_eq_error(capsys):
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=1)
    m.con = pyo.Constraint(expr=(2 * m.x == 1))
    m.obj = pyo.Objective(expr=0)

    dh = DegeneracyHunter(m)
    with pytest.raises(
        ValueError,
        match="Model needs at least 2 equality constraints to " "perform svd_analysis.",
    ):
        dh.svd_analysis()

    dh.check_rank_equality_constraints()
    captured = capsys.readouterr()
    assert captured.out == (
        "\nChecking rank of Jacobian of equality constraints...\n"
        "Model contains 1 equality constraints and 1 variables.\n"
        "Only singular value: 2.0\n"
    )


# This was from
# @pytest.fixture()
def problem1():
    m = pyo.ConcreteModel()

    m.I = pyo.Set(initialize=[i for i in range(5)])

    m.x = pyo.Var(m.I, bounds=(-10, 10), initialize=1.0)

    m.con1 = pyo.Constraint(expr=m.x[0] + m.x[1] - m.x[3] >= 10)
    m.con2 = pyo.Constraint(expr=m.x[0] * m.x[3] + m.x[1] >= 0)
    m.con3 = pyo.Constraint(expr=m.x[4] * m.x[3] + m.x[0] * m.x[3] - m.x[4] == 0)

    m.obj = pyo.Objective(expr=sum(m.x[i] ** 2 for i in m.I))

    return m


def example2(with_degenerate_constraint=True):
    """Create the Pyomo model for Example 2

    Arguments:
        with_degenerate_constraint: Boolean, if True, include the redundant linear constraint

    Returns:
        m2: Pyomo model
    """

    m2 = pyo.ConcreteModel()

    m2.I = pyo.Set(initialize=[i for i in range(1, 4)])

    m2.x = pyo.Var(m2.I, bounds=(0, 5), initialize=1.0)

    m2.con1 = pyo.Constraint(expr=m2.x[1] + m2.x[2] >= 1)
    m2.con2 = pyo.Constraint(expr=m2.x[1] + m2.x[2] + m2.x[3] == 1)
    m2.con3 = pyo.Constraint(expr=m2.x[2] - 2 * m2.x[3] <= 1)
    m2.con4 = pyo.Constraint(expr=m2.x[1] + m2.x[3] >= 1)

    if with_degenerate_constraint:
        m2.con5 = pyo.Constraint(expr=m2.x[1] + m2.x[2] + m2.x[3] == 1)

    m2.obj = pyo.Objective(expr=sum(m2.x[i] for i in m2.I))

    return m2


def extract_constraint_names(cs):
    """Get constraint names from ComponentSet

    Arguments:
        cs: ComponentSet object

    Return:
        constraint_names: list of constraint names (strings)

    """

    constraint_names = []
    for i in cs:
        constraint_names.append(i.name)
    return constraint_names


# Problem 1
@pytest.mark.skipif(not pyo.SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem1():
    # Create test problem
    m = problem1()

    # Specify Ipopt as the solver
    opt = pyo.SolverFactory("ipopt")

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options["max_iter"] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m, tee=True)

    # Create Degeneracy Hunter object
    dh = DegeneracyHunter(m)

    # Find constraints with residuals > 0.1
    initial_point_constraints = dh.check_residuals(tol=0.1)

    # Check there are 2 constraints with large residuals
    assert len(initial_point_constraints) == 2

    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)

    # Check first constraint
    assert initial_point_constraint_names[0] == "con1"

    # Check second constraint
    assert initial_point_constraint_names[1] == "con3"

    opt.options["max_iter"] = 50

    # Solve
    opt.solve(m, tee=True)

    # Find constraints with residuals > 0.1
    solution_constraints = dh.check_residuals(tol=1e-6)

    # Check at the solution no constraints are violated
    assert len(solution_constraints) == 0

    # Check no constraints are near their bounds
    solution_bounds = dh.check_variable_bounds(tol=0.1)

    # Check at the solution no constraints are violated
    assert len(solution_bounds) == 0


# Problem 2 without degenerate constraint
@pytest.mark.skipif(not pyo.SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem2_without_degenerate_constraint():

    # Create test problem instance
    m2 = example2(with_degenerate_constraint=False)

    # Specify Ipopt as the solver
    opt = pyo.SolverFactory("ipopt")

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options["max_iter"] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m2, tee=True)

    # Create Degeneracy Hunter object
    dh2 = DegeneracyHunter(m2)

    # Check for violated constraints at the initial point
    initial_point_constraints = dh2.check_residuals(tol=0.1)

    # Check there are 1 constraints with large residuals
    assert len(initial_point_constraints) == 1

    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)

    # Check first constraint
    assert initial_point_constraint_names[0] == "con2"

    # Resolve
    opt.options["max_iter"] = 500
    opt.solve(m2, tee=True)

    # Check solution
    x_sln = []

    for i in m2.I:
        x_sln.append(m2.x[i]())

    assert pytest.approx(x_sln[0], abs=1e-6) == 1.0
    assert pytest.approx(x_sln[1], abs=1e-6) == 0.0
    assert pytest.approx(x_sln[2], abs=1e-6) == 0.0


# Problem 2 with degenerate constraint
@pytest.mark.skipif(not pyo.SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem2_with_degenerate_constraint():

    # Create test problem instance
    m2 = example2(with_degenerate_constraint=True)

    # Specify Ipopt as the solver
    opt = pyo.SolverFactory("ipopt")

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options["max_iter"] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m2, tee=True)

    # Create Degeneracy Hunter object
    dh2 = DegeneracyHunter(m2)

    # Check for violated constraints at the initial point
    initial_point_constraints = dh2.check_residuals(tol=0.1)

    # Check there are 2 constraints with large residuals
    assert len(initial_point_constraints) == 2

    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)

    # Check first constraint
    assert initial_point_constraint_names[0] == "con2"

    # Check first constraint
    assert initial_point_constraint_names[1] == "con5"

    # Resolve
    opt.options["max_iter"] = 500
    opt.solve(m2, tee=True)

    # Check solution
    x_sln = []

    for i in m2.I:
        x_sln.append(m2.x[i]())

    assert pytest.approx(x_sln[0], abs=1e-6) == 1.0
    assert pytest.approx(x_sln[1], abs=1e-6) == 0.0
    assert pytest.approx(x_sln[2], abs=1e-6) == 0.0

    # Check the rank
    n_rank_deficient = dh2.check_rank_equality_constraints()

    # Test DH with SVD

    assert n_rank_deficient == 1

    # TODO: Add MILP solver to idaes get-extensions and add more tests
