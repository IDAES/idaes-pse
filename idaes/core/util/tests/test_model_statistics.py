#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains miscellaneous utility functions for use in IDAES models.
"""

import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    Set,
    Var,
    TransformationFactory,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.collections import ComponentSet

from idaes.core.util.model_statistics import *
from idaes.core.util.model_statistics import _iter_indexed_block_data_objects


@pytest.mark.unit
def test_iter_indexed_block_data_objects_scalar_block():
    m = ConcreteModel()

    m.s = Set(initialize=["a", "b"])
    m.b = Block()

    m.b.v = Var(m.s, initialize=1, bounds=(0, 10))

    chk_lst = list(_iter_indexed_block_data_objects(m.b, Var, True, True))
    assert len(chk_lst) == 2
    for v in chk_lst:
        assert v.name in ["b.v[a]", "b.v[b]"]


@pytest.mark.unit
def test_iter_indexed_block_data_objects_indexed_block():
    m = ConcreteModel()

    m.s = Set(initialize=["a", "b"])
    m.b = Block(m.s)

    m.b["a"].v = Var(m.s, initialize=1, bounds=(0, 10))
    m.b["b"].v = Var(m.s, initialize=1, bounds=(0, 10))

    chk_lst = list(_iter_indexed_block_data_objects(m.b, Var, True, True))
    assert len(chk_lst) == 4
    for v in chk_lst:
        assert v.name in ["b[a].v[a]", "b[a].v[b]", "b[b].v[a]", "b[b].v[b]"]


# Author: Andrew Lee
@pytest.fixture()
def m():
    m = ConcreteModel()

    m.s = Set(initialize=["a", "b"])
    m.cs = ContinuousSet(bounds=(0, 1))

    m.v = Var(m.cs, initialize=1, bounds=(0, 10))
    m.dv = DerivativeVar(m.v)

    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=10, wrt=m.cs, scheme="BACKWARD")

    m.e = Expression(expr=m.v[1])

    m.b1 = Block()
    m.b1.v1 = Var(initialize=1)
    m.b1.v2 = Var(m.s, initialize=1)
    m.b1.v1.fix(1)
    m.b1.c1 = Constraint(expr=1 == m.b1.v1)
    m.b1.c2 = Constraint(expr=1 <= m.b1.v1)

    m.b1.sb = Block()
    m.b1.sb.v1 = Var(initialize=1)
    m.b1.sb.v2 = Var(m.s, initialize=1)
    m.b1.sb.v1.fix(1)
    m.b1.sb.e1 = Expression(expr=m.b1.sb.v1)
    m.b1.sb.o1 = Objective(expr=m.b1.sb.v1)
    m.b1.sb.o2 = Objective(expr=m.b1.sb.v1)
    m.b1.sb.o2.deactivate()
    m.b1.sb.c1 = Constraint(expr=1 == m.b1.sb.v1)
    m.b1.sb.c2 = Constraint(expr=1 >= m.b1.sb.v1)

    m.b1.deactivate()

    m.b2 = Block(m.s)
    for i in m.s:
        m.b2[i].v1 = Var(initialize=1, bounds=(1, 10))
        m.b2[i].v2 = Var(m.s, initialize=1, bounds=(0, 1))
        m.b2[i].v1.fix(1)
        m.b2[i].e1 = Expression(expr=m.b2[i].v1)
        m.b2[i].c1 = Constraint(expr=2 == m.b2[i].v1)
        m.b2[i].c2 = Constraint(expr=2 <= m.b2[i].v1)

        if i == "a":
            m.b2[i].o1 = Objective(expr=m.b2[i].v1)
            m.b2[i].o2 = Objective(expr=m.b2[i].v1)
            m.b2[i].o2.deactivate()
            m.b2[i].c1.deactivate()
            m.b2[i].c2.deactivate()

    return m


# -------------------------------------------------------------------------
# Block methods
@pytest.mark.unit
def test_total_blocks_set(m):
    assert len(total_blocks_set(m)) == 5
    assert len(total_blocks_set(m.b2)) == 1


@pytest.mark.unit
def test_number_total_blocks(m):
    assert number_total_blocks(m) == 5
    assert number_total_blocks(m.b2) == 1


@pytest.mark.unit
def test_activated_blocks_set(m):
    assert len(activated_blocks_set(m)) == 3
    assert len(activated_blocks_set(m.b2)) == 1


@pytest.mark.unit
def test_number_activated_blocks(m):
    assert number_activated_blocks(m) == 3
    assert number_activated_blocks(m.b2) == 1


@pytest.mark.unit
def test_deactivated_blocks_set(m):
    assert len(deactivated_blocks_set(m)) == 2
    assert len(deactivated_blocks_set(m.b2)) == 0


@pytest.mark.unit
def test_number_deactivated_blocks(m):
    assert number_deactivated_blocks(m) == 2
    assert number_deactivated_blocks(m.b2) == 0


# -------------------------------------------------------------------------
# Basic Constraint methods
@pytest.mark.unit
def test_total_constraints_set(m):
    assert len(total_constraints_set(m)) == 14
    assert len(total_constraints_set(m.b2)) == 4


@pytest.mark.unit
def test_number_total_constraints(m):
    assert number_total_constraints(m) == 14
    assert number_total_constraints(m.b2) == 4


@pytest.mark.unit
def test_activated_constraints_set(m):
    assert len(activated_constraints_set(m)) == 12
    assert len(activated_constraints_set(m.b2)) == 2


@pytest.mark.unit
def test_number_activated_constraints(m):
    assert number_activated_constraints(m) == 12
    assert number_activated_constraints(m.b2) == 2


@pytest.mark.unit
def test_deactivated_constraints_set(m):
    assert len(deactivated_constraints_set(m)) == 2
    assert len(deactivated_constraints_set(m.b2)) == 2


@pytest.mark.unit
def test_number_deactivated_constraints(m):
    assert number_deactivated_constraints(m) == 2
    assert number_deactivated_constraints(m.b2) == 2


# -------------------------------------------------------------------------
# Equality Constraints
@pytest.mark.unit
def test_total_equalities_set(m):
    assert len(total_equalities_set(m)) == 12
    assert len(total_equalities_set(m.b2)) == 2


@pytest.mark.unit
def test_number_total_equalities(m):
    assert number_total_equalities(m) == 12
    assert number_total_equalities(m.b2) == 2


@pytest.mark.unit
def test_activated_equalities_set(m):
    assert len(activated_equalities_set(m)) == 11
    assert len(activated_equalities_set(m.b2)) == 1


@pytest.mark.unit
def test_number_activated_equalities(m):
    assert number_activated_equalities(m) == 11
    assert number_activated_equalities(m.b2) == 1


@pytest.mark.unit
def test_deactivated_equalities_set(m):
    assert len(deactivated_equalities_set(m)) == 1
    assert len(deactivated_equalities_set(m.b2)) == 1


@pytest.mark.unit
def test_number_deactivated_equalities(m):
    assert number_deactivated_equalities(m) == 1
    assert number_deactivated_equalities(m.b2) == 1


# -------------------------------------------------------------------------
# Inequality Constraints
@pytest.mark.unit
def test_total_inequalities_set(m):
    assert len(total_inequalities_set(m)) == 2
    assert len(total_inequalities_set(m.b2)) == 2


@pytest.mark.unit
def test_number_total_inequalities(m):
    assert number_total_inequalities(m) == 2
    assert number_total_inequalities(m.b2) == 2


@pytest.mark.unit
def test_activated_inequalities_set(m):
    assert len(activated_inequalities_set(m)) == 1
    assert len(activated_inequalities_set(m.b2)) == 1


@pytest.mark.unit
def test_number_activated_inequalities(m):
    assert number_activated_inequalities(m) == 1
    assert number_activated_inequalities(m.b2) == 1


@pytest.mark.unit
def test_deactivated_inequalities_set(m):
    assert len(deactivated_inequalities_set(m)) == 1
    assert len(deactivated_inequalities_set(m.b2)) == 1


@pytest.mark.unit
def test_number_deactivated_inequalities(m):
    assert number_deactivated_inequalities(m) == 1
    assert number_deactivated_inequalities(m.b2) == 1


# -------------------------------------------------------------------------
# Basic Variable Methods
# Always use ComponentSets for Vars to avoid duplication of References
# i.e. number methods should always use the ComponentSet, not a generator
@pytest.mark.unit
def test_variables_set(m):
    assert len(variables_set(m)) == 28
    assert len(variables_set(m.b2)) == 6


@pytest.mark.unit
def test_number_variables(m):
    assert number_variables(m) == 28
    assert number_variables(m.b2) == 6


@pytest.mark.unit
def test_fixed_variables_set(m):
    assert len(fixed_variables_set(m)) == 2
    assert len(fixed_variables_set(m.b2)) == 2


@pytest.mark.unit
def test_number_fixed_variables(m):
    assert number_fixed_variables(m) == 2
    assert number_fixed_variables(m.b2) == 2


@pytest.mark.unit
def test_unfixed_variables_set(m):
    assert len(unfixed_variables_set(m)) == 26
    assert len(unfixed_variables_set(m.b2)) == 4


@pytest.mark.unit
def test_number_unfixed_variables(m):
    assert number_unfixed_variables(m) == 26
    assert number_unfixed_variables(m.b2) == 4


@pytest.mark.unit
def test_variables_near_bounds_tol_deprecation(m, caplog):
    variables_near_bounds_set(m, tol=1e-3)

    msg = (
        "DEPRECATED: variables_near_bounds_generator has deprecated the tol argument. "
        "Please set abs_tol and rel_tol arguments instead.  (deprecated in "
        "2.2.0, will be removed in (or after) 3.0.0)"
    )
    assert msg.replace(" ", "") in caplog.records[0].message.replace("\n", "").replace(
        " ", ""
    )


@pytest.mark.unit
def test_variables_near_bounds_relative_deprecation(m, caplog):
    variables_near_bounds_set(m, relative=False)

    msg = (
        "DEPRECATED: variables_near_bounds_generator has deprecated the relative argument. "
        "Please set abs_tol and rel_tol arguments instead.  (deprecated in "
        "2.2.0, will be removed in (or after) 3.0.0)"
    )
    assert msg.replace(" ", "") in caplog.records[0].message.replace("\n", "").replace(
        " ", ""
    )


@pytest.mark.unit
def test_variables_near_bounds_set():
    m = ConcreteModel()
    m.v = Var(initialize=0.5, bounds=(0, 1))

    # Small value, both bounds
    # Away from bounds
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Near lower bound, relative
    m.v.set_value(1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near lower bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Near upper bound, relative
    m.v.set_value(1 - 1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near upper bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Small value, lower bound
    # Away from bounds
    m.v.set_value(0.5)
    m.v.setub(None)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Near lower bound, relative
    m.v.set_value(1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    # Lower bound of 0 means relative tolerance is 0
    assert len(vset) == 0

    # Near lower bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Small value, upper bound
    # Away from bounds
    m.v.set_value(0.5)
    m.v.setub(1)
    m.v.setlb(None)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Near upper bound, relative
    m.v.set_value(1 - 1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near upper bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Large value, both bounds
    # Relative tolerance based on magnitude of 100
    m.v.setlb(450)
    m.v.setub(550)
    m.v.set_value(500)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Near lower bound, relative
    m.v.set_value(451)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near lower bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1

    # Near upper bound, relative
    m.v.set_value(549)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near upper bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1

    # Large value, lower bound
    # Relative tolerance based on magnitude of 450
    m.v.setlb(450)
    m.v.setub(None)
    m.v.set_value(500)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Near lower bound, relative
    m.v.set_value(451)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near lower bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1

    # Large value, upper bound
    # Relative tolerance based on magnitude of 550
    m.v.setlb(None)
    m.v.setub(550)
    m.v.set_value(500)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Near upper bound, relative
    m.v.set_value(549)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near upper bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1


@pytest.mark.unit
def test_number_variables_near_bounds(m):
    assert number_variables_near_bounds(m) == 6
    assert number_variables_near_bounds(m.b2) == 6


# -------------------------------------------------------------------------
# Variables in Constraints
@pytest.mark.unit
def test_variables_in_activated_constraints_set(m):
    assert len(variables_in_activated_constraints_set(m)) == 22
    assert len(variables_in_activated_constraints_set(m.b2)) == 1


@pytest.mark.unit
def test_number_variables_in_activated_constraints(m):
    assert number_variables_in_activated_constraints(m) == 22
    assert number_variables_in_activated_constraints(m.b2) == 1


@pytest.mark.unit
def test_variables_not_in_activated_constraints_set(m):
    assert len(variables_not_in_activated_constraints_set(m)) == 6
    assert len(variables_not_in_activated_constraints_set(m.b2)) == 5


@pytest.mark.unit
def test_number_variables_not_in_activated_constraints(m):
    assert number_variables_not_in_activated_constraints(m) == 6
    assert number_variables_not_in_activated_constraints(m.b2) == 5


@pytest.mark.unit
def test_variables_in_activated_equalities_set(m):
    assert len(variables_in_activated_equalities_set(m)) == 22
    assert len(variables_in_activated_equalities_set(m.b2)) == 1


@pytest.mark.unit
def test_number_variables_in_activated_equalities(m):
    assert number_variables_in_activated_equalities(m) == 22
    assert number_variables_in_activated_equalities(m.b2) == 1


@pytest.mark.unit
def test_variables_in_activated_inequalities_set(m):
    assert len(variables_in_activated_inequalities_set(m)) == 1
    assert len(variables_in_activated_inequalities_set(m.b2)) == 1


@pytest.mark.unit
def test_number_variables_in_activated_inequalities(m):
    assert number_variables_in_activated_inequalities(m) == 1
    assert number_variables_in_activated_inequalities(m.b2) == 1


@pytest.mark.unit
def test_variables_only_in_inequalities(m):
    m.v3 = Var(initialize=1)
    m.c3 = Constraint(expr=m.v3 >= 1)
    assert len(variables_only_in_inequalities(m)) == 1
    assert len(variables_only_in_inequalities(m.b2)) == 0


@pytest.mark.unit
def test_number_variables_only_in_inequalities(m):
    m.v3 = Var(initialize=1)
    m.c3 = Constraint(expr=m.v3 >= 1)
    assert number_variables_only_in_inequalities(m) == 1
    assert number_variables_only_in_inequalities(m.b2) == 0


# -------------------------------------------------------------------------
# Fixed Variables in Constraints
@pytest.mark.unit
def test_fixed_variables_in_activated_equalities_set(m):
    assert len(fixed_variables_in_activated_equalities_set(m)) == 1
    assert len(fixed_variables_in_activated_equalities_set(m.b2)) == 1


@pytest.mark.unit
def test_number_fixed_variables_in_activated_equalities(m):
    assert number_fixed_variables_in_activated_equalities(m) == 1
    assert number_fixed_variables_in_activated_equalities(m.b2) == 1


@pytest.mark.unit
def test_fixed_variables_only_in_inequalities(m):
    m.v3 = Var(initialize=1)
    m.v3.fix(1)
    m.c3 = Constraint(expr=m.v3 >= 1)
    assert len(fixed_variables_only_in_inequalities(m)) == 1
    assert len(fixed_variables_only_in_inequalities(m.b2)) == 0


@pytest.mark.unit
def test_number_fixed_variables_only_in_inequalities(m):
    m.v3 = Var(initialize=1)
    m.v3.fix(1)
    m.c3 = Constraint(expr=m.v3 >= 1)
    assert number_fixed_variables_only_in_inequalities(m) == 1
    assert number_fixed_variables_only_in_inequalities(m.b2) == 0


# -------------------------------------------------------------------------
# Unused and un-Transformed Variables
@pytest.mark.unit
def test_unused_variables_set(m):
    assert len(unused_variables_set(m)) == 6
    assert len(unused_variables_set(m.b2)) == 5


@pytest.mark.unit
def test_number_unused_variables(m):
    assert number_unused_variables(m) == 6
    assert number_unused_variables(m.b2) == 5


@pytest.mark.unit
def test_fixed_unused_variables_set(m):
    assert len(fixed_unused_variables_set(m)) == 1
    assert len(fixed_unused_variables_set(m.b2)) == 1


@pytest.mark.unit
def test_number_fixed_unused_variables(m):
    assert number_fixed_unused_variables(m) == 1
    assert number_fixed_unused_variables(m.b2) == 1


@pytest.mark.unit
def test_derivative_variables_set():
    m = ConcreteModel()

    m.cs = ContinuousSet(bounds=(0, 1))

    m.v = Var(m.cs, initialize=1)
    m.dv = DerivativeVar(m.v)

    assert len(derivative_variables_set(m)) == 2

    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=10, wrt=m.cs, scheme="BACKWARD")

    assert len(derivative_variables_set(m)) == 0


@pytest.mark.unit
def test_number_derivative_variables():
    m = ConcreteModel()

    m.cs = ContinuousSet(bounds=(0, 1))

    m.v = Var(m.cs, initialize=1)
    m.dv = DerivativeVar(m.v)

    assert number_derivative_variables(m) == 2

    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=10, wrt=m.cs, scheme="BACKWARD")

    assert number_derivative_variables(m) == 0


# -------------------------------------------------------------------------
# Objective methods
@pytest.mark.unit
def test_total_objectives_set(m):
    assert len(total_objectives_set(m)) == 2
    assert len(total_objectives_set(m.b2)) == 2


@pytest.mark.unit
def test_number_total_objectives(m):
    assert number_total_objectives(m) == 2
    assert number_total_objectives(m.b2) == 2


@pytest.mark.unit
def test_activated_objectives_set(m):
    assert len(activated_objectives_set(m)) == 1


@pytest.mark.unit
def test_number_activated_objectives(m):
    assert number_activated_objectives(m) == 1
    assert number_activated_objectives(m.b2) == 1


@pytest.mark.unit
def test_deactivated_objectives_set(m):
    assert len(deactivated_objectives_set(m)) == 1
    assert len(deactivated_objectives_set(m.b2)) == 1


@pytest.mark.unit
def test_number_deactivated_objectives(m):
    assert number_deactivated_objectives(m) == 1
    assert number_deactivated_objectives(m.b2) == 1


# -------------------------------------------------------------------------
# Expression methods
@pytest.mark.unit
def test_expressions_set(m):
    assert len(expressions_set(m)) == 3
    assert len(expressions_set(m.b2)) == 2


@pytest.mark.unit
def test_number_expressions(m):
    assert number_expressions(m) == 3
    assert number_expressions(m.b2) == 2


# -------------------------------------------------------------------------
# Other model statistics
@pytest.mark.unit
def test_degrees_of_freedom(m):
    assert degrees_of_freedom(m) == 10
    assert degrees_of_freedom(m.b2) == -1


@pytest.mark.unit
def test_large_residuals_set(m):
    # Initialize derivative var values so no errors occur
    for v in m.dv.keys():
        m.dv[v] = 0
    assert len(large_residuals_set(m)) == 2


@pytest.mark.unit
def test_number_large_residuals(m):
    # Initialize derivative var values so no errors occur
    for v in m.dv.keys():
        m.dv[v] = 0
    assert number_large_residuals(m) == 2


@pytest.mark.unit
def test_active_variables_in_deactivated_blocks_set(m):
    assert len(active_variables_in_deactivated_blocks_set(m)) == 0
    assert (
        len(active_variables_in_deactivated_blocks_set(m.b2)) == 1
    )  # TODO: Not sure why?

    m.c = Constraint(expr=m.b1.v1 >= 2)

    assert len(active_variables_in_deactivated_blocks_set(m)) == 1
    assert (
        len(active_variables_in_deactivated_blocks_set(m.b2)) == 1
    )  # TODO: Not sure why?


@pytest.mark.unit
def test_number_active_variables_in_deactivated_blocks(m):
    assert number_active_variables_in_deactivated_blocks(m) == 0
    assert (
        number_active_variables_in_deactivated_blocks(m.b2) == 1
    )  # TODO: Not sure why?

    m.c = Constraint(expr=m.b1.v1 >= 2)

    assert number_active_variables_in_deactivated_blocks(m) == 1
    assert (
        number_active_variables_in_deactivated_blocks(m.b2) == 1
    )  # TODO: Not sure why?


# -------------------------------------------------------------------------
# Reporting methods
@pytest.mark.unit
def test_report_statistics(m):
    report_statistics(m)
