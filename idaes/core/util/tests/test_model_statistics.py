#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
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
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)

from idaes.core.util.model_statistics import *
from idaes.core.util.model_statistics import _iter_indexed_block_data_objects
from idaes.core.scaling.util import set_scaling_factor


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
        "2.2.0, will be removed in (or after) 2.11.0)"
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
        "2.2.0, will be removed in (or after) 2.11.0)"
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
def test_variables_near_bounds_set_scaled():
    m = ConcreteModel()
    m.v = Var(initialize=0.5, bounds=(0, 1))

    # Small value, both bounds
    # Away from bounds
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Small value, large scaling factor, both bounds
    set_scaling_factor(m.v, 1e6)
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Small value, tiny scaling factor, both bounds
    set_scaling_factor(m.v, 1e-6, overwrite=True)
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=0)
    assert len(vset) == 0

    # Near lower bound, relative
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.set_value(1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near lower bound, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=0)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near lower bound, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=0)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near lower bound, absolute
    set_scaling_factor(m.v, 1, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Near upper bound, relative
    m.v.set_value(1 - 1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near upper bound, relative, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=0)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near upper bound, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=0)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near upper bound, absolute
    set_scaling_factor(m.v, 1, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Small value, lower bound
    # Away from bounds
    m.v.set_value(0.5)
    m.v.setub(None)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Small value, lower bound
    # Away from bounds
    # Enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    m.v.setub(None)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 0

    # Small value, lower bound
    # Away from bounds
    # Tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    m.v.setub(None)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 1

    # Near lower bound, relative
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.set_value(1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    # Lower bound of 0 means relative tolerance is 0
    assert len(vset) == 0

    # Near lower bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Near lower bound, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-4)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=0)
    assert len(vset) == 1

    # Near lower bound, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=0)
    assert len(vset) == 0

    # Small value, upper bound
    # Away from bounds
    m.v.set_value(0.5)
    m.v.setub(1)
    m.v.setlb(None)
    set_scaling_factor(m.v, 1, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Small value, upper bound
    # Away from bounds
    # Tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 1

    # Small value, upper bound
    # Away from bounds
    # Enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 0

    # Near upper bound, relative
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.set_value(1 - 1e-6)
    vset = variables_near_bounds_set(m, abs_tol=1e-8, rel_tol=1e-4)
    assert len(vset) == 1

    # Near upper bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1

    # Near upper bound, absolute, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=0)
    assert len(vset) == 1

    # Near upper bound, absolute, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=1e-8)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1e-4, rel_tol=0)
    assert len(vset) == 0

    # Large value, both bounds
    # Relative tolerance based on magnitude of 100
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.setlb(450)
    m.v.setub(550)
    m.v.set_value(500)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Large value, both bounds
    # Tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 1

    # Large value, both bounds
    # Enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 0

    # Near lower bound, relative
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.set_value(451)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near lower bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1

    # Near lower bound, absolute, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1, rel_tol=0)
    assert len(vset) == 1

    # Near lower bound, absolute, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1, rel_tol=0)
    assert len(vset) == 0

    # Near upper bound, relative
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.set_value(549)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near upper bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1

    # Near upper bound, absolute, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-6)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-6, rel_tol=0)
    assert len(vset) == 1

    # Near upper bound, absolute, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1, rel_tol=0)
    assert len(vset) == 0

    # Large value, lower bound
    # Relative tolerance based on magnitude of 450
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.setlb(450)
    m.v.setub(None)
    m.v.set_value(500)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Large value, lower bound, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=0)
    assert len(vset) == 0

    # Large value, lower bound, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, rel_tol=0)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=0)
    assert len(vset) == 0

    # Near lower bound, relative
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.set_value(451)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near lower bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1

    # Near lower bound, absolute, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-6)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-6, rel_tol=0)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=0)
    assert len(vset) == 0

    # Near lower bound, absolute, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-6)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1e-6, rel_tol=0)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=0, rel_tol=1e-2)
    assert len(vset) == 1

    # Large value, upper bound
    # Relative tolerance based on magnitude of 550
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.setlb(None)
    m.v.setub(550)
    m.v.set_value(500)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0

    # Large value, upper bound
    # Tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, rel_tol=0, abs_tol=1e-6)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=0)
    assert len(vset) == 0

    # Large value, upper bound
    # Enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, rel_tol=0, abs_tol=1e-6)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=0)
    assert len(vset) == 0

    # Near upper bound, relative
    set_scaling_factor(m.v, 1, overwrite=True)
    m.v.set_value(549)
    vset = variables_near_bounds_set(m, rel_tol=1e-2)
    assert len(vset) == 1

    # Near upper bound, absolute
    vset = variables_near_bounds_set(m, abs_tol=1)
    assert len(vset) == 1

    # Near upper bound, absolute, tiny scaling factor
    set_scaling_factor(m.v, 1e-8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-6)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-6, rel_tol=0)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=0)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=0, rel_tol=1e-2)
    assert len(vset) == 1

    # Near upper bound, absolute, enormous scaling factor
    set_scaling_factor(m.v, 1e8, overwrite=True)
    vset = variables_near_bounds_set(m, abs_tol=1e-6)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=1e-6, rel_tol=1e-2)
    assert len(vset) == 1
    vset = variables_near_bounds_set(m, abs_tol=1e-6, rel_tol=0)
    assert len(vset) == 0
    vset = variables_near_bounds_set(m, abs_tol=0)
    assert len(vset) == 0


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


@pytest.mark.unit
def test_uninitialized_variables_in_activated_constraints():
    m = ConcreteModel()
    m.u = Var()
    m.w = Var(initialize=1)
    m.x = Var(initialize=1)
    m.y = Var(range(3), initialize=[1, None, 2])
    m.z = Var()

    m.con_w = Constraint(expr=m.w == 0)
    m.con_x = Constraint(expr=m.x == 1)

    def rule_con_y(b, i):
        return b.y[i] == 3

    m.con_y = Constraint(range(3), rule=rule_con_y)
    m.con_z = Constraint(expr=m.x + m.z == -1)

    m.block1 = Block()
    m.block1.a = Var()
    m.block1.b = Var(initialize=7)
    m.block1.c = Var(initialize=-5)
    m.block1.d = Var()

    m.block1.con_a = Constraint(expr=m.block1.a**2 + m.block1.b == 1)
    m.block1.con_b = Constraint(expr=m.block1.b + m.block1.c == 3)
    m.block1.con_c = Constraint(expr=m.block1.c + m.block1.a == -4)

    active_uninit_set = ComponentSet([m.y[1], m.z, m.block1.a])

    assert variables_with_none_value_in_activated_equalities_set(m) == active_uninit_set
    assert number_variables_with_none_value_in_activated_equalities(m) == len(
        active_uninit_set
    )

    m.block1.deactivate()

    active_uninit_set = ComponentSet([m.y[1], m.z])

    assert variables_with_none_value_in_activated_equalities_set(m) == active_uninit_set
    assert number_variables_with_none_value_in_activated_equalities(m) == len(
        active_uninit_set
    )

    m.block1.activate()
    m.con_z.deactivate()

    active_uninit_set = ComponentSet([m.y[1], m.block1.a])

    assert variables_with_none_value_in_activated_equalities_set(m) == active_uninit_set
    assert number_variables_with_none_value_in_activated_equalities(m) == len(
        active_uninit_set
    )


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
def test_degrees_of_freedom_with_graybox():
    """non functional graybox model added to m fixture, to test DOFs

    GreyBoxModel has 3 inputs and 2 outputs calculated an unknown function,
    input a1 and a2 are bound by equality constraint through internal graybox model"""

    class BasicGrayBox(ExternalGreyBoxModel):
        def input_names(self):
            return ["a1", "a2", "a3"]

        def output_names(self):
            return ["o1", "o2"]

        def equality_constraint_names(self):
            return ["a_sum"]

        def evaluate_equality_constraints(self):
            a1 = self._input_values[0]
            a2 = self._input_values[1]
            return [a1 * 0.5 + a2]

    m = ConcreteModel()

    m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
    m.gb_inactive = ExternalGreyBoxBlock(external_model=BasicGrayBox())
    m.gb_inactive.deactivate()
    # test counting functions
    assert number_greybox_blocks(m) == 2
    assert number_deactivated_greybox_block(m) == 1
    assert number_activated_greybox_blocks(m) == 1
    assert number_of_greybox_variables(m) == 5
    assert number_of_unfixed_greybox_variables(m) == 5
    assert number_activated_greybox_equalities(m) == 3
    assert number_variables_in_activated_constraints(m) == 5
    # verify DOFS works on stand alone greybox
    assert degrees_of_freedom(m) == 2
    m.gb.inputs.fix()
    m.gb.inputs["a1"].unfix()
    assert number_of_unfixed_greybox_variables(m) == 3
    assert degrees_of_freedom(m) == 0
    m.gb.outputs.fix()
    assert degrees_of_freedom(m) == -2
    m.gb.outputs.unfix()

    # verify DOFs works on greybox connected to other vars on a model via constraints
    m.a1 = Var(initialize=1)
    m.a1.fix()
    m.gb.inputs["a2"].unfix()
    m.a1_eq = Constraint(expr=m.a1 == m.gb.inputs["a1"])
    assert degrees_of_freedom(m) == 0
    m.o1 = Var(initialize=1)
    m.o1_eq = Constraint(expr=m.o1 == m.gb.outputs["o1"])
    m.o1.fix()
    assert degrees_of_freedom(m) == -1
    assert number_variables_in_activated_constraints(m) == 7
    assert number_total_constraints(m) == 5
    assert number_total_equalities(m) == 5
    assert number_deactivated_equalities(m) == 3
    assert number_deactivated_constraints(m) == 3


@pytest.mark.unit
def test_large_residuals_set(m):
    # Initialize derivative var values so no errors occur
    for v in m.dv.keys():
        m.dv[v] = 0

    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b2["b"].c1 in lrs
    assert m.b2["b"].c2 in lrs

    m.b2["b"].c1.deactivate()
    lrs = large_residuals_set(m)
    assert len(lrs) == 1
    assert m.b2["b"].c2 in lrs

    m.b1.sb.v1.fix(0)
    # m.b1.sb.c1 now is not satisfied, but because
    # b1 is deactivated it should not be counted
    lrs = large_residuals_set(m)
    assert len(lrs) == 1
    assert m.b2["b"].c2 in lrs

    m.b1.activate()
    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b1.sb.c1 in lrs
    assert m.b2["b"].c2 in lrs

    m.b1.sb.v1.unfix()
    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b1.sb.c1 in lrs
    assert m.b2["b"].c2 in lrs

    # m.b1.sb.c2 now is not satisfied
    m.b1.sb.v1.set_value(2)
    lrs = large_residuals_set(m)
    assert len(lrs) == 3
    assert m.b1.sb.c1 in lrs
    assert m.b1.sb.c2 in lrs
    assert m.b2["b"].c2 in lrs

    # The tiny scaling factor on c1 means we can
    # get away with extremely large errors
    set_scaling_factor(m.b1.sb.c1, 1e-10)
    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b1.sb.c2 in lrs
    assert m.b2["b"].c2 in lrs

    # m.b.sb.c2 is not satisfied, but the residual
    # is less than the tolerance
    # Testing above upper bound
    m.b1.sb.v1.set_value(1 + 1e-6)
    lrs = large_residuals_set(m)
    assert len(lrs) == 1
    assert m.b2["b"].c2 in lrs

    # By setting a huge scaling factor on c2,
    # the constraint is no longer satisfied
    set_scaling_factor(m.b1.sb.c2, 1e10)
    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b1.sb.c2 in lrs
    assert m.b2["b"].c2 in lrs

    # Testing below upper bound
    m.b1.sb.v1.set_value(1 - 1e-6)
    lrs = large_residuals_set(m)
    assert len(lrs) == 1
    assert m.b2["b"].c2 in lrs

    # Upset this variable to make m.b1.c1 and m.b1.c2
    # not satisfied.
    m.b1.v1.fix(1 - 0.1)
    lrs = large_residuals_set(m)
    assert len(lrs) == 3
    assert m.b1.c1 in lrs
    assert m.b1.c2 in lrs
    assert m.b2["b"].c2 in lrs

    # Variable scaling shouldn't affect anything
    set_scaling_factor(m.b1.v1, 1e-12, overwrite=True)
    lrs = large_residuals_set(m)
    assert len(lrs) == 3
    assert m.b1.c1 in lrs
    assert m.b1.c2 in lrs
    assert m.b2["b"].c2 in lrs

    set_scaling_factor(m.b1.v1, 1e12, overwrite=True)
    lrs = large_residuals_set(m)
    assert len(lrs) == 3
    assert m.b1.c1 in lrs
    assert m.b1.c2 in lrs
    assert m.b2["b"].c2 in lrs

    # Constraint scaling lets us ignore the residual
    set_scaling_factor(m.b1.c1, 1e-6, overwrite=True)
    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b1.c2 in lrs
    assert m.b2["b"].c2 in lrs

    assert len(large_residuals_set(m, tol=1)) == 0
    set_scaling_factor(m.b1.c1, 1e5, overwrite=True)
    lrs = large_residuals_set(m, tol=1)
    assert len(lrs) == 1
    assert m.b1.c1 in lrs

    set_scaling_factor(m.b1.c1, 1, overwrite=True)

    # Now try the inequality constraint
    set_scaling_factor(m.b1.c2, 1e-5, overwrite=True)
    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b1.c1 in lrs
    assert m.b2["b"].c2 in lrs

    assert len(large_residuals_set(m, tol=1)) == 0
    set_scaling_factor(m.b1.c2, 1e5, overwrite=True)
    lrs = large_residuals_set(m, tol=1)
    assert len(lrs) == 1
    assert m.b1.c2 in lrs

    # Put us on the other side of the inequality constraint
    m.b1.v1.fix(1 + 0.1)
    lrs = large_residuals_set(m)
    assert len(lrs) == 2
    assert m.b1.c1 in lrs
    assert m.b2["b"].c2 in lrs


@pytest.mark.unit
def test_large_residuals_set_none(m):
    # Without initializing the derivative variables,
    # all the discretization equations get added
    # to the large residuals set.
    # Contains b2[b].c1, b2[b].c2, and dv_disc_eq[j]
    # from j=0.1 to j=1.
    lrs = large_residuals_set(m)
    assert len(lrs) == 12
    assert m.b2["b"].c1 in lrs
    assert m.b2["b"].c2 in lrs
    assert len(m.dv_disc_eq) == 10
    for condata in m.dv_disc_eq.values():
        assert condata in lrs


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
