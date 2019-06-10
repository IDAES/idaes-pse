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
This module contains miscalaneous utility functions for use in IDAES models.
"""

import pytest

from pyomo.environ import (Block,
                           ConcreteModel,
                           Constraint,
                           Expression,
                           Objective,
                           Set,
                           Var,
                           TransformationFactory)
from pyomo.dae import ContinuousSet, DerivativeVar

from idaes.core.util.model_statistics import *


# Author: Andrew Lee
@pytest.fixture()
def m():
    m = ConcreteModel()

    m.s = Set(initialize=["a", "b"])
    m.cs = ContinuousSet(bounds=(0, 1))

    m.v = Var(m.cs, initialize=1)
    m.dv = DerivativeVar(m.v)

    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m,
                           nfe=10,
                           wrt=m.cs,
                           scheme="BACKWARD")

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
    m.b1.sb.o1 = Objective(expr=1 <= m.b1.sb.v1)
    m.b1.sb.o2 = Objective(expr=1 <= m.b1.sb.v1)
    m.b1.sb.o2.deactivate()
    m.b1.sb.c1 = Constraint(expr=1 == m.b1.sb.v1)
    m.b1.sb.c2 = Constraint(expr=1 >= m.b1.sb.v1)

    m.b1.deactivate()

    m.b2 = Block(m.s)
    for i in m.s:
        m.b2[i].v1 = Var(initialize=1)
        m.b2[i].v2 = Var(m.s, initialize=1)
        m.b2[i].v1.fix(1)
        m.b2[i].e1 = Expression(expr=m.b2[i].v1)
        m.b2[i].c1 = Constraint(expr=2 == m.b2[i].v1)
        m.b2[i].c2 = Constraint(expr=2 <= m.b2[i].v1)

        if i == "a":
            m.b2[i].o1 = Objective(expr=1 <= m.b2[i].v1)
            m.b2[i].o2 = Objective(expr=1 <= m.b2[i].v1)
            m.b2[i].o2.deactivate()
            m.b2[i].c1.deactivate()
            m.b2[i].c2.deactivate()

    return m


# -------------------------------------------------------------------------
# Block methods
def test_total_blocks_set(m):
    assert len(total_blocks_set(m)) == 5


def test_number_total_blocks(m):
    assert number_total_blocks(m) == 5


def test_activated_blocks_set(m):
    assert len(activated_blocks_set(m)) == 3


def test_number_activated_blocks(m):
    assert number_activated_blocks(m) == 3


def test_deactivated_blocks_set(m):
    assert len(deactivated_blocks_set(m)) == 2


def test_number_deactivated_blocks(m):
    assert number_deactivated_blocks(m) == 2


# -------------------------------------------------------------------------
# Basic Constraint methods
def test_total_constraints_set(m):
    assert len(total_constraints_set(m)) == 14


def test_number_total_constraints(m):
    assert number_total_constraints(m) == 14


def test_activated_constraints_set(m):
    assert len(activated_constraints_set(m)) == 12


def test_number_activated_constraints(m):
    assert number_activated_constraints(m) == 12


def test_deactivated_constraints_set(m):
    assert len(deactivated_constraints_set(m)) == 2


def test_number_deactivated_constraints(m):
    assert number_deactivated_constraints(m) == 2


# -------------------------------------------------------------------------
# Equality Constraints
def test_total_equalities_set(m):
    assert len(total_equalities_set(m)) == 12


def test_number_total_equalities(m):
    assert number_total_equalities(m) == 12


def test_activated_equalities_set(m):
    assert len(activated_equalities_set(m)) == 11


def test_number_activated_equalities(m):
    assert number_activated_equalities(m) == 11


def test_deactivated_equalities_set(m):
    assert len(deactivated_equalities_set(m)) == 1


def test_number_deactivated_equalities(m):
    assert number_deactivated_equalities(m) == 1


# -------------------------------------------------------------------------
# Inequality Constraints
def test_total_inequalities_set(m):
    assert len(total_inequalities_set(m)) == 2


def test_number_total_inequalities(m):
    assert number_total_inequalities(m) == 2


def test_activated_inequalities_set(m):
    assert len(activated_inequalities_set(m)) == 1


def test_number_activated_inequalities(m):
    assert number_activated_inequalities(m) == 1


def test_deactivated_inequalities_set(m):
    assert len(deactivated_inequalities_set(m)) == 1


def test_number_deactivated_inequalities(m):
    assert number_deactivated_inequalities(m) == 1


# -------------------------------------------------------------------------
# Basic Variable Methods
# Always use ComponentSets for Vars to avoid duplication of References
# i.e. number methods should alwys use the ComponentSet, not a generator
def test_variables_set(m):
    assert len(variables_set(m)) == 28


def test_number_variables(m):
    assert number_variables(m) == 28


def test_fixed_variables_set(m):
    assert len(fixed_variables_set(m)) == 2


def test_number_fixed_variables(m):
    assert number_fixed_variables(m) == 2


## -------------------------------------------------------------------------
## Variables in Constraints
#
#def variables_in_activated_constraints_set(block):
#    var_set = ComponentSet()
#    for c in block.component_data_objects(
#            ctype=Constraint, active=True, descend_into=True):
#        for v in identify_variables(c.body):
#            var_set.add(v)
#    return var_set
#
#
#def number_variables_in_activated_constraints(block):
#    return len(variables_in_activated_constraints_set(block))
#
#
#def variables_in_activated_equalities_set(block):
#    var_set = ComponentSet()
#    for c in activated_equalities_generator(block):
#        for v in identify_variables(c.body):
#            var_set.add(v)
#    return var_set
#
#
#def number_variables_in_activated_equalities(block):
#    return len(variables_in_activated_equalities_set(block))
#
#
#def variables_in_activated_inequalities_set(block):
#    var_set = ComponentSet()
#    for c in activated_inequalities_generator(block):
#        for v in identify_variables(c.body):
#            var_set.add(v)
#    return var_set
#
#
#def number_variables_in_activated_inequalities(block):
#    return len(variables_in_activated_inequalities_set(block))
#
#
#def variables_only_in_inequalities(block):
#    return (variables_in_activated_inequalities_set(block) -
#            variables_in_activated_equalities_set(block))
#
#
#def number_variables_only_in_inequalities(block):
#    return len(variables_only_in_inequalities(block))
#
#
## -------------------------------------------------------------------------
## Fixed Variables in Constraints
#def fixed_variables_in_activated_equalities_set(block):
#    var_set = ComponentSet()
#    for v in variables_in_activated_equalities_set(block):
#        if v.fixed:
#            var_set.add(v)
#    return var_set
#
#
#def number_fixed_variables_in_activated_equalities(block):
#    return len(fixed_variables_in_activated_equalities_set(block))
#
#
#def fixed_variables_only_in_inequalities(block):
#    var_set = ComponentSet()
#    for v in variables_only_in_inequalities(block):
#        if v.fixed:
#            var_set.add(v)
#    return var_set
#
#
#def number_fixed_variables_only_in_inequalities(block):
#    return len(fixed_variables_only_in_inequalities(block))
#
#
## -------------------------------------------------------------------------
## Unused and un-Transformed Variables
#def unused_variables_set(block):
#    return variables_set(block) - variables_in_activated_constraints_set(block)
#
#
#def number_unused_variables(block):
#    return len(unused_variables_set(block))
#
#
#def fixed_unused_variables_set(block):
#    var_set = ComponentSet()
#    for v in unused_variables_set(block):
#        if v.fixed:
#            var_set.add(v)
#    return var_set
#
#
#def number_fixed_unused_variables(block):
#    return len(fixed_unused_variables_set(block))
#
#
def test_derivative_variables_set():
    m = ConcreteModel()

    m.cs = ContinuousSet(bounds=(0, 1))

    m.v = Var(m.cs, initialize=1)
    m.dv = DerivativeVar(m.v)

    assert len(derivative_variables_set(m)) == 2

    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m,
                           nfe=10,
                           wrt=m.cs,
                           scheme="BACKWARD")

    assert len(derivative_variables_set(m)) == 0


def test_number_derivative_variables():
    m = ConcreteModel()

    m.cs = ContinuousSet(bounds=(0, 1))

    m.v = Var(m.cs, initialize=1)
    m.dv = DerivativeVar(m.v)

    assert number_derivative_variables(m) == 2

    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m,
                           nfe=10,
                           wrt=m.cs,
                           scheme="BACKWARD")

    assert number_derivative_variables(m) == 0


# -------------------------------------------------------------------------
# Objective methods
def test_total_objectives_set(m):
    assert len(total_objectives_set(m)) == 2


def test_number_total_objectives(m):
    assert number_total_objectives(m) == 2


def test_activated_objectives_set(m):
    assert len(activated_objectives_set(m)) == 1


def test_number_activated_objectives(m):
    assert number_activated_objectives(m) == 1


def test_deactivated_objectives_set(m):
    assert len(deactivated_objectives_set(m)) == 1


def test_number_deactivated_objectives(m):
    assert number_deactivated_objectives(m) == 1


# -------------------------------------------------------------------------
# Expression methods
def test_expressions_set(m):
    assert len(expressions_set(m)) == 3


def test_number_expressions(m):
    assert number_expressions(m) == 3


## -------------------------------------------------------------------------
## Other model statistics
#def degrees_of_freedom(block):
#    return (number_variables_in_activated_equalities(block) -
#            number_fixed_variables_in_activated_equalities(block) -
#            number_activated_equalities(block))
#
#
#def large_residuals_set(block, tol=1e-5):
#    large_residuals_set = ComponentSet()
#    for c in block.component_data_objects(
#            ctype=Constraint, active=True, descend_into=True):
#        if c.active and value(c.lower - c.body()) > tol:
#            large_residuals_set.add(c)
#        elif c.active and value(c.body() - c.upper) > tol:
#            large_residuals_set.add(c)
#    return large_residuals_set
#
#
#def number_large_residuals(block, tol=1e-5):
#    lr = 0
#    for c in block.component_data_objects(
#            ctype=Constraint, active=True, descend_into=True):
#        if c.active and value(c.lower - c.body()) > tol:
#            lr += 1
#        elif c.active and value(c.body() - c.upper) > tol:
#            lr += 1
#    return lr
#
#
#def active_variables_in_deactivated_blocks_set(block):
#    var_set = ComponentSet()
#    block_set = activated_blocks_set(block)
#    for v in variables_in_activated_constraints_set(block):
#        if v.parent_block() not in block_set:
#            var_set.add(v)
#    return var_set
#
#
#def number_active_variables_in_deactivated_blocks(block):
#    return len(active_variables_in_deactivated_blocks_set(block))


# -------------------------------------------------------------------------
# Reporting methods
def test_report_statistics(m):
    report_statistics(m)
