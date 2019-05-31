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

from idaes.core.util.model_statistics import (
        block_set,
        constraint_set,
        derivative_variables_set,
        expression_set,
        objective_set,
        variable_set,
        equality_constraint_set,
        inequality_constraints_set,
        activated_component_set,
        variables_in_constraints_set,
        fixed_variable_set,
        unfixed_variable_set,
        calculate_degrees_of_freedom,
        report_model_statistics,
        active_variables_in_deactived_blocks_set,
        large_residual_set)
from idaes.core.util.exceptions import ConfigurationError


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

    return m


def test_block_set(m):
    assert len(block_set(m, active=None)) == 4

    assert len(block_set(m, active=None, descend_into=False)) == 3

    assert len(block_set(m, active=True)) == 2

    assert len(block_set(m, active=True, descend_into=False)) == 2

    m.b1.sb.activate()
    assert len(block_set(m, active=True)) == 2

    assert len(block_set(m, active=True, descend_into=False)) == 2

    m.b1.activate()
    m.b1.sb.deactivate()
    assert len(block_set(m, active=True)) == 3

    assert len(block_set(m, active=True, descend_into=False)) == 3


def test_variable_set(m):
    assert len(variable_set(m, active=True, descend_into=True)) == 28

    assert len(variable_set(m, active=None, descend_into=True)) == 34

    assert len(variable_set(m, active=True, descend_into=False)) == 22

    assert len(variable_set(m, active=None, descend_into=False)) == 22


def test_derivative_variables_set(m):
    assert len(derivative_variables_set(m)) == 0

    m.v2 = Var(m.cs)
    m.dv2 = DerivativeVar(m.v2)
    assert len(derivative_variables_set(m)) == 11


def test_expression_set(m):
    assert len(expression_set(m, active_blocks=True, descend_into=True)) == 3

    assert len(expression_set(m, active_blocks=None, descend_into=True)) == 4

    assert len(expression_set(m, active_blocks=True, descend_into=False)) == 1

    assert len(expression_set(m, active_blocks=None, descend_into=False)) == 1


def test_objective_set(m):
    assert len(objective_set(m, active_blocks=True, descend_into=True)) == 2

    assert len(objective_set(m, active_blocks=None, descend_into=True)) == 4


def test_constraint_set(m):
    assert len(constraint_set(m, active_blocks=True, descend_into=True)) == 14

    assert len(constraint_set(m, active_blocks=None, descend_into=True)) == 18

    assert len(constraint_set(m, active_blocks=True, descend_into=False)) == 10


def test_equality_constraint_set(m):
    assert len(equality_constraint_set(m,
                                       active_blocks=True,
                                       descend_into=True)) == 12

    assert len(equality_constraint_set(m,
                                       active_blocks=None,
                                       descend_into=True)) == 14

    assert len(equality_constraint_set(m,
                                       active_blocks=True,
                                       descend_into=False)) == 10

    c = constraint_set(m, active_blocks=True, descend_into=True)
    assert len(equality_constraint_set(c,
                                       active_blocks=True,
                                       descend_into=False)) == 12
    assert len(equality_constraint_set(c,
                                       active_blocks=None,
                                       descend_into=False)) == 12
    assert len(equality_constraint_set(c,
                                       active_blocks=True,
                                       descend_into=True)) == 12

    c = constraint_set(m, active_blocks=None, descend_into=True)
    assert len(equality_constraint_set(c)) == 14

    c = constraint_set(m, active_blocks=True, descend_into=False)
    assert len(equality_constraint_set(c)) == 10


def test_inequality_constraints_set(m):
    assert len(inequality_constraints_set(m,
                                          active_blocks=True,
                                          descend_into=True)) == 2

    assert len(inequality_constraints_set(m,
                                          active_blocks=None,
                                          descend_into=True)) == 4

    assert len(inequality_constraints_set(m,
                                          active_blocks=True,
                                          descend_into=False)) == 0

    c = constraint_set(m, active_blocks=True, descend_into=True)
    assert len(inequality_constraints_set(c,
                                          active_blocks=True,
                                          descend_into=False)) == 2
    assert len(inequality_constraints_set(c,
                                          active_blocks=None,
                                          descend_into=False)) == 2
    assert len(inequality_constraints_set(c,
                                          active_blocks=True,
                                          descend_into=True)) == 2

    c = constraint_set(m, active_blocks=None, descend_into=True)
    assert len(inequality_constraints_set(c)) == 4

    c = constraint_set(m, active_blocks=True, descend_into=False)
    assert len(inequality_constraints_set(c)) == 0


def test_activated_component_set(m):
    c = constraint_set(m, active_blocks=True, descend_into=True)
    assert len(activated_component_set(c)) == 14

    m.b2["a"].c1.deactivate()
    assert len(activated_component_set(c)) == 13


def test_variables_in_constraints_set(m):
    c = constraint_set(m, active_blocks=True, descend_into=True)
    assert len(variables_in_constraints_set(c)) == 23

    assert len(variables_in_constraints_set(m,
                                            active_blocks=True,
                                            descend_into=True)) == 23

    c = constraint_set(m, active_blocks=None, descend_into=True)
    assert len(variables_in_constraints_set(c)) == 25

    c = constraint_set(m, active_blocks=True, descend_into=False)
    assert len(variables_in_constraints_set(c)) == 21


def test_fixed_variable_set(m):
    assert len(fixed_variable_set(m)) == 2

    v = variable_set(m, active=None)
    assert len(fixed_variable_set(v)) == 4


def test_unfixed_variable_set(m):
    assert len(unfixed_variable_set(m)) == 26

    v = variable_set(m, active=None)
    assert len(unfixed_variable_set(v)) == 30


def test_calculate_degrees_of_freedom(m):
    assert calculate_degrees_of_freedom(m) == 9


def test_large_residual_set(m):
    # Initialize derivative var values so no errors occur
    for v in m.dv.keys():
        m.dv[v] = 0
    assert len(large_residual_set(m)) == 4

    m.b2["a"].deactivate()
    assert len(large_residual_set(m)) == 2


def test_report_model_statistics(m):
    # Run method to ensure no exceptions occur
    report_model_statistics(m)


def test_act_vars_deact_blocks(m):
    assert len(active_variables_in_deactived_blocks_set(m)) == 0

    m.c1 = Constraint(expr=10 == m.b1.v1)

    assert len(active_variables_in_deactived_blocks_set(m)) == 1

    with pytest.raises(ConfigurationError):
        active_variables_in_deactived_blocks_set(m.b1)
