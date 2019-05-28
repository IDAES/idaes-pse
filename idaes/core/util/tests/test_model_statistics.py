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
        enumerate_blocks,
        enumerate_constraints,
        enumerate_derivative_variables,
        enumerate_expressions,
        enumerate_objectives,
        enumerate_variables,
        enumerate_equality_constraints,
        enumerate_inequality_constraints,
        enumerate_activated_components,
        enumerate_variables_in_constraints,
        enumerate_fixed_variables,
        calculate_degrees_of_freedom,
        report_model_statistics,
        enumerate_active_varaibles_in_deactived_blocks)
from idaes.core.util.exceptions import ConfigurationError


# Author: Andrew Lee
@pytest.fixture()
def m():
    m = ConcreteModel()

    m.s = Set(initialize=["a", "b"])
    m.cs = ContinuousSet(bounds=(0, 1))

    m.v = Var(m.cs)
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


def test_enumerate_blocks(m):
    assert len(enumerate_blocks(m, deactivated_blocks=True)) == 4

    assert len(enumerate_blocks(m,
                                deactivated_blocks=True,
                                descend_into=False)) == 3

    assert len(enumerate_blocks(m, deactivated_blocks=False)) == 2

    assert len(enumerate_blocks(m,
                                deactivated_blocks=False,
                                descend_into=False)) == 2

    m.b1.sb.activate()
    assert len(enumerate_blocks(m, deactivated_blocks=False)) == 2

    assert len(enumerate_blocks(m,
                                deactivated_blocks=False,
                                descend_into=False)) == 2

    m.b1.activate()
    m.b1.sb.deactivate()
    assert len(enumerate_blocks(m, deactivated_blocks=False)) == 3

    assert len(enumerate_blocks(m,
                                deactivated_blocks=False,
                                descend_into=False)) == 3


def test_enumerate_variables(m):
    assert len(enumerate_variables(m,
                                   deactivated_blocks=False,
                                   descend_into=True)) == 28

    assert len(enumerate_variables(m,
                                   deactivated_blocks=True,
                                   descend_into=True)) == 34

    assert len(enumerate_variables(m,
                                   deactivated_blocks=False,
                                   descend_into=False)) == 22

    assert len(enumerate_variables(m,
                                   deactivated_blocks=True,
                                   descend_into=False)) == 22


def test_enumerate_derivative_variables(m):
    assert len(enumerate_derivative_variables(m)) == 0

    m.v2 = Var(m.cs)
    m.dv2 = DerivativeVar(m.v2)
    assert len(enumerate_derivative_variables(m)) == 11


def test_enumerate_expressions(m):
    assert len(enumerate_expressions(m,
                                     deactivated_blocks=False,
                                     descend_into=True)) == 3

    assert len(enumerate_expressions(m,
                                     deactivated_blocks=True,
                                     descend_into=True)) == 4

    assert len(enumerate_expressions(m,
                                     deactivated_blocks=False,
                                     descend_into=False)) == 1

    assert len(enumerate_expressions(m,
                                     deactivated_blocks=True,
                                     descend_into=False)) == 1


def test_enumerate_objectives(m):
    assert len(enumerate_objectives(m,
                                    deactivated_blocks=False,
                                    descend_into=True)) == 2

    assert len(enumerate_objectives(m,
                                    deactivated_blocks=True,
                                    descend_into=True)) == 4


def test_enumerate_constraints(m):
    assert len(enumerate_constraints(m,
                                     deactivated_blocks=False,
                                     descend_into=True)) == 14

    assert len(enumerate_constraints(m,
                                     deactivated_blocks=True,
                                     descend_into=True)) == 18

    assert len(enumerate_constraints(m,
                                     deactivated_blocks=False,
                                     descend_into=False)) == 10


def test_enumerate_equality_constraints(m):
    assert len(enumerate_equality_constraints(m,
                                              deactivated_blocks=False,
                                              descend_into=True)) == 12

    assert len(enumerate_equality_constraints(m,
                                              deactivated_blocks=True,
                                              descend_into=True)) == 14

    assert len(enumerate_equality_constraints(m,
                                              deactivated_blocks=False,
                                              descend_into=False)) == 10

    c = enumerate_constraints(m, deactivated_blocks=False, descend_into=True)
    assert len(enumerate_equality_constraints(c,
                                              deactivated_blocks=False,
                                              descend_into=False)) == 12
    assert len(enumerate_equality_constraints(c,
                                              deactivated_blocks=True,
                                              descend_into=False)) == 12
    assert len(enumerate_equality_constraints(c,
                                              deactivated_blocks=False,
                                              descend_into=True)) == 12

    c = enumerate_constraints(m, deactivated_blocks=True, descend_into=True)
    assert len(enumerate_equality_constraints(c)) == 14

    c = enumerate_constraints(m, deactivated_blocks=False, descend_into=False)
    assert len(enumerate_equality_constraints(c)) == 10


def test_enumerate_inequality_constraints(m):
    assert len(enumerate_inequality_constraints(m,
                                                deactivated_blocks=False,
                                                descend_into=True)) == 2

    assert len(enumerate_inequality_constraints(m,
                                                deactivated_blocks=True,
                                                descend_into=True)) == 4

    assert len(enumerate_inequality_constraints(m,
                                                deactivated_blocks=False,
                                                descend_into=False)) == 0

    c = enumerate_constraints(m, deactivated_blocks=False, descend_into=True)
    assert len(enumerate_inequality_constraints(c,
                                                deactivated_blocks=False,
                                                descend_into=False)) == 2
    assert len(enumerate_inequality_constraints(c,
                                                deactivated_blocks=True,
                                                descend_into=False)) == 2
    assert len(enumerate_inequality_constraints(c,
                                                deactivated_blocks=False,
                                                descend_into=True)) == 2

    c = enumerate_constraints(m, deactivated_blocks=True, descend_into=True)
    assert len(enumerate_inequality_constraints(c)) == 4

    c = enumerate_constraints(m, deactivated_blocks=False, descend_into=False)
    assert len(enumerate_inequality_constraints(c)) == 0


def test_enumerate_activated_components(m):
    c = enumerate_constraints(m, deactivated_blocks=False, descend_into=True)
    assert len(enumerate_activated_components(c)) == 14

    m.b2["a"].c1.deactivate()
    assert len(enumerate_activated_components(c)) == 13


def test_enumerate_variables_in_constraints(m):
    c = enumerate_constraints(m, deactivated_blocks=False, descend_into=True)
    assert len(enumerate_variables_in_constraints(c)) == 23

    assert len(enumerate_variables_in_constraints(m,
                                                  deactivated_blocks=False,
                                                  descend_into=True)) == 23

    c = enumerate_constraints(m, deactivated_blocks=True, descend_into=True)
    assert len(enumerate_variables_in_constraints(c)) == 25

    c = enumerate_constraints(m, deactivated_blocks=False, descend_into=False)
    assert len(enumerate_variables_in_constraints(c)) == 21


def test_enumerate_fixed_variables(m):
    assert len(enumerate_fixed_variables(m)) == 2

    v = enumerate_variables(m, deactivated_blocks=True)
    assert len(enumerate_fixed_variables(v)) == 4


def test_calculate_degrees_of_freedom(m):
    assert calculate_degrees_of_freedom(m) == 9


def test_report_model_statistics(m):
    # Run method to ensure no exceptions occur
    report_model_statistics(m)


def test_act_vars_deact_blocks(m):
    assert len(enumerate_active_varaibles_in_deactived_blocks(m)) == 0

    m.c1 = Constraint(expr=10 == m.b1.v1)

    assert len(enumerate_active_varaibles_in_deactived_blocks(m)) == 1

    with pytest.raises(ConfigurationError):
        enumerate_active_varaibles_in_deactived_blocks(m.b1)
