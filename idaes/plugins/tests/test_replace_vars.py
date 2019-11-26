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
import pyomo.environ as pyo
import idaes

__author__ = "John Eslick"

def test_1():
    # Test scalar variables
    rp = pyo.TransformationFactory("replace_variables")
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=2)
    m.y = pyo.Var(initialize=3)
    m.z = pyo.Var(initialize=0)
    m.c1 = pyo.Constraint(expr=m.z==m.x+m.y)
    m.e1 = pyo.Expression(expr=m.x**m.y)

    assert(m.c1.body() == -5)
    assert(pyo.value(m.e1) == 8)
    rp.apply_to(m, substitute=[(m.y, m.x)])
    assert(m.c1.body() == -4)
    assert(pyo.value(m.e1) == 4)

def test_2():
    # Test vector variables and sums
    rp = pyo.TransformationFactory("replace_variables")
    m = pyo.ConcreteModel()
    m.x = pyo.Var(["a", "b", "c"], initialize=2)
    m.y = pyo.Var(initialize=3)
    m.z = pyo.Var(initialize=0)
    m.c1 = pyo.Constraint(expr=m.z==m.x["a"] + m.x["b"] + m.x["c"])
    m.e1 = pyo.Expression(expr=sum(m.x[i] for i in m.x))

    assert(m.c1.body() == -6)
    assert(pyo.value(m.e1) == 6)
    rp.apply_to(m, substitute=[(m.x["c"], m.y)])
    assert(m.c1.body() == -7)
    assert(pyo.value(m.e1) == 7)
