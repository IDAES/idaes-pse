##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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

from pyomo.environ import ConcreteModel, Expression, Set, Block, Var
from pyomo.network import Port

from idaes.core.util.misc import (add_object_reference, copy_port_values,
                                  TagReference, VarLikeExpression)


# Author: Andrew Lee
@pytest.mark.unit
def test_add_object_reference():
    m = ConcreteModel()

    m.s = Set(initialize=[1, 2, 3])

    add_object_reference(m, "test_ref", m.s)

    assert hasattr(m, "test_ref")
    assert m.test_ref == m.s


# Author: Andrew Lee
@pytest.mark.unit
def test_add_object_reference_fail():
    m = ConcreteModel()

    with pytest.raises(AttributeError):
        add_object_reference(m, "test_ref", m.s)


# Author: John Eslick
@pytest.mark.unit
def test_port_copy():
    m = ConcreteModel()
    m.b1 = Block()
    m.b2 = Block()
    m.b1.x = Var(initialize=3)
    m.b1.y = Var([0, 1], initialize={0: 4, 1: 5})
    m.b1.z = Var([0, 1], ["A", "B"], initialize={
        (0, "A"): 6, (0, "B"): 7, (1, "A"): 8, (1, "B"): 9})
    m.b2.x = Var(initialize=1)
    m.b2.y = Var([0, 1], initialize=1)
    m.b2.z = Var([0, 1], ["A", "B"], initialize=1)
    m.b1.port = Port()
    m.b2.port = Port()
    m.b1.port.add(m.b1.x, "x")
    m.b1.port.add(m.b1.y, "y")
    m.b1.port.add(m.b1.z, "z")
    m.b2.port.add(m.b2.x, "x")
    m.b2.port.add(m.b2.y, "y")
    m.b2.port.add(m.b2.z, "z")
    copy_port_values(m.b2.port, m.b1.port)
    assert(m.b2.x == 3)
    assert(m.b2.y[0] == 4)
    assert(m.b2.y[1] == 5)
    assert(m.b2.z[0, "A"] == 6)
    assert(m.b2.z[0, "B"] == 7)
    assert(m.b2.z[1, "A"] == 8)
    assert(m.b2.z[1, "B"] == 9)


# Author: John Eslick
@pytest.mark.unit
def test_tag_reference():
    m = ConcreteModel()
    m.z = Var([0, 1], ["A", "B"], initialize={
        (0, "A"): 6, (0, "B"): 7, (1, "A"): 8, (1, "B"): 9})
    test_tag = {}
    test_tag["MyTag34&@!e.5"] = TagReference(m.z[:, "A"], description="z tag")
    assert(len(test_tag["MyTag34&@!e.5"]) == 2)
    assert(test_tag["MyTag34&@!e.5"][0].value == 6)
    assert(test_tag["MyTag34&@!e.5"][1].value == 8)
    assert(test_tag["MyTag34&@!e.5"].description == "z tag")
    m.b = Block([0, 1])
    m.b[0].y = Var(initialize=1)
    m.b[1].y = Var(initialize=2)
    test_tag = TagReference(m.b[:].y, description="y tag")
    assert(test_tag[0] == 1)
    assert(test_tag[1] == 2)
    assert(test_tag.description == "y tag")


@pytest.mark.unit
def test_SimpleVarLikeExpression():
    m = ConcreteModel()

    # Need a Var to use in the Expression to avoid being able to set the value
    # of a float
    m.v = Var(initialize=42)

    m.e = VarLikeExpression(expr=m.v)

    assert m.e.type() is Expression
    assert not m.e.is_indexed()
    assert m.e.value == 42

    with pytest.raises(TypeError,
                       match="e is an Expression and does not have a value "
                       "which can be set."):
        m.e.value = 10

    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use an inequality Constraint instead."):
        m.e.setub(10)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use an inequality Constraint instead."):
        m.e.setlb(0)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be fixed. "
                       "Use an equality Constraint instead."):
        m.e.fix(8)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be unfixed."):
        m.e.unfix()


@pytest.mark.unit
def test_IndexedVarLikeExpression():
    m = ConcreteModel()

    # Need a Var to use in the Expression to avoid being able to set the value
    # of a float
    m.v = Var(initialize=42)

    m.e = VarLikeExpression([1, 2, 3, 4], expr=m.v)

    assert m.e.type() is Expression
    assert m.e.is_indexed()

    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use inequality Constraints instead."):
        m.e.setub(10)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not have bounds. "
                       "Use inequality Constraints instead."):
        m.e.setlb(0)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be fixed. "
                       "Use equality Constraints instead."):
        m.e.fix(8)
    with pytest.raises(TypeError,
                       match="e is an Expression and can not be unfixed."):
        m.e.unfix()

    for i in m.e:
        assert m.e[i].value == 42

        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and does not have "
                           "a value which can be set.".format(i)):
            m.e[i].value = 10

        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not have "
                           "bounds. Use an inequality Constraint instead."
                           .format(i)):
            m.e[i].setub(10)
        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not have "
                           "bounds. Use an inequality Constraint instead."
                           .format(i)):
            m.e[i].setlb(0)
        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not be "
                           "fixed. Use an equality Constraint instead."
                           .format(i)):
            m.e[i].fix(8)
        with pytest.raises(TypeError,
                           match="e\[{}\] is an Expression and can not be "
                           "unfixed.".format(i)):
            m.e[i].unfix()
