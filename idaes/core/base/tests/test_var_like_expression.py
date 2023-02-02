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
This module tests for the VarLikeExpression component.
"""

import pytest

from pyomo.environ import ConcreteModel, Expression, Var

from idaes.core.base.var_like_expression import VarLikeExpression


@pytest.mark.unit
def test_SimpleVarLikeExpression():
    m = ConcreteModel()

    # Need a Var to use in the Expression to avoid being able to set the value
    # of a float
    m.v = Var(initialize=42)

    m.e = VarLikeExpression(expr=m.v)

    assert m.e.ctype is Expression
    assert not m.e.is_indexed()

    with pytest.raises(
        TypeError,
        match="e is an Expression and does not have a value "
        "attribute. Use the 'value\(\)' method instead.",
    ):
        assert m.e.value == 42

    with pytest.raises(
        TypeError,
        match="e is an Expression and does not have a value " "which can be set.",
    ):
        m.e.set_value(10)

    with pytest.raises(
        TypeError,
        match="e is an Expression and does not have a value " "which can be set.",
    ):
        m.e.value = 10

    with pytest.raises(
        TypeError,
        match="e is an Expression and can not have bounds. "
        "Use an inequality Constraint instead.",
    ):
        m.e.setub(10)
    with pytest.raises(
        TypeError,
        match="e is an Expression and can not have bounds. "
        "Use an inequality Constraint instead.",
    ):
        m.e.setlb(0)
    with pytest.raises(
        TypeError,
        match="e is an Expression and can not be fixed. "
        "Use an equality Constraint instead.",
    ):
        m.e.fix(8)
    with pytest.raises(TypeError, match="e is an Expression and can not be unfixed."):
        m.e.unfix()

    m.e.set_value(10, force=True)
    assert m.e._expr == 10


@pytest.mark.unit
def test_IndexedVarLikeExpression():
    m = ConcreteModel()

    # Need a Var to use in the Expression to avoid being able to set the value
    # of a float
    m.v = Var(initialize=42)

    m.e = VarLikeExpression([1, 2, 3, 4], expr=m.v)

    assert m.e.ctype is Expression
    assert m.e.is_indexed()

    with pytest.raises(
        TypeError,
        match="e is an Expression and can not have bounds. "
        "Use inequality Constraints instead.",
    ):
        m.e.setub(10)
    with pytest.raises(
        TypeError,
        match="e is an Expression and can not have bounds. "
        "Use inequality Constraints instead.",
    ):
        m.e.setlb(0)
    with pytest.raises(
        TypeError,
        match="e is an Expression and can not be fixed. "
        "Use equality Constraints instead.",
    ):
        m.e.fix(8)
    with pytest.raises(TypeError, match="e is an Expression and can not be unfixed."):
        m.e.unfix()

    for i in m.e:
        with pytest.raises(
            TypeError,
            match=f"e\[{i}\] is an Expression and does not have"
            f" a value attribute. Use the 'value\(\)' method "
            "instead",
        ):
            assert m.e[i].value == 42

        with pytest.raises(
            TypeError,
            match=f"e\[{i}\] is an Expression and does not "
            f"have a value which can be set.",
        ):
            m.e[i].set_value(10)

        with pytest.raises(
            TypeError,
            match="e\[{}\] is an Expression and does not have "
            "a value which can be set.".format(i),
        ):
            m.e[i].value = 10

        with pytest.raises(
            TypeError,
            match="e\[{}\] is an Expression and can not have "
            "bounds. Use an inequality Constraint instead.".format(i),
        ):
            m.e[i].setub(10)
        with pytest.raises(
            TypeError,
            match="e\[{}\] is an Expression and can not have "
            "bounds. Use an inequality Constraint instead.".format(i),
        ):
            m.e[i].setlb(0)
        with pytest.raises(
            TypeError,
            match="e\[{}\] is an Expression and can not be "
            "fixed. Use an equality Constraint instead.".format(i),
        ):
            m.e[i].fix(8)
        with pytest.raises(
            TypeError,
            match="e\[{}\] is an Expression and can not be " "unfixed.".format(i),
        ):
            m.e[i].unfix()

        m.e[i].set_value(i, force=True)
        assert m.e[i]._expr == i
