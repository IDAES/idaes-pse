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
Tests for dynamic utility methods.
"""

import pytest
from pyomo.environ import (
    ConcreteModel,
    Block,
    Constraint,
    Var,
    Set,
    TransformationFactory,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.collections import ComponentSet
import idaes.logger as idaeslog
from idaes.core.util.dyn_utils import *

__author__ = "Robert Parker"


@pytest.mark.unit
def test_fix_and_deactivate():
    m = ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))
    m.space = ContinuousSet(bounds=(0, 5))
    m.set1 = Set(initialize=["a", "b", "c"])
    m.set2 = Set(initialize=["d", "e", "f"])
    m.fs = Block()

    m.fs.v0 = Var(m.space, initialize=1)

    @m.fs.Block()
    def b1(b):
        b.v = Var(m.time, m.space, initialize=1)
        b.dv = DerivativeVar(b.v, wrt=m.time)

        b.con = Constraint(
            m.time, m.space, rule=lambda b, t, x: b.dv[t, x] == 7 - b.v[t, x]
        )

        @b.Block(m.time)
        def b2(b, t):
            b.v = Var(initialize=2)

    @m.fs.Block(m.time, m.space)
    def b2(b, t, x):
        b.v = Var(m.set1, initialize=2)

        @b.Block(m.set1)
        def b3(b, c):
            b.v = Var(m.set2, initialize=3)

            @b.Constraint(m.set2)
            def con(b, s):
                return 5 * b.v[s] == m.fs.b2[m.time.first(), m.space.first()].v[c]

    @m.fs.Constraint(m.time)
    def con1(fs, t):
        return fs.b1.v[t, m.space.last()] == 5

    @m.fs.Constraint(m.space)
    def con2(fs, x):
        return fs.b1.v[m.time.first(), x] == fs.v0[x]

    disc = TransformationFactory("dae.collocation")
    disc.apply_to(m, wrt=m.time, nfe=5, ncp=2, scheme="LAGRANGE-RADAU")
    disc.apply_to(m, wrt=m.space, nfe=5, ncp=2, scheme="LAGRANGE-RADAU")

    for t in m.time:
        m.fs.b1.v[t, m.space.first()].fix()

    active_dict = get_activity_dict(m.fs)
    for comp in m.fs.component_data_objects(Constraint, Block):
        assert active_dict[id(comp)] is True

    deactivate_model_at(m, m.time, m.time.at(2))
    assert m.fs.con1[m.time.at(1)].active
    assert not m.fs.con1[m.time.at(2)].active
    assert m.fs.con2[m.space.at(1)].active
    assert not m.fs.b1.con[m.time.at(2), m.space.at(1)].active
    assert not m.fs.b2[m.time.at(2), m.space.last()].active
    assert m.fs.b2[m.time.at(2), m.space.last()].b3["a"].con["e"].active

    deactivate_model_at(m, m.time, [m.time.at(1), m.time.at(3)], outlvl=idaeslog.ERROR)
    # Higher outlvl threshold as will encounter warning trying to deactivate
    # disc equations at time.first()
    assert not m.fs.con1[m.time.at(1)].active
    assert not m.fs.con1[m.time.at(3)].active
    assert not m.fs.b1.con[m.time.at(1), m.space.at(1)].active
    assert not m.fs.b1.con[m.time.at(3), m.space.at(1)].active

    init_derivs = get_derivatives_at(m, m.time, m.time.first())[m.time.first()]
    init_deriv_names = [var.name for var in init_derivs]

    assert m.fs.b1.dv[m.time.first(), m.space.first()] in init_derivs
    assert m.fs.b1.dv[m.time.at(1), m.space.at(1)].name in init_deriv_names

    deriv_dict = get_derivatives_at(m, m.time, [m.time.first(), m.time.last()])
    deriv_name_dict = {t: [d.name for d in deriv_dict[t]] for t in deriv_dict.keys()}
    assert m.time.first() in deriv_name_dict.keys()
    assert m.time.last() in deriv_name_dict.keys()
    assert (
        m.fs.b1.dv[m.time.last(), m.space.at(1)].name in deriv_name_dict[m.time.last()]
    )
    assert (
        m.fs.b1.dv[m.time.last(), m.space.at(1)].name
        not in deriv_name_dict[m.time.first()]
    )

    vars_unindexed = fix_vars_unindexed_by(m, m.time)
    cons_unindexed = deactivate_constraints_unindexed_by(m, m.time)

    unindexed_vars = ComponentSet(vars_unindexed)
    assert m.fs.v0[m.space.at(1)] in unindexed_vars
    assert m.fs.b1.b2[m.time.at(1)].v not in unindexed_vars
    assert m.fs.con2[m.space.at(2)] in cons_unindexed
    assert m.fs.con1[m.time.at(1)] not in cons_unindexed
    assert not m.fs.con2[m.space.at(1)].active
    assert m.fs.v0[m.space.at(1)].fixed

    path = path_from_block(
        m.fs.b2[m.time.at(1), m.space.at(1)].b3["a"].v, m, include_comp=False
    )
    assert path == [("fs", None), ("b2", (m.time.at(1), m.space.at(1))), ("b3", "a")]
    path = path_from_block(
        m.fs.b2[m.time.at(1), m.space.at(1)].b3["a"].v, m, include_comp=True
    )
    assert path == [
        ("fs", None),
        ("b2", (m.time.at(1), m.space.at(1))),
        ("b3", "a"),
        ("v", None),
    ]
    path = path_from_block(
        m.fs.b2[m.time.at(1), m.space.at(1)].b3["a"].v["f"], m, include_comp=True
    )
    assert path == [
        ("fs", None),
        ("b2", (m.time.at(1), m.space.at(1))),
        ("b3", "a"),
        ("v", "f"),
    ]
    path = path_from_block(
        m.fs.b2[m.time.at(1), m.space.at(1)].b3["a"].v["f"],
        m.fs.b2[m.time.at(1), m.space.at(1)],
        include_comp=True,
    )
    assert path == [("b3", "a"), ("v", "f")]
    path = path_from_block(m.fs.b1.con[m.time.at(1), m.space.at(1)], m.fs)
    assert path == [("b1", None)]

    m.fs.b1.b2[m.time.at(1)].v.set_value(-1)
    for x in m.space:
        m.fs.b1.v[m.time.at(1), x].set_value(-1)
        m.fs.b1.dv[m.time.at(1), x].set_value(-1)
        for c1 in m.set1:
            m.fs.b2[m.time.at(1), x].v[c1].set_value(-1)
            for c2 in m.set2:
                m.fs.b2[m.time.at(1), x].b3[c1].v[c2].set_value(-1)

    copy_values_at_time(m, m, m.time.last(), m.time.at(1), copy_fixed=False)
    assert m.fs.b1.b2[m.time.last()].v.value == -1
    for x in m.space:
        if x != m.space.first():
            assert m.fs.b1.v[m.time.last(), x].value == -1
        else:
            assert m.fs.b1.v[m.time.last(), x].value == 1
            assert m.fs.b1.v[m.time.last(), x].fixed
        assert m.fs.b1.dv[m.time.last(), x].value == -1
        for c1 in m.set1:
            assert m.fs.b2[m.time.last(), x].v[c1].value == -1
            for c2 in m.set2:
                assert m.fs.b2[m.time.at(1), x].b3[c1].v[c2].value == -1


@pytest.mark.unit
def test_copy_non_time_indexed_values():
    m1 = ConcreteModel()
    m1.time = Set(initialize=[1, 2, 3, 4, 5])
    m1.v1 = Var(m1.time, initialize=1)
    m1.v2 = Var(initialize=1)

    @m1.Block(["a", "b"])
    def b1(b, i):
        b.v3 = Var(initialize=1)

        @b.Block(m1.time)
        def b2(b, t):
            b.v4 = Var(initialize=1)

        @b.Block()
        def b4(b):
            b.v6 = Var(initialize=1)

    @m1.Block(m1.time)
    def b3(b, t):
        b.v5 = Var(initialize=1)

    m2 = ConcreteModel()
    m2.time = Set(initialize=[1, 2, 3, 4, 5])
    m2.v1 = Var(m2.time, initialize=2)
    m2.v2 = Var(initialize=2)

    @m2.Block(["a", "b"])
    def b1(b):
        b.v3 = Var(initialize=2)

        @b.Block(m2.time)
        def b2(b, t):
            b.v4 = Var(initialize=2)

        @b.Block()
        def b4(b):
            b.v6 = Var(initialize=2)

    @m2.Block(m1.time)
    def b3(b, t):
        b.v5 = Var(initialize=2)

    copy_non_time_indexed_values(m1, m2)
    assert m1.v1[1].value != m2.v1[1].value
    assert m1.v2.value == m2.v2.value == 2
    assert m1.b1["a"].v3.value == m2.b1["a"].v3.value == 2
    assert m1.b1["b"].b4.v6.value == m2.b1["b"].b4.v6.value == 2
    assert m1.b3[3].v5.value != m2.b3[3].v5.value


@pytest.mark.unit
def test_find_comp_in_block():
    m1 = ConcreteModel()

    @m1.Block([1, 2, 3])
    def b1(b):
        b.v = Var([1, 2, 3])

    m2 = ConcreteModel()

    @m2.Block([1, 2, 3])
    def b1(b):
        b.v = Var([1, 2, 3, 4])

    @m2.Block([1, 2, 3])
    def b2(b):
        b.v = Var([1, 2, 3])

    v1 = m1.b1[1].v[1]

    assert find_comp_in_block(m2, m1, v1) is m2.b1[1].v[1]

    v2 = m2.b2[1].v[1]
    v3 = m2.b1[3].v[4]

    # These should result in Attribute/KeyErrors

    with pytest.raises(AttributeError, match=r".*has no attribute.*"):
        find_comp_in_block(m1, m2, v2)
    with pytest.raises(KeyError, match=r".*is not a valid index.*"):
        find_comp_in_block(m1, m2, v3)
    assert find_comp_in_block(m1, m2, v2, allow_miss=True) is None
    assert find_comp_in_block(m1, m2, v3, allow_miss=True) is None


@pytest.mark.unit
def test_find_comp_in_block_at_time():
    m1 = ConcreteModel()
    m1.time = Set(initialize=[1, 2, 3])

    @m1.Block(m1.time)
    def b1(b):
        b.v = Var(m1.time)

    @m1.Block([1, 2, 3])
    def b(bl):
        bl.v = Var(m1.time)
        bl.v2 = Var(m1.time, ["a", "b", "c"])

    m2 = ConcreteModel()
    m2.time = Set(initialize=[1, 2, 3, 4, 5, 6])

    @m2.Block(m2.time)
    def b1(b):
        b.v = Var(m2.time)

    @m2.Block([1, 2, 3, 4])
    def b(bl):
        bl.v = Var(m2.time)
        bl.v2 = Var(m2.time, ["a", "b", "c"])

    @m2.Block([1, 2, 3])
    def b2(b):
        b.v = Var(m2.time)

    v1 = m1.b1[1].v[1]
    v3 = m1.b[3].v[1]

    assert find_comp_in_block_at_time(m2, m1, v1, m2.time, 4) is m2.b1[4].v[4]
    assert find_comp_in_block_at_time(m2, m1, v3, m2.time, 5) is m2.b[3].v[5]

    v = m1.b[1].v2[1, "a"]
    assert find_comp_in_block_at_time(m2, m1, v, m2.time, 6) is m2.b[1].v2[6, "a"]

    v2 = m2.b2[1].v[1]
    v4 = m2.b[4].v[1]

    # Should result in exceptions:
    with pytest.raises(AttributeError, match=r".*has no attribute.*"):
        find_comp_in_block_at_time(m1, m2, v2, m1.time, 3)
    with pytest.raises(KeyError, match=r".*is not a valid index.*"):
        find_comp_in_block_at_time(m1, m2, v4, m1.time, 3)

    assert find_comp_in_block_at_time(m1, m2, v2, m1.time, 3, allow_miss=True) is None
    assert find_comp_in_block_at_time(m1, m2, v4, m1.time, 3, allow_miss=True) is None


@pytest.mark.unit
def test_get_location_of_coordinate_set():
    m = ConcreteModel()
    m.s1 = Set(initialize=[1, 2, 3])
    m.s2 = Set(initialize=[("a", 1), ("b", 2)])
    m.s3 = Set(initialize=[("a", 1, 0), ("b", 2, 1)])
    m.v1 = Var(m.s1)
    m.v2 = Var(m.s1, m.s2)
    m.v121 = Var(m.s1, m.s2, m.s1)
    m.v3 = Var(m.s3, m.s1, m.s2)

    assert get_location_of_coordinate_set(m.v1.index_set(), m.s1) == 0
    assert get_location_of_coordinate_set(m.v2.index_set(), m.s1) == 0
    assert get_location_of_coordinate_set(m.v3.index_set(), m.s1) == 3

    # These should each raise a ValueError
    with pytest.raises(ValueError):
        get_location_of_coordinate_set(m.v1.index_set(), m.s2)
    with pytest.raises(ValueError):
        get_location_of_coordinate_set(m.v121.index_set(), m.s1)


@pytest.mark.unit
def test_get_index_of_set():
    m = ConcreteModel()
    m.s1 = Set(initialize=[1, 2, 3])
    m.s2 = Set(initialize=[("a", 1), ("b", 2)])
    m.s3 = Set(initialize=[("a", 1, 0), ("b", 2, 1)])
    m.v0 = Var()
    m.v1 = Var(m.s1)
    m.v2 = Var(m.s1, m.s2)
    m.v121 = Var(m.s1, m.s2, m.s1)
    m.v3 = Var(m.s3, m.s1, m.s2)

    assert get_index_of_set(m.v1[2], m.s1) == 2
    assert get_index_of_set(m.v2[2, "a", 1], m.s1) == 2
    assert get_index_of_set(m.v3["b", 2, 1, 3, "b", 2], m.s1) == 3

    with pytest.raises(ValueError) as exc_test:
        get_index_of_set(m.v0, m.s1)
    with pytest.raises(ValueError) as exc_test:
        get_index_of_set(m.v2[1, "a", 1], m.s2)
    with pytest.raises(ValueError) as exc_test:
        get_index_of_set(m.v2[1, "b", 2], m.s3)


@pytest.mark.unit
def test_get_implicit_index_of_set():
    m = ConcreteModel()
    m.s1 = Set(initialize=[1, 2, 3])
    m.s2 = Set(initialize=["a", "b", "c"])
    m.s3 = Set(initialize=[("d", 4), ("e", 5)])

    @m.Block()
    def b1(b):
        @b.Block(m.s3, m.s1)
        def b2(b):
            @b.Block()
            def b3(b):
                b.v1 = Var(m.s2)
                b.v2 = Var(m.s1)

    assert get_implicit_index_of_set(m.b1.b2["d", 4, 1].b3.v1["a"], m.s1) == 1
    assert get_implicit_index_of_set(m.b1.b2["d", 4, 1].b3.v1["a"], m.s2) == "a"

    with pytest.raises(ValueError) as exc_test:
        get_implicit_index_of_set(m.b1.b2["e", 5, 2].b3.v2[1], m.s1)
