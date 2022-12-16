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

__Author__ = "John Eslick"


import pytest
import pyomo.environ as pyo
from pyomo.environ import ConcreteModel, Set, Block, Var
from pyomo.network import Port
from idaes.core.util import ModelTag, ModelTagGroup
from idaes.core.util.tags import svg_tag


@pytest.fixture()
def model():
    m = pyo.ConcreteModel()
    m.w = pyo.Var([1, 2, 3], ["a", "b"], initialize=4, units=pyo.units.kg)
    m.x = pyo.Var([1, 2, 3], initialize=5, units=pyo.units.kg)
    m.y = pyo.Var(initialize=6, units=pyo.units.s)
    m.z = pyo.Var(initialize=7)
    m.e = pyo.Expression(expr=m.w[1, "a"] / m.x[1])
    m.f = pyo.Expression(expr=m.x[1] / m.y)

    @m.Expression([1, 2, 3], ["a", "b"])
    def g(b, i, j):
        return m.w[i, j] / m.x[i] * 100

    return m


@pytest.mark.unit
def test_tag_display(model):
    m = model
    tw = ModelTag(expr=m.w, format_string="{:.3f}", doc="Tag for w")
    tf = ModelTag(expr=m.f, format_string="{:.3f}", doc="Tag for f")
    tg = ModelTag(expr=m.g, format_string="{:.1f}", doc="Tag for g", display_units="%")
    tz = ModelTag(expr=m.z, format_string="{:.1f}", doc="Tag for z")

    m.w[1, "a"] = 4.0
    m.x[1] = 8.0
    m.y = 4.0
    assert tw.doc == "Tag for w"
    assert str(tw[1, "a"]) == "4.000 kg"
    assert tw[1, "a"](units=False) == "4.000"
    assert tw[1, "a"].display(units=False) == "4.000"
    assert tw.display(units=False, index=(1, "a")) == "4.000"
    assert str(tf) == "2.000 kg/s"
    tf.str_include_units = False
    assert str(tf) == "2.000"
    assert str(tg[1, "a"]) == "50.0%"
    assert str(tg[1, "a"].display(format_string="{:.2f}")) == "50.00%"

    m.z = 7
    assert str(tz) == "7.0"


@pytest.mark.unit
def test_basic_str_float_int():
    ti = ModelTag(expr=1, format_string="{:.3f}", doc="Tag for just an int")
    tf = ModelTag(expr=2.0, format_string="{:.3f}", doc="Tag for just a float")
    ts = ModelTag(expr="hi", format_string="{:.3f}", doc="Tag for just a str")
    tn = ModelTag(expr=None, format_string="{:.3f}", doc="Tag for just None")

    assert ti.value == 1
    assert tf.value == 2.0
    assert ts.value == "hi"
    assert tn.value == None

    ti.set("new")
    tf.set("new")
    ts.set("new")
    tn.set("new")

    assert ti.value == "new"
    assert tf.value == "new"
    assert ts.value == "new"
    assert tn.value == "new"


@pytest.mark.unit
def test_tag_dict_like(model):
    m = model
    tw = ModelTag(expr=m.w, format_string="{:.3f}", doc="Tag for w")
    ty = ModelTag(expr=m.y, format_string="{:.3f}", doc="Tag for y")

    assert len(tw) == 6
    assert len(ty) == 1

    for i in ty.keys():
        assert i == None
        ty[i].set(1)

    check_w_keys = [(1, "a"), (1, "b"), (2, "a"), (2, "b"), (3, "a"), (3, "b")]
    for i, k in enumerate(tw.keys()):
        assert check_w_keys[i] == k
        tw[k].set(i)

    for i, v in enumerate(tw.values()):
        assert pyo.value(v.expression) == i

    for i, (k, v) in enumerate(tw.items()):
        assert pyo.value(v.expression) == i
        assert check_w_keys[i] == k

    assert tw[1, "a"].value == 0
    assert tw[1, "b"].value == 1
    assert tw.value[1, "b"] == 1
    assert ty.value == 1


@pytest.mark.unit
def test_tag_conditional_formatting(model):
    m = model

    tx = ModelTag(
        expr=m.x,
        format_string=lambda x: "{:,.0f}" if x >= 100 else "{:.2f}",
        doc="Tag for x",
        display_units=pyo.units.g,
    )

    tx.set(1 * pyo.units.g)
    assert str(tx[1]) == "1.00 g"
    tx.set(1 * pyo.units.kg)
    assert str(tx[1]) == "1,000 g"


@pytest.mark.unit
def test_tag_errors(model):
    m = model

    tx = ModelTag(
        expr=m.x,
        format_string=lambda x: "{:,.0f}" if x >= 100 else "{:.2f}",
        display_units=pyo.units.g,
    )
    tg = ModelTag(expr=m.g, format_string="{:,.0f}", display_units=None)
    m.x.fix(0)
    assert str(tg[1, "a"]) == "ZeroDivisionError"
    m.x.fix(None)
    assert str(tg[1, "a"]) == "None"


@pytest.mark.unit
def test_tag_display_convert(model):
    m = model
    tw = ModelTag(
        expr=m.w, format_string="{:.3f}", doc="Tag for w", display_units=pyo.units.g
    )
    tx = ModelTag(
        expr=m.x, format_string="{:.3f}", doc="Tag for x", display_units=pyo.units.g
    )
    ty = ModelTag(
        expr=m.y, format_string="{:.1f}", doc="Tag for y", display_units=pyo.units.hr
    )
    tf = ModelTag(
        expr=m.f,
        format_string="{:.1f}",
        doc="Tag for f",
        display_units=pyo.units.g / pyo.units.hr,
    )

    assert str(tf) == "3000000.0 " + str(pyo.units.g / pyo.units.hr)
    m.x[1].value = 4
    assert str(tf) == "2400000.0 " + str(pyo.units.g / pyo.units.hr)
    assert str(tx[1]) == "4000.000 " + str(pyo.units.g)
    m.x[1].value = 3
    assert str(tx[1]) == "3000.000 " + str(pyo.units.g)
    assert str(tw[1, "a"]) == "4000.000 " + str(pyo.units.g)
    assert tw[1, "a"]._index == (1, "a")
    assert tw._cache_display_value[1, "a"] == pytest.approx(4000.0)
    assert tw._cache_validation_value[1, "a"] == 4
    m.w[1, "a"].value = 1
    m.w[2, "a"].value = 2
    m.w[3, "a"].value = 3
    assert str(tw[1, "a"]) == "1000.000 " + str(pyo.units.g)
    assert str(tw[2, "a"]) == "2000.000 " + str(pyo.units.g)
    assert str(tw[3, "a"]) == "3000.000 " + str(pyo.units.g)
    assert tw._cache_display_value[1, "a"] == pytest.approx(1000.0)
    assert tw._cache_validation_value[1, "a"] == 1


@pytest.mark.unit
def test_tag_input(model):
    m = model
    tw = ModelTag(expr=m.w, format_string="{:.3f}", doc="Tag for w")
    ty = ModelTag(expr=m.y, format_string="{:.2f}", doc="Tag for y")

    tw.fix(1)
    for v in m.w.values():
        assert v.value == 1
        assert v.fixed
    tw.unfix()
    tw[:, "a"].fix(2)
    for k, v in m.w.items():
        if k[1] == "a":
            assert v.value == 2
            assert v.fixed
        else:
            assert v.value == 1
            assert not v.fixed
    tw[:, "b"].fix(4)
    for k, v in m.w.items():
        if k[1] == "b":
            assert v.value == 4

    tw.var[1, "a"] = 3
    assert issubclass(tw.var[1, "a"].ctype, pyo.Var)
    assert pyo.value(tw[1, "a"].var) == 3

    ty.fix(2)
    assert m.y.fixed
    assert pyo.value(m.y) == 2
    ty.set(4)
    assert pyo.value(m.y) == 4


@pytest.mark.unit
def test_tag_input_convert(model):
    m = model
    tw = ModelTag(
        expr=m.w, format_string="{:.3f}", doc="Tag for w", display_units=pyo.units.g
    )
    ty = ModelTag(expr=m.y, format_string="{:.2f}", doc="Tag for y")

    ty.set(1 * pyo.units.hr)
    assert str(ty) == "3600.00 s"

    ty.fix(2 * pyo.units.hr)
    assert str(ty) == "7200.00 s"
    assert m.y.fixed

    tw[:, "a"].set(3000, in_display_units=True)
    tw[:, "b"].set(2 * pyo.units.kg, in_display_units=True)
    assert pyo.value(m.w[1, "a"]) == pytest.approx(3.0)
    assert pyo.value(m.w[2, "a"]) == pytest.approx(3.0)
    assert pyo.value(m.w[3, "a"]) == pytest.approx(3.0)
    assert pyo.value(m.w[1, "b"]) == pytest.approx(2.0)
    assert pyo.value(m.w[2, "b"]) == pytest.approx(2.0)
    assert pyo.value(m.w[3, "b"]) == pytest.approx(2.0)

    tw.set_in_display_units = True
    tw[:, "a"].set(1000)
    tw[:, "b"].set(4 * pyo.units.kg)
    assert pyo.value(m.w[1, "a"]) == pytest.approx(1.0)
    assert pyo.value(m.w[2, "a"]) == pytest.approx(1.0)
    assert pyo.value(m.w[3, "a"]) == pytest.approx(1.0)
    assert pyo.value(m.w[1, "b"]) == pytest.approx(4.0)
    assert pyo.value(m.w[2, "b"]) == pytest.approx(4.0)
    assert pyo.value(m.w[3, "b"]) == pytest.approx(4.0)

    tw.set_in_display_units = False
    tw[:, "a"].set(2)
    assert pyo.value(m.w[1, "a"]) == pytest.approx(2.0)
    assert pyo.value(m.w[2, "a"]) == pytest.approx(2.0)
    assert pyo.value(m.w[3, "a"]) == pytest.approx(2.0)


@pytest.mark.unit
def test_tag_ref(model):
    m = model
    m.rw = pyo.Reference(m.w[:, "a"])
    m.ry = pyo.Reference(m.y)

    rw = ModelTag(expr=m.rw, format_string="{:.3f}", doc="Tag for rw")
    ry = ModelTag(expr=m.ry, format_string="{:.2f}", doc="Tag for ry")

    m.w[1, "a"] = 3
    assert str(rw[1]) == "3.000 kg"
    assert rw[1].is_var

    ry[:].set(2)
    assert str(ry[None]) == "2.00 s"


@pytest.mark.unit
def test_tag_data_object(model):
    m = model

    rw = ModelTag(expr=m.w[1, "a"], format_string="{:.3f}", doc="Tag for rw")

    m.w[1, "a"] = 3
    assert str(rw) == "3.000 kg"
    assert rw.is_var

    rw.set(5)
    assert str(rw) == "5.000 kg"

    rw = ModelTag(expr=m.w[1, "a"], format_string="{:.3f}", display_units=pyo.units.g)

    assert str(rw) == "5000.000 g"
    rw.set_in_display_units = True
    rw.set(5)
    assert str(rw) == "5.000 g"


@pytest.mark.unit
def test_tag_group(model):
    m = model
    g = ModelTagGroup()
    m.rw = pyo.Reference(m.w[:, "a"])
    g.add(
        "w",
        expr=m.rw,
        format_string="{:.3f}",
        doc="make sure this works",
        display_units=pyo.units.g,
    )
    g.add("x", expr=m.x, format_string="{:.3f}")
    g.add("y", expr=m.y, format_string="{:.3f}")
    g.add("z", expr=m.z, format_string="{:.3f}")
    g.add("e", expr=m.e, format_string="{:.3f}")
    g.add("f", ModelTag(expr=m.f, format_string="{:.3f}"))
    g.add("g", expr=m.g, format_string="{:.1f}", display_units="%")

    assert g["w"].doc == "make sure this works"
    g.str_include_units = True

    assert g.str_include_units
    assert g["w"].str_include_units
    assert not g.set_in_display_units
    assert not g["w"].set_in_display_units
    g["w"].fix(2)
    g["x"].fix(1)
    g["y"].fix(3)
    assert str(g["w"][1]) == "2000.000 g"
    assert str(g["x"][1]) == "1.000 kg"
    assert str(g["y"]) == "3.000 s"

    g.str_include_units = False
    assert str(g["w"][1]) == "2000.000"
    assert str(g["x"][1]) == "1.000"
    assert str(g["y"]) == "3.000"

    g.set_in_display_units = True
    g["w"][:].set(1000)
    assert str(g["w"][1]) == "1000.000"
    assert str(g["w"][2]) == "1000.000"
    assert str(g["w"][2]) == "1000.000"

    g.set_in_display_units = True
    g["w"][:].set(2 * pyo.units.kg)
    assert str(g["w"][1]) == "2000.000"
    assert str(g["w"][2]) == "2000.000"
    assert str(g["w"][2]) == "2000.000"

    g["w"][:].fix(3 * pyo.units.kg)
    assert str(g["w"][1]) == "3000.000"
    assert str(g["w"][2]) == "3000.000"
    assert str(g["w"][2]) == "3000.000"

    g["w"][:].fix(4000)
    assert str(g["w"][1]) == "4000.000"
    assert str(g["w"][2]) == "4000.000"
    assert str(g["w"][2]) == "4000.000"


@pytest.mark.unit
def test_tabulate_runs(model):
    m = model
    g = ModelTagGroup()
    g.add("w", expr=m.w, format_string="{:.3f}", display_units=pyo.units.g)
    g.add("x", expr=m.x, format_string="{:.3f}")
    g.add("y", expr=m.y, format_string="{:.3f}")
    g.add("z", expr=m.z, format_string="{:.3f}")
    g.add("e", expr=m.e, format_string="{:.3f}")
    g.add("f", ModelTag(expr=m.f, format_string="{:.3f}"))
    g.add("g", expr=m.g, format_string="{:.1f}", display_units="%")

    columns = (
        ["w", (1, "a")],
        ["w", (2, "a")],
        ["y", None],
    )

    head = g.table_heading(tags=columns, units=True)

    assert head[0] == "w[(1, 'a')] (g)"
    assert head[1] == "w[(2, 'a')] (g)"
    assert head[2] == "y (s)"

    g["w"][1, "a"].set(1 * pyo.units.g)
    g["w"][2, "a"].set(2 * pyo.units.g)
    g["y"].set(1 * pyo.units.s)

    row = g.table_row(tags=columns, units=False)
    assert row[0] == "1.000"
    assert row[1] == "2.000"
    assert row[2] == "1.000"

    row = g.table_row(tags=columns, units=True)
    assert row[0] == "1.000 g"
    assert row[1] == "2.000 g"
    assert row[2] == "1.000 s"

    row = g.table_row(tags=columns, numeric=True)
    assert row[0] == pytest.approx(1.000)
    assert row[1] == pytest.approx(2.000)
    assert row[2] == pytest.approx(1.000)

    g["w"][1, "a"].set(1 * pyo.units.g)
    g["w"][1, "b"].set(2 * pyo.units.g)
    g["w"][2, "a"].set(3 * pyo.units.g)
    g["w"][2, "b"].set(4 * pyo.units.g)
    g["w"][3, "a"].set(5 * pyo.units.g)
    g["w"][3, "b"].set(6 * pyo.units.g)

    head = g.table_heading(tags=("w", "y"), units=True)
    assert head[0] == "w[(1, 'a')] (g)"
    assert head[1] == "w[(1, 'b')] (g)"
    assert head[2] == "w[(2, 'a')] (g)"
    assert head[3] == "w[(2, 'b')] (g)"
    assert head[4] == "w[(3, 'a')] (g)"
    assert head[5] == "w[(3, 'b')] (g)"
    assert head[6] == "y (s)"

    row = g.table_row(tags=("w", "y"), units=True)
    assert row[0] == "1.000 g"
    assert row[1] == "2.000 g"
    assert row[2] == "3.000 g"
    assert row[3] == "4.000 g"
    assert row[4] == "5.000 g"
    assert row[5] == "6.000 g"
    assert row[6] == "1.000 s"


@pytest.mark.unit
def test_doc_example_and_bound(model):
    m = model
    g = ModelTagGroup()

    g["w"] = ModelTag(expr=m.w, format_string="{:.3f}")
    g["x"] = ModelTag(expr=m.x, format_string="{:.3f}", display_units=pyo.units.g)
    g["y"] = ModelTag(expr=m.y, format_string="{:.3f}")
    g["z"] = ModelTag(expr=m.z, format_string="{:.3f}")
    g["e"] = ModelTag(expr=m.e, format_string="{:.3f}")
    g["f"] = ModelTag(expr=m.f, format_string="{:.3f}")
    g["g"] = ModelTag(expr=m.g, format_string="{:.3f}")

    g.set_in_display_units = True
    g.str_include_units = False

    g["x"].set(2)
    g["x"].setlb(1)
    g["x"].setub(3)

    assert str(g["x"][1]) == "2.000"
    assert abs(g["x"][1].expression.lb - 0.001) < 1e-5  # x is in kg
    assert abs(g["x"][1].expression.ub - 0.003) < 1e-5  # x is in kg


@pytest.mark.unit
def test_tag_slice(model):
    m = model

    tw = ModelTag(
        expr=m.w[1, :],
        format_string=lambda x: "{:,.0f}" if x >= 100 else "{:.2f}",
        doc="Tag for x",
        display_units=pyo.units.g,
    )

    tw.set(1 * pyo.units.g)
    assert str(tw["a"]) == "1.00 g"
    tw.set(1 * pyo.units.kg)
    assert str(tw["b"]) == "1,000 g"
