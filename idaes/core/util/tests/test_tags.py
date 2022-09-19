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
from idaes.core.util.misc import TagReference
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


"""
DEPRECATED this module provides tests for functions that have been deprecated,
and will be removed.

This module contains miscalaneous utility functions for use in IDAES models.
"""

# simple svg from drawn in incscape, with text ids:
#   TAGME_4.x, TAGME_4.y, and TAGME_4.f
svg_test_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   width="210mm"
   height="297mm"
   viewBox="0 0 210 297"
   version="1.1"
   id="svg2045"
   inkscape:version="0.92.4 (5da689c313, 2019-01-14)"
   sodipodi:docname="test.svg">
  <defs
     id="defs2039" />
  <sodipodi:namedview
     id="base"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:pageopacity="0.0"
     inkscape:pageshadow="2"
     inkscape:zoom="0.98994949"
     inkscape:cx="174.80933"
     inkscape:cy="316.53691"
     inkscape:document-units="mm"
     inkscape:current-layer="layer1"
     showgrid="false"
     inkscape:window-width="3840"
     inkscape:window-height="2035"
     inkscape:window-x="-13"
     inkscape:window-y="-13"
     inkscape:window-maximized="1" />
  <metadata
     id="metadata2042">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <g
     inkscape:label="Layer 1"
     inkscape:groupmode="layer"
     id="layer1">
    <ellipse
       style="stroke-width:0.26499999;stroke-miterlimit:4;stroke-dasharray:none"
       id="path2590"
       cx="113.18864"
       cy="130.22382"
       rx="30.869631"
       ry="31.003265" />
    <path
       style="fill:none;stroke:#000000;stroke-width:0.26458332px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       d="m 45.70309,75.433559 37.685,33.943231"
       id="path2592"
       inkscape:connector-curvature="0" />
    <path
       style="fill:none;stroke:#000000;stroke-width:0.26458332px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
       d="m 136.04019,155.61442 49.17759,45.16855"
       id="path2594"
       inkscape:connector-curvature="0" />
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="63.610146"
       y="68.751816"
       id="text2598"><tspan
         sodipodi:role="line"
         id="tspan2596"
         x="63.610146"
         y="68.751816"
         style="stroke-width:0.26458332">x:</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="66.282845"
       y="86.391609"
       id="text2602"><tspan
         sodipodi:role="line"
         x="66.282845"
         y="86.391609"
         style="stroke-width:0.26458332"
         id="tspan2604">y:</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="166.77618"
       y="171.11606"
       id="text2610"><tspan
         sodipodi:role="line"
         id="tspan2608"
         x="166.77618"
         y="171.11606"
         style="stroke-width:0.26458332">f:</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="77.723328"
       y="68.57209"
       id="TAGME_4.x"
       inkscape:label="#text2598-9"><tspan
         sodipodi:role="line"
         id="tspan2596-1"
         x="77.723328"
         y="68.57209"
         style="stroke-width:0.26458332">1</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="80.396027"
       y="86.211876"
       id="TAGME_4.y"
       inkscape:label="#text2602-0"><tspan
         sodipodi:role="line"
         x="80.396027"
         y="86.211876"
         style="stroke-width:0.26458332"
         id="tspan2604-1">1</tspan></text>
    <text
       xml:space="preserve"
       style="font-style:normal;font-weight:normal;font-size:10.58333302px;line-height:1.25;font-family:sans-serif;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;stroke-width:0.26458332"
       x="176.65611"
       y="171.99467"
       id="TAGME_4.f"
       inkscape:label="#text2610-0"><tspan
         sodipodi:role="line"
         id="tspan2608-1"
         x="176.65611"
         y="171.99467"
         style="stroke-width:0.26458332"></tspan></text>
    <flowRoot
       xml:space="preserve"
       id="flowRoot2642"
       style="fill:black;fill-opacity:1;stroke:none;font-family:sans-serif;font-style:normal;font-weight:normal;font-size:40px;line-height:1.25;letter-spacing:0px;word-spacing:0px"><flowRegion
         id="flowRegion2644"><rect
           id="rect2646"
           width="533.36053"
           height="577.80725"
           x="714.17786"
           y="232.57529" /></flowRegion><flowPara
         id="flowPara2648"></flowPara></flowRoot>  </g>
</svg>"""


@pytest.mark.unit
def test_tag_reference():
    m = ConcreteModel()
    m.x = Var([0, 1], initialize={0: 2.22, 1: 3.33})
    m.y = Var([0, 1], initialize={0: 4.44, 1: 5.55})
    m.f = Var([0, 1], initialize={0: 6.66, 1: 7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag

    xml_str = svg_tag(m.tag, svg_test_str, idx=0)
    # lazy testing
    assert "2.2200e" in xml_str
    assert "4.4400e" in xml_str
    assert "6.6600e" in xml_str

    xml_str = svg_tag(m.tag, svg_test_str, idx=1)
    # lazy testing
    assert "3.3300e" in xml_str
    assert "5.5500e" in xml_str
    assert "7.7700e" in xml_str

    xml_str = svg_tag(m.tag, svg_test_str, show_tags=True)
    # lazy testing
    assert "TAGME@4.x" in xml_str
    assert "TAGME@4.y" in xml_str
    assert "TAGME@4.f" in xml_str

    tag_data_like = {}
    tag_data_like["TAGME@4.x"] = 1.1212
    tag_data_like["TAGME@4.y"] = 2.1212
    tag_data_like["TAGME@4.f"] = "3.1212 Hello"
    xml_str = svg_tag(tag_data_like, svg_test_str, idx=None)
    assert "1.1212" in xml_str
    assert "2.1212" in xml_str
    assert "3.1212 Hello" in xml_str


@pytest.mark.unit
def test_tag_reference_default_format():
    m = ConcreteModel()
    m.x = Var([0, 1], initialize={0: 2.22, 1: 3.33})
    m.y = Var([0, 1], initialize={0: 4.44, 1: 5.55})
    m.f = Var([0, 1], initialize={0: 6.66, 1: 7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag

    xml_str = svg_tag(m.tag, svg_test_str, idx=0, tag_format_default="{:.2f}X")
    # lazy testing
    assert "2.22X" in xml_str
    assert "4.44X" in xml_str
    assert "6.66X" in xml_str


@pytest.mark.unit
def test_tag_reference_tag_format():
    m = ConcreteModel()
    m.x = Var([0, 1], initialize={0: 2.22, 1: 3.33})
    m.y = Var([0, 1], initialize={0: 4.44, 1: 5.55})
    m.f = Var([0, 1], initialize={0: 6.66, 1: 7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag

    xml_str = svg_tag(m.tag, svg_test_str, idx=0, tag_format={"TAGME@4.y": "{:.5e}"})
    # lazy testing
    assert "2.2200e" in xml_str
    assert "4.44000e" in xml_str
    assert "6.6600e" in xml_str


@pytest.mark.unit
def test_tag_reference_tag_format_conditional():
    m = ConcreteModel()
    m.x = Var([0, 1], initialize={0: 10000, 1: 3.33})
    m.y = Var([0, 1], initialize={0: 1.1, 1: 5.55})
    m.f = Var([0, 1], initialize={0: 4.5, 1: 7.77})
    test_tag = {}
    test_tag["TAGME@4.x"] = TagReference(m.x[:], description="x tag")
    test_tag["TAGME@4.y"] = TagReference(m.y[:], description="y tag")
    test_tag["TAGME@4.f"] = TagReference(m.f[:], description="z tag")
    m.tag = test_tag

    xml_str = svg_tag(
        m.tag,
        svg_test_str,
        idx=0,
        tag_format_default=lambda x: "{:,.0f} kPa" if x >= 100 else "{:.2f} kPa",
    )
    # lazy testing
    assert "10,000 kPa" in xml_str
    assert "1.10 kPa" in xml_str
    assert "4.50 kPa" in xml_str


if __name__ == "__main__":
    # Check deprication warnings
    test_tag_reference_tag_format_conditional()
