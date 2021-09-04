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
from idaes.core.util import ModelTag, ModelTagGroup

@pytest.fixture()
def model():
    m = pyo.ConcreteModel()
    m.w = pyo.Var([1,2,3], ["a", "b"], initialize=4, units=pyo.units.kg)
    m.x = pyo.Var([1,2,3], initialize=5, units=pyo.units.kg)
    m.y = pyo.Var(initialize=6, units=pyo.units.s)
    m.z = pyo.Var(initialize=7)
    m.e = pyo.Expression(expr=m.w[1, "a"]/m.x[1])
    m.f = pyo.Expression(expr=m.x[1]/m.y)
    @m.Expression([1, 2, 3], ["a", "b"])
    def g(b, i, j):
        return m.w[i, j]/m.x[i]*100
    return m

@pytest.mark.unit
def test_tag_display(model):
    m = model
    tw = ModelTag(expr=m.w, format="{:.3f}", doc="Tag for w")
    tf = ModelTag(expr=m.f, format="{:.3f}", doc="Tag for f")
    tg = ModelTag(expr=m.g, format="{:.1f}", doc="Tag for g", units="%")
    tz = ModelTag(expr=m.z, format="{:.1f}", doc="Tag for z")

    m.w[1, "a"] = 4.0
    m.x[1] = 8.0
    m.y = 4.0
    assert tw.doc == "Tag for w"
    assert str(tw[1, "a"]) == "4.000 kg"
    assert tw[1, "a"](units=False) == "4.000"
    assert tw[1, "a"].display(units=False) == "4.000"
    assert tw.display(units=False, index=(1, "a")) == "4.000"
    assert str(tf) == "2.000 kg/s"
    tf.str_include_units(False)
    assert str(tf) == "2.000"
    assert str(tg[1, "a"]) == "50.0%"
    assert str(tg[1, "a"].display(format="{:.2f}")) == "50.00%"

    m.z = 7
    assert str(tz) == "7.0"

@pytest.mark.unit
def test_tag_input(model):
    m = model
    tw = ModelTag(expr=m.w, format="{:.3f}", doc="Tag for w")
    ty = ModelTag(expr=m.y, format="{:.2f}", doc="Tag for y")

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
def test_tag_ref(model):
    m = model
    m.rw = pyo.Reference(m.w[:,"a"])
    m.ry = pyo.Reference(m.y)

    rw = ModelTag(expr=m.rw, format="{:.3f}", doc="Tag for rw")
    ry = ModelTag(expr=m.ry, format="{:.2f}", doc="Tag for ry")

    m.w[1, "a"] = 3
    assert str(rw[1]) == "3.000 kg"
    assert rw[1].is_var

    ry[:].set(2)
    assert str(ry[None]) == "2.00 s"

@pytest.mark.unit
def test_tag_group(model):
    m = model
    g = ModelTagGroup()
    m.rw = pyo.Reference(m.w[:,"a"])
    g.add("w", expr=m.rw, format="{:.3f}", doc="make sure this works")
    g.add("x", expr=m.x, format="{:.3f}")
    g.add("y", expr=m.y, format="{:.3f}")
    g.add("z", expr=m.z, format="{:.3f}")
    g.add("e", expr=m.e, format="{:.3f}")
    g.add("f", ModelTag(expr=m.f, format="{:.3f}"))
    g.add("g", expr=m.g, format="{:.1f}", units="%")


    assert g["w"].doc == "make sure this works"
    g.str_index(1)
    g.str_include_units(True)
    g['w'].fix(2)

    g['x'].fix(1)
    g['y'].fix(3)
    assert str(g["w"]) == "2.000 kg"
    assert str(g["x"]) == "1.000 kg"
    assert str(g["y"]) == "3.000 s"

    g.str_include_units(False)
    assert str(g["w"]) == "2.000"
    assert str(g["x"]) == "1.000"
    assert str(g["y"]) == "3.000"
