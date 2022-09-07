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
""" Tests for subclasses of Var
"""

import pyomo.environ as pyo
from idaes.apps.caprese.nmpc_var import (
    NmpcVar,
    _NmpcVector,
    DiffVar,
    DerivVar,
    AlgVar,
    InputVar,
    FixedVar,
    MeasuredVar,
)
import pytest


@pytest.mark.unit
def test_NmpcVar():
    m = pyo.ConcreteModel()

    with pytest.raises(NotImplementedError, match=r".*component must be indexed.*"):
        m.v = NmpcVar()

    m.s1 = pyo.Set(initialize=[0, 1, 2, 3])
    m.s2 = pyo.Set(initialize=[2, 4, 6, 8])

    m.v0 = NmpcVar(m.s1)
    assert m.v0.setpoint is None
    assert m.v0.weight is None
    assert m.v0.nominal is None
    assert m.v0.variance is None

    assert m.v0.ctype is NmpcVar
    for v in m.component_objects(NmpcVar):
        assert v is m.v0

    m.v1 = NmpcVar(
        m.s1,
        setpoint=1.0,
        weight=2.0,
        variance=3.0,
        nominal=4.0,
    )
    assert m.v1.setpoint == 1.0
    assert m.v1.weight == 2.0
    assert m.v1.variance == 3.0
    assert m.v1.nominal == 4.0

    m.v2 = NmpcVar(m.s1, m.s2)
    for i, j in m.s1 * m.s2:
        assert (i, j) in m.v2


@pytest.mark.unit
def test_custom_vars():
    m = pyo.ConcreteModel()
    m.s = pyo.Set(initialize=[0, 1, 2])

    m.diff = DiffVar(m.s)
    assert m.diff.ctype == DiffVar
    assert m.diff._attr == "differential"

    m.alg = AlgVar(m.s)
    assert m.alg.ctype == AlgVar
    assert m.alg._attr == "algebraic"

    m.inp = InputVar(m.s)
    assert m.inp.ctype == InputVar
    assert m.inp._attr == "input"

    m.deriv = DerivVar(m.s)
    assert m.deriv.ctype == DerivVar
    assert m.deriv._attr == "derivative"

    m.meas = MeasuredVar(m.s)
    assert m.meas.ctype == MeasuredVar
    assert m.meas._attr == "measurement"


@pytest.mark.unit
def test_NmpcVector():
    m = pyo.ConcreteModel()
    m.coords = pyo.Set(initialize=[0, 1, 2, 3])
    m.time = pyo.Set(initialize=[0.0, 0.5, 1.0, 1.5, 2.0])

    @m.Block(m.coords)
    def b(b, i):
        b.var = NmpcVar(m.time)

    m.vector = pyo.Reference(m.b[:].var[:], ctype=_NmpcVector)

    assert type(m.vector) is _NmpcVector

    # Test that `vector` is a proper reference
    for i, t in m.coords * m.time:
        assert m.vector[i, t] is m.b[i].var[t]

    # Test that we can generate the underlying NmpcVars
    for v, i in zip(m.vector._generate_referenced_vars(), m.coords):
        assert v is m.b[i].var

    # `set_setpoint`
    m.vector.set_setpoint(3.14)
    for i in m.coords:
        assert m.b[i].var.setpoint == 3.14

    setpoint = tuple(i / 10.0 for i in m.coords)
    m.vector.set_setpoint(setpoint)
    for i, sp in zip(m.coords, setpoint):
        assert m.b[i].var.setpoint == sp

    # `get_setpoint`
    sp_get = tuple(m.vector.get_setpoint())
    assert setpoint == sp_get

    # `set_values`
    val = -1.0
    m.vector.values = -1.0
    for i, t in m.coords * m.time:
        assert m.b[i].var[t].value == val

    newvals = tuple(i * 1.1 for i in m.coords)
    m.vector.values = newvals
    for i, val in zip(m.coords, newvals):
        for t in m.time:
            assert m.b[i].var[t].value == val

    # `get_values`
    val_lil = m.vector.values
    for var_values, target_val in zip(val_lil, newvals):
        for var_val in var_values:
            assert var_val == target_val
