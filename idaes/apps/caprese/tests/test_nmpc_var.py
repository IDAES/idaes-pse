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

def test_NmpcVar():
    m = pyo.ConcreteModel()

    with pytest.raises(NotImplementedError,
            match=r".*component must be indexed.*"):
        m.v = NmpcVar()

    m.s1 = pyo.Set(initialize=[0,1,2,3])
    m.s2 = pyo.Set(initialize=[2,4,6,8])

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
            setpoint=1.,
            weight=2.,
            variance=3.,
            nominal=4.,
            )
    assert m.v1.setpoint is 1.
    assert m.v1.weight is 2.
    assert m.v1.variance is 3.
    assert m.v1.nominal is 4.

    m.v2 = NmpcVar(m.s1, m.s2)
    for i, j in m.s1*m.s2:
        assert (i,j in m.v2)

def test_custom_vars():
    m = pyo.ConcreteModel()
    m.s = pyo.Set(initialize=[0,1,2])

    m.diff = DiffVar(m.s)
    assert m.diff.ctype == DiffVar
    assert m.diff._attr == 'differential'

    m.alg = AlgVar(m.s)
    assert m.alg.ctype == AlgVar
    assert m.alg._attr == 'algebraic'

    m.inp = InputVar(m.s)
    assert m.inp.ctype == InputVar
    assert m.inp._attr == 'input'

    m.deriv = DerivVar(m.s)
    assert m.deriv.ctype == DerivVar
    assert m.deriv._attr == 'derivative'

    m.meas = MeasuredVar(m.s)
    assert m.meas.ctype == MeasuredVar
    assert m.meas._attr == 'measurement'

