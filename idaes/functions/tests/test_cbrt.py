from idaes.functions import functions_lib, functions_available
import pyomo.environ as pyo
import pytest

@pytest.mark.skipif(not functions_available(), reason="functions.so not available")
def test_cbrt_values():
    m = pyo.ConcreteModel()
    flib = functions_lib()
    m.cbrt = pyo.ExternalFunction(library=flib, function="cbrt")
    assert(abs(pyo.value(m.cbrt(-27.0)) + 3.0) < 0.00001)
    assert(abs(pyo.value(m.cbrt(0.0))) < 0.00001)
    assert(abs(pyo.value(m.cbrt(27.0)) - 3.0) < 0.00001)

def test_cbrt_derivs():
    m = pyo.ConcreteModel()
    flib = functions_lib()
    m.cbrt = pyo.ExternalFunction(library=flib, function="cbrt")
    h = 1e-6
    tol = 1e-5
    for v in [-27.0, -0.5, 0.5, 27]:
        f1, g1, h1 = m.cbrt.evaluate_fgh(args=(v,))
        f2, g2, h2 = m.cbrt.evaluate_fgh(args=(v + h,))
        gfd = (f2 - f1)/h
        assert(abs(g1[0] - gfd) < tol)

def test_cbrt_hes():
    m = pyo.ConcreteModel()
    flib = functions_lib()
    m.cbrt = pyo.ExternalFunction(library=flib, function="cbrt")
    h = 1e-6
    tol = 1e-5
    for v in [-27.0, -0.5, 0.5, 27]:
        f1, g1, h1 = m.cbrt.evaluate_fgh(args=(v,))
        f2, g2, h2 = m.cbrt.evaluate_fgh(args=(v + h,))
        hfd = (g2[0] - g1[0])/h
        assert(abs(h1[0] - hfd) < tol)
