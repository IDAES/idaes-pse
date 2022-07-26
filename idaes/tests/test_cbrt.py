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
from idaes.core.util.functions import functions_lib, functions_available
import pyomo.environ as pyo
import pytest


@pytest.mark.skipif(not functions_available(), reason="functions.so not available")
@pytest.mark.unit
def test_cbrt_values():
    m = pyo.ConcreteModel()
    flib = functions_lib()
    m.cbrt = pyo.ExternalFunction(library=flib, function="cbrt")
    assert abs(pyo.value(m.cbrt(-27.0)) + 3.0) < 0.00001
    assert abs(pyo.value(m.cbrt(0.0))) < 0.00001
    assert abs(pyo.value(m.cbrt(27.0)) - 3.0) < 0.00001


@pytest.mark.skipif(not functions_available(), reason="functions.so not available")
@pytest.mark.unit
def test_cbrt_derivs():
    m = pyo.ConcreteModel()
    flib = functions_lib()
    m.cbrt = pyo.ExternalFunction(library=flib, function="cbrt")
    h = 1e-6
    tol = 1e-5
    for v in [-27.0, -0.5, 0.5, 27]:
        f1, g1, h1 = m.cbrt.evaluate_fgh(args=(v,))
        f2, g2, h2 = m.cbrt.evaluate_fgh(args=(v + h,))
        gfd = (f2 - f1) / h
        assert abs(g1[0] - gfd) < tol


@pytest.mark.skipif(not functions_available(), reason="functions.so not available")
@pytest.mark.unit
def test_cbrt_hes():
    m = pyo.ConcreteModel()
    flib = functions_lib()
    m.cbrt = pyo.ExternalFunction(library=flib, function="cbrt")
    h = 1e-6
    tol = 1e-5
    for v in [-27.0, -0.5, 0.5, 27]:
        f1, g1, h1 = m.cbrt.evaluate_fgh(args=(v,))
        f2, g2, h2 = m.cbrt.evaluate_fgh(args=(v + h,))
        hfd = (g2[0] - g1[0]) / h
        assert abs(h1[0] - hfd) < tol
