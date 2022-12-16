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
Basic test to make sure the pynumero library is available and working
"""

import pyomo.environ as pyo
import pytest

try:
    from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
except ImportError:
    PyomoNLP = None


@pytest.fixture
def model():
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2, 3], initialize=4.0)
    m.c = pyo.Constraint(expr=m.x[3] ** 2 + m.x[1] == 25)
    m.d = pyo.Constraint(expr=m.x[2] ** 2 + m.x[1] <= 18.0)
    m.o = pyo.Objective(expr=m.x[1] ** 4 - 3 * m.x[1] * m.x[2] ** 3 + m.x[3] ** 2 - 8.0)
    m.x[1].setlb(0.0)
    m.x[2].setlb(0.0)
    return m


@pytest.mark.unit
def test_import():
    assert PyomoNLP is not None


@pytest.mark.unit
def test_have_pynumero(model):
    assert PyomoNLP is not None
    nlp = PyomoNLP(model)
    f = nlp.evaluate_objective()
    assert f == pytest.approx(-504.0)

    jac = nlp.evaluate_jacobian().toarray()
    assert jac[0][0] == pytest.approx(0)
    assert jac[0][1] == pytest.approx(8)
    assert jac[0][2] == pytest.approx(1)
    assert jac[1][0] == pytest.approx(8)
    assert jac[1][1] == pytest.approx(0)
    assert jac[1][2] == pytest.approx(1)
