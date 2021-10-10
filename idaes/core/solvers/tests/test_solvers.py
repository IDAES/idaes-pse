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
from pyomo.environ import SolverFactory
import pyomo.environ as pyo
import pytest
import idaes.core.plugins
from idaes.core.solvers.interogate import (
    lp, milp, nlp, minlp, ipopt_has_linear_solver
)

@pytest.mark.unit
def test_couenne_available():
    if not SolverFactory('couenne').available():
        raise Exception("Could not find couenne.")

@pytest.mark.unit
def test_couenne_available():
    if not SolverFactory('bonmin').available():
        raise Exception("Could not find bonmin.")

@pytest.mark.unit
def test_sipopt_available():
    if not SolverFactory('ipopt_sens').available():
        raise Exception("Could not find ipopt_sens.")

@pytest.mark.unit
def test_ipopt_idaes_available():
    """
    Tries to set-up the IPOPT with the IDAES SolverFactory wrapper
    """
    if not pyo.SolverFactory('ipopt').available():
        raise Exception(
            "Could not find IPOPT. Users are strongly encouraged to have a "
            "version of IPOPT available, as it is the default solver assumed "
            "by many IDAES examples and tests. See the IDAES install "
            "documentation for instructions on how to get IPOPT.")

@pytest.mark.unit
def test_cbc_available():
    if not SolverFactory('cbc').available():
        raise Exception("Could not find cbc.")

@pytest.mark.unit
def test_sipopt_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = nlp()
    solver = SolverFactory('ipopt_sens')
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)

@pytest.mark.unit
def test_ipopt_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = nlp()
    solver = SolverFactory('ipopt')
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)

@pytest.mark.unit
def test_ipopt_has_ma27():
    if not ipopt_has_linear_solver("ma27"):
        raise Exception(
            "The ma27 linear solver is not available to Ipopt. Models may solve"
            " more reliably with HSL linear solvers see https://www.hsl.rl.ac.uk/,"
            " or use solvers distributed by the IDAES project. See IDAES install"
            " guide.")

@pytest.mark.unit
def test_bonmin_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = minlp()
    solver = SolverFactory('bonmin')
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)

@pytest.mark.unit
def test_couenne_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = minlp()
    solver = SolverFactory('couenne')
    solver.solve(m)
    assert pytest.approx(x) == pyo.value(m.x)

@pytest.mark.unit
def test_cbc_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    m, x = milp()
    solver = SolverFactory('cbc')
    solver.solve(m)
    m.display()
    assert pytest.approx(x) == pyo.value(m.x)
