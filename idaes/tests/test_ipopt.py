##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
from pyomo.environ import SolverFactory
import pyomo.environ as pyo
import pytest
import idaes.core.plugins


@pytest.mark.unit
def test_ipopt_available():
    """
    Tries to set-up the IPOPT and returns exception if not available
    """
    if not SolverFactory('ipopt').available():
        raise Exception(
            "Could not find IPOPT. Users are strongly encouraged to have a "
            "version of IPOPT available, as it is the default solver assumed "
            "by many IDAES examples and tests. See the IDAES install "
            "documentation for instructions on how to get IPOPT.")

@pytest.mark.unit
def test_ipopt_idaes_available():
    """
    Tries to set-up the IPOPT and returns exception if not available
    """
    if not SolverFactory('ipopt-idaes').available():
        raise Exception(
            "Could not find IPOPT. Users are strongly encouraged to have a "
            "version of IPOPT available, as it is the default solver assumed "
            "by many IDAES examples and tests. See the IDAES install "
            "documentation for instructions on how to get IPOPT.")

@pytest.mark.skipif(not SolverFactory('ipopt-idaes').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_ipopt_idaes_config():
    """
    Test that the default solver options are set
    """
    solver = SolverFactory('ipopt-idaes')
    for k, v in idaes.cfg["ipopt-idaes"]["options"].items():
        solver.options[k] = v

@pytest.mark.skipif(not SolverFactory('ipopt-idaes').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_ipopt_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    solver = SolverFactory('ipopt-idaes')
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=-0.1)
    m.y = pyo.Var(initialize=1)
    m.c1 = pyo.Constraint(expr=1==m.x**3)
    m.c2 = pyo.Constraint(expr=1==m.y)
    solver.solve(m)
    assert pytest.approx(1) == pyo.value(m.x)
