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
import pyomo.environ as pyo
import pytest
import idaes
import idaes.core.solvers as isolve

isolve.use_idaes_solver_configuration_deafults()


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

@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_ipopt_idaes_config():
    """
    Test that the default solver options are set
    """
    orig = idaes.cfg["ipopt"]["options"]["nlp_scaling_method"]
    idaes.cfg["ipopt"]["options"]["nlp_scaling_method"] = "gradient-based"
    solver = pyo.SolverFactory('ipopt')
    assert solver.options["nlp_scaling_method"] == "gradient-based"
    idaes.cfg["ipopt"]["options"]["nlp_scaling_method"] = orig
    solver = pyo.SolverFactory('ipopt', options={"tol":1})
    assert solver.options["tol"] == 1
    isolve.use_idaes_solver_configuration_deafults(False)
    solver = pyo.SolverFactory('ipopt')
    assert "nlp_scaling_method" not in solver.options
    isolve.use_idaes_solver_configuration_deafults()


@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_ipopt_idaes_solve():
    """
    Make sure there is no issue with the solver class or default settings that
    break the solver object.  Passing a bad solver option will result in failure
    """
    solver = pyo.SolverFactory('ipopt')
    m = pyo.ConcreteModel()
    m.x = pyo.Var(initialize=-0.1)
    m.y = pyo.Var(initialize=1)
    m.c1 = pyo.Constraint(expr=1==m.x**3)
    m.c2 = pyo.Constraint(expr=1==m.y)
    solver.solve(m)
    assert pytest.approx(1) == pyo.value(m.x)

@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
@pytest.mark.unit
def test_default_solver():
    """Test that default solver returns the correct solver type
    """
    assert type(pyo.SolverFactory("default")) \
        == type(pyo.SolverFactory(idaes.cfg.default_solver))
