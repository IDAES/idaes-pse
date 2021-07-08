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
import pyomo.environ as pyo
import pytest
import idaes
import idaes.core.solvers as isolve

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
    with idaes.temporary_config_ctx():
        # in this context use idaes default solver options
        isolve.use_idaes_solver_configuration_defaults()
        idaes.cfg.ipopt.options.nlp_scaling_method = "toast-based"
        solver = pyo.SolverFactory('ipopt')
        assert solver.options["nlp_scaling_method"] == "toast-based"
        solver = pyo.SolverFactory('ipopt', options={"tol":1})
        assert solver.options["tol"] == 1
        idaes.cfg.ipopt.options.tol = 1
        solver = pyo.SolverFactory('ipopt')
        assert solver.options["tol"] == 1
    #
    # TODO(JCE): Bring back the rest of this test after the tests
    # untangled in the GT PR
    #
    # back to the original config, don't use idaes default option settings
    #solver = pyo.SolverFactory('ipopt')
    # should not be using idaes defaults here so options should be empty
    # if this fails there is a very good chance another test was set to use
    # idaes defaults.
    #assert "nlp_scaling_method" not in solver.options

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
