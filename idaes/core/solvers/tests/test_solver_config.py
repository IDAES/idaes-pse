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


@pytest.mark.skipif(not pyo.SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_ipopt_idaes_config():
    """
    Test that the default solver options are set
    """
    with idaes.temporary_config_ctx():
        # in this context use idaes default solver options
        isolve.use_idaes_solver_configuration_defaults()
        idaes.cfg.ipopt.options.nlp_scaling_method = "toast-based"
        solver = pyo.SolverFactory("ipopt")
        assert solver.options["nlp_scaling_method"] == "toast-based"
        solver = pyo.SolverFactory("ipopt", options={"tol": 1})
        assert solver.options["tol"] == 1
        idaes.cfg.ipopt.options.tol = 1
        solver = pyo.SolverFactory("ipopt")
        assert solver.options["tol"] == 1


@pytest.mark.skipif(not pyo.SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_default_solver():
    """Test that default solver returns the correct solver type"""
    assert type(pyo.SolverFactory("default")) == type(
        pyo.SolverFactory(idaes.cfg.default_solver)
    )
