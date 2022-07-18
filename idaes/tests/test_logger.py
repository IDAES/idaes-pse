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
import pytest
import idaes.logger as idaeslog
import logging
import pyomo.environ as pyo

__author__ = "John Eslick"


@pytest.mark.unit
def test_get_idaes_logger(caplog):
    caplog.set_level(logging.DEBUG)
    log = idaeslog.getLogger("My Test Logger 1")
    log.setLevel(logging.DEBUG)
    assert log.name == "idaes.My Test Logger 1"
    log.info_high("Hello!")
    log.info("Hello!")
    log.info_low("Hello!")
    assert caplog.records[0].levelname == "INFO"
    assert caplog.records[1].levelname == "INFO"
    assert caplog.records[2].levelname == "INFO"

    log = idaeslog.getLogger("idaes.My Test Logger 2")
    assert log.name == "idaes.My Test Logger 2"


@pytest.mark.unit
def test_get_model_logger(caplog):
    log = idaeslog.getModelLogger("My Model 1")
    assert isinstance(log, logging.LoggerAdapter)
    assert log.name == "idaes.model.My Model 1"
    caplog.set_level(idaeslog.INFO)
    log.info_low("Hello! from INFO_LOW")
    log.info("Hello! from INFO")
    log.info_high("Hello! from INFO_HIGH")
    for record in caplog.records:
        assert record.message in ["Hello! from INFO", "Hello! from INFO_LOW"]
    log = idaeslog.getModelLogger("idaes.My Model 2")
    assert log.name == "idaes.model.My Model 2"


@pytest.mark.unit
def test_get_init_logger():
    log = idaeslog.getInitLogger("My Init 1")
    assert log.name == "idaes.init.My Init 1"


@pytest.mark.unit
def test_solver_condition():
    # test the results that can be tested without a solver
    assert idaeslog.condition(None) == "Error, no result"
    assert idaeslog.condition("something else") == "something else"


@pytest.mark.unit
def test_tags(caplog):
    idaeslog.remove_log_tag("model")

    def a(tag):
        caplog.clear()
        caplog.set_level(logging.DEBUG)
        log = idaeslog.getLogger("Unit", tag=tag)
        log.setLevel(logging.DEBUG)
        log.info_high("Hello!")
        log.info("Hello!")
        log.info_low("Hello!")
        if tag not in idaeslog.log_tags() and tag is not None:
            assert len(caplog.records) == 0
        else:
            assert caplog.records[0].levelname == "INFO"
            assert caplog.records[1].levelname == "INFO"
            assert caplog.records[2].levelname == "INFO"

    for m in idaeslog.valid_log_tags():
        a(m)


@pytest.mark.unit
def test_add_remove_tags():
    assert "framework" in idaeslog.valid_log_tags()
    assert "framework" in idaeslog.log_tags()
    idaeslog.remove_log_tag("framework")
    assert "framework" not in idaeslog.log_tags()
    idaeslog.add_log_tag("framework")
    assert "framework" in idaeslog.log_tags()
    idaeslog.add_valid_log_tag("model2")
    assert "model2" in idaeslog.valid_log_tags()
    idaeslog.set_log_tags(idaeslog.valid_log_tags())
    assert "model2" in idaeslog.log_tags()


@pytest.mark.skipif(not pyo.SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_solver_condition2():
    solver = pyo.SolverFactory("ipopt")
    model = pyo.ConcreteModel("Solver Result Test Model")
    model.x = pyo.Var([1, 2])
    model.y = pyo.Var(initialize=5)
    model.x.fix(2)
    model.y.unfix()
    model.c = pyo.Constraint(expr=model.x[1] + model.x[2] == model.y)
    res = solver.solve(model)
    assert idaeslog.condition(res).startswith("optimal")  # better solve
    model.c2 = pyo.Constraint(expr=model.x[1] == model.y)
    res = solver.solve(model)
    assert idaeslog.condition(res).startswith("other")  # too few degrees of freedom


@pytest.mark.skipif(not pyo.SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_solver_log(caplog):
    solver = pyo.SolverFactory("ipopt")
    model = pyo.ConcreteModel("Solver Result Test Model")
    model.x = pyo.Var([1, 2])
    model.y = pyo.Var(initialize=5)
    model.x.fix(2)
    model.y.unfix()
    model.c = pyo.Constraint(expr=model.x[1] + model.x[2] == model.y)

    log = idaeslog.getLogger("solver")
    caplog.set_level(idaeslog.DEBUG)
    log.setLevel(idaeslog.DEBUG)

    idaeslog.solver_capture_on()
    with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
        res = solver.solve(model, tee=True)
    assert not slc.thread.is_alive()  # make sure logging thread is down
    s = ""
    for record in caplog.records:
        s += record.message
    assert "Optimal" in s

    # test that an exception still results in the thread terminating
    try:
        with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
            res = solver.solve(modelf, tee=True)
    except NameError:
        pass  # expect name error
    assert not slc.thread.is_alive()  # make sure logging thread is down
