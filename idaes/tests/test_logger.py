import pytest
import idaes.logger as idaeslog
import logging
import pyomo.environ as pyo

__author__ = "John Eslick"

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

def test_get_model_logger(caplog):
    log = idaeslog.getModelLogger("My Model 1")
    assert isinstance(log, logging.LoggerAdapter)
    assert log.name == "idaes.model.My Model 1"
    caplog.set_level(idaeslog.INFO)
    log.info_low("Hello! from unit")
    log.info("Hello! from info")
    log.info_high("Hello! from unit high")
    for record in caplog.records:
        assert record.message in ["Hello! INFO", "Hello! from INFO_LOW"]
    log = idaeslog.getModelLogger("idaes.My Model 2")
    assert log.name == "idaes.model.My Model 2"

def test_get_init_logger():
    log = idaeslog.getInitLogger("My Init 1")
    assert log.name == "idaes.init.My Init 1"

def test_solver_condition():
    # test the results that can be tested without a solver
    assert idaeslog.condition(None) == "Error, no result"
    assert idaeslog.condition("something else") == "something else"

def test_tags(caplog):
    def a(tag):
        caplog.set_level(logging.DEBUG)
        log = idaeslog.getLogger("Unit", tag=tag)
        log.setLevel(logging.DEBUG)
        log.info_high("Hello!")
        log.info("Hello!")
        log.info_low("Hello!")
        if tag not in idaeslog.log_tags():
             assert len(caplog.records) == 0
        else:
            assert caplog.records[0].levelname == "INFO"
            assert caplog.records[1].levelname == "INFO"
            assert caplog.records[2].levelname == "INFO"
    for m in idaeslog.valid_log_tags():
        a(m)

@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
def test_solver_condition2():
    solver = pyo.SolverFactory('ipopt')
    model = pyo.ConcreteModel("Solver Result Test Model")
    model.x = pyo.Var([1,2])
    model.y = pyo.Var(initialize=5)
    model.x.fix(2)
    model.y.unfix()
    model.c = pyo.Constraint(expr=model.x[1] + model.x[2]==model.y)
    res = solver.solve(model)
    assert idaeslog.condition(res).startswith("optimal") # better solve
    model.c2 = pyo.Constraint(expr=model.x[1]==model.y)
    res = solver.solve(model)
    assert idaeslog.condition(res).startswith("other") # too few degrees of freedom

@pytest.mark.skipif(not pyo.SolverFactory('ipopt').available(False), reason="no Ipopt")
def test_solver_log(caplog):
    solver = pyo.SolverFactory('ipopt')
    model = pyo.ConcreteModel("Solver Result Test Model")
    model.x = pyo.Var([1,2])
    model.y = pyo.Var(initialize=5)
    model.x.fix(2)
    model.y.unfix()
    model.c = pyo.Constraint(expr=model.x[1] + model.x[2]==model.y)

    log = idaeslog.getLogger("solver")
    caplog.set_level(idaeslog.DEBUG)
    log.setLevel(idaeslog.DEBUG)

    idaeslog.solver_capture_on()
    with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
        res = solver.solve(model, tee=True)
    assert(not slc.thread.is_alive()) # make sure logging thread is down
    s = ""
    for record in caplog.records:
        s += record.message
    assert "Optimal" in s

    # test that an exception still results in the thread terminating
    try:
        with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
            res = solver.solve(modelf, tee=True)
    except NameError:
        pass # expect name error
    assert(not slc.thread.is_alive()) # make sure logging thread is down
