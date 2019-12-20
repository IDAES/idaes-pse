
import idaes.logger as idaeslog
import logging

__author__ = "John Eslick"

def test_get_idaes_logger(caplog):
    caplog.set_level(logging.DEBUG)
    log = idaeslog.getLogger("My Test Logger 1")
    assert log.name == "idaes.My Test Logger 1"
    log.info_least("Hello!")
    log.info_less("Hello!")
    log.info("Hello!")
    log.info_more("Hello!")
    log.info_most("Hello!")
    for record in caplog.records:
        assert record.levelname == "INFO"
    log = idaeslog.getLogger("idaes.My Test Logger 2")
    assert log.name == "idaes.My Test Logger 2"

def test_get_model_logger(caplog):
    log = idaeslog.getModelLogger("My Model 1")
    assert isinstance(log, logging.Logger)
    assert log.name == "idaes.model.My Model 1"
    caplog.set_level(idaeslog.INFO_LESS)
    log.info_least("Hello! from least")
    log.info_less("Hello! from less")
    log.info("Hello! from info")
    log.info_more("Hello! from more")
    log.info_most("Hello! from most")
    for record in caplog.records:
        assert record.message in ["Hello! from least", "Hello! from less"]
    log = idaeslog.getModelLogger("idaes.My Model 2")
    assert log.name == "idaes.model.My Model 2"

def test_get_init_logger(caplog):
    log = idaeslog.getInitLogger("My Init 1")
    assert log.name == "idaes.init.My Init 1"

if __name__ == "__main__":
    test_get_model_logger()
