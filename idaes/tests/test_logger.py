
import idaes
import idaes.logger as ill
import logging

__author__ = "John Eslick"

def test_custom_levels():
    log = ill.getModelLogger("My Model")

    #assert logging.getLevelName(22) == 'INFO_LEAST'
    #assert logging.getLevelName(21) == 'INFO_LESS'
    #assert logging.getLevelName(19) == 'INFO_MORE'
    #assert logging.getLevelName(18) == 'INFO_MOST'

def test_get_model_logger():
    log = ill.getModelLogger("My Model")
    assert isinstance(log, logging.Logger)
    assert log.name == "idaes.model.My Model"

    log = ill.getModelLogger("My Model")
    log.info_least("hello")
    #log.setLevel(logging.INFO_MORE)
    #log.info_more("Hello, nice to meet you!")

if __name__ == "__main__":
    test_get_model_logger()
