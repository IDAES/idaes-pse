import logging
import bisect
import threading
from collections.abc import Iterable

from contextlib import contextmanager
from pyutilib.misc import capture_output

# There isn't really any reason to neforce this other that to catch typing errors
# and enforce standards.
_valid_modules = set(
    [
        None,
        "model",
        "flowsheet",
        "unit",
        "control_volume",
        "properties",
        "reactions",
    ]
)


_config = {
    "solver_capture":True,
    "modules": set(["model", "flowsheet", "unit"]),
}


# Throw the standard levels in here, just let you access it all in one place
CRITICAL = logging.CRITICAL # 50
ERROR = logging.ERROR # 40
WARNING = logging.WARNING # 30
INFO_LOW = 21 #Most important info
INFO = logging.INFO # 20  #Medium info (default)
INFO_HIGH = 19 #Less improtant important info
DEBUG = logging.DEBUG # 10
NOTSET = logging.NOTSET # 0


levelname = { # the level name of all our extra info levels is "INFO"
    INFO_HIGH: "INFO",
    INFO_LOW: "INFO",
}


class _ModuleFilter(logging.Filter):
    """Filter applied to IDAES loggers returned by this modulue."""
    def filter(record):
        """Add in the custom level name and let the record through"""
        if record.levelno in levelname:
            record.levelname = levelname[record.levelno]
        if record.levelno >= WARNING:
            return True
        try:
            if record.idaesmod is None or record.idaesmod in _config["modules"]:
                return True
        except AttributeError:
            return True
        return False


def __info_low(self, *args, **kwargs):
    self.log(INFO_LOW, *args, **kwargs)


def __info_high(self, *args, **kwargs):
    self.log(INFO_HIGH, *args, **kwargs)


def __add_methods(log, module=None):
    log.addFilter(_ModuleFilter)
    log = logging.LoggerAdapter(log, {"idaesmod": module})
    log.info_high = __info_high.__get__(log)
    log.info_low = __info_low.__get__(log)
    # hopfully adding this multiple times is not a problem
    return log


def _getLogger(name, logger_name="idaes", level=None, module=None):
    assert module in _valid_modules
    if name.startswith("idaes."):
        name = name[6:]
    name = ".".join([logger_name, name])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return __add_methods(logging.getLogger(name), module)


def getIdaesLogger(name, level=None, module=None):
    """ Return an idaes logger.

    Args:
        name: usually __name__
        level: standard IDAES logging level (default use IDAES config)

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes", level=level, module=module)


getLogger = getIdaesLogger

def getSolveLogger(name, level=None, module=None):
    """ Get a model initialization logger

    Args:
        name: logger name is "idaes.solve." + name (if name starts with "idaes."
            it is removed before creating the logger name)
        level: Log level

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.solve", level=level, module=module)

def getInitLogger(name, level=None, module=None):
    """ Get a model initialization logger

    Args:
        name: Object name (usually Pyomo Component name)
        level: Log level

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.init", level=level, module=module)


def getModelLogger(name, level=None, module=None):
    """ Get a logger for an IDAES model. This function helps users keep their
    loggers in a standard location and using the IDAES logging config.

    Args:
        name: Name (usually __name__).  Any starting 'idaes.' is stripped off, so
            if a model is part of the idaes package, idaes won't be repeated.
        level: Standard Python logging level (default use IDAES config)

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.model", level=level, module=module)

def condition(res):
    """Get the solver termination condition to log.  This isn't a specifc value
    that you can really depend on, just a message to pass on from the solver for
    the user's benefit. Somtimes the solve is in a try-except, so we'll handle
    None and str for those cases, where you don't have a real result."""

    if res is None:
        return "Error, no result"
    elif isinstance(res, str):
        return res
    else:
        s = str(res.solver.termination_condition)

    try:
        return "{} - {}".format(s, str(res.solver.message))
    except:
        return s

def solver_capture_on():
    """This function turns on the solver capture for the solver_log context
    manager. If this is on, solver output within the solver_log contex is captured
    and sent to the logger.
    """
    _config["solver_capture"] = True

def solver_capture_off():
    """This function turns off the solver capture for the solver_log context
    manager. If this is off solver output within the solver_log contex is just
    sent to stdout like normal.
    """
    _config["solver_capture"] = False

def solver_capture():
    """Return True if solver capture is on or False otherwise."""
    return _config["solver_capture"]

def log_modules():
    return _config["modules"]

def set_log_modules(modules):
    for m in modules:
        if m not in _valid_modules:
            raise ValueError("{} is not a valid logging module".format(m))
    _config["modules"] = set(modules)

def add_log_module(module):
    if m not in _valid_modules:
        raise ValueError("{} is not a valid logging module".format(m))
    _config["modules"].add(module)

def remove_log_module(module):
    try:
        _config["modules"].remove(module)
    except ValueError:
        pass

def valid_log_modules():
    return _config["modules"]

def add_valid_log_module(module):
    _config["modules"].add(module)

def remove_valid_log_module(module):
    try:
        __config["modules"].remove(module)
    except ValueError:
        pass


class IOToLogTread(threading.Thread):
    """This is a Thread class that can log solver messages and show them as
    they ar produced, while the main thread is waiting on the solver to finish"""

    def __init__(self, stream, logger, sleep=1.0, level=logging.ERROR):
        super().__init__(daemon=True)
        self.log = logger
        self.level = level
        self.stream = stream
        self.sleep = sleep
        self.stop = threading.Event()
        self.pos=0

    def log_value(self):
        try:
            v = self.stream.getvalue()[self.pos:]
        except ValueError:
            self.stop.set()
            return
        self.pos += len(v)
        for l in v.split("\n"):
            if l:
                self.log.log(self.level, l.strip())

    def run(self):
        while True:
            self.log_value()
            self.stop.wait(self.sleep)
            if self.stop.isSet():
                self.log_value()
                self.pos=0
                return


class SolverLogInfo(object):
    def __init__(self, tee=True, thread=None):
        self.tee = tee
        self.thread = thread


@contextmanager
def solver_log(logger, level=logging.ERROR):
    """Context manager to send solver output to a logger.  This uses a seperate
    thread to log solver output while the solver is running"""
    # wait 3 seconds to  join thread.  Should be plenty of time.  In case
    # something goes horibly wrong though don't want to hang.  The logging
    # thread is daemonic, so it will shut down with the main process even if it
    # stays around for some mysterious reason while the model is running.
    join_timeout = 3
    tee = logger.isEnabledFor(level)
    if not solver_capture():
        yield SolverLogInfo(tee=tee)
    else:
        with capture_output() as s:
            lt = IOToLogTread(s, logger=logger, level=level)
            lt.start()
            try:
                yield SolverLogInfo(tee=tee, thread=lt)
            except:
                lt.stop.set()
                lt.join(timeout=join_timeout)
                raise
        # thread should end when s is closed, but the setting stop makes sure
        # the last of the output gets logged before closing s
        lt.stop.set()
        lt.join(timeout=join_timeout)
