import logging
import bisect
import threading

from contextlib import contextmanager
from pyutilib.misc import capture_output

_config = {
    "solver_capture":True,
}

# Throw the standard levels in here, just let you access it all in one place
CRITICAL = logging.CRITICAL # 50
ERROR = logging.ERROR # 40
WARNING = logging.WARNING # 30
INFO_LEAST = 22 # Level for info you almost always want
INFO_LESS = 21
INFO = logging.INFO # 20
INFO_MORE = 19
INFO_MOST = 18 # Level for info you usually don't want
SOLVER = 17 # see solver output and info between INFO and DEBUG
DEBUG = logging.DEBUG # 10
NOTSET = logging.NOTSET # 0

levelname = { # the level name of all our extra info levels is "INFO"
    INFO_LEAST: "INFO",
    INFO_LESS: "INFO",
    INFO_MORE: "INFO",
    INFO_MOST: "INFO",
    SOLVER: "SOLVER",
}

# Unfortunatly it seens to be a common practice in initialization routines
# to call sub model initialization functions with an output level that's slightly
# less.  This list makes it easy to shift up or down one NAMED level.
_defined_levels = (
    DEBUG,
    SOLVER,
    INFO_MOST,
    INFO_MORE,
    INFO,
    INFO_LESS,
    INFO_LEAST,
    WARNING,
    ERROR,
    CRITICAL,
)

class _levelNamesFilter(logging.Filter):
    """Filter applied to IDAES loggers returned by this modulue."""
    def filter(record):
        """Add in the custom level name and let the record through"""
        if record.levelno in levelname:
            record.levelname = levelname[record.levelno]
        return True


def __info_least(self, *args, **kwargs):
    self.log(INFO_LEAST, *args, **kwargs)


def __info_less(self, *args, **kwargs):
    self.log(INFO_LESS, *args, **kwargs)


def __info_more(self, *args, **kwargs):
    self.log(INFO_MORE, *args, **kwargs)


def __info_most(self, *args, **kwargs):
    self.log(INFO_MOST, *args, **kwargs)


def __solver(self, *args, **kwargs):
    self.log(SOLVER, *args, **kwargs)


def __add_methods(log):
    log.info_least = __info_least.__get__(log)
    log.info_less = __info_less.__get__(log)
    log.info_more = __info_more.__get__(log)
    log.info_most = __info_most.__get__(log)
    log.solver = __solver.__get__(log)
    # hopfully adding this multiple times is not a problem
    log.addFilter(_levelNamesFilter)
    return log


def _getLogger(name, logger_name="idaes", level=None):
    if name.startswith("idaes."):
        name = name[6:]
    name = ".".join([logger_name, name])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return __add_methods(logging.getLogger(name))


def getIdaesLogger(name, level=None):
    """ Return an idaes logger.

    Args:
        name: usually __name__
        level: standard IDAES logging level (default use IDAES config)

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes", level=level)


getLogger = getIdaesLogger

def getSolveLogger(name, level=None):
    """ Get a model initialization logger

    Args:
        name: logger name is "idaes.solve." + name (if name starts with "idaes."
            it is removed before creating the logger name)
        level: Log level

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.solve", level=level)

def getInitLogger(name, level=None):
    """ Get a model initialization logger

    Args:
        name: Object name (usually Pyomo Component name)
        level: Log level

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.init", level=level)


def getModelLogger(name, level=None):
    """ Get a logger for an IDAES model. This function helps users keep their
    loggers in a standard location and using the IDAES logging config.

    Args:
        name: Name (usually __name__).  Any starting 'idaes.' is stripped off, so
            if a model is part of the idaes package, idaes won't be repeated.
        level: Standard Python logging level (default use IDAES config)

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.model", level=level)


def increased_output(logger):
    """Get the a logging level that produces one level more output than logger.

    Args:
        logger (logging.Logger|int): Logger to read the level from or level

    Returns:
        (int): Number for the next named level with more output.  At most DEBUG.
    """
    if isinstance(logger, logging.Logger):
        lvl = logger.getEffectiveLevel()
    else:
        lvl = logger
    i = bisect.bisect_left(_defined_levels, lvl) - 1
    if i < 0:
        i = 0
    return _defined_levels[i]


def decreased_output(logger):
    """Get the a logging level that produces one level less output than loggers

    Args:
        logger (logging.Logger|int): Logger to read the level from or level

    Returns:
        (int): Number for the next named level with less output.  At least CRITICAL.
    """
    if isinstance(logger, logging.Logger):
        lvl = logger.getEffectiveLevel()
    else:
        lvl = logger

    i = bisect.bisect_left(_defined_levels, lvl) + 1
    if i >= len(_defined_levels):
        i = -1
    return _defined_levels[i]


def solver_tee(logger, tee_level=SOLVER):
    """Function to produce solver output based on the logging level of a specific
    logger. This function just helps standardize the level for solver output to
    appear and make code a bit cleaner.

    Args:
        logger: logger to get output level from
        tee_level: Level at which to show solver output, usually use default

    Returns
        (bool)
    """
    return logger.isEnabledFor(tee_level)


def init_tee(logger, tee_level=2):
    """Function to use in initialization to determine at a given output level
    whether to use the sovler tee option to print solver output. This function
    just helps standardize the level for solver output to appear and make the
    initialization routine code a bit cleaner.

    Args:
        logger: logger to get output level from
        tee_level: Level at which to show solver output, usually use default

    Returns
        (bool)
    """
    logging.getLogger(__name__).critical("WARNING THIS WILL BE REMOVED")
    return logger.getEffectiveLevel() <= tee_level


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


@contextmanager
def solver_log(logger, level=logging.ERROR):
    """Context manager to send solver output to a logger.  This uses a seperate
    thread to log solver output while the solver is running"""
    # wait 3 seconds to  join thread.  Should be plenty of time.  In case
    # something goes horibly wrong though don't want to hang.  The logging
    # trhead is daemonic, so it will shut down with the main process even if it
    # stays around for some mysterious reason.
    join_timeout = 3
    if not solver_capture():
        yield
    else:
        with capture_output() as s:
            lt = IOToLogTread(s, logger=logger, level=level)
            lt.start()
            try:
                yield lt
            except:
                lt.stop.set()
                lt.join(timeout=join_timeout)
                raise
        # thread should end when s is closed, but the setting stop makes sure
        # the last of the output gets logged before closing s
        lt.stop.set()
        lt.join(timeout=join_timeout)
