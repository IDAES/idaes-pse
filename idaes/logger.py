import logging
import bisect

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


def getIdaesLogger(name, level=None):
    """ Return an idaes logger.

    Args:
        name: usually __name__
        level: standard IDAES logging level (default use IDAES config)

    Returns:
        logger
    """
    # this function is fairly useless right now, but it helps to standardize

    # This probably does nothing, but if you ask for an idaes logger outside
    # the idaes package it makes sure you get one.
    name = ".".join(["idaes", name.lstrip("idaes.")])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return __add_methods(logging.getLogger(name))


getLogger = getIdaesLogger


def getInitLogger(name, level=None):
    """ Get a model initialization logger

    Args:
        name: Object name (usually Pyomo Component name)
        level: Log level

    Returns:
        logger
    """
    name = ".".join(["idaes.init", name])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return __add_methods(l)


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
    name = ".".join(["idaes.model", name.lstrip("idaes.")])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return __add_methods(logging.getLogger(name))


def increased_output(logger):
    i = bisect.bisect_left(_defined_levels, logger.getEffectiveLevel()) - 1
    if i < 0:
        i = 0
    return _defined_levels[i]


def decreased_output(logger):
    i = bisect.bisect_left(_defined_levels, logger.getEffectiveLevel()) + 1
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
