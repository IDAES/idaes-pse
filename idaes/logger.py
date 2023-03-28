#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
# TODO: Missing doc strings
# pylint: disable=missing-module-docstring
# pylint: disable=missing-class-docstring

import logging
from contextlib import contextmanager

from pyomo.common.tee import capture_output
from pyomo.common.log import LogStream

import idaes

# Throw the standard levels in here, just let you access it all in one place
CRITICAL = logging.CRITICAL  # 50
ERROR = logging.ERROR  # 40
WARNING = logging.WARNING  # 30
INFO_LOW = 21  # Most important info
INFO = logging.INFO  # 20  #Medium info (default)
INFO_HIGH = 19  # Less important important info
DEBUG = logging.DEBUG  # 10
NOTSET = logging.NOTSET  # 0


levelname = {  # the level name of all our extra info levels is "INFO"
    INFO_HIGH: "INFO",
    INFO_LOW: "INFO",
}


class _TagFilter(logging.Filter):
    """Filter applied to IDAES loggers returned by this module."""

    @staticmethod
    def filter(record):
        """Add in the custom level name and let the record through"""
        if record.levelno in levelname:
            record.levelname = levelname[record.levelno]
        if record.levelno >= WARNING:
            return True
        try:
            if record.tag is None or record.tag in idaes.cfg.logger_tags:
                return True
        except AttributeError:
            return True
        return False


def __info_low(self, *args, **kwargs):
    self.log(INFO_LOW, *args, **kwargs)


def __info_high(self, *args, **kwargs):
    self.log(INFO_HIGH, *args, **kwargs)


def __add_methods(log, tag=None):
    # pylint: disable=assignment-from-no-return,no-value-for-parameter
    log.addFilter(_TagFilter)
    log = logging.LoggerAdapter(log, {"tag": tag})
    log.info_high = __info_high.__get__(log)
    log.info_low = __info_low.__get__(log)
    # hopefully adding this multiple times is not a problem
    return log


def _getLogger(name, logger_name="idaes", level=None, tag=None):
    assert tag in idaes.cfg.valid_logger_tags.union({None})
    if name.startswith("idaes."):
        name = name[6:]
    name = ".".join([logger_name, name])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return __add_methods(logging.getLogger(name), tag)


def getIdaesLogger(name, level=None, tag=None):
    """Return an idaes logger.

    Args:
        name: usually __name__
        level: standard IDAES logging level (default use IDAES config)
        tag: logger tag for filtering, see valid_log_tags()

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes", level=level, tag=tag)


getLogger = getIdaesLogger


def getSolveLogger(name, level=None, tag=None):
    """Get a solver logger

    Args:
        name: logger name is "idaes.solve." + name (if name starts with "idaes."
            it is removed before creating the logger name)
        level: Log level
        tag: logger tag for filtering, see valid_log_tags()

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.solve", level=level, tag=tag)


def getInitLogger(name, level=None, tag=None):
    """Get a model initialization logger

    Args:
        name: Object name (usually Pyomo Component name)
        level: Log level
        tag: logger tag for filtering, see valid_log_tags()

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.init", level=level, tag=tag)


def getModelLogger(name, level=None, tag=None):
    """Get a logger for an IDAES model. This function helps users keep their
    loggers in a standard location and use the IDAES logging config.

    Args:
        name: Name (usually __name__).  Any starting 'idaes.' is stripped off, so
            if a model is part of the idaes package, 'idaes' won't be repeated.
        level: Standard Python logging level (default use IDAES config)
        tag: logger tag for filtering, see valid_log_tags()

    Returns:
        logger
    """
    return _getLogger(name=name, logger_name="idaes.model", level=level, tag=tag)


def condition(res):
    """Get the solver termination condition to log.  This isn't a specific value
    that you can really depend on, just a message to pass on from the solver for
    the user's benefit. Sometimes the solve is in a try-except, so we'll handle
    None and str for those cases, where you don't have a real result."""

    if res is None:
        return "Error, no result"
    elif isinstance(res, str):
        return res
    else:
        s = str(res.solver.termination_condition)

    try:
        if "ipopt" in str(res.solver.message).lower():
            solver_message = " ".join(str(res.solver.message).split(" ")[2:])
            return f"{s} - {solver_message}"
        else:
            return f"{s} - {str(res.solver.message)}"
    except:  # pylint: disable=bare-except
        return s


def solver_capture_on():
    """This function turns on the solver capture for the solver_log context
    manager. If this is on, solver output within the solver_log context
    is captured and sent to the logger.

    """
    idaes.cfg.logger_capture_solver = True


def solver_capture_off():
    """This function turns off the solver capture for the solver_log context
    manager. If this is off, solver output within the solver_log context
    is just sent to stdout like normal.

    """
    idaes.cfg.logger_capture_solver = False


def solver_capture():
    """Return True if solver capture is on or False otherwise."""
    return idaes.cfg.logger_capture_solver


def log_tags():
    """Returns a set of logging tags to be logged.

    Returns:
        (set) tags to be logged
    """
    return idaes.cfg.logger_tags


def set_log_tags(tags):
    """Specify a set of tags to be logged

    Args:
        tags(iterable of str): Tags to log

    Returns:
        None
    """
    for m in tags:
        if m not in idaes.cfg.valid_logger_tags.union({None}):
            raise ValueError(f"{m} is not a valid logging tag")
    idaes.cfg.logger_tags = set(tags)


def add_log_tag(tag):
    """Add a tag to the list of tags to log.

    Args:
        tag(str): Tag to log

    Returns:
        None
    """
    if tag not in idaes.cfg.valid_logger_tags.union({None}):
        raise ValueError(f"{tag} is not a valid logging tag")
    idaes.cfg.logger_tags.add(tag)


def remove_log_tag(tag):
    """Remove a tag from the list of tags to log.

    Args:
        tag(str): Tag to no longer log

    Returns:
        None
    """
    try:
        idaes.cfg.logger_tags.remove(tag)
    except ValueError:
        pass


def valid_log_tags():
    """Returns a set of valid logging tag names.

    Returns:
        (set) valid tag names
    """
    return idaes.cfg.valid_logger_tags.union({None})


def add_valid_log_tag(tag):
    """Add a tag name to the list of valid names.

    Args:
        tag(str): A tag name

    Returns:
        None
    """
    assert isinstance(tag, str)
    idaes.cfg.valid_logger_tags.add(tag)


class SolverLogInfo(object):
    def __init__(self, tee=True):
        self.tee = tee


@contextmanager
def solver_log(logger, level=logging.ERROR):
    """Context manager to send solver output to a logger."""
    tee = logger.isEnabledFor(level)
    if not solver_capture():
        yield SolverLogInfo(tee=tee)
    else:
        with capture_output(LogStream(level, logger)):
            yield SolverLogInfo(tee=tee)
