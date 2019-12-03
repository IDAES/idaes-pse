import logging

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
    return logging.getLogger(name)

def solver_tee(logger, tee_level=logging.DEBUG):
    """Function to produce solver output based on the logging level of a specific
    logger. This function just helps standardize the level for solver output to
    appear and make code a bit cleaner.

    Args:
        logger: logger to get output level from
        tee_level: Level at which to show solver output, usually use default

    Returns
        (bool)
    """
    return logger.getEffectiveLevel() <= tee_level

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
    return logger.getEffectiveLevel() <= tee_level

def condition(res):
    """Get the solver termination condition.  Since it seems to be common to
    have an if block to check for None if the solver call raised a handled
    exception"""

    if res is None:
        return "Error, no result"
    elif isinstance(res, str):
        return res
    else:
        return str(res.solver.termination_condition)

def getInitLogger(name, level=None):
    """ Get a model initialization logger

    Args:
        name: Object name (usually Pyomo Component name)
        level: Logging detail level (for initialization routines 1 to 6)
             * 0 = Use default idaes.init logger setting
             * 1 = Maximum output
             * 2 = Include solver output
             * 3 = Return solver state for each step in subroutines
             * 4 = Return solver state for each step in routine
             * 5 = Final initialization status and exceptions
             * 6 = No output

    Returns:
        logger
    """
    name = ".".join(["idaes.init", name])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return l

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
    return logging.getLogger(name)
