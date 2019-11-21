import logging

def getIdaesLogger(name, level=None):
    """ Return an idaes logger.

    Args:
        name: usually __name__
        level: standard IDAES logging level (default use IDAES config)

    Returns:
        logger
    """
    # this function is fairly usless right now, but it helps to standardize

    # This probably does nothing, but if you ask for an idaes logger outside
    # the idaes package it makes sure you get one.
    name = ".".join(["idaes", name.lstrip("idaes.")])
    l = logging.getLogger(name)
    if level is not None:
        l.setLevel(level)
    return logging.getLogger(name)

def solver_tee(logger, tee_level=2):
    """Function to use in initialization to determine at a given output level
    whether to use the sovler tee option to print solver output. This function
    just helps standadize the level for solver output to appear and make the
    initilization routine code a bit cleaner.

    Args:
        logger: logger to get output level from
        tee_level: Level at which to show solver output, usually use default

    Returns
        (bool)
    """
    return logger. getEffectiveLevel() <= tee_level

def getInitLogger(name, level=None):
    """ Get a model initilization logger

    Args:
        name: Object name (usually Pyomo Component name)
        level: Logging detail level (for initilization routines 1 to 6)
             * 0 = Use default idaes.init logger setting
             * 1 = Maximum output
             * 2 = Include solver output
             * 3 = return solver state for each step in subroutines
             * 4 = return solver state for each step in routine
             * 5 = Indicate finial initliation status
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
    """ Get a logger for an IDAES model. This function helps users keep thier
    loggers in a stndard location and using the IDAES logging config.

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
