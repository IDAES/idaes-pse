Logging
=======

.. module:: idaes.logger
  :noindex:

IDAES provides some logging extensions to provide finer control over information
logging and to allow solver output to be logged. Logging can be a useful tool for debugging.

Getting Loggers
---------------

There are four main roots of IDAES loggers (``idaes``, ``idaes.model``,
``idaes.init``, ``idaes.solve``).  All of these loggers are standard Python
loggers, and can be used as such.  The main differences between using the IDAES
logging functions to get the loggers and plain Python methods are that
the IDAES functions make it a little easier to get loggers that fit into IDAES's
standard logging hierarchy, and the IDAES loggers have a few additional named
logging levels, which allow for finer control over the information displayed.
Logging levels are described in detail later.

A tag can also be specified and used to filter logging records.  By default the
tag is None and log records won't be filtered. Valid tags are in the set ``{None,
"framework", "model", "flowsheet", "unit", "control_volume", "properties",
"reactions"}``.  Users may add to the set of valid names. To see how
to control which logging tags are logged, see section "Tags" below. To avoid
filtering out import warning and error messages, records logged at the WARNING
level and above are not filtered out regardless of tag.

idaes Logger
~~~~~~~~~~~~

Loggers descending from ``idaes`` (other than ``idaes.init``, ``idaes.model``, or
``idaes.solve``) are used for general IDAES framework logging. Typically the
module name ``__name__`` is used for the logger name. Modules in the ``idaes``
package already start with ``idaes``, but if an IDAES logger is requested for a
module outside of the ``idaes`` package ``idaes.`` is prepended to the name.

.. autofunction:: getLogger
  :noindex:

*Example*

.. testcode::

  import idaes.logger as idaeslog

  _log = idaeslog.getLogger(__name__, tag="framework")


idaes.init Loggers
~~~~~~~~~~~~~~~~~~

The init logger will always descend from "idaes.init". This logger is used in
IDAES model initialization methods, and can be used in user models as well.
Initialization methods are usually attached to a Pyomo Block. Blocks have a
``name`` attribute. So the logger name is usually given as the block name, and
the ``getInitLogger()`` function prepends ``idaes.init.``. The advantage of using
the block name over the module name is that users can see exactly which model
instance the initialization log messages are coming from.

.. autofunction:: getInitLogger
  :noindex:

*Example*

.. testcode::

  import idaes.logger as idaeslog

  class DummyBlock(object):
    """A dummy block for demonstration purposes"""
    def __init__(name):
      self.name = name

    def initialize(outlvl=idaeslog.INFO):
      init_log = idaeslog.getInitLogger(self.name, level=outlvl, tag="unit")

idaes.model Loggers
~~~~~~~~~~~~~~~~~~~

The model logger is used to provide a standard way to produce log messages from
user models that are not part of the ``idaes`` package. The logger name has
``idaes.model`` prepended to the name provided by the user. This is convenient
because it provides a way to use a standard configuration system for user model
loggers. The user can choose any name they like for these loggers.

.. autofunction:: getModelLogger
  :noindex:

*Example*

.. testcode::

  import idaes.logger as idaeslog

  _log = idaeslog.getModelLogger("my_model", level=idaeslog.DEBUG, tag="model")

idaes.solve Loggers
~~~~~~~~~~~~~~~~~~~

The solve logger will always descend from "idaes.solve". This logger is
used to log solver output. Since solvers may produce a lot of output,
it can be useful to specify different handlers for the solve logger to
direct it to a separate file.

.. autofunction:: getSolveLogger
  :noindex:

Tags
----

Logger tags are provided to allow control over what types of log records
to display. The logger tag is just a string that gets attached to a
logger, which specifies that a logger generates records of a certain
type. You can then specify what tags you want to see information from.
A filter removes any tags that are not in the list of tags to display at
levels below WARNING.

The set of tags to display information from is a global setting in the
idaes.logger module. When getting a logger, you can set its tag by
providing the ``tag`` argument, see "Getting Loggers" above.

The following functions can be used to specify which logging tags to
display:

.. autofunction:: log_tags
  :noindex:

.. autofunction:: set_log_tags
  :noindex:

.. autofunction:: add_log_tag
  :noindex:

.. autofunction:: remove_log_tag
  :noindex:

The tags are validated against a list of valid tags to provide error checking
for typos and to enforce some standard tag names. To provide more flexibility,
users can add to the list of valid tag names, but cannot remove names.

.. autofunction:: valid_log_tags
  :noindex:

.. autofunction:: add_valid_log_tag
  :noindex:

Levels
------

Several logging level constants are defined in the ``idaes.logger`` module. These
include the standard Python Levels. The following levels are provided for IDAES
loggers. The additional levels of info provide finer control over the amount of
logging information produced by IDAES loggers.

===================== ====== ============ ============================
Constant Name         Value  Name         Log Method
===================== ====== ============ ============================
CRITICAL              50     CRITICAL     ``critical()``
ERROR                 40     ERROR        ``error()``, ``exception()``
WARNING               30     WARNING      ``warning()``
INFO_LOW              21     INFO         ``info_low()``
INFO                  20     INFO         ``info()``
INFO_HIGH             19     INFO         ``info_high()``
DEBUG                 10     DEBUG        ``debug()``
NOTSET                0      NOTSET       --
===================== ====== ============ ============================


Utility Functions
-----------------

There are some additional utility functions to perform logging tasks that are
common in the IDAES framework.

.. autofunction:: condition
  :noindex:


Logging Solver Output
---------------------

The solver output can be captured directly as part of the ``solve``
command by passing in a ``LogStream`` object and specifying the
preferred logger and logging level. This can be combined with other
supported ``tee`` values. Note that this **ONLY** works with the newest
version of the Pyomo solvers (also known as ``v2`` solvers).

*Example*

.. testcode::

  import sys, logging
  import pyomo.environ as pyo
  from pyomo.common.log import LogStream

  logger = logging.getLogger("logging.demo")

  # Only available with v2 solvers
  solver = pyo.SolverFactory("ipopt_v2")

  model = pyo.ConcreteModel()
  model.x = pyo.Var()
  model.y = pyo.Var()
  model.x.fix(3)
  model.c = pyo.Constraint(expr=model.y==model.x**2)

  # Direct all output to a logger
  res = solver.solve(model, tee=LogStream(level=logging.INFO, logger=logger))
  # Direct output to both stdout and a logger
  res = solver.solve(model, tee=[sys.stdout, LogStream(level=logging.INFO, logger=logger)])

If you **MUST** use version 1 solvers, solver output can be captured using the
``idaes.logger.solver_log(logger, level)`` context manager. The ``logger``
argument is the logger to log to, and the ``level`` argument is the
level at which records are sent to the logger. The output is logged by a
separate logging thread, so output can be logged as it is produced
instead of after the solve completes.  If the ``solver_log()`` context
manager is used, it can be turned on and off by using the
``idaes.logger.solver_capture_on()`` and
``idaes.logger.solver_capture_off()`` functions.  If the capture is off
solver output won't be logged and it will go to standard output as
usual.

The ``solver_log`` context yields an object with the ``tee``
attribute, which should be passed to the ``tee``
argument of the ``solve`` method. Tee tells the Pyomo solver to
display solver output. The solver log context can provide this argument
by determining if the solver output would be logged at the given level.

*Example*

.. testcode::

  import idaes.logger as idaeslog
  import pyomo.environ as pyo

  solver = pyo.SolverFactory("ipopt")

  model = pyo.ConcreteModel()
  model.x = pyo.Var()
  model.y = pyo.Var()
  model.x.fix(3)
  model.c = pyo.Constraint(expr=model.y==model.x**2)

  log = idaeslog.getSolveLogger("solver.demo")
  log.setLevel(idaeslog.DEBUG)

  with idaeslog.solver_log(log, idaeslog.DEBUG) as slc:
    res = solver.solve(model, tee=slc.tee)
