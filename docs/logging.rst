Logging
=======

.. module:: idaes.logger

IDAES provides the ``idaes.logger`` module to assist with a few common logging
tasks, and the make it to get a logger that descends from one of the IDAES standard
loggers (``idaes``, ``idaes.model``, ``idaes.init``, ``idaes.solve``).  Extra
logging levels are also provided for ``idaes`` loggers to allow finer control over
information output.


Getting Loggers
---------------

There are three main roots of IDAES loggers.  All these loggers are just
standard Python loggers, and can be used as such.  The main differences between
using the IDAES logging functions to the loggers and plain Python methods, is that
the IDAES functions make it a little easier to get loggers that fit into IDAES's
standard logging hierarchy, and the IDAES loggers have a few additional named
logging levels, which allow for a bit finer control over information output.
Logging levels are described in detail later.

idaes Logger
~~~~~~~~~~~~

Loggers descending from ``idaes`` (other that ``idaes.init`` or ``idaes.model``)
are used for general IDAES framework logging. Typically the module name
``__name__`` is used for the logger name. Modules in the ``idaes`` package already
start with ``idaes``, but if an idaes logger is requested for a module outside the
``idaes`` package "idaes." is prepended to the name.

.. autofunction:: getLogger

*Example*

.. testcode::

  import idaes.logger as idaeslog

  _log = idaeslog.getLogger(__name__)


idaes.init Logger
~~~~~~~~~~~~~~~~~

The init logger will always descend from "idaes.init".  This logger is used in
IDAES model initialization methods, and can be used in user models as well.
Initialization methods are usually attached to a Pyomo Block. Blocks have a
``name`` attribute.  So the logger name is usually given as the block name, and
the ``getInitLogger()`` function prepends ``idaes.init.`` The advantage of using
the block name over the module name is that users can see exactly which model
instance the initialization log messages are coming from.

.. autofunction:: getInitLogger

*Example*

.. testcode::

  import idaes.logger as idaeslog

  class DummyBlock(object):
    """A dummy block for demonstration purposes"""
    def __init__(name):
      self.name = name

    def initialize():
      init_log = idaeslog.getInitLogger(self.name)

idaes.model Logger
~~~~~~~~~~~~~~~~~~

The model logger is used to provide a standard way to produce log messages from
user models that are not part of the ``idaes`` package. The logger name has
"idaes.model" prepended to the name provided by the user.  This is convenient
because it provides a way to use a standard configuration system for user model
loggers.  The name the user can choose any name they like for these loggers.

.. autofunction:: getModelLogger

*Example*

.. testcode::
  import idaes.logger as idaeslog

  _log = idaeslog.getModelLogger("my_model")

  idaes.solve Logger
  ~~~~~~~~~~~~~~~~~

  The solve logger will always descend from "idaes.solve".  This logger is used in
  to log solver output.  Since solvers may produce a lot of output, it can be useful
  to specify different handlers for the solve logger to direct it to a separate file.

  .. autofunction:: getSolveLogger


Levels
------

Several logging level constants are defined in the ``idaes.logger`` module these
include the standard Python Levels.  The following levels are provided for IDAES
loggers.  The additional levels of info provide finer control over the amount of
logging information produced by IDAES loggers.

===================== ====== ============ ============================
Constant Name         Value  Name         Log Method
===================== ====== ============ ============================
CRITICAL              50     CRITICAL     ``critial()``
ERROR                 40     ERROR        ``error()``, ``exception()``
WARNING               30     WARNING      ``warning()``
FLOWSHEET             23     FLOWSHEET    ``flowsheet()
UNIT                  22     UNIT         ``unit()``
INFO                  20     INFO         ``info()``
UNIT_HIGH             19     UNIT         ``unit_high()``
CV                    18     CV           ``cv()``
PROP                  17     PROP         ``prop()``
DEBUG                 10     DEBUG        ``debug()``
NOTSET                0      NOTSET       --
===================== ====== ============ ============================


Utility Functions
-----------------

There are some additional utility functions to perform logging tasks that are
common in the IDAES framework.

.. autofunction:: condition


Logging Solver Output
---------------------

The solver output can be captured and directed to a logger using the
``idaes.logger.solver_log(logger, level)`` context manager, which uses
``pyutilib.misc.capture_output()`` to temporarily redirect ``sys.stdout`` and
``sys.stderr`` to a string buffer.  The logger argument is the logger to log to,
and the level argument is the level to sent records to the logger at. The output
is logged by a separate logging thread, so output can be logged as it is produced
instead of after the solve completes.  If the  ``solver_log()`` context manager
is used, it can be turned on and off by using the ``idaes.logger.solver_capture_on()``
and ``idaes.logger.solver_capture_off()`` functions.  If the capture is off solver
output won't be logged and it will go to standard output as usual.

The ``solver_log`` context yields and object with ``tee`` and ``thread`` attributes.
Thread is the logging thread, which is not needed for most uses. Tee is the ``solve``
``tee`` argument.  Tee tells the Pyomo solvers to produce solver output.  The solver
log context can provide this argument by determining if the solver output would be
logged at the given level.

*Example*

.. testcode::

  import idaes.logging as idaeslog
  import pyomo.environ as pyo

  solver = pyo.SolverFactory("ipopt")

  model.x = pyo.Var()
  model.y = pyo.Var()
  model.x.fix(3)
  model.c = pyo.Constraint(expr=model.y==model.x**2)

  log = idaes.log("solver.demo")
  log.setLevel(idaeslog.DEBUG)

  with solver_log(log, idaeslog.DEBUG) as slc:
      res = solver.solve(model, tee=slc.tee)
