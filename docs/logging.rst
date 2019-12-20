Logging
=======

.. module:: idaes.logger

IDAES provides the ``idaes.logger`` module to assist with a few common logging
tasks, and the make it to get a logger that descends from one of the IDAES standard
loggers (``idaes``, ``idaes.model``, or ``idaes.init``).  Extra logging levels
are also provided for ``idaes`` loggers to allow finer control over information
output.


Getting Loggers
---------------

There are three main roots of IDAES loggers.  All these loggers are just
standard Python loggers, and can be used as such.  The main differences between
using the IDAES logging functions and plain Python methods, is that the IDAES
functions make it a little easier to get loggers that fit into IDAES's standard
logging hierarchy, and the IDAES loggers have a few additional named logging
levels, which allow for a bit finer control over information output. Logging
levels are described in detail later.

idaes
~~~~~

The usually __name__ should be used when creating this logger.  Regardless of
what the module name is, though the root of an IDAES logger will always be "idaes."
Since this is usually used withing the IDAES framework the logger name is usually
``__name__``, but if this is call outside the ``idaes`` package, "idaes." will be
prepended to the name.

.. autofunction:: getLogger

*Example*

.. testcode::

  import idaes.logger as idaeslog

  _log = idaeslog.getLogger(__name__)


idaes.init
~~~~~~~~~~

The init logger will always descend from "idaes.init".  This logger is used in
IDAES model initialization methods. These models are instances of Pyomo's Block
and have a name attribute.  When a used sees initialization log output it is
useful to see which object instance produced the output.  The root "idaes.init"
logger can be configured differently than "idaes" to provided different default
behavior.

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

idaes.model
~~~~~~~~~~~

The model logger is used to provide a standard way to produce log messages from
user models that are not part of the ``idaes`` package. These logger names
prepend "idaes.model" to the name provided by the user.  This is convenient
because it provides a way to use a standard configuration system for user model
loggers.

.. autofunction:: getModelLogger

*Example*


.. testcode::
  import idaes.logger as idaeslog

  _log = idaeslog.getLogger("my_model")

Levels
------

Several logging level constants are defined in the ``idaes.logger`` module these
include the standard Python Levels.  The following levels are provided for IDAES
loggers.  The additional levels of info provide finer control over the amount of
logging information produced by IDAES loggers.

================== ====== ============ ============================
Constant Name      Value  Name         Log Method
================== ====== ============ ============================
CRITICAL           50     CRITICAL     ``critial()``
ERROR              40     ERROR        ``error()``, ``exception()``
WARNING            30     WARNING      ``warning()``
INFO_LEAST         22     INFO         ``info_least()``
INFO_LESS          21     INFO         ``info_less()``
INFO               20     INFO         ``info()``
INFO_MORE          19     INFO         ``info_more()``
INFO_MOST          18     INFO         ``info_most()``
SOLVER             17     SOLVER       ``solver()``
DEBUG              10     DEBUG        ``debug()``
NOTSET             0      NOTSET       --
================== ====== ============ ============================


Utility Functions
-----------------

There are some additional utility functions to perform logging tasks that are
common in the IDAES framework.

.. autofunction:: solver_tee

.. autofunction:: condition

.. autofunction:: increased_output

.. autofunction:: decreased_output
