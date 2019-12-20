Logging
=======

IDAES provides the ``idaes.logger`` module to assist with a few common logging
tasks, and the make it to get a logger that descends from one of the IDAES standard
loggers (``idaes``, ``idaes.model``, or ``idaes.init``).


Getting Loggers
---------------

There are three main roots of IDAES loggers.

``idaes``
~~~~~~~~~

.. testcode::

  import idaes.logger as idaeslog

  _log = idaeslog.getLogger(__name__)

The 
