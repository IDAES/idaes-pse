idaes get-extensions: Get solvers and libraries
===============================================

This page lists the options for the idaes "get-extensions" subcommand.
This is invoked like::

    idaes [general options] get-extensions [subcommand options]


.. program:: idaes

general options
---------------
The following general options from the `idaes` base command
affect the get-extensions subcommand. They should be placed *before* the
"get-extensions" subcommand, on the command-line.

* -v/--verbose
* -q/--quiet

See :ref:`idaes-base-command` for details.

idaes get-extensions
--------------------

This subcommand gets the compiled solvers and libraries
from a remote repository, and installs them locally.

.. program:: idaes get-extensions

options
^^^^^^^

.. option:: --help

    Show the help message and exit.

.. option:: --url

    URL from which to download the solvers/libraries.
