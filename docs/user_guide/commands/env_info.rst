idaes environment-info: Get information about the IDAES environment
===================================================================

The idaes "environment-info" command return information about IDAES, Pyomo,
dependencies, and solvers. This information is particularly useful for
debugging.

This page lists the options for the idaes "environment-info" subcommand.
This is invoked like::

    idaes [general options] environment-info [subcommand options]


.. program:: idaes

general options
---------------
The following general options from the `idaes` base command
affect the get-extensions subcommand. They should be placed *before* the
"get-extensions" subcommand, on the command-line.

* -v/--verbose
* -q/--quiet

See :ref:`idaes-base-command` for details.

idaes environment-info
----------------------

This subcommand gets the compiled solvers and libraries
from a remote repository, and installs them locally.

.. program:: environment-info

options
^^^^^^^

.. option:: --help

    Show the help message and exit.

.. option:: --solver <Pyomo solver name>

    The environment-info command returns version for a set of known solvers. You
    can use the "--solver" option to specify additional solvers that are not
    part of the IDAES binary distribution.  The solver option can be supplied
    multiple times to add multiple solvers.

.. option:: --json <file path>

    Write output to the specified json file, instead of to the screen.
