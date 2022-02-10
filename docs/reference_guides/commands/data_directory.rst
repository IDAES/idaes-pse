idaes data-directory: Show IDAES data directory
===============================================

This page lists the options for the idaes "data-directory" subcommand.
This is invoked like::

    idaes [general options] data-directory [subcommand options]


.. program:: idaes

general options
---------------
The following general options from the `idaes` base command
affect the data-directory subcommand. They should be placed *before* the
"data-directory" subcommand, on the command-line.

* -v/--verbose
* -q/--quiet

See the :ref:`idaes-base-command` for details.

idaes data-directory
--------------------

This subcommand shows the IDAES data directory.

.. program:: idaes data-directory

options
^^^^^^^

.. option:: --help

    Show the help message and exit.

.. option::  --exists

  Show if the directory exists.

.. option::  --create

    Create the directory.
