idaes lib-directory: Show IDAES library file directory
======================================================

This page lists the options for the idaes "lib-directory" subcommand.
This is invoked like::

    idaes [general options] lib-directory [subcommand options]


.. program:: idaes

general options
---------------
The following general options from the `idaes` base command
affect the lib-directory subcommand. They should be placed *before* the
"lib-directory" subcommand, on the command-line.

* -v/--verbose
* -q/--quiet

See the :ref:`idaes-base-command` for details.

idaes lib-directory
-------------------

This subcommand shows the IDAES library file directory.

.. program:: idaes lib-directory

options
^^^^^^^

.. option:: --help

    Show the help message and exit.

.. option::  --exists

  Show if the directory exists.

.. option::  --create

    Create the directory.
