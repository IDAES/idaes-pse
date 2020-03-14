Command-line interface
======================
The IDAES PSE Toolkit includes a command-line tool that can be invoked
by typing `idaes` in a UNIX or Mac OSX shell, or Windows Powershell,
that is in an installed IDAES environment. For the most part, this means
that wherever you installed IDAES will have this command available.

This section of the documentation describes the capabilities of this
command-line program.

.. _idaes-base-command:

idaes command
-------------
The base `idaes` command does not do anything by itself, besides set some
shared configuration values. All the real work is done by one of the sub-commands,
each of which is described on a separate page below.


.. toctree::
    :maxdepth: 1

    bin_directory
    copyright
    data_directory
    get_examples
    get_extensions
    lib_directory

shared configuration
^^^^^^^^^^^^^^^^^^^^

.. option:: -v
.. option:: --verbose

Increase verbosity. Show warnings if given once, then info, and then
debugging messages.

.. option:: -q
.. option:: --quiet

Increase quietness. If given once, only show critical messages.
If given twice, show no messages.
