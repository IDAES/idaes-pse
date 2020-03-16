idaes get-examples: Fetch example scripts and Jupyter Notebooks
===============================================================

This page lists the options for the idaes "get-examples" subcommand.
This is invoked like::

    idaes [general options] get-examples [subcommand options]


.. program:: idaes

general options
---------------
The following general options from the `idaes` base command
affect the get-examples subcommand. They should be placed *before* the
"get-examples" subcommand, on the command-line.

* -v/--verbose
* -q/--quiet

See the :ref:`idaes-base-command` for details.

idaes get-examples
------------------

This subcommand fetches example scripts and Jupyter Notebooks from
a given release in `Github <https://github.com/IDAES/examples-pse/releases>`_.
and puts them in a directory of the users' choosing. If the user does not
specify a directory, the default is `examples`.

.. program:: idaes get-examples

options
^^^^^^^

.. option:: --help

    Show the help message and exit.

.. option:: -d,--dir TEXT

Select the installation target directory. See :ref:`get-examples-usage` for
several examples of this option.

.. _idaes-get-examples-no-install:

.. option:: -I, --no-install

Do *not* install examples into 'idaes_examples' package.
Examples are installed by default so they can be imported directly
from Python. Not installing them might cause some tests, which import
the examples, to fail.

.. _idaes-get-examples-list-releases:

.. option:: -l, --list-releases

List all available released versions, and stop.
This lets people browse the releases and select one. By default,
the release that matches the version of the currently installed "idaes"
package is used. See also the :ref:`--unstable <idaes-get-examples-unstable>` option.

.. option:: -N, --no-download

Do not download anything. If the :ref:`--no-install <idaes-get-examples-no-install>` option
is also given, this means the command will essentially do nothing. Or, looked
at another way, this option means that only  action will be the installation
of the "idaes_examples" package from the selected directory.

.. _idaes-get-examples-unstable:

.. option:: -U, --unstable

Allow and list unstable/pre-release versions. This applies to both download
and the :ref:`--list-releases <idaes-get-examples-list-releases>` option.
Unstable releases are marked with "rcN" or similar suffixes.

.. option::  -V, --version TEXT

Version of examples to download. The default version, which is shown for the
`--help` option, is the same as the version of the IDAES PSE toolkit for which
the `idaes` command is installed. If the version to install is unstable
(ends with "rcN") then you will need to add the :ref:`--unstable <idaes-get-examples-unstable>`
option to avoid errors.

.. _get-examples-usage:

example usage
^^^^^^^^^^^^^

idaes get-examples
    Download examples from release matching release for the ``idaes`` command,
    install them in the `examples` subdirectory of this directory, and
    install the modules found under `examples/src` as a package named `idaes_examples`.
    The `examples` directory must not exist, i.e. the program will refuse to
    overwrite the contents of an existing directory.

idaes get-examples -d /tmp/examples
    Same as above, but put downloaded code in `/tmp/examples` instead.

idaes get-examples -d /tmp/examples -I
    Download to `/tmp/examples`, but do not install.

idaes get-examples -d /tmp/examples -N
    Install the examples found under `/tmp/examples`.

idaes get-examples --version 1.4.2-pre
    Download examples from release `1.4.2-pre`,
    install them in the `examples` subdirectory of this directory, and
    install the modules found under `examples/src` as a package named `idaes_examples`.

idaes get-examples --list-releases
    List available releases of the examples in Github repository, `idaes/examples-pse`.
    Do not attempt to download or install anything.

idaes get-examples --list-releases --unstable
    Same as above, but include non-stable releases.

