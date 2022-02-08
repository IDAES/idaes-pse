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

.. option:: --license

    Display the license info for the currently installed extensions, doesn't
    install anything.

.. option:: --release <release number e.g. 2.4.4>

    Specify an official binary release version number to download.  If this is
    not specified the default release will be used.

.. option:: --url <url>

    URL from which to download the solvers/libraries.  If this is not
    provided the URL of the default release will be used.

.. option:: --insecure

    Don't verify download location

.. option:: --cacert <ca>

    Specify certificate file to verify download location

.. option:: --nochecksum

    Don't verify the file checksum

.. option:: --library-only

    Only install shared physical property function libraries, and any specified
    extras not solvers.

.. option:: --no-download

    Don't download anything, but report what would be done

.. option:: --show-current-version

    Just show the version information if any for the currently installed solvers
    and libraries.

.. option:: --show-platforms

    Just show the platform options

.. option:: --show-extras

    Just show list of binary extras

.. option:: --extra <extra>

    Add an extra binary package to the things to install. You can specify the
    extra option multiple times for multiple extras.

.. option:: --to <path>

    Put extensions in a alternate location.  This can be used to just download
    and extract the binaries. It lets you download the files without putting
    them in IDAES's bin directory.

.. option:: --verbose

    Show details
