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

.. option:: --release <release number e.g. 4.0.0>

    Specify an official binary release version number to download.  If this is
    not specified the default release will be used. Releases 4.x and newer use
    a different set of files and install layout than 3.x; see
    :ref:`Binary Packages <tutorials/getting_started/binaries:Binary Packages>`.

.. option:: --url <url>

    URL from which to download the solvers/libraries.  If this is not
    provided the URL of the default release will be used. This cannot be
    combined with ``--release``, so the files at the URL must use the same
    layout (4.x or 3.x) as the default release.

.. option:: --distro <distro>

    OS or Linux distribution to download binaries for, e.g. ``ubuntu2404``
    (default: detect automatically). See
    :ref:`Supported Systems <tutorials/getting_started/binaries:Supported Systems>`.

.. option:: --insecure

    Don't verify download location

.. option:: --cacert <ca>

    Specify certificate file to verify download location

.. option:: --nochecksum

    Don't verify the file checksum

.. option:: --library-only

    Only install shared physical property function libraries (and, in 3.x, any
    specified extras), not solvers.

.. option:: --no-download

    Don't download anything, but report what would be done, including the
    release and platform that would be used.

.. option:: --info

    Show information about available binary builds for the release (the
    default release, or the one given with ``--release``). Lists all platforms
    and architectures and shows which (if any) matches the current system.
    If given, all other options are ignored.

.. option:: --extra <extra>

    Add an extra binary package to the things to install. You can specify the
    extra option multiple times for multiple extras. Extras only exist in 3.x
    (``petsc``); in 4.x, PETSc is included in the main package and this option
    is ignored.

.. option:: --extras-only

    Only install the extras given with ``--extra`` (3.x only).

.. option:: --to <path>

    Put extensions in a alternate location.  This can be used to just download
    and extract the binaries. It lets you download the files without putting
    them in IDAES's data directory.

.. option:: --verbose

    Show details

Related commands
^^^^^^^^^^^^^^^^

- ``idaes extensions-version``: show the version of the installed extensions.
- ``idaes extensions-license``: show the license of the installed extensions.
- ``idaes bin-platform``: show the binary build that matches the current system.
