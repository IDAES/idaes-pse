Binary Packages
===============

This section describes the IDAES binary distribution, where it is
installed, how to use it, and alternative installation methods. The IDAES binary
distribution contains mainly third-party solvers compiled for user convenience and
function libraries used for some IDAES physical property packages.

Installation is generally done through the
:doc:`idaes get-extensions command<../../reference_guides/commands/get_extensions>`.

There are two major versions of the IDAES binary distribution (the "IDAES
extensions"):

- **4.x**, installed by default by newer versions of idaes-pse. This page
  describes 4.x unless noted otherwise.
- **3.x**, installed by default by older versions of idaes-pse. Differences
  are noted below.

To see which release your version of idaes-pse installs by default, run
``idaes get-extensions --no-download``.

What is Included
----------------

The National Energy Technology Laboratory (NETL) has obtained a distribution
license for HSL linear solvers that allows inclusion in IDAES's compiled version
of Ipopt.

All technical papers, sales and publicity material resulting from use of
the HSL codes within Ipopt must contain the following acknowledgement: HSL, a
collection of Fortran codes for large-scale scientific computation. See
http://www.hsl.rl.ac.uk.

.. note:: On macOS (Apple Silicon), HSL is built without METIS due to
    compatibility issues, so METIS-based orderings are not available in the HSL
    linear solvers on that platform.

    In 3.x, the macOS x86_64 (Intel) binaries do not include HSL, so HSL
    linear solvers in Ipopt and k_aug and dot_sens are not available on that
    platform.

Solvers
~~~~~~~

- `ipopt <https://coin-or.github.io/Ipopt/>`_ with `HSL linear solvers <http://www.hsl.rl.ac.uk>`_.
- `ipopt_sens <https://github.com/coin-or/Ipopt>`_, sIpopt solver with `HSL linear solvers <http://www.hsl.rl.ac.uk>`_.
- shared ipopt library intended for `cyipopt <https://cyipopt.readthedocs.io/en/stable/>`_ with `HSL linear solvers <http://www.hsl.rl.ac.uk>`_.
- `dot_sens <https://github.com/IDAES/k_aug>`_
- `k_aug <https://github.com/IDAES/k_aug>`_
- `pynumero <https://pyomo.readthedocs.io/en/stable/contributed_packages/pynumero/index.html>`_
- `clp <https://github.com/coin-or/Clp>`_
- `cbc <https://github.com/coin-or/Cbc>`_
- `bonmin <https://github.com/coin-or/Bonmin>`_ (uses Ipopt with HSL)
- `couenne <https://github.com/coin-or/Couenne>`_ (uses Ipopt with HSL)
- petsc, `AMPL solver library <https://ampl.com/REFS/hooking2.pdf>`_ `interface wrapper <https://github.com/IDAES/idaes-ext/tree/main/petsc>`_ for `PETSc's <https://petsc.org/release/>`_ SNES and TS solvers.

In 3.x, PETSc is an optional extra installed with
``idaes get-extensions --extra petsc``. 3.x also includes
`ipopt_l1 <https://github.com/IDAES/Ipopt/tree/restoration_mod>`_, Ipopt with a
modified restoration phase algorithm, which is not included in 4.x.

Function Libraries
~~~~~~~~~~~~~~~~~~

- `general_helmholtz_external <https://github.com/IDAES/idaes-ext/tree/main/src/general_helmholtz>`_, functions for Helmholtz equations of state
- `cubic_roots <https://github.com/IDAES/idaes-ext/tree/main/src/cubic>`_, functions for finding roots of a cubic equation
- `functions <https://github.com/IDAES/idaes-ext/tree/main/src/functions>`_, miscellaneous simple math functions, examples, and test code

Release Files
~~~~~~~~~~~~~

Each release is published on
`the binary releases page <https://github.com/IDAES/idaes-ext/releases>`_ as
a set of tarballs for each platform:

- 4.x: ``idaes-solvers-{platform}.tar.gz`` (solvers, including PETSc) and
  ``idaes-functions-{platform}.tar.gz`` (function libraries).
- 3.x: ``idaes-solvers-{platform}.tar.gz`` (solvers),
  ``idaes-lib-{platform}.tar.gz`` (function libraries), and the optional
  ``idaes-petsc-{platform}.tar.gz``.

Supported Systems
-----------------

Builds are available for Windows, macOS, and several Linux distributions.
Although the binaries are compiled on a smaller number of platforms, one of the
available builds should work on most recent Linux distributions. The table
below shows distributions that should work automatically, and which build is
used for each release series ("--" means no build is available). If you have a
distribution not on the list, you can try to specify a similar distribution
with ``idaes get-extensions --distro <ID>``. Run ``idaes get-extensions --info``
to see which build matches your system.

.. note:: macOS Intel (x86_64) is not supported as of 4.0.0. The last release
    with these binaries is 3.4.2, which can be installed with
    ``idaes get-extensions --release 3.4.2``.

+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Distro/OS                  | ID                                   | 4.x build (arch)             | 3.x build (arch)             |
+============================+======================================+==============================+==============================+
| AlmaLinux 8                | almalinux8                           | el8 (x86_64, aarch64)        | el8 (x86_64, aarch64)        |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| AlmaLinux 9                | almalinux9                           | el9 (x86_64, aarch64)        | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| CentOS 7                   | centos7                              | --                           | el7 (x86_64)                 |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| CentOS 8                   | centos8                              | el8 (x86_64, aarch64)        | el8 (x86_64, aarch64)        |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| CentOS Stream 9            | centos9                              | el9 (x86_64, aarch64)        | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Debian 9                   | debian9                              | --                           | el7 (x86_64)                 |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Debian 10                  | debian10                             | el8 (x86_64, aarch64)        | el8 (x86_64, aarch64)        |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Debian 11                  | debian11                             | ubuntu2004 (x86_64, aarch64) | ubuntu2004 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Debian 12                  | debian12                             | ubuntu2204 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Debian 13                  | debian13                             | ubuntu2404 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Linux Mint 20              | linuxmint20                          | ubuntu2004 (x86_64, aarch64) | ubuntu2004 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Linux Mint 21              | linuxmint21                          | ubuntu2204 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Linux Mint 22              | linuxmint22                          | ubuntu2404 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| macOS                      | darwin                               | darwin (arm64)               | darwin (x86_64, aarch64)     |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Pop!_OS 22.04              | pop22                                | ubuntu2204 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Red Hat Enterprise Linux 7 | rhel7                                | --                           | el7 (x86_64)                 |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Red Hat Enterprise Linux 8 | rhel8                                | el8 (x86_64, aarch64)        | el8 (x86_64, aarch64)        |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Red Hat Enterprise Linux 9 | rhel9                                | el9 (x86_64, aarch64)        | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Rocky Linux 8              | rocky8                               | el8 (x86_64, aarch64)        | el8 (x86_64, aarch64)        |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Rocky Linux 9              | rocky9                               | el9 (x86_64, aarch64)        | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Scientific Linux 7         | scientific7                          | --                           | el7 (x86_64)                 |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Ubuntu 18.04               | ubuntu1804, kubuntu1804, xubuntu1804 | --                           | ubuntu1804 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Ubuntu 20.04               | ubuntu2004, kubuntu2004, xubuntu2004 | ubuntu2004 (x86_64, aarch64) | ubuntu2004 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Ubuntu 22.04               | ubuntu2204, kubuntu2204, xubuntu2204 | ubuntu2204 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Ubuntu 24.04               | ubuntu2404, kubuntu2404, xubuntu2404 | ubuntu2404 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Ubuntu 26.04               | ubuntu2604, kubuntu2604, xubuntu2604 | ubuntu2604 (x86_64, aarch64) | ubuntu2204 (x86_64, aarch64) |
+----------------------------+--------------------------------------+------------------------------+------------------------------+
| Windows 10, 11             | windows                              | windows (x86_64)             | windows (x86_64)             |
+----------------------------+--------------------------------------+------------------------------+------------------------------+

Install Location
----------------

The binaries are installed in the IDAES data directory, which can be found with
the command ``idaes data-directory``. On Windows, by default it is
``%LOCALAPPDATA%\idaes``. On other operating systems it is ``$HOME/.idaes`` by
default. The location of the IDAES data directory can be changed by setting the
``IDAES_DATA`` environment variable. The environment variable method can
accommodate nonstandard setups or multiple IDAES installations.

In 4.x, the data directory uses a standard install layout:

- ``bin/``: solver executables (``idaes bin-directory`` shows this location)
- ``lib/``: shared libraries, Helmholtz data (``lib/helm_data/``), PETSc
  Python utilities (``lib/petscpy/``), and pkg-config files
  (``lib/pkgconfig/``)
- ``include/``: headers (e.g., Ipopt headers under ``include/coin-or/``)
- ``share/``: documentation from the bundled solvers
- version and license files

In 3.x, everything is installed directly in the ``bin/`` directory.

.. note:: When upgrading from 3.x to 4.x, remove the old 3.x files from the
    ``bin/`` directory first. ``idaes get-extensions`` logs an error if it finds
    them, because old libraries left in ``bin/`` can conflict with the new ones.

Manual Installation
-------------------

The install location is not configurable, but the installation step can easily be
done manually. This is occasionally necessary when, for example, a firewall
blocks downloading the binary file from GitHub.

Download the release files for your platform from
`the binary releases page <https://github.com/IDAES/idaes-ext/releases>`_ (see
`Release Files`_), then:

- 4.x: extract ``idaes-solvers-{platform}.tar.gz`` and
  ``idaes-functions-{platform}.tar.gz`` in the IDAES data directory (shown by
  ``idaes data-directory``).
- 3.x: extract ``idaes-solvers-{platform}.tar.gz``,
  ``idaes-lib-{platform}.tar.gz``, and any optional packages in the IDAES
  binary directory (shown by ``idaes bin-directory``).

Using Solvers Without IDAES
---------------------------

Generally, the environment variables are set up to use the IDAES solvers when any
IDAES module is imported in Python.  If you would like to use the solvers in
Pyomo without importing IDAES, you will need to add the install locations to your
executable and dynamic library search path environment variables:

- 4.x: on Linux and macOS, add ``bin/`` to ``PATH`` and ``lib/`` to
  ``LD_LIBRARY_PATH`` (Linux) or ``DYLD_LIBRARY_PATH`` (macOS). On Windows, add
  both ``bin/`` and ``lib/`` to ``PATH``. To build against the bundled Ipopt
  (e.g., for cyipopt), add ``lib/pkgconfig/`` to ``PKG_CONFIG_PATH``.
- 3.x: add the ``bin/`` directory to the same variables.

If you would like to use the IDAES binary distribution with Pyomo and have IDAES
installed the simplest way to set the appropriate paths is just to
``import idaes``.
