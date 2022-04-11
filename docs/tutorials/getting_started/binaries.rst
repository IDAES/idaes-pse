Binary Packages
===============

This section describes what is in the IDAES binary distribution, where it is
installed, how to use it, and alternative installation methods.  The IDAES binary
distribution contains mainly third-party solvers compiled for user convenience and
function libraries used for some IDAES physical property packages.

Installation is generally done through the
:doc:`idaes get-extensions command<../../reference_guides/commands/get_extensions>`.

What is Included
----------------

The National Energy Technology Laboratory (NETL) has obtained a distribution
license for HSL linear solvers that allows inclusion in IDAES's compiled version
of Ipopt.

All technical papers, sales and publicity material resulting from use of
the HSL codes within IPOPT must contain the following acknowledgement: HSL, a
collection of Fortran codes for large-scale scientific computation. See
http://www.hsl.rl.ac.uk.

Main Package
~~~~~~~~~~~~

- `ipopt <https://coin-or.github.io/Ipopt/>`_ with `HSL linear solvers <http://www.hsl.rl.ac.uk>`_.
- `ipopt_sens <https://projects.coin-or.org/Ipopt/wiki/sIpopt>`_, sIpopt solver with `HSL linear solvers <http://www.hsl.rl.ac.uk>`_.
- `ipopt_l1 <https://github.com/IDAES/Ipopt/tree/restoration_mod>`_, Ipopt with a modified restoration phase algorithm with `HSL linear solvers <http://www.hsl.rl.ac.uk>`_.
- shared ipopt library intended for `cyipopt <https://cyipopt.readthedocs.io/en/stable/>`_ with `HSL linear solvers <http://www.hsl.rl.ac.uk>`_.
- `dot_sens <https://github.com/dthierry/k_aug>`_
- `k_aug <https://github.com/dthierry/k_aug>`_
- `pynumero <https://pyomo.readthedocs.io/en/stable/contributed_packages/pynumero/index.html>`_
- `clp <https://projects.coin-or.org/Clp>`_
- `cbc <https://projects.coin-or.org/Cbc>`_
- `bonmin, <https://projects.coin-or.org/Bonmin>`_ (uses Ipopt with HSL)
- `couenne <https://projects.coin-or.org/Couenne/>`_ (uses Ipopt with HSL)
- `iapws95_external <https://github.com/IDAES/idaes-ext/tree/main/src/helmholtz>`_, function for the IAPWS95 EoS water
- `swco2_external <https://github.com/IDAES/idaes-ext/tree/main/src/helmholtz>`_, functions for the Span-Wagner EoS for CO2
- `cubic_roots <https://github.com/IDAES/idaes-ext/tree/main/src/cubic>`_, function for finding root of a cubic equation
- `functions <https://github.com/IDAES/idaes-ext/tree/main/src/functions>`_, miscellaneous simple math functions, examples, and test code

Extras
~~~~~~

- petsc, `AMPL solver library <https://ampl.com/REFS/hooking2.pdf>`_ `interface wrapper <https://github.com/IDAES/idaes-ext/tree/main/petsc>`_ for `PETSc's <https://petsc.org/release/>`_ SNES and TS solvers.

Supported Systems
-----------------

Currently builds are available for Windows and several Linux distributions.
Although the binaries are complied on a small number of platforms, one of the
available sets should work on most recent Linux distributions. The table below
show distributions that should work automatically. If you have a distribution
not on the list, you can try to specify a similar distribution.

+-----------------------------+---------------+--------------------+-------------------+
| Distro/OS                   | ID            | Arch               | Build Distro/OS   |
+=============================+===============+====================+===================+
| AlmaLinux 8                 | almalinux7    | x86_64             | el8               |
+-----------------------------+---------------+--------------------+-------------------+
| CentOS 7                    | centos7       | x86_64             | el7               |
+-----------------------------+---------------+--------------------+-------------------+
| CentOS 8                    | centos8       | x86_64             | el8               |
+-----------------------------+---------------+--------------------+-------------------+
| Debian 9                    | debian9       | x86_64             | el7               |
+-----------------------------+---------------+--------------------+-------------------+
| Debian 10                   | debian10      | x86_64             | el8               |
+-----------------------------+---------------+--------------------+-------------------+
| Debian 11                   | debian11      | x86_64, aarch64    | ubuntu2004        |
+-----------------------------+---------------+--------------------+-------------------+
| Kubuntu 18.04               | kubuntu1804   | x86_64             | ubuntu1804        |
+-----------------------------+---------------+--------------------+-------------------+
| Kubuntu 20.04               | kubuntu2004   | x86_64, aarch64    | ubuntu2004        |
+-----------------------------+---------------+--------------------+-------------------+
| Kubuntu 22.04               | kubuntu2204   | x86_64             | ubuntu2204        |
+-----------------------------+---------------+--------------------+-------------------+
| Linux Mint 20               | ubuntu2004    | x86_64, aarch64    | ubuntu2004        |
+-----------------------------+---------------+--------------------+-------------------+
| Red Hat Enterprise Linux 7  | rhel7         | x86_64             | el7               |
+-----------------------------+---------------+--------------------+-------------------+
| Red Hat Enterprise Linux 8  | rhel8         | x86_64             | el8               |
+-----------------------------+---------------+--------------------+-------------------+
| Rocky Linux 8               | rocky8        | x86_64             | el8               |
+-----------------------------+---------------+--------------------+-------------------+
| Scientific Linux 7          | scientific7   | x86_64             | el7               |
+-----------------------------+---------------+--------------------+-------------------+
| Ubuntu 18.04                | ubuntu1804    | x86_64             | ubuntu1804        |
+-----------------------------+---------------+--------------------+-------------------+
| Ubuntu 20.04                | ubuntu2004    | x86_64, aarch64    | ubuntu2004        |
+-----------------------------+---------------+--------------------+-------------------+
| Ubuntu 22.04                | ubuntu2204    | x86_64             | ubuntu2204        |
+-----------------------------+---------------+--------------------+-------------------+
| Windows 10, 11              | windows       | x86_64             | windows           |
+-----------------------------+---------------+--------------------+-------------------+
| Xubuntu 18.04               | xubuntu1804   | x86_64             | ubuntu1804        |
+-----------------------------+---------------+--------------------+-------------------+
| Xubuntu 20.04               | xubuntu2004   | x86_64, aarch64    | ubuntu2004        |
+-----------------------------+---------------+--------------------+-------------------+
| Xubuntu 22.04               | xubuntu2204   | x86_64             | ubuntu2204        |
+-----------------------------+---------------+--------------------+-------------------+

Install Location
----------------

The location of the binary file installation is not currently configurable and
it can be found with the command line command ``idaes bin-directory``.

Manual Installation
-------------------

The install location is not configurable, but the installation step can easily be
done manually. This is occasionally necessary when, for example, a firewall
blocks downloading the binary file from GitHub.

The first step for a manual install is to determine the location of the IDAES
binary directory, which can be done with the command ``idaes bin-directory``.

Download the IDAES binary release files from
`the binary releases page <https://github.com/IDAES/idaes-ext/releases>`_.
You will need ``idaes-lib-{platform}.tar.gz``, ``idaes-solvers-{platform}.tar.gz``,
and ``idaes-{optional package}-{platform}.tar.gz``. Extract the tar files in the
IDAES binary directory.

Using Solvers Without IDAES
---------------------------

Generally, the environment variables are set up to use the IDAES solvers when any
idaes module is imported in Python.  If you would like to use the solvers in
Pyomo without importing IDAES, you will need to add the install location to your
executable and dynamic library search path environment variables as appropriate
for your operating system (e.g. add the install location to $PATH in Linux).

If you would like to use the IDAES binary distribution with Pyomo and have IDAES
installed the simplest way to set the appropriate paths is just to
``import idaes``.
