Binary Packages
===============

This section describes what is in the IDAES binary distribution, where it is
installed, how to use it, and alternative installation methods.  The IDAES binary
distribution contains mainly third-party solvers compiled for user convenience and
function libraries used for some IDAES physical property packages.

Installation is generally done through the
:doc:`idaes get-extensions command<../user_guide/commands/get_extensions>`.

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
- `bonmin, <https://petsc.org/release/>`_ (uses Ipopt with HSL)
- `couenne <https://projects.coin-or.org/Couenne/>`_ (uses Ipopt with HSL)
- iapws95_external, function for the IAPWS95 EoS water
- swco2_external, functions for the Span-Wagner EoS for CO2
- cubic_roots, function for finding root of a cubic equation
- functions, miscellaneous simple math functions, examples, and test code

Extras
~~~~~~

- petsc, `AMPL solver library <https://ampl.com/REFS/hooking2.pdf>`_ interface wrapper for `PETSc's <https://petsc.org/release/>`_ SNES and TS solvers

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

Using the WSL
-------------

The Linux binary builds can be run on Windows using the
`WSL <https://docs.microsoft.com/en-us/windows/wsl/about>`_.  This is useful
when solvers depend on POSIX standards.  Currently the optional PETSc solver is
the only one which does not run on Windows.

In the WSL environment, download and extract the desired Linux release of
the solver (i.e. PETSc) you would like to run following the instructions under manual
installation above. Install required libraries libfortran, libblas, liblapack i.e.
``apt-get install liblapack3 libblas3 ligbomp1 libgfortran5``.

In the WSL environment you can check that PETSc is installed correctly. In your 
desired Linux release navigate to the folder containint the compiled PETSc binary
and run the command ``./petsc --version``. If setup correctly details of the PETSc
release will be printed on the screen.

In the Windows IDAES binary directory, create a
batch file with the format ``{solver}_wsl.bat`` for example for PETSc
``petsc_wsl.bat``. Assuming you have put the binary file in the ``$HOME/local/bin``
and your WSL user name is ``john``, and the distribution installed is ``Ubuntu-20.04``
the contents of the batch file to run the petsc are below:

.. code-block ::

  @echo off
  idaes solver-wsl --distribution Ubuntu-20.04 --user john --executable ~/local/bin/petsc %*

If you are unsure of the name of your installed distribution, navigate to your linux
release and run the command ``wsl --list --all``. This will list all your linux
distributions (if you have more than one) with their names.

There are a few options to check that your ``{solver}_wsl.bat`` is installed
correctly. In a power shell window run the command
``$Env:LOCALAPPDATA\idaes\bin\petsc_wsl.bat --version``. If setup correctly details
of the solver release version will be printed on the screen. For users of git bash,
open a git bash window and navigate to the environment (if you used environments)
where the solver i.e. PETSc was installed. Run the command
``$LOCALAPPDATA/idaes/bin/petsc_wsl.bat --version``.  If setup correctly details of
the PETSc release version will be printed on the screen.