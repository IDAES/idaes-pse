Solvers
=======

This section provides an overview of using and configuring solvers for IDAES.
In general, standard Pyomo solver interfaces and features are used in IDAES,
but IDAES provides a few extensions to make working with solvers slightly easier.
Some IDAES solver features are documented in other sections, so references are
provided as appropriate.

Default Solver Config
---------------------

The global solver settings can be set via the
:ref:`IDAES configuration system<reference_guides/configuration:Configuration>`.
This feature is handy in IDAES where multiple solver objects are used for
initialization before finally solving a problem. Since IDAES default solver
settings differ from Pyomo, users must explicitly enable the IDAES solver
configuration system with the ``use_idaes_solver_configuration_defaults()``
function.

.. autofunction:: idaes.core.solvers.config.use_idaes_solver_configuration_defaults

Getting a Solver
----------------

Typically users can use the standard Pyomo SoverFactory to get a solver.  If a
solver is needed in a general model or utility, a
:ref:`utility function <reference_guides/core/util/misc:get_solver>` (``idaes.core.util.misc.get_solver``)
provides a default or user configured solver at runtime. This is used by IDAES
core models and tests.

Solver Logging
--------------

A logger for solver-related log messages can be obtained from the
``idaes.logger.getSolveLogger()`` function
(:ref:`documented here<reference_guides/logging:idaes.solve Loggers>`).
IDAES also has features for redirecting solver output to a log (see
:ref:`Logging Solver Output<reference_guides/logging:Logging Solver Output>`).


Solver Feature Checking
-----------------------

There are some functions available to check what features are available to
solvers and to help with basic solver testing.

.. autofunction:: idaes.core.solvers.ipopt_has_linear_solver


PETSc Utilities
---------------

IDAES provides an AMPL solver interface for the PETSc solver suite,
`(see the PETSc website) <https://petsc.org/release/overview/>`_.  PETSc provides
nonlinear equation (NLE) and differential algebraic equation (DAE) solvers. Both
NLE and DAE solvers are capable of solving simulation problems with zero degrees
of freedom. These solvers may be useful for initial model development,
initialization, and running simulation cases without optimization.

PETSc includes optimization solvers, but they are not currently supported by the
IDAES AMPL solver wrapper. Optimization support will likely be added in the
future.

DAE Terminology
~~~~~~~~~~~~~~~

For the following discussion regarding the PETSc solver interface, the following
terminology is used.

* Derivative variable: a time derivative
* Differential variable: a variable that is differentiated with respect to time
* Algebraic variable: a variable with no explicit time derivative appearing in the problem
* State variables: the set of algebraic and differential variables
* Time variable: a variable representing time

DAE problems do not need to include a time variable, but, if they do, there can
only be one. Differential variables do not need to explicitly appear in
constraints, but their time derivatives do.  DAE problems must have zero degrees
of freedom, which means the number of constraints must equal the number of state
variables.

Installing PETSc
~~~~~~~~~~~~~~~~

The PETSc solver is an extra binary package, and not installed by default.  If
you are using a supported Linux distribution, you can use the command
``idaes get-extensions --extra petsc`` to install it.

There is no precompiled PETSc solver for Windows, but here are two options for
Windows installation. The easiest option is to run the available precompiled
Linux version via the WSL
(see :ref:`binary installation <tutorials/getting_started/binaries:Using the WSL>`
for details). Expert users may wish to compile their own solver. Source code is
available in the `idaes-ext repo <https://github.com/IDAES/idaes-ext/tree/main/petsc>`_.
If you can compile PETSc for Windows, compiling the interface is trivial (see
`PETSc's windows installation documentation <https://petsc.org/main/install/windows/>`_).

The IDAES PETSc package also includes Python modules for reading binary data
written by the PETSc solver.  On Windows, some manual installation of the Python
modules is required.  If you are using the WSL method to run PETSc, copy the
``petscpy`` directory from the Linux package you are using to the IDAES binary
directory.  You can find the IDAES binary directory by running the command
``idaes bin-directory`` in the OS command shell (e.g. Bash, Windows CMD,
PowerShell).  If IDAES is installed in a Python environment, the environment
must be active. The primary use for these Python modules is to read trajectory
files saved by the TS solver.

Registered Solvers
~~~~~~~~~~~~~~~~~~

Importing ``idaes.core.solvers.petsc`` registers two new solvers "petsc_snes"
and "petsc_ts."  The "petsc_snes" solver provides nonlinear equation solvers.
The SNES (Scalable Nonlinear Equation Solvers) solvers are strictly nonlinear
equations solvers, so they cannot directly handle optimization problems and the
problem must have zero degrees of freedom. The TS (time-stepping) solvers require
specialized suffixes to designate derivative, differential, algebraic, and time
variables and the associations between derivative and differential variables.
Both the SNES and TS solvers accept the standard scaling factor suffixes, but
for TS solvers, derivatives and differential variables cannot be scaled
independently, so the differential variable scale is used.  Currently, time
cannot be scaled for TS solvers.

Standard PETSc command line options are available to the solvers except that for
compatibility reasons, they are specified with a double dash instead of single.
Command line options can be used to set up the SNES and TS solvers. Currently
only implicit TS solvers are supported.  Commonly used TS types are:

* "beuler", implicit Euler,
* "cn", Crank-Nicolson, and
* "alpha", generalized-alpha method.

To get started, important command line options are for SNES solvers are described
`<https://petsc.org/release/docs/manualpages/SNES/SNESSetFromOptions.html>`_
and TS options are described
`<https://petsc.org/release/docs/manualpages/TS/TSSetFromOptions.html>`_.
Remember that options specified through the IDAES AMPL interface use a double
dash rather than the single dash shown in the PETSc documentation. Users can
also set linear solver and preconditioner options, and are encouraged to
read the PETSc documentation if needed.

Utilities for DAEs with Pyomo.DAE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to use the "petsc_ts" solver is to use the utility method that
converts a standard Pyomo.DAE to the form used by the solver.

Discretization
""""""""""""""

The utility for solving Pyomo.DAE problems uses the PETSc TS solvers to integrate
between each time point in the Pyomo.DAE discretization.  The results are stored
for each time point in the Pyomo model. This can be used to initialize and verify
the results of the full time-discretized model. For example this could be used to
determine if the time steps used in the discretization are too big by comparing
the integrator solution to the fully discretized model solution.

To quickly run a DAE model, a time discretization with one element can be made.
In this case, the PETSc TS solver will integrate from the initial condition to
the end point. Results for intermediate times can be read from PETSc's stored
trajectory data if the proper solver options are specified.

Time Variable
"""""""""""""

Although it is probably not typical of Pyomo.DAE models, a time variable can be
specified.  Constraints can be written as explicit functions of time.  For
example, some model input could be ramped up or down as a function of time.

Limitations
"""""""""""

The integrator approach does not support some constraints that can be solved
using the full discretized model.  For example, you can have constraints to
calculate initial conditions, but cannot have constraints that specify final
or intermediate conditions.  Optimization is not directly possible, but future
implementation of optimization solvers in combination with adjoint sensitivity
calculations may enable optimization.

Non-time-indexed variables and constraints should usually be solved with the
initial conditions in the first step. Non-time-indexed variables can optionally
be detected and added to the equations solved for the initial conditions, or
explicitly specified by the user.  Users will have to take care not to include
non-time indexed constraints that contain time-indexed variables at times other
than the initial time.  If such constraints exist for the fully discretized
model users should deactivate them as appropriate.

Solving
"""""""

The following function can be used to solve the DAE.

.. autofunction:: idaes.core.solvers.petsc.petsc_dae_by_time_element

Reading Trajectory Data
"""""""""""""""""""""""

Usually if you want to read the trajectory data from the solver, you will want to
solve the whole time domain at once, so you will want to specify one time element
in the Pyomo.DAE discretization. By specifying the ``--ts_save_trajectory=1`` and
``--ts_trajectory_type=visualization`` options the trajectory information will
be saved.  Supplying the ``vars_stub`` argument to the ``petsc_dae_by_time_element()``
function will copy the `*.col` and `*.typ` files needed to interpret the trajectory
data to the current working directory.

The ``PetscTrajectory`` class has methods to read in trajectory data and
interpolate time points as needed.

.. autoclass:: idaes.core.solvers.petsc.PetscTrajectory
    :members:

Using TS Solvers without Pyomo.DAE
""""""""""""""""""""""""""""""""""

Most IDAES models use Pyomo.DAE and that is probably the easiest way to set up
a DAE problem, however you may directly construct a DAE problem.

There are two suffixes that need to be specified to use the PETSc TS solvers from
Pyomo.  The first is an integer suffix ``dae_suffix``, which specifies the variable
types.  The algebraic variables do not need to be included, but 0 specifies
algebraic variables, 1 specifies differential variables, 2 specifies derivative
variables, and 3 specifies the time variable.  A variable for time is optional,
and only one time variable can be specified. The other suffix is an integer
suffix ``dae_link`` which contains differential and derivative variables. The
integer in the suffix links the derivative to it's differential variable, by
specifying an integer greater than 0 that is unique to the pair.

If there are differential variables that do not appear in the constraints,
they can be supplied to the ``export_nonlinear_variables`` argument of solve.
For the trajectory data, you will also want to use ``symbolic_solver_labels``.

To solve the problem, start with the initial conditions in the Pyomo model.
After the solve the final conditions will be in the Pyomo model.  To get
intermediate results, you will need to store the solver trajectory as described
previously.

Planned Future PETSc Support
""""""""""""""""""""""""""""

This section provides PETSc features that are planned to be supported in the
future, but are not currently supported.

* Enable parallel methods
* Enable IMEX methods for TS solvers
* Enable TAO optimization solvers
* Provide PETSc Python functions for reading trajectory data (rather than
  requiring users to get them manually).

Test Models
-----------

The ``idaes.core.solvers.features`` module provides functions to return simple
models of various types.  These models can be used to test if solvers are
available and functioning properly.  They can also be used to test that various
optional solver features are available.  These functions all return a tuple where
the first element is a model of the specified type, and the remaining elements are
the correct solved values for select variables.

.. autofunction:: idaes.core.solvers.features.lp

.. autofunction:: idaes.core.solvers.features.milp

.. autofunction:: idaes.core.solvers.features.nlp

.. autofunction:: idaes.core.solvers.features.minlp

.. autofunction:: idaes.core.solvers.features.nle

.. autofunction:: idaes.core.solvers.features.dae
