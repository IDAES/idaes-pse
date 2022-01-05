Solvers
=======

This section provides information on using and configuring solvers for IDAES.
Some IDAES solver features are documented in other sections, so references are
provided as appropriate. In general, standard Pyomo solver interfaces and
features can be used in IDAES, but IDAES provides a few extensions to make
working with solvers slightly easier.

Default Solver Config
---------------------

The default solver settings can be set via the IDAES :ref:`IDAES configuration system<user_guide/configuration:Configuration>`.  By
default the solver settings are not taken from the IDAES configuration system
to avoid confusion where IDAES defaults differ from Pyomo defaults.  To use the
IDAES solver defaults a user must explicitly enable them with the
``use_idaes_solver_configuration_defaults()`` function.

.. autofunction:: idaes.core.solvers.config.use_idaes_solver_configuration_defaults

Getting a Solver
----------------

A :ref:`utility function <technical_specs/core/util/misc:get_solver>` (``idaes.core.util.misc.get_solver``) provides
functionality to get a default or user specified solver at run time.  This is
usually used by IDAES core models and tests.  If no solver is specified, the
IDAES default (usually Ipopt) is used.  This function is handy for methods such
as model initialization, where a default solver and options is usually used, but
a user can specify a different solver or options.

Solver Logging
--------------

A logger for solver-related log messages can be obtained from the
:ref:`IDAES logging system<user_guide/logging:idaes.solve Loggers>`.

IDAES has features for redirecting solver output to a log.  This is documented
in the :ref:`logging section<user_guide/logging:Logging Solver Output>`.


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

Installing PETSc
~~~~~~~~~~~~~~~~

The PETSc solver is an extra binary package, and not installed by default.  If
you are using a supported Linux distribution, you can use the command
``idaes get-extensions --extra petsc`` to install it.

There is no precompiled PETSc solver for Windows, but here are two options for
Windows installation. The easiest option is to run the available precompiled
Linux version via the WSL
(see :ref:`binary installation <getting_started/binaries:Using the WSL>`
for details). Expert users may wish to compile their own solver. Source code is
available in the `idaes-ext repo <https://github.com/IDAES/idaes-ext/tree/main/petsc>`_
If you can compile PETSc for Windows, compiling the interface is trivial (see
`PETSc's windows installation documentation <https://petsc.org/main/install/windows/>`_).

Registered Solvers
~~~~~~~~~~~~~~~~~~

Importing ``idaes.core.solvers.petsc`` registers two new solvers "petsc_snes"
and "petsc_ts."  The "petsc_snes" solver provides nonlinear equation solvers.
The SNES (Scalable Non-Linear Equation Solvers) solvers are strictly nonlinear
equations solvers, so they cannot directly handle optimization problems and the
problem must have zero degrees of freedom. The TS (time stepping) solvers require
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

To get started, import command line options are for SNES solvers are described
here and TS options are described here. Remember that options specified through
the IDAES AMPL interface use a double dash rather than the single dash shown in
the PETSc documentation.

Utilities for DAEs with Pyomo.DAE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to use the "petsc_ts" solver is to use the utility method that
converts a standard Pyomo.DAE to the form used by the solver.

Discretization
""""""""""""""

The utility for solving Pyomo.DAE problems uses the PETSc TS solvers to integrate
between each time point in the Pyomo.DAE discretization.  The results are stored
for each time point. This can be used to initialize and verify the results of the
full time-discretized model.

To quickly run a DAE model, a time discretization with one element can be made.
In this case the PETSc TS solver will integrate from the initial condition to the
end point. Results for intermediate times can be read from PETSc's stored
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
calculate initial conditions, but cannot have constraints that specify a final
or intermediate conditions.  Optimization is not directly possible, but future
implementation of optimization solvers in combination with adjoint sensitivity
calculations may enable optimization.

Non-time-indexed variables and constraints should usually be solved with the
initial conditions in the first step.  Currently these must be explicitly listed
by the user, but in the future additional utilities may be provided to automate
identification of such constraints.

Solving
"""""""

The following function can be used to solve the DAE.

.. autofunction:: idaes.core.solvers.petsc.petsc_dae_by_time_element

Reading Trajectory Data
"""""""""""""""""""""""

Usually if you want to read the trajectory data from the solver, you will want to
solve the whole time domain, so you will want to specify one time element in the
Pyomo.DAE discretization.  By specifying the ``--ts_save_trajectory=1`` and
``--ts_trajectory_type=visualization`` options the trajectory information will
be saved.  Supplying the ``vars_stub`` argument to the ``petsc_dae_by_time_element()``
function will copy the `*.col` and `*.typ` files needed to interpret the trajectory
data to the current working directory.

The following example shows how to read and plot some data. A more streamlined
method of providing the PETSc Python functions is forthcoming.

Using TS Solvers without Pyomo.DAE
""""""""""""""""""""""""""""""""""

Most IDAES models use Pyomo.DAE and that is probably the easiest way to set up
a DAE problem, however you may directly construct a DAE problem.

There are two suffixes that need to be specified to use the PETSc TS solvers from
Pyomo.  The first is an integer suffix ``dae_suffix``, which specifies the variables
types.  The algebraic variables do not need to be included, but 0 specifies
algebraic variables, 1 specifies differential variables, 2 specifies derivative
variables, and 3 specifies the time variable.  A variable for time does not need
to be included.  The other suffix is an integer suffix ``dae_link`` which contains
differential and derivative variables.  The integer in the suffix links the
derivative to it's differential variable, by specifying an integer that is unique
to the pair.

To solve the problem, start with the initial conditions in the variable values.
After the solve the finial conditions will be in the Pyomo model.  To get
intermediate results, you will need to store the solver trajectory as described
previously.

.. code-block:: python

  """
  Visualize trajectory, assuming problem was run with these options:

      --ts_save_trajectory=1
      --ts_trajectory_type=visualization

  That should put results from each time step in the directory Visualization-data

  To run this you need either:

  Option 1) $PETSC_DIR/lib/petsc/bin/ in PYTHONPATH

  Option 2) Get these files and put them in PYTHONPATH
    * https://github.com/petsc/petsc/blob/main/lib/petsc/bin/petsc_conf.py
    * https://github.com/petsc/petsc/blob/main/lib/petsc/bin/PetscBinaryIO.py
    * https://github.com/petsc/petsc/blob/main/lib/petsc/bin/PetscBinaryIOTrajectory.py
  """
  import sys
  import PetscBinaryIOTrajectory as pbt

  def main():
      try:
          stub = sys.argv[1]
      except:
          print("Specify a file stub")
          return
      (t,v,names) = pbt.ReadTrajectory("Visualization-data")
      with open(f'{stub}.col') as f:
          names = list(map(str.strip, f.readlines()))
      with open(f'{stub}.typ') as f:
          typ = list(map(int,f.readlines()))
      names = [name for i, name in enumerate(names) if typ[i] in [0,1]]
      pbt.PlotTrajectories(t,v,names,["my_variable_name"])

  if __name__ == '__main__':
      main()

Planned Future PETSc Support
""""""""""""""""""""""""""""

This section provides PETSc features that are planned to be supported in the
future, but are not currently supported.

* Enable parallel methods
* Enable IMEX methods for TS solvers
* Enable TAO optimization solvers

Test Models
-----------

The ``idaes.core.solvers.features`` module provides functions to return simple
models of various types.  These models can be used to test if solvers are
available and functioning properly.  They can also be used to test that various
optional solver features are available.  These functions all return a tuple where
the first element is a model of the specified type, and the remaining element is
the correct solved value for select variables.

.. autofunction:: idaes.core.solvers.features.lp

.. autofunction:: idaes.core.solvers.features.milp

.. autofunction:: idaes.core.solvers.features.nlp

.. autofunction:: idaes.core.solvers.features.minlp

.. autofunction:: idaes.core.solvers.features.nle

.. autofunction:: idaes.core.solvers.features.dae

.. autofunction:: idaes.core.solvers.features.dae_with_non_time_indexed_constraint
