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
usually used by IDAES core models and tests.  If no solver is specified the
IDAES default (usually Ipopt) is used.  This is handy for method such as problem
initialization, where a default solver and options is usually used, but a user
can also choose a different solver or different options.

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

Test Models
-----------

The ``idaes.core.solvers.features`` provides functions to return simple models
of various types.  These models can be used to test if solvers are available and
functioning properly.  They can also be used to test that various option solver
features are available.  These functions all return a tuple where the first
element is a model of the specified type, and the second element is the correct
value for the x variable after the problem is solved.

.. autofunction:: idaes.core.solvers.features.lp

.. autofunction:: idaes.core.solvers.features.milp

.. autofunction:: idaes.core.solvers.features.nlp

.. autofunction:: idaes.core.solvers.features.minlp
