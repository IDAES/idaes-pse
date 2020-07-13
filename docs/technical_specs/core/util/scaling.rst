Scaling Methods
===============

.. module:: idaes.core.util.scaling

This section describes scaling utility functions and methods.

Context
-------
Creating well scaled models is important for increasing the efficiency and
reliability of solvers. Since the standard units for IDAES are SI, oftentimes
variables and constraints are considered badly scaled (values less than
:math:`10^{-3}` and greater than :math:`10^3`).

Scaling factors can be specified for any variable or constraint. Pyomo and many
solvers support the ``scaling_factor`` suffix. To eliminate the possibility of
defining conflicting scaling factors in various places in the model, the IDAES
standard is to define the ``scaling_factor`` suffixes in the same block as the
variable or constraint that they are scaling.  This ensures that each scale
factor is defined in only one place, and is organized based on the model block
structure.

Scaling factors in IDAES (and Pyomo) are multiplied by the variable or constraint
they scale.  For example, a Pressure variable in Pa units may be expected to have
a magnitude of around :math:`10^6` for a specific process.  To scale the
variable to a more reasonable magnitude, the scale factor for the variable could
be defined to be :math:`1 \times 10^{-5}`.

While many scaling factors are given good default values in the property packages,
some (e.g. flow rates or material holdups) must be given scale factors by the
user for a specific process model. Still other scale factors can be calculated
from supplied scale factors, for example, mass balance scale factors could be
determined from flow rate scale factors. To calculate scale factors, models may
have a standard ``calculate_scaling_factors()`` method.  For more specific scaling
information, see the model documentation.

Specifying Scaling
------------------
Suffixes are used to specify scaling factors for IDAES models. These suffixes
are created when needed by calling the ``set_scaling_factor()`` function. Using
the ``set_scaling_factor()``, ``get_scaling_factor()``, and
``unset_scaling_factor()`` eliminates the need to deal directly with scaling
suffixes, and ensures that scaling factors are stored in the IDAES standard
location.

.. autofunction:: set_scaling_factor

.. autofunction:: get_scaling_factor

.. autofunction:: unset_scaling_factor

Calculation in Model
~~~~~~~~~~~~~~~~~~~~

Some scaling factors may also be calculated by a call to a model's
``calculate_scaling_factors()`` method.  For more information see specific model
documentation.

Constraint Auto-Scaling
~~~~~~~~~~~~~~~~~~~~~~~

Constraints can be scaled to automatically reduce very large entries in the Jacobian
matrix with the ``constraint_autoscale_large_jac()`` function.

.. autofunction:: constraint_autoscale_large_jac


Constraint Transformation
-------------------------
Depending on the solver various tolerances may be applied to either the scaled
or unscaled problem.  It can also be desirable to perform constraint scaling in
two steps.  First scale the constraint so that a constraint violation tolerance
of :math:`1 \times 10^{-8}` is a reasonable criteria for convergence. Then
transform the constraint.  This will ensure that various solver tolerances make
sense, and ensure that things like selection of units of measure don't affect
the problem solution. Other constraint scaling methods such as Jacobian based
scaling can be applied to the transformed constraints.

To transform constraints, the following functions can be used.  These functions
are often used in general unit models to scale constraints such as mass an energy
balances where proper constraint scaling must be determined based on process size
and units of measure.  Scaling factors supplied are usually based on some known
variable scale factors.

.. autofunction:: constraint_scaling_transform

.. autofunction:: constraint_scaling_transform_undo

.. autofunction:: get_constarint_tranform_applied_scaling_factor


Inspect Scaling
---------------

Models can be large, so it is often difficult to identify where scaling is needed
and where the problem may be poorly scaled.  The functions below may be helpful
in inspecting a models scaling. Additionally ``constraint_autoscale_large_jac()``
described above can provide Jacobian information at the current variable values.


.. autofunction:: badly_scaled_var_generator

.. autofunction:: unscaled_variables_generator

.. autofunction:: unscaled_constraints_generator


Applying Scaling
----------------

Scale factor suffixes can be passed directly to a solver.  How the scale factors
are used may vary by solver. Pyomo also contains tools to transform a problem to
a scaled version.

Ipopt is the standard solver in IDAES.  To use scale factors the
``nlp_scaling_method`` option should be set to ``user-scaling``.  Be aware that
this deactivates any NLP automatic scaling.
