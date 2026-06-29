Scaling Methods
===============

.. note::

  In v2.0, IDAES is beginning to deploy a new suite of scaling tools. This documentation refers to the older scaling tools.
  For documentation of the new scaling toolbox, :ref:`see here<reference_guides/scaling/scaling:Scaling Toolbox>`.

.. module:: idaes.core.util.scaling
  :noindex:

This section describes scaling utility functions and methods.

Context
-------
Creating well scaled models is important for increasing the efficiency and
reliability of solvers. Depending on property package units of measure and
process scale, variables and constraints are often badly scaled.

Scaling factors can be specified for any variable or constraint. Pyomo and many
solvers support the ``scaling_factor`` suffix. To eliminate the possibility of
defining conflicting scaling factors in various places in the model, the IDAES
standard is to define the ``scaling_factor`` suffixes in the same block as the
variable or constraint that they are scaling. This ensures that each scale
factor is defined in only one place, and is organized based on the model block
structure.

Scaling factors in IDAES (and Pyomo) are multiplied by the variable or constraint
they scale.  For example, a Pressure variable in Pa units may be expected to have
a magnitude of around :math:`10^6` for a specific process.  To scale the
variable to a more reasonable magnitude, the scale factor for the variable could
be defined to be :math:`1 \times 10^{-5}`.

While many scaling factors should be give good default values in the property
packages, some (e.g. flow rates or material holdups) must be given scale factors
by the user for a specific process model. Still other scale factors can be
calculated from supplied scale factors, for example, mass balance scale factors
could be determined from flow rate scale factors. To calculate scale factors,
models may have a standard ``calculate_scaling_factors()`` method.  For more
specific scaling information, see the model documentation.

For much of the core IDAES framework, model constraints are automatically scaled
via a simple transformation where both sides of the constraint are multiplied by
a scale factor determined based on supplied variable and expression scaling
factors. The goal of this is to ensure that solver tolerances are meaningful for
each constraint.  A constraint violation of :math:`1 \times 10^{-8}` should be
acceptable, but not too tight to achieve given machine precision limits.  IDAES
model constraints should conform approximately to this guideline after the
``calculate_scaling_factors()`` method is executed.  Users should follow this
guideline for constraints they write.  The scaling of constraints for reasonable
residual tolerances is done as a constraint transformation independent of the
scaling factor suffix.  Scaling factors for constraints can still be set based
on other methods such as reducing very large Jacobian matrix entries.

Specifying Scaling
------------------
Suffixes are used to specify scaling factors for IDAES models. These suffixes
are created when needed by calling the ``set_scaling_factor()`` function. Using
the ``set_scaling_factor()``, ``get_scaling_factor()``, and
``unset_scaling_factor()`` eliminates the need to deal directly with scaling
suffixes, and ensures that scaling factors are stored in the IDAES standard
location.

.. autofunction:: set_scaling_factor
  :noindex:

.. autofunction:: get_scaling_factor
  :noindex:

.. autofunction:: unset_scaling_factor
  :noindex:


Constraint Transformation
-------------------------
As mentioned previously, constraints in the IDAES framework are transformed such
that :math:`1 \times 10^{-8}` is a reasonable criteria for convergence before any
other scaling factors are applied. There are a few utility functions for scaling
transformation of constraints. When transforming constraints with these functions,
the scaling applies to the original constraint, not combined with any previous
transformation.

.. autofunction:: constraint_scaling_transform

.. autofunction:: constraint_scaling_transform_undo

.. autofunction:: get_constraint_transform_applied_scaling_factor


Calculation in Model
~~~~~~~~~~~~~~~~~~~~

Some scaling factors may also be calculated by a call to a model's
``calculate_scaling_factors()`` method.  For more information see specific model
documentation.

Sometimes a scaling factor may be set on an indexed component and prorogated to
it's data objects later can be useful for example in models that use the DAE
transformation, not all data objects exist until after the transformation.

.. autofunction:: propagate_indexed_component_scaling_factors

Constraint Auto-Scaling
~~~~~~~~~~~~~~~~~~~~~~~

Constraints can be scaled to automatically reduce very large entries in the Jacobian
matrix with the ``constraint_autoscale_large_jac()`` function.

.. autofunction:: constraint_autoscale_large_jac


Inspect Scaling
---------------

Models can be large, so it is often difficult to identify where scaling is needed
and where the problem may be poorly scaled.  The functions below may be helpful
in inspecting a models scaling. Additionally ``constraint_autoscale_large_jac()``
described above can provide Jacobian information at the current variable values.


.. autofunction:: extreme_jacobian_columns
  :noindex:

.. autofunction:: extreme_jacobian_rows
  :noindex:

.. autofunction:: badly_scaled_var_generator
  :noindex:

.. autofunction:: unscaled_variables_generator
  :noindex:

.. autofunction:: unscaled_constraints_generator
  :noindex:

.. autofunction:: map_scaling_factor

.. autofunction:: min_scaling_factor

.. autofunction:: get_jacobian

.. autofunction:: jacobian_cond
  :noindex:

Applying Scaling
----------------

Scale factor suffixes can be passed directly to a solver.  How the scale factors
are used may vary by solver. Pyomo also contains tools to transform a problem to
a scaled version.

Ipopt is the standard solver in IDAES.  To use scale factors with Ipopt, the
``nlp_scaling_method`` option should be set to ``user-scaling``.  Be aware that
this deactivates any NLP automatic scaling.
