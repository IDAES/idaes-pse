Scaling
=======

.. module:: idaes.core.util.scaling

This section describes scaling utility functions and methods.

Context
-------
Creating well scaled models is important for increasing the efficiency and 
reliability of solvers. Since SI units are the standard for IDAES, several 
variables and constraints would be considered badly scaled (values less than 
:math:`10^-3` and greater than :math:`10^3`).

Scaling factors can be specified for any variable or constraint. Pyomo and many
solvers support the ``scaling_factor`` suffix. IDAES, as described below, also
supports the ``scaling_expression`` suffix which can be used to calculate
``scaling_factor`` values (e.g. based on the scaling factors of other variables).

To eliminate the possibility of defining conflicting scaling factors in various
places in the model, the IDAES standard is to define the ``scaling_factor`` and
``scaling_expression`` suffixes in the same block as the variable or constraint
that they are scaling.  This ensures that each scale factor is defined in only
one place, and is organized based on the model block structure.

Scaling factors in IDAES (and Pyomo) are multiplied by the variable or constraint
they scale.  For example, a Pressure variable in Pa units may be expected to have
a magnitude of around :math:`1\times10^6` for a specific process.  To scale the
variable to a more reasonable magnitude the scale factor for the variable could
be defined to be :math:`1\times10^-5`.

Specifying Scaling
------------------
Suffixes are used to specify scaling factors for IDAES models.

To supply variable and constraint scaling factors, an export suffix called 
``scaling_factor``should be created in the same block as the variable or constraint. 
Additionally, if the scaling factor for variables or constraints will be based on
other variables, a local suffix called ``scaling_expression`` should be created. 

Implementing Scaling
--------------------
While some solvers, such as Ipopt, support scaling factors, Pyomo also supplies scaling
transformations for models when solver scaling is not supported.

Neither Pyomo or other solvers use the IDAES ``scaling_expression``, so the scaling 
expressions must be converted to scaling factors with the 
``calculate_scaling_factors(m, basis)`` function. This function replaces the variables in 
the scaling expression with the specified basis values, calculates the scaling factors,
and puts the scaling factor in the ``scaling factor`` suffix.

Specifically for Ipopt, it is recommended to use the ``scale_constraints(m, basis)`` function 
to scale the constraints before sending the model to the solver. For the variables, the user
must set the ``nlp_scaling_method`` option to "user-scaling".

Example
-------
.. testcode::

    from pyomo.environ import Suffix, ConcreteModel, Var, Constraint
    from idaes.core.util.scaling import (
        ScalingBasis,
        calculate_scaling_factors,
    )

    m = ConcreteModel()
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)
    m.scaling_expression = Suffix(direction=Suffix.LOCAL)

    m.x = Var(initialize=1e6)
    m.y = Var(initialize=1e6)
    m.z = Var(initialize=1e12)

    m.scaling_factor[m.x] = 1e-5
    m.scaling_factor[m.y] = 1e-5
    m.scaling_expression[m.z] = 1/(m.x*m.y)

    m.c = Constraint(expr=m.z == m.x*m.y)
    m.scaling_expression[m.c] = 1/(m.x*m.y)

    calculate_scaling_factors(m, basis=ScalingBasis.InverseVarScale)

    # Show that the constraint scaling factor is 1/((1/1e-5)*(1/1e-5))
    assert(m.scaling_factor[m.c] - 1e-10 < 1e-12)
    # Show that the z variable scaling factor is 1/((1/1e-5)*(1/1e-5))
    assert(m.scaling_factor[m.z] - 1e-10 < 1e-12)


Scaling expression basis
------------------------
The general guideline for scaling expressions is that a scaling factor is
calculated based on the expected magnitude of the variables. The magnitude
could be estimated in different ways, but the IDAES standard is the inverse 
variable scale. The list below shows variable scaling bases that are provided.

ScalingBasis.InverseVarScale:
  Use the inverse variable scaling factors in scaling expressions.
ScalingBasis.Value:
  Use the current variable values in scaling expressions.
ScalingBasis.Mid:
  Use the mid-point between the upper and lower bounds in scaling expressions.
ScalingBasis.Lower:
  Use the lower bound of variables in scaling expressions.
ScalingBasis.Upper:
  Use the lower bound of variables in scaling expressions.
ScalingBasis.VarScale:
  This is less common, but it uses the variable scales directly. This can be
  used if you are using alternative scaling methods with divide by the scaling
  factor.


Scaling Utility Functions
-------------------------

IDAES includes some utility functions to help evaluate model scaling and to
auto-scale constraints.

.. autofunction:: calculate_scaling_factors

.. autofunction:: scale_constraints

.. autofunction:: badly_scaled_var_generator

.. autofunction:: grad_fd

.. autofunction:: constraint_fd_autoscale

.. autofunction:: set_scaling_factor


