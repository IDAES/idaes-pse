Scaling
=======

.. module:: idaes.core.util.scaling

This section describes scaling utility functions and methods.

Standard Practice
-----------------

Scaling factors can be specified for any variable or constraint. Pyomo and many
solvers support the ``scaling_factor`` suffix. IDAES, as described below, also
supports the ``scaling_expression`` suffix which can be used to calculate
``scaling_factor`` values (e.g. based on state block units of measure).

To eliminate the possibility of defining, conflicting scaling factors in various
places in the model the IDAES standard is to define the ``scaling_factor`` and
``scaling_expression`` suffixes in the same block as the variable or constraint
that they are scaling.  This ensures that each scale factor is defined in only
one place, and is organized based on the model block structure.

Scaling factors in IDAES (and Pyomo) are multiplied by the variable or constraint
they scale.  For example, a Pressure variable in Pa units may be expected to have
a magnitude of around :math:`1\times10^6` for a specific process.  To scale the
variable to a more reasonable magnitude the scale factor for the variable could
be defined to be :math:`1\times10^-5`.


Specifying Variable Scaling
---------------------------

Suffixes are used to specify scaling factors for IDAES models.  Some solvers, such
as Ipopt, support supplying scale factors. Pyomo also supplies scaling
transformations for models when solver scaling is not supported.

To supply variable and constraint scaling factors, a suffix called ``scaling_factor``
should be created in the same block as the variable or constraint. For example:

.. testcode::

    from pyomo.environ import Suffix, ConcreteModel, Var

    m = ConcreteModel()
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)
    m.P = Var(initialize=1e6, doc="Pressure [Pa]")
    m.conc = Var(["Na+", "Cl-"], initialize=1e-4)
    m.scaling_factor[m.P] = 1e-5
    m.scaling_factor[m.conc["Na+"]] = 1e3
    m.scaling_factor[m.conc["Cl-"]] = 1e3

Variable scaling in state blocks is provided by the developer of a state
block and can be used as a basis for scaling other model variables and
constraints. Scaling factors can be modified by users to better
represent the process they are modeling.

Specifying Scaling Factor Expressions
-------------------------------------

Scaling factors for variables and constraints can be calculated based on variable
scaling factors, bounds, or values that have been provided.  The calculation for
a scaling factor can be provided as a python expression using model variables in
the ``scaling_expression`` suffix. For variables, generally the expression should
only depend on variables where scaling factors have been defined.

The ``calculate_scaling_factors(m, basis)`` function replaces the variables in the scaling
expression with the specified basis values, calculates the scaling factors, and
puts the scaling factor in the ``scaling_factor`` suffix.


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

In the scaling expression the general guideline is that a scaling factor is being
calculated based on the expected magnitude of the variable values.  The magnitude
could be estimated in different ways, but the most common way should be the inverse
variable scale.  The list below shows variable scaling bases that are provided.

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


.. autofunction:: calculate_scaling_factors


Scaling with Ipopt
------------------

To use the supplied scaling factors with Ipopt the ``nlp_scaling_method`` solver
option should be set to "user-scaling."
