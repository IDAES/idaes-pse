Scaling
=======
.. module:: idaes.core.util.scaling

This section describes scaling utility functions and methods.

Context
-------
Creating well scaled models is important for increasing the efficiency and 
reliability of solvers. Since the standard units for IDAES are SI, oftentimes 
variables and constraints are considered badly scaled (values less than 
:math:`10^{-3}` and greater than :math:`10^3`).

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
a magnitude of around :math:`10^6` for a specific process.  To scale the
variable to a more reasonable magnitude the scale factor for the variable could
be defined to be :math:`10^{-5}`.

Specifying Scaling
------------------
Suffixes are used to specify scaling factors for IDAES models.

To supply variable and constraint scaling factors, an export suffix called 
``scaling_factor`` should be created in the same block as the variable or constraint. 
Additionally, if the scaling factor for variables or constraints will be based on
other variables, a local suffix called ``scaling_expression`` should be created. 

Implementing Scaling
--------------------
While some solvers, such as Ipopt, support scaling factors, Pyomo also supplies scaling
transformations for models when solver scaling is not supported.

Neither Pyomo or other solvers use the IDAES ``scaling_expression``, so the scaling 
expressions must be converted to scaling factors with the 
``calculate_scaling_factors(m, basis)`` function. This function replaces the variables in 
the scaling expression with the specified basis value, calculates the scaling factors, 
and puts the scaling factor in the ``scaling_factor`` suffix.

Specifically, for Ipopt, it is recommended to use the ``scale_constraints(m, basis)`` function 
to scale the constraints before sending the model to the solver. To scale the variables, the user
must set the ``nlp_scaling_method`` option to "user-scaling".

Example
-------
.. testcode::

    from pyomo.environ import Suffix, ConcreteModel, Var, NonNegativeReals, \
        Constraint, Objective, SolverFactory
    from idaes.core.util.scaling import (scale_constraints, ScalingBasis,
                                     calculate_scaling_factors)
    from math import isclose

    # create badly scaled model
    var_value_magnitude = 1e-12
    m = ConcreteModel()

    m.x = Var(initialize=var_value_magnitude, domain=NonNegativeReals)
    m.y = Var(initialize=var_value_magnitude, domain=NonNegativeReals)
    m.z = Var(initialize=var_value_magnitude ** 2, domain=NonNegativeReals)

    m.c_1 = Constraint(expr=m.x + m.y <= var_value_magnitude)
    m.c_2 = Constraint(expr=m.z == (m.x * m.y))
    m.obj = Objective(expr=-m.z)

    # create and specify scaling factors and expressions
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)
    m.scaling_expression = Suffix(direction=Suffix.LOCAL)
    # NOTE: the direction of the scaling_expression is LOCAL
    m.scaling_factor[m.x] = 1 / var_value_magnitude
    m.scaling_factor[m.y] = 1 / var_value_magnitude
    m.scaling_factor[m.c_1] = 1 / var_value_magnitude
    m.scaling_expression[m.z] = 1 / (m.x * m.y)
    m.scaling_expression[m.c_2] = 1 / (m.x * m.y)
    m.scaling_expression[m.obj] = 1 / (m.x * m.y)

    # calculate scaling factors from scaling expression
    calculate_scaling_factors(m, basis=ScalingBasis.InverseVarScale)
    # scale constraints
    scale_constraints(m)
    # NOTE: After the constraints are scaled, their scaling factor and expression
    # are set to 1.

    # call solver with user-scaling option
    solver = SolverFactory('ipopt')
    solver.options = {'nlp_scaling_method': 'user-scaling'}
    results = solver.solve(m, tee=False)
    assert (isclose(m.z.value, (0.5 * var_value_magnitude) ** 2, rel_tol=1e-3))


Scaling Expression Basis
------------------------
The general guideline for calculating scaling factors from scaling expressions
is to use the expected magnitude of the variables. The magnitude
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



