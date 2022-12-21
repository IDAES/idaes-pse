Scaling Toolbox
===============

.. note::

  In v2.0, IDAES is beginning to deploy a new suite of scaling tools. This documentation refers to the new tools.
  For documentation of the older scaling tools, :ref:`see here<reference_guides/core/util/scaling:Scaling Methods>`.

.. module:: idaes.core.util.scaling

.. contents:: Contents
    :depth: 2

Importance of Model Scaling
---------------------------

Creating well scaled models is important for increasing the efficiency and reliability of solvers. However, depending on units of measure and process scale, variables and constraints for process applications are often badly scaled unless efforts are taken to rescale the problem.

The IDAES-IP takes the approach of writing the problem in the most natural form, and then applying model transformations to convert the problem to a more tractable form for the solver. This is achieved through leveraging the `Pyomo Scaling Transformation <https://pyomo.readthedocs.io/en/stable/model_transformations/scaling.html>`_ and model ``Suffixes``. In this way, users may assign scaling factors to any variable, constraint or objective function in the problem in order to convert it from its natural form into a better scaled equivalent. The advantage of using ``Suffixes`` and transformations over more traditional approaches is that this allows scaling to be adjusted as required without needing to rework the model code.

As a general rule of thumb, all variables, component and objective functions in a problem should be scaled to have a magnitude between 1 and 10 in order to ensure maximum robustness of the problem (although this can vary between different solvers). As quantities within process models can vary significantly (even within the same process), achieving good scaling requires a lot of input from the user to determine the expected magnitude of each variable and constraint and thus the appropriate scaling factor, To assist users with this, the IDAES-IP provides a number of utility functions for setting scaling factors and calculating initial guesses for scaling factors from available information. These tools can be found in the ``idaes.core.util.scaling`` module.

Setting Scaling Factors
-----------------------

Suffixes are used to specify scaling factors for any component in an IDAES model. These suffixes are created when needed by calling the ``set_scaling_factor()`` function. Using the ``set_scaling_factor()``, ``get_scaling_factor()``, and ``unset_scaling_factor()`` eliminates the need for users to deal directly with scaling suffixes, and ensures that scaling factors are stored in the correct location.

.. autofunction:: set_scaling_factor

.. autofunction:: get_scaling_factor

.. autofunction:: unset_scaling_factor

For variables, the ``set_variable_scaling_from_current_value`` can be used to automatically set the scaling factor for a given variable (or all variables in a given model) based on the current value of the variable.

.. autofunction:: set_variable_scaling_from_current_value

.. note ::

  If this function is used on a variable which has a current value of 0 or no current value (i.e., ``var.value == None``) then a warning will be logged and no scaling factor will be set for the variable.

Default Scaling Factors
-----------------------

Process models are generally large and contain a large number of components, thus it is often infeasible for a user to manually set scaling factors for all components individually. Additionally, getting good initial values in order to use the ``set_variable_scaling_from_current_value`` function requires solving the model which in turn requires initial scaling of the model. In order to provide a starting point for initialization, all IDAES models contain a ``default_scaling_factors`` dict which allows developers and users to assign default scaling factors for any component in a given model.

The following methods are available on all IDAES models which can be used to access and manipulate the default scaling factors:

.. module:: idaes.core.base.process_base
  :noindex:
  
.. autoclass:: ProcessBlockData
  :members: set_default_scaling, get_default_scaling, unset_default_scaling
  :noindex:

These default scaling factors can be applied to a model using the ``set_scaling_from_default`` utility function.

.. module:: idaes.core.util.scaling
  :noindex:

.. autofunction:: set_scaling_from_default

.. note::

  Default scaling factors are NOT automatically applied to a model. If ``set_scaling_from_default`` is not called, then the default scaling factors will not be used.

Calculating Constraint and Objective Scaling Factors
----------------------------------------------------

If all variables in a problem have been assigned scaling factors, it is possible to automatically evaluate the terms in all expression (i.e., constraints and objective functions) using the inverse of the variable scaling factors as a nominal value for each variable. This information can then be used to estimate scaling factors for each expression. The IDAES-IP provides a number of utility functions for automatically calculating scaling factors for constraints and objective functions based on different approaches.

.. autofunction:: set_constraint_scaling_harmonic_magnitude

.. autofunction:: set_constraint_scaling_max_magnitude

.. autofunction:: set_constraint_scaling_min_magnitude

These functions make use of the ``NominalValueExtractionVisitor`` class which automatically walks the entire expression tree and determines the nominal value (expected magnitude and sign) for each additive term in the expression. Given an expression of the form :math:`f(x) = A(x) + B(x) + C(x)`, this class will return a list of the nominal values of :math:`A(x)`, :math:`B(x)` and :math:`C(x)` based on the scaling factors assigned to the variables in each sub-expression. These values can then be used to determine the best scaling factor for the overall expression.

.. autoclass:: NominalValueExtractionVisitor
  :members:

Identifying Scaling Issues
-----------------------------

A number of utility functions are available to help identify potential scaling issues in a problem.

.. autofunction:: report_scaling_issues

.. autofunction:: unscaled_variables_generator

.. autofunction:: list_unscaled_variables

.. autofunction:: badly_scaled_var_generator

.. autofunction:: list_badly_scaled_variables

.. autofunction:: unscaled_constraints_generator

.. autofunction:: list_unscaled_constraints

.. autofunction:: extreme_jacobian_entries

.. autofunction:: extreme_jacobian_rows

.. autofunction:: extreme_jacobian_columns

.. autofunction:: jacobian_cond
