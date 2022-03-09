Miscellaneous Utility Methods
=============================

.. contents:: Contents
    :depth: 2

.. module:: idaes.core.util.misc

get_solver
----------

The `get_solver` method is primarily intended for use when a solver object is required in a general model. The `get_solver` method takes two optional arguments which allow the user to specify a specific solver and/or solver options if required, or to use the default solver specified in the IDAES Configuration if no arguments are provided.

.. autofunction:: get_solver

Variable-Like Expressions
-------------------------

There are a number of cases within IDAES where a modeler may wish to use an `Expression` in place of a `Var` to reduce the complexity of their model. A common example of this is in the ideal `Separator` unit where the outlet `Ports` use `Expressions` for the state variable in order to reduce the number of variables (and thus constraints) in the model.

In these cases, it is possible that a user might mistake the `Expression` for a `Var` and attempt to use methods such as `fix()` on it. In order to provide the user with a useful error message informing them that this will not work, IDAES has created a derived `VarLikeExpression` component for these situations. This component derives directly from Pyomo’s `Expression` component and implements common methods associated with `Vars` which will return an error message informing the user that the component is an `Expression`, and a suggestion on how to proceed.


.. autoclass:: VarLikeExpression
  :members:

.. autoclass:: SimpleVarLikeExpression
  :members:

.. autoclass:: IndexedVarLikeExpression
  :members:

.. autoclass:: _GeneralVarLikeExpressionData
  :members:

Other Utility Methods
---------------------

This library also contains a number of other utility methods that do not fall under another category.

.. autofunction:: add_object_reference

.. autofunction:: extract_data

.. autofunction:: set_param_from_config
