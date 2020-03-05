Variable Replacement
====================

.. module:: idaes.core.plugins.variable_replace

There are a number of cases where it can be convenient to replace one variable
for another. IDAES offers a convenient variable replacement transformation. This
transformation is not reversible and can significantly alter the model structure.

An example use of this transformation, is a parameter estimation problem where
a model contains several instances of a particular sub-model and each model
contains a variable (:math:`\beta`) for a model parameter to be estimated. In
many cases :math:`\beta` should be the same across all sub-models. One approach
to this problem would be to add equality constraints to equate all the
:math:`\beta`'s.  Another approach would be to use the variable replacement
transformation to replace the individual :math:`\beta`'s with a single global
:math:`\beta` variable.

Example
-------

The following example demonstrates the basic usage of the transformation.

.. testcode::

  import idaes.core.plugins # Load IDAES plugins
  import pyomo.environ as pyo

  # Use Pyomo's transformation factory to create the transformation object
  rp = pyo.TransformationFactory("replace_variables")

  # Create an example model
  m = pyo.ConcreteModel()
  m.x = pyo.Var({1,2,3}, initialize=2)
  m.new_x = pyo.Var({1,2,3}, initialize=3)
  m.e1 = pyo.Expression(expr=sum(m.x[i] for i in m.x))

  # Apply the transformation to the model, the substitute argument contains a list
  # of replacements, each element is a list-like object where the first element is
  # a variable to be replaced by the second element.
  rp.apply_to(m, substitute=[(m.x, m.new_x)])

  # See that the variable was replaced
  print(pyo.value(m.e1)) # since new_x has a value of 3 the expression value is 9

Output:

.. testoutput::

    9

Usage
-----

There are three basic steps to using the variable replacement transformation.

  1. Import anything from the ``idaes`` package; this will cause the IDAES
     plugins to be loaded.
  2. Use Pyomo's transformation factory to create a variable replacement
     transformation object (e.g.
     ``rp = TransformationFactory("replace_variables")``.
  3. Call the transformation object's ``apply_to()`` method to apply the
     transformation.

The ``apply_to(instance, substitute)`` method takes two arguments ``instance`` and
``substitute``. The instance argument is a model or block to apply the transformation
to. The substitute argument is a list-like object with substitutions.  Each element
is a two-element list-like object where the first element is a Pyomo Var, IndexedVar
element or Reference to the variable to replace and the second element is a Pyomo
Var, IndexedVar element or Reference to replace the first element with.

Indexed variables are allowed. The index set of the variable to replace must be
a subset of the index set of the variable to replace it with. It can also be
useful to use a Pyomo Reference to emulate an indexed variable, so this is also
supported.

ReplaceVariables Class
----------------------

The transformation object class is ReplaceVariables.

.. autoclass:: ReplaceVariables
  :members:
