Model Statistics Methods
========================

The IDAES toolset contains a number of utility functions which are useful for quantifying model statistics such as the number of variable and constraints, and calculating the available degrees of freedom in a model. These methods can be found in ``idaes.core.util.model_statistics``.

The most commonly used methods are ``degrees_of_freedom`` and ``report_statistics``, which are described below.

Degrees of Freedom Method
-------------------------

The ``degrees_of_freedom`` method calculates the number of degrees of freedom available in a given model. The calcuation is based on the number of unfixed variables which appear in active constraints, minus the number of active equality constraints in the model. Users should note that this method does not consider inequality or deactivated constraints, or variables which do not appear in active equality constraints.

.. autofunction:: idaes.core.util.model_statistics.degrees_of_freedom

Report Statistics Method
------------------------

The ``report_statistics`` method provides the user with a summary of the contents of their model, including the degrees of freedom and a break down of the different ``Variables``, ``Constraints``, ``Objectives``, ``Blocks`` and ``Expressions``. This method also includes numbers of deactivated components for the user to use in debugging complex models.

.. note::

    This method only considers Pyomo components in activated ``Blocks``. The number of deactivated ``Blocks`` is reported, but any components within these ``Blocks`` are not included.

.. admonition:: Example Output

    Model Statistics  

    Degrees of Freedom: 0

    Total No. Variables: 52

        No. Fixed Variables: 12

        No. Unused Variables: 0 (Fixed: 0)

        No. Variables only in Inequalities: 0 (Fixed: 0)

    Total No. Constraints: 40

        No. Equality Constraints: 40 (Deactivated: 0)

        No. Inequality Constraints: 0 (Deactivated: 0)

    No. Objectives: 0 (Deactivated: 0)

    No. Blocks: 14 (Deactivated: 0)

    No. Expressions: 2

.. autofunction:: idaes.core.util.model_statistics.report_statistics

Other Statistics Methods
------------------------

In addition to the methods discussed above, the ``model_statistics`` module also contains a number of methods for quantifying model statistics which may be of use to the user in debugging models. These methods come in three types:

* Number methods (start with ``number_``) return the number of components which meet a given criteria, and are useful for quickly quantifying different types of components within a model for determining where problems may exist.
* Set methods (end with ``_set``) return a Pyomo ``ComponentSet`` containing all components which meet a given criteria. These methods are useful for determining where a problem may exist, as the ``ComponentSet`` indicates which components may be causing a problem.
* Generator methods (end with ``_generator``) contain Python ``generators`` which return all components which meet a given criteria.

Available Methods
^^^^^^^^^^^^^^^^^

.. automodule:: idaes.core.util.model_statistics
    :exclude-members: degrees_of_freedom, report_statistics
    :members:

