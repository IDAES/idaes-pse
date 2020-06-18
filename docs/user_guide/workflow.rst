Workflow
========

The section describes the recommended workflow for the IDAES framework. Users are encouraged 
to follow these best practices in order to efficiently create reliable and accurate models.

This workflow is demonstrated in the tutorials and examples on the
`IDAES examples site <https://examples-pse.readthedocs.io/en/stable/>`_.

1. Importing modules
--------------------
To build a model on the IDAES framework, the following should be imported (as demonstrated in 
the IDAES tutorials):

* Pyomo classes and functions (e.g. ConcreteModel, SolverFactory, Arc, Var, Constraint)
* IDAES flowsheet block (`from idaes.core import FlowsheetBlock`)
* Unit models (from the model libraries or user defined custom models)
* Property models (from the model libraries or user defined custom models)
* Utility tools (common tools include degrees of freedom and scaling, full list is provided here)


2. Building flowsheets
----------------------

3. Scaling
----------

4. Specifying
-------------

5. Initializing
---------------

6. Optimizing
-------------

7. Analyzing
------------
