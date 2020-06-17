Utility Methods
===============

.. warning:: This section is currently being developed

IDAES provides several utility methods to help in formulating and analyzing models and flowsheets.
Here, we list and described some commonly used utility methods, while a full list of IDAES Utility 
Methods is found in the 
:ref:`technical specifications<technical_specs/core/util/index:Utility Methods>`.

Degrees of freedom
------------------
Knowing and manipulating the degrees of freedom of a model is a core component of the recommended 
:ref:`workflow<user_guide/workflow:Workflow>` for IDAES. The function 
``degrees_of_freedom(m)`` in the ``idaes.core.util.model_statistics`` module readily determines the
degrees of freedome for a model. Additional details are provided
:ref:`here<technical_specs/core/util/model_statistics:Degrees of Freedom Method>`.

Scaling
-------
Creating well scaled models is important for increasing the efficiency and 
reliability of solvers. IDAES provides several functions that help scale the model in the
``idaes.core.util.scaling`` module. Additional details are provided
:ref:`here<technical_specs/core/util/scaling:Scaling Methods>`.