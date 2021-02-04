=======================
uncertainty_propagation
=======================
The uncertainty_propagation module provides tools for calculating error propagation of the objective function and constraints that include uncertain parameters to be estimated for a given model. 
For a given model with a dataset, it first estimates parameters. The estimated parameters are passed to the optimization model. Finally, the variance of objective function and constraints with respect to the estimated parameters at the optimal solution are calculated.


**Dependencies**

`k_aug <https://github.com/dthierry/k_aug>`_ is required to use uncertainty_propagation module. If you have the IDAES framework installed, you can get k_aug with the IDAES extensions by running the following command:

.. code-block:: python

   idaes get-extensions

Basic Usage
------------
The ``uncertainty_propagation``'s main function is **quantify_propagate_uncertainty**. The quantify_propagate_uncertainty python function can have eight arguments and retunrs five outputs. For example

.. code-block:: python

    obj, theta, cov, propagation_f, propagation_c = quantify_propagate_uncertainty(model_function, model_uncertain,  data, theta_names, obj_function)

, where the **model_function** is a python function that generates an instance of the Pyomo model using **data** as the input argument. The **model_function** should have parameters with the name, **theta_names** to be estimated. The **model_function** is used to estimate parameters by minimizing a given **obj_function**. An example of **obj_function** is the SSE (the sum of squares error). Some outputs are generated in this step. The **obj** is an objective function value for the given **obj_function**, the **theta** is estimated parameters, and the **cov** is the covariance matrix of **theta**.

The **model_uncertain** is a python function that generates an instance of the Pyomo model or Pyomo ConcreteModel. The objective functions and constraints of **model_uncertain** can include variables with the name, **theta_names**. Note that **model_uncertain** requires an objective function in the model. If the model is not an optimization problem, the objection fucntion can be an arbitrary number, e.g 0. The variable **theta_names** are automatically fixed with the estimated **theta**. Then, the Pyomo model is solved with the ipopt solver.
The gradients vector of the objective function and the gradients vectors constraints with respect to **theta_names**  are calculated at the optimal solution with the k_aug.  Let the :math:`g_f` is the gradients vector of the objective function, and :math:`g_{c1},g_{c2},\ldots,g_{cn}` be the gradients vectors constraints, where :math:`n` constraints includes some variables in **theta_names**. Finally, the error propagation of the objective function (**propagation_f**) and constraints (**propagation_c**) are calculated as 


   propagation_f :math:`= g_{f}\times cov \times g_{f}^T`

   propagation_c :math:`=\Big\{constraints \ 1: g_{c1}\times cov \times g_{c1}^T, \ constraints \ 2: g_{c2}\times cov \times g_{c2}^T,\ \ldots \ , \ constraints \ n: g_{cn}\times cov \times g_{cn}^T\}`

, where the superscript :math:`T` means transpose.

Available Methods
------------------

.. automethod:: idaes.apps.uncertainty_propagation.uncertainties.quantify_propagate_uncertainty
.. automethod:: idaes.apps.uncertainty_propagation.uncertainties.propagate_uncertainty
.. automethod:: idaes.apps.uncertainty_propagation.uncertainties.get_sensitivity
.. automethod:: idaes.apps.uncertainty_propagation.uncertainties.clean_variable_name

Examples
--------
Two `examples
<https://github.com/IDAES/idaes-pse/tree/main/idaes/apps/uncertainty_propagation/examples>`_ are provided to illustrate the detailed usage of uncertainty_propagation. In each case, a Jupyter notebook with explanations as well as an equivalent Python script is provided.

