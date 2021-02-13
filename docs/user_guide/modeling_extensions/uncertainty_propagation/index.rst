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
The ``uncertainty_propagation``'s main function is **quantify_propagate_uncertainty**. The quantify_propagate_uncertainty python function can have eight arguments and retunrs a namedtuple output. 
The following example shows a usage of **quantify_propagate_uncertainty** with a Rooney and Biegler model [1].

.. code:: python

   # Required imports
   >>> from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_uncertainty
   >>> import pyomo.environ as *
   >>> import pandas as pd

   # Define rooney_biegler_model
   >>> def rooney_biegler_model(data):
           model = ConcreteModel()
           model.asymptote = Var(initialize = 15)
           model.rate_constant = Var(initialize = 0.5)
           
           def response_rule(m, h):
               expr = m.asymptote * (1 - exp(-m.rate_constant * h))
               return expr
           model.response_function = Expression(data.hour, rule = response_rule)
    
           def SSE_rule(m):
               return sum((data.y[i] - m.response_function[data.hour[i]])**2 for i in data.index)
           model.SSE = Objective(rule = SSE_rule, sense=minimize)
           return model

   # Define rooney_biegler_model_opt
   >>> def rooney_biegler_model_opt(data):
           model = ConcreteModel()
           model.asymptote = Var(initialize = 15)
           model.rate_constant = Var(initialize = 0.5)
           model.obj = Objective(expr=model.asymptote*(1-exp(-model.rate_constant*10))
                                 , sense=minimize)
           return model
   
   # Define data
   >>> data = pd.DataFrame(data=[[1,8.3],
                          [2,10.3],
                          [3,19.0],
                          [4,16.0],
                          [5,15.6],
                          [7,19.8]],
                    columns=['hour', 'y'])

   # Define variable_name
   >>> variable_name = ['asymptote', 'rate_constant']

   # Define obj_function
   >>> def obj_function(model, data):
           expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
           return expr

   # Run quantify_propagate_uncertainty
   >>> results = quantify_propagate_uncertainty(rooney_biegler_model, rooney_biegler_model_opt,  
                                                data, variable_name, obj_function)  
    


, where the **rooney_biegler_model** is a python function that generates an instance of the Pyomo model using **data** as the input argument. The **rooney_biegler_model** should have parameters with the name, **variable_name** to be estimated. The **rooney_biegler_model** is used to estimate parameters by minimizing a given **obj_function**. An example of **obj_function** is the SSE (the sum of squares error). Some outputs are generated in this step. The **results.obj** is an objective function value for the given **obj_function**, the **theta** is estimated parameters, and the **results.cov** is the covariance matrix of **results.theta**.

The **rooney_biegler_model_opt** is a python function that generates an instance of the Pyomo model or Pyomo ConcreteModel. The objective functions and constraints of **rooney_biegler_model_opt** can include variables with the name, **variable_name**. Note that **rooney_biegler_model_opt** requires an objective function in the model. If the model is not an optimization problem, the objection fucntion can be an arbitrary number, e.g 0. The variable **results.theta_names** are automatically fixed with the estimated **results.theta**. Then, the Pyomo model is solved with the ipopt solver.
The gradients vector of the objective function and the gradients vectors constraints with respect to **variable_name**  are calculated at the optimal solution with the k_aug.  Let the :math:`g_f` is the gradients vector of the objective function, and :math:`g_{c1},g_{c2},\ldots,g_{cn}` be the gradients vectors constraints, where :math:`n` constraints includes some variables in **variable_name**. Finally, the error propagation of the objective function (**results.propagation_f**) and constraints (**results.propagation_c**) are calculated as 


   results.propagation_f :math:`= g_{f}\times cov \times g_{f}^T`

   results.propagation_c :math:`=\Big\{constraints \ 1: g_{c1}\times cov \times g_{c1}^T, \ constraints \ 2: g_{c2}\times cov \times g_{c2}^T,\ \ldots \ , \ constraints \ n: g_{cn}\times cov \times g_{cn}^T\}`

, where the superscript :math:`T` means transpose.

Available Functions
-------------------

.. currentmodule:: idaes.apps.uncertainty_propagation.uncertainties

.. autofunction:: quantify_propagate_uncertainty

.. autofunction:: propagate_uncertainty

.. autofunction:: get_sensitivity

.. autofunction:: clean_variable_name

Examples
--------
Two `examples
<https://github.com/IDAES/idaes-pse/tree/main/idaes/apps/uncertainty_propagation/examples>`_ are provided to illustrate the detailed usage of uncertainty_propagation. In each case, a Jupyter notebook with explanations as well as an equivalent Python script is provided.

References
----------
[1] Rooney, W. C. and Biegler, L. T. (2001). Design for model parameter uncertainty using nonlinear confidence regions. *AIChE Journal*,  47(8), 1794-1804.