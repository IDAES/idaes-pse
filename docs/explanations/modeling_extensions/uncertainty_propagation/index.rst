===============================
Uncertainty Propagation Toolbox
===============================
The uncertainty_propagation module quantifies and propagates parametric uncertainties through optimization and simulation problem based on IDAES models. The module has two core features:

1. Given a parameter estimation model and data, calculate the least squares best fit parameter values and estimate their covariance.

2. Given a process model and the covariance for its parameters, estimate the variance of the optimal solution and the objective function.

Consider the optimization problem:

.. math::
    :nowrap:
    
    \begin{align*}
        \mathrm{minimize} ~~ & f(x,p) \\
        \mathrm{s.t.} \quad ~ & c(x,p) = 0 \\
        & x_{lb} \leq x \leq x_{ub}
    \end{align*}

Here :math:`x \in \mathbb{R}^{n\ \times\ 1}` are the decision variables, :math:`p \in \mathbb{R}^{m\ \times\ 1}` are the parameters, :math:`f(x,p):\ \mathbb{R}^{n\ \times\ 1}\ \times \mathbb{R}^{m\ \times\ 1} \rightarrow \mathbb{R}` is the objective function, :math:`c(x,p) = \{c_1(x,p), \ldots, c_k(x,p)\}\ :\ \mathbb{R}^{n\ \times\ 1}\ \times \mathbb{R}^{m\ \times\ 1} \rightarrow \mathbb{R}^{k\ \times\ 1}` are the constraints, and :math:`x_{lb}` and :math:`x_{ub}` are the lower and upper bounds, respectively.

Let :math:`x^*` represent the optimal solution given parameters :math:`p^*`. In many process systems engineering problems, :math:`p^*` is estimated from data and has some uncertainty represented with covariance matrix :math:`\Sigma_p`. This toolbox estimates the uncertainty in the optimal solution :math:`x^*` and objective function value :math:`f(x^*, p^*)` induced by uncertainty in :math:`p^*`.

Based on a first-order error propagation formula, the variance of the objective function is:

.. math::

    \text{Var}[f(x^*,p^*)] = \left(\frac{\partial f}{\partial p} + \frac{\partial f}{\partial x}\frac{\partial x}{\partial p} \right) \Sigma_p  \left(\frac{\partial f}{\partial p} + \frac{\partial f}{\partial x}\frac{\partial x}{\partial p} \right)^T 


Likewise, the variance in constraint :math:`k` is:

.. math::

    \text{Var}[c_k(x^*,p^*)] = \left(\frac{\partial c_k}{\partial p} + \frac{\partial c_k}{\partial x}\frac{\partial x}{\partial p} \right)  \Sigma_p \left(\frac{\partial c_k}{\partial p} + \frac{\partial c_k}{\partial x}\frac{\partial x}{\partial p} \right)^T 
    
Note that :math:`\text{Var}[c_k(x^*,p^*)] \approx 0` because the constraints remain feasible for a small perturbation in :math:`p`, i.e.,  :math:`c(x,p) = 0`.

All gradients are calculated with `k_aug <https://github.com/dthierry/k_aug>`_  [1]. More specifically, :math:`\frac{\partial f}{\partial p}, \frac{\partial f}{\partial x}, \frac{\partial c_1}{\partial p}, \frac{\partial c_1}{\partial x}, \ldots, \frac{\partial c_k}{\partial p}, \frac{\partial c_k}{\partial x}` evaluated at :math:`(x^*, p)` are computed via automatic differentiation whereas :math:`\frac{\partial x}{\partial p}` are computed via nonlinear programming sensitivity theory.

The covariance matrix :math:`\Sigma_p` is either user supplied or obtained via regression (with ``Pyomo.Parmest``). 

**Dependencies**

`k_aug <https://github.com/dthierry/k_aug>`_ [1] is required to use uncertainty_propagation module. The ``k_aug`` solver executable is easily installed via the ``idaes get-extensions`` command.

Basic Usage
------------

This toolbox has two core functions:

1. **propagate_uncertainty**: Given an IDAES (Pyomo) process model with parameters :math:`p` and covariance :math:`\Sigma_p`, estimate :math:`\text{Var}[f(x^*,p)]`.

2. **quantify_propagate_uncertainty**: Given an IDAES (Pyomo) regression model and data, first estimate parameters :math:`p` and covariance :math:`\Sigma_p`. Then given a second IDAES (Pyomo) process model, estimate :math:`\text{Var}[f(x^*,p)]`.

The following example shows the usage of **quantify_propagate_uncertainty** for a reaction kinetic problem from Rooney and Biegler [2]. Consider the mathematical model:

.. math::
    y_{i}(\theta_1, \theta_2, t_i) = \theta_1 (1 - e^{-\theta_2 t_i}), \quad i = 1, ..., n

Here :math:`y_i` is the concentration of the chemical product at time :math:`t_i` and :math:`p = (\theta_1, \theta_2)` are the fitted model parameters. 

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
   >>> data = pd.DataFrame(data=[[1,8.3],[2,10.3],[3,19.0],[4,16.0],[5,15.6],[7,19.8]],
                           columns=['hour', 'y'])

   # Define variable_name
   >>> variable_name = ['asymptote', 'rate_constant']

   # Define SSE_obj_function
   >>> def SSE_obj_function(model, data):
           expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
           return expr

   # Run quantify_propagate_uncertainty
   >>> results = quantify_propagate_uncertainty(rooney_biegler_model, rooney_biegler_model_opt,  
                                                data, variable_name, obj_function)  
    


The Python function **rooney_biegler_model** generates a Pyomo regression model using the input Pandas DataFrame **data**. This model contains the attributes ``asymptote`` and ``rate_constant`` which are the parameters :math:`p` to be estimated by minimizing the sum of squared errors (SSE). The list **variable_name** contains strings with these attribute names.

Similarly, the Python function **rooney_biegler_model_opt** returns a concrete Pyomo model which is the process optimization problem. This specific example has no degrees of freedom once :math:`p` is specified; the objective function computes the product concentration at time :math:`t=10`.

The function **quantify_propagate_uncertainty** returns the object **results** which contains several important attributes:

* **results.theta_names** contains the names of parameters :math:`p`
* **results.theta** contains the estimate values for parameters :math:`p`
* **results.gradient_f** contains the gradient :math:`\frac{\partial f}{\partial x},~ \frac{\partial f}{\partial p}`
* **results.gradient_c** contains the Jacobians :math:`\frac{\partial c}{\partial x},~ \frac{\partial c}{\partial p}`
* **results.dsdp** contains the Jacobians :math:`\frac{\partial x}{\partial p},~\frac{\partial p}{\partial p}``
* **results.propagation_f** contains the estimate variance of the objective function 

**Important Notes:**

1. The uncertain parameters :math:`p` should be declared as ``Var`` in Pyomo.

2. The uncertain parameters :math:`p` should not be fixed in Pyomo. Instead, set the upper and lower bounds equal. If they are fixed, ``k_aug`` will give an error message that the optimization problem has too few degrees of freedom.

Available Functions
-------------------

.. currentmodule:: idaes.apps.uncertainty_propagation.uncertainties

.. autofunction:: quantify_propagate_uncertainty

.. autofunction:: propagate_uncertainty

.. autofunction:: clean_variable_name

Examples
--------
Two `examples
<https://github.com/IDAES/idaes-pse/tree/main/idaes/apps/uncertainty_propagation/examples>`_ are provided to illustrate the detailed usage of uncertainty_propagation. In each case, a Jupyter notebook with explanations as well as an equivalent Python script is provided.

References
----------
[1] David Thierry (2020). k_aug, https://github.com/dthierry/k_aug

[2] Rooney, W. C. and Biegler, L. T. (2001). Design for model parameter uncertainty using nonlinear confidence regions. *AIChE Journal*,  47(8), 1794-1804.
