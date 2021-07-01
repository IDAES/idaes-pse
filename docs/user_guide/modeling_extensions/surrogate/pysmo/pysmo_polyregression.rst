Generating Polynomial Models with PySMO
===========================================

The *pysmo.polynomial_regression* method learns polynomial models from data. Presented with a small
number of samples generated from experiments or computer simulations, the approach determines the most
accurate polynomial approximation by comparing the accuracy and performance of polynomials of different
orders and basis function forms.

*pysmo.polynomial_regression* considers three types of basis functions

* univariate polynomials,
* second-degree bivariate polynomials, and
* user-specified basis functions.

Thus, for a problem with :math:`m` sample points and :math:`n` input variables, the resulting polynomial is of the form

.. math::
    \begin{equation}
    y_{k}={\displaystyle \sum_{i=1}^{n}\beta_{i}x_{ik}^{\alpha}}+\sum_{i,j>i}^{n}\beta_{ij}x_{ik}x_{jk}+\beta_{\Phi}\Phi\left(x_{ik}\right)\qquad i,j=1,\ldots,n;i\neq j;k=1,\ldots,m;\alpha \leq 10\qquad\quad\label{eq:poly_eq}
    \end{equation}

Basic Usage
------------
To generate a polynomial model with PySMO, the  *pysmo.polynomial_regression* class is first initialized,
and then the method ``training`` is called on the initialized object:

.. code:: python

   # Required imports
   >>> from idaes.surrogate.pysmo import polynomial_regression
   >>> import pandas as pd

   # Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   # Initialize the PolynomialRegression class, extract the list of features and train the model
   >>> pr_init = polynomial_regression.PolynomialRegression(xy_data, xy_data, maximum_polynomial_order=3, *kwargs)
   >>> features = pr_init.get_feature_vector()
   >>> pr_init.training()

* **xy_data** is a two-dimensional python data structure containing the input and output training data. The output values **MUST** be in the last column.
* **maximum_polynomial_order** refers to the maximum polynomial order to be considered when training the surrogate.

**Optional Arguments**

* **multinomials** - boolean option which determines whether bivariate terms are considered in polynomial generation.
* **training_split** - option which determines fraction of training data to be used for training (the rest will be for testing). Default is 0.8.
* **number_of_crossvalidations** - Number of cross-validations during training. Default number is 3.

*pysmo.polynomial_regression* Output
---------------------------------------
The result of the ``pysmo.polynomial_regression`` method is a python object containing information
about the problem set-up, the final optimal polynomial order, the polynomial coefficients and different error and quality-of-fit metrics such as
the mean-squared-error (MSE) and the :math:`R^{2}` coefficient-of-fit. A Pyomo expression can be generated from the
object simply passing a list of variables into the function *generate_expression*:

.. code:: python

   # Create a python list from the headers of the dataset supplied for training
   >>> list_vars = []
   >>> for i in features.keys():
   >>>     list_vars.append(features[i])
   # Pass list to generate_expression function to obtain a Pyomo expression as output
   >>> print(pr_init.generate_expression(list_vars))

Prediction with *pysmo.polynomial_regression* models
-----------------------------------------------------
Once a polynomial model has been trained, predictions for values at previously unsampled points :math:*x_unsampled* can be evaluated by calling the
``predict_output()`` method on the unsampled points:

.. code:: python

   # Create a python list from the headers of the dataset supplied for training
   >>> y_unsampled = pr_init.predict_output(x_unsampled)
   
The confidence intervals for the regression paramaters may be viewed using the method ``confint_regression``.

Flowsheet Integration
----------------------

The result of the polynomial training process can be passed directly into a process flowsheet as an objective or a constraint.
The following code snippet demonstrates how an output polynomial model may be integrated directly into a Pyomo flowsheet as
an objective:

.. code:: python

   # Required imports
   >>> import pyomo.environ as pyo
   >>> from idaes.surrogate.pysmo import polynomial_regression
   >>> import pandas as pd

   # Create a Pyomo model
   >>> m = pyo.ConcreteModel()
   >>> i = pyo.Set(initialize=[1, 2])

   # Create a Pyomo variable with indexed by the 2D-set i with initial values {0, 0}
   >>> init_x = {1: 0, 2: 0}
   >>> def x_init(m, i):
   >>>     return (init_x[i])
   >>> m.x = pyo.Var(i, initialize=x_init)

   # Train a simple polynomial model on data available in csv format, resulting in the Python object polyfit
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)
   >>> pr_init = polynomial_regression.PolynomialRegression(xy_data, xy_data, maximum_polynomial_order=3)
   >>> features = pr_init.get_feature_vector()
   >>> polyfit = pr_init.training()

   # Use the resulting polynomial as an objective, passing in the Pyomo variable x
   >>> m.obj = pyo.Objective(expr=polyfit.generate_expression([m.x[1], m.x[2]]))

   # Solve the model
   >>> instance = m
   >>> opt = pyo.SolverFactory("ipopt")
   >>> result = opt.solve(instance, tee=True)

Further details about *pysmo.polynomial_regression* may be found by consulting the examples or reading the paper [...]

.. module:: idaes.surrogate.pysmo.polynomial_regression

Available Methods
------------------

.. autoclass:: idaes.surrogate.pysmo.polynomial_regression.FeatureScaling
    :members:

.. autoclass:: idaes.surrogate.pysmo.polynomial_regression.PolynomialRegression
    :members: __init__, get_feature_vector, set_additional_terms, training, predict_output, generate_expression, confint_regression 
	
References:
----------------
[1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide", https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801