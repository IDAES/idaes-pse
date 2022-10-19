Generating Polynomial Models with PySMO
===========================================

.. note::
   The IDAES surrogate API is a wrapper around the original PySMO surrogate interface.


The *pysmo_surrogate.PysmoPolyTrainer* method learns polynomial models from data. Presented with a small number of samples generated from experiments or computer simulations, the approach determines the most accurate polynomial approximation by comparing the accuracy and performance of polynomials of different orders and basis function forms.

*pysmo_surrogate.PysmoPolyTrainer* considers three types of basis functions

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
To generate a polynomial model with PySMO, the  *pysmo_surrogate.PysmoPolyTrainer* trainer is called with the desired optional arguments:

.. code:: python

   # Required imports
   >>> from idaes.core.surrogate.pysmo_surrogate import PysmoPolyTrainer
   >>> import pandas as pd

   # Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   # Define the input and output labels
   input_labels = ['X1', 'X2']
   output_labels = ['Y1','Y2']

   # Create the PySMO regression trainer object
   >>> pr_trainer = PysmoPolyTrainer(input_labels=input_labels, output_labels=output_labels, training_dataframe = data_training)

   # Set PySMO options
   >>> pr_trainer.config.maximum_polynomial_order = 4

   # Train the model
   >>> poly_train = pr_trainer.train_surrogate()


.. note::
   **maximum_polynomial_order** refers to the maximum polynomial order to be considered when training the surrogate and must be specified.


Configuration Options
----------------------
In addition to **maximum_polynomial_order**, ``PysmoPolyTrainer`` takes the following optional arguments that help define/improve the final regression model:

.. list-table:: PySMO regression options
   :widths: 20 20 60
   :header-rows: 1

   * - **Option**
     - Configuration argument
     - Description
   * - **multinomials**
     - *pysmo_poly_trainer.config.multinomials*
     - Boolean option which determines whether bivariate terms are considered in polynomial generation.
   * - **solution_method**
     - *pysmo_poly_trainer.config.solution_method*
     - | Method used to solve the regression option:
       | BFGS ('BFGS'), maximum likelihood ('MLE') or Pyomo least squares minimization ('pyomo'). 
       | Default is 'mle'.
   * - **training_split** 
     - *pysmo_poly_trainer.config.training_split*
     - | Option which determines fraction of training data to be used for regression training (the rest will be for testing). 
       | Note that this is different from the model training/test validation process. 
       | PySMO's regression tool internally validates the quality of the model before it is returns to the user. 
       | Default is 0.8.
   * - **extra_features** 
     - *pysmo_poly_trainer.config.extra_features*
     - | Option for defining additional desired non-regular regression terms. 
       | See section on custom basis functions for more details.

Custom Basis Functions
----------------------

PySMO has advanced capabilities that allows the user to supply custom basis functions to be parts of the regression process. Custom basis functions can be added to the built-in functions to expand the functional forms available. This custom basis function feature allow the user to incorporate apriori knowledge or physics-based information about the system into the regression model. 

The custom basis functions may be simple multivariate relationships (e.g. :math:`x_{1}/x_{2}`), trigonometric terms (e.g. :math:`sin(x_{1})`) or logarithmic terms (e.g. :math:`log(x_{2})`).

To use this advanced capability in PySMO, the *config.extra_features* optional argument is supplied. Note it is necessary to use the xlabels assigned to the input parameters.

.. code-block:: python
  
  pr_trainer.config.custom_basis_functions = ["x1/ x2", "sin(x2)", "...", "..." ...]


Output
-------
The result of the ``pysmo_surrogate.PysmoPolyTrainer`` method is a python object containing information about the problem set-up, the final optimal polynomial order, the polynomial coefficients and different error and quality-of-fit metrics such as the mean-squared-error (MSE) and the :math:`R^{2}` coefficient-of-fit. 


Visualization
-------------

.. rubric:: Visualizing Surrogate Model Results

For visualizing PySMO-trained surrogates via parity and residual plots, see :ref:`Visualizing Surrogate Model Results<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>`.


Building an IDAES Surrogate Object
------------------------------------
To add the model to an IDAES flowsheet, the SurrogateTrainer object needs to be transformed into an IDAES SurrogateObject object. This is done by calling PySMOSurrogate and passing the generated surrogate expressions, along with variable labels and optionally the bounds:

.. code:: python

   >>> surr = PysmoSurrogate(poly_train, input_labels, output_labels, input_bounds)


Prediction with *pysmo_surrogate.PysmoPolyTrainer* models
-----------------------------------------------------------
Once a polynomial model has been trained and the SurrogateObject object created, predictions for values at previously unsampled points :math:`x_{unsampled}` (a Pandas dataframe) can be evaluated by calling the ``evaluate_surrogate()`` method on the unsampled points:

.. code:: python

   >>> y_unsampled = surr.evaluate_surrogate(x_unsampled)
 

Confidence intervals for *pysmo_surrogate.PysmoPolyTrainer* models
--------------------------------------------------------------------
[Needs to be moved over to new intweface ---coming soon]
The confidence intervals for the regression paramaters may be viewed using the method ``confint_regression``.


Saving and Loading Surrogates
------------------------------


Flowsheet Integration
----------------------

The result of the polynomial training process can be passed directly into a process flowsheet using the IDAES SurrogateBlock option. 
The following code snippet demonstrates how an output polynomial model may be integrated directly into a Pyomo flowsheet as
an objective:

.. code:: python

   # Required imports
   >>> import pyomo.environ as pyo
   >>> from idaes.core.surrogate.pysmo import polynomial_regression
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



PySMO Examples
----------------

For an example of optimizing a flowsheet containing an PySMO-trained surrogate model, see the `Autothermal reformer flowsheet optimization example <https://github.com/IDAES/examples-pse/blob/main/src/Examples/SurrMod/FlowsheetOptimization/PySMO_flowsheet_optimization.ipynb>`_.


References:
----------------
[1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide", https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801