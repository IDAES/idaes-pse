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
To generate a polynomial model with PySMO, the  *pysmo_surrogate.PysmoPolyTrainer* trainer is instantiated and initialized with the desired configuration arguments, and the training function ``train_surrogate`` is called:

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
     - *PysmoPolyTrainer.config.multinomials*
     - Boolean option which determines whether bivariate terms are considered in polynomial generation.
   * - **solution_method**
     - *PysmoPolyTrainer.config.solution_method*
     - | Method used to solve the regression option:
       | BFGS ('BFGS'), maximum likelihood ('MLE') or Pyomo least squares minimization ('pyomo'). 
       | Default is 'mle'.
   * - **training_split** 
     - *PysmoPolyTrainer.config.training_split*
     - | Option which determines fraction of training data to be used for regression training (the rest will be for testing). 
       | Note that this is different from the model training/test validation process. 
       | PySMO's regression tool internally validates the quality of the model before it is returns to the user. 
       | Default is 0.8.
   * - **extra_features** 
     - *PysmoPolyTrainer.config.extra_features*
     - | Option for defining additional desired non-regular regression terms. 
       | See section on custom basis functions for more details.

Custom Basis Functions
----------------------

PySMO has advanced capabilities that allows the user to supply custom basis functions to be parts of the regression process. Custom basis functions can be added to the built-in functions to expand the functional forms available. This custom basis function feature allow the user to incorporate *a priori* knowledge or physics-based information about the system into the regression model. 

The custom basis functions may be simple multivariate relationships (e.g. :math:`x_{1}/x_{2}`), trigonometric terms (e.g. :math:`sin(x_{1})`) or logarithmic terms (e.g. :math:`log(x_{2})`).

To use this advanced capability in PySMO, the *config.extra_features* optional argument is supplied. Note it is necessary to use the xlabels assigned to the input parameters.

.. code-block:: python
  
  pr_trainer.config.custom_basis_functions = ["x1/ x2", "sin(x2)", "...", "..." ...]


Output
-------
The result of the ``pysmo_surrogate.PysmoPolyTrainer`` method is a python object containing information about the problem set-up, the final optimal polynomial order, the polynomial coefficients and different error and quality-of-fit metrics such as the mean-squared-error (MSE) and the :math:`R^{2}` coefficient-of-fit. 


Confidence intervals for *PysmoPolyTrainer* models
--------------------------------------------------------------------

PySMO provides the user with the capability to compute confidence intervals for the regression parameters using the ``get_confidence_intervals`` method. This can be done by passing the result of the model training and the confidence interval value (*default=0.95*) into the ``pysmo_surrogate.PysmoPolyTrainer`` object:

.. code-block:: python
  
  >>> conf_int = pr_trainer.get_confidence_intervals(poly_train, confidence=0.99)

The resulting object ```conf_int`` is a dictionary containing upper and lower confidence bounds as well as the estimated standard errors for all of the regression parameters of the trained models in ``poly_train``. The dictionary keys correspond to the output variable list supplied in ``output labels``. 

Surrogate Visualization
------------------------
For visualizing PySMO-trained surrogates via parity and residual plots, see :ref:`Visualizing Surrogate Model Results<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>`.


Building an IDAES Surrogate Object
------------------------------------
To add the model to an IDAES flowsheet or generate model predictions, the SurrogateTrainer object needs to be transformed into an IDAES SurrogateObject object. This is done by calling ``PySMOSurrogate`` and passing the generated surrogate expressions, along with variable labels and optionally the bounds:

.. code:: python

   >>> surr = PysmoSurrogate(poly_train, input_labels, output_labels, input_bounds)

The resulting ``PysmoSurrogate`` object may be saved to (and reloaded from) a JSON file; for details, see :ref:`the PySMO main page<explanations/modeling_extensions/surrogate/api/pysmo/index:PySMO: Python-based Surrogate Modeling Objects>`.


Prediction with *PysmoPolyTrainer* models
-----------------------------------------------------------
Once a polynomial model has been trained and the SurrogateObject object created, predictions for values at previously unsampled points :math:`x_{unsampled}` (a Pandas dataframe) can be evaluated by calling the ``evaluate_surrogate()`` method on the unsampled points:

.. code:: python

   >>> y_unsampled = surr.evaluate_surrogate(x_unsampled)


Flowsheet Integration
----------------------
The result of the polynomial training process can be passed directly into a process flowsheet using the IDAES ``SurrogateBlock`` option. 
The following code snippet demonstrates how a saved polynomial model may be integrated directly into an IDAES flowsheet:

.. code:: python

   # Required imports
   >>> from pyomo.environ import Var, ConcreteModel, Constraint, SolverFactory, Objective, minimize
   >>> from idaes.core import FlowsheetBlock
   >>> from idaes.core.surrogate.pysmo_surrogate import PysmoSurrogate
   >>> from idaes.core.surrogate.surrogate_block import SurrogateBlock

   # Create a Pyomo model
   >>> m = pyo.ConcreteModel()
   >>> m.fs = FlowsheetBlock(default={"dynamic": False})

   # create input and output variables
   >>> m.fs.X1 = Var(initialize=0, bounds=(0, 5)) 
   >>> m.fs.X2 = Var(initialize=0, bounds=(0, 5)) 
   >>> m.fs.Y1 = Var(initialize=0) 
   >>> m.fs.Y2 = Var(initialize=0) 

   # create list of surrogate inputs and outputs for flowsheet
   >>> inputs = [m.fs.X1, m.fs.X2]
   >>> outputs = [m.fs.Y1, m.fs.Y2]

   # create the Pyomo/IDAES block that corresponds to the surrogate
   >>> m.fs.surrogate = SurrogateBlock(concrete=True)
   >>> surrogates_obj =PysmoSurrogate.load_from_file('poly_surrogate.json') # poly_surrogate.json is an existing surrogate JSON file
   >>> m.fs.surrogate.build_model(surrogates_obj, input_vars=inputs, output_vars=outputs)
   >>> m.fs.surrogate.pprint()

   # Set the variable Y1 as the model objective
   >>> m.fs.obj = Objective(expr=m.fs.Y1, sense=minimize)

   # Solve the model
   >>> solver = SolverFactory('ipopt')
   >>> res = solver.solve(m, tee=True)
   >>> m.fs.display()


For an example of optimizing a flowsheet containing a PySMO-trained polynomial regression surrogate model, see the `Autothermal reformer flowsheet optimization example <https://github.com/IDAES/examples-pse/blob/main/src/Examples/SurrMod/FlowsheetOptimization/PySMO_flowsheet_optimization.ipynb>`_.


References:
----------------
[1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide", https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801