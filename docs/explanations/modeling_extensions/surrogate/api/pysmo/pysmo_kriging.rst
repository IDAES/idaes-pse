Generating Kriging Models with PySMO
===================================================

.. note::
   The IDAES surrogate API is a wrapper around the original PySMO surrogate interface.


*pysmo_surrogate.PysmoKrigingTrainer* trains ordinary Kriging models. Interpolating kriging models assume that the outputs :math:`\hat{y}\in\mathbb{R}^{m\times1}` are correlated and may be treated as a normally distributed stochastic process. For a set of input measurements :math:`X=\left\{ x_{1},x_{2},\ldots,x_{m}\right\} ;x_{i}\in\mathbb{R}^{n}`, the output :math:`\hat{y}` is modeled
as the sum of a mean :math:`\left(\mu\right)`$` and a Gaussian process error,

.. math::
    \begin{equation}
    \hat{y_{k}}=\mu+\epsilon\left(x_{k}\right)\qquad k=1,\ldots,m \qquad\quad
    \end{equation}

Kriging models assume that the errors in the outputs :math:`\epsilon` are correlated proportionally to the distance
between corresponding points,

.. math::

    \begin{equation}
    \text{cor}\left[\epsilon\left(x_{j}\right),\epsilon\left(x_{k}\right)\right]=\exp\left(-\sum_{i=1}^{n}\theta_{i}\mid x_{ij}-x_{ik}\mid^{\tau_{i}}\right)\qquad j,k=1,\ldots,m;\:\tau_{i}\in\left[1,2\right];\:\theta_{i}\geq0\qquad\quad\label{eq:corr-function}
    \end{equation}

The hyperparameters of the Kriging model :math:`\left(\mu,\sigma^{2},\theta_{1},\ldots,\theta_{n},\tau_{1},\ldots,\tau_{n}\right)`
are selected such that the concentrated log likelihood function is maximized.

Basic Usage
------------
To generate a Kriging model with PySMO, the  *pysmo_surrogate.PysmoKrigingTrainer* trainer is instantiated and initialized with the desired optional arguments, and the training function ``train_surrogate`` is called:

.. code:: python

   # Required imports
   >>> from idaes.core.surrogate.pysmo_surrogate import PysmoKrigingTrainer
   >>> import pandas as pd

   # Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   # Define the input and output labels
   input_labels = ['X1', 'X2']
   output_labels = ['Y1','Y2']

   # Create the Kriging trainer object
   >>> krg_trainer = PysmoKrigingTrainer(input_labels=input_labels, output_labels=output_labels, training_dataframe = data_training)

   # Set desired PySMO kriging options
   >>> krg_trainer.config.numerical_gradients = False

   # Train the model
   >>> krg_train = krg_trainer.train_surrogate()


Configuration Options
----------------------
``PysmoKrigingTrainer`` takes the following optional arguments:

.. list-table:: PySMO Kriging options
   :widths: 20 20 60
   :header-rows: 1

   * - **Option**
     - Configuration argument
     - Description
   * - **numerical_gradients**
     - *PysmoKrigingTrainer.config.numerical_gradients*
     - | Whether or not numerical gradients should be used in training. This choice determines the algorithm used to solve the problem.
       |    - True: The problem is solved with BFGS using central differencing with :math:`\Delta=10^{-6}` to evaluate numerical gradients.
       |    - False: The problem is solved with Basinhopping, a stochastic optimization algorithm.
   * - **regularization**
     - *PysmoKrigingTrainer.config.regularization*
     - | Boolean argument which determines whether model regularization is applied during training.
       |    - When regularization is turned on, the resulting model is a regressing kriging model.
       |    - When regularization is turned off, the resulting model is an interpolating kriging model.
       | Default is True.

Output
-------
The result of the ``pysmo_surrogate.PysmoKrigingTrainer`` method is a python object containing information about the optimal Kriging hyperparameters :math:`\left(\mu,\sigma^{2},\theta_{1},\ldots,\theta_{n}\right)` and different error and quality-of-fit metrics such as the mean-squared-error (MSE) and the :math:`R^{2}` coefficient-of-fit.

Surrogate Visualization
------------------------
For visualizing PySMO-trained surrogates via parity and residual plots, see :ref:`Visualizing Surrogate Model Results<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>`.

Building the IDAES Surrogate Object
------------------------------------
To add the model to an IDAES flowsheet or generate model predictions, the SurrogateTrainer object needs to be transformed into an IDAES SurrogateObject object. This is done by calling ``PySMOSurrogate`` and passing the generated surrogate expressions, along with variable labels and optionally the bounds:

.. code:: python

   >>> surr = PysmoSurrogate(krg_train, input_labels, output_labels, input_bounds)

The resulting ``PysmoSurrogate`` object may be saved to (and reloaded from) a JSON file; for details, see :ref:`the PySMO main page<explanations/modeling_extensions/surrogate/api/pysmo/index:PySMO: Python-based Surrogate Modeling Objects>`.

Prediction with *PysmoKrigingTrainer* models
----------------------------------------------------------
Once the Kriging model has been trained and the SurrogateObject object created, predictions for values at previously unsampled points *x_unsampled* can be evaluated by calling SurrogateObject's ``evaluate_surrogate()`` function on the unsampled points:

.. code:: python

   >>> y_unsampled = surr.evaluate_surrogate(x_unsampled)



Flowsheet Integration
----------------------
The final Kriging model can be passed into a process flowsheet using the IDAES ``SurrogateBlock`` option. The following code snippet demonstrates how a saved Kriging model may be integrated directly into an IDAES flowsheet:

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
   >>> surrogates_obj =PysmoSurrogate.load_from_file('krg_surrogate.json') # krg_surrogate.json is an existing surrogate JSON file containing a kriging model
   >>> m.fs.surrogate.build_model(surrogates_obj, input_vars=inputs, output_vars=outputs)
   >>> m.fs.surrogate.pprint()

   # Set the variable Y1 as the model objective
   >>> m.fs.obj = Objective(expr=m.fs.Y1, sense=minimize)

   # Solve the model
   >>> solver = SolverFactory('ipopt')
   >>> res = solver.solve(m, tee=True)
   >>> m.fs.display()


For an example of optimizing a flowsheet containing a PySMO-trained Kriging surrogate model, see the `Autothermal reformer flowsheet optimization example <https://github.com/IDAES/examples-pse/blob/main/src/Examples/SurrMod/FlowsheetOptimization/PySMO_flowsheet_optimization.ipynb>`_.


References:
----------------
[1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide", https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801

[2] D. R. Jones, A taxonomy of global optimization methods based on response surfaces, Journal of Global Optimization, https://link.springer.com/article/10.1023%2FA%3A1012771025575