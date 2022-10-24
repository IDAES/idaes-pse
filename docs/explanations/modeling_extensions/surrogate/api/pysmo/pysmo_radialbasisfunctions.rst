Generating Radial Basis Function (RBF) models with PySMO
==========================================================

.. note::
   The IDAES surrogate API is a wrapper around the original PySMO surrogate interface.


The *pysmo_surrogate.PysmoRBFTrainer* method has the capability to generate different types of RBF surrogates from data
based on the basis function selected. RBFs models are usually of the form
where

.. math::
    \begin{equation}
    y_{k}=\sum_{j=1}^{\Omega}w_{j}\psi\left(\Vert x_{k}-z_{j}\Vert\right)\qquad k=1,\ldots,m\quad\label{eq:RBF-expression}
    \end{equation}

where :math:`z_{j}` are basis function centers (in this case, the training data points), :math:`w_{j}` are the radial
weights associated with each center :math:`z_{j}`,  and  :math:`\psi` is a basis function transformation of the
Euclidean distances.

PySMO offers a range of basis function transformations :math:`\psi`, as shown in the table below.

.. list-table:: List of available RBF basis transformations, :math:`d = \parallel x_{k}-z_{j}\parallel`
   :widths: 25 25 50
   :header-rows: 1

   * - Transformation type
     - PySMO option name
     - :math:`\psi(d)`
   * - Linear
     - 'linear'
     -  :math:`d`
   * - Cubic
     - 'cubic'
     - :math:`d^{3}`
   * - Thin-plate spline
     - 'spline'
     - :math:`d^{2}\ln(d)`
   * - Gaussian
     - 'gaussian'
     - :math:`e^{\left(-d^{2}\sigma^{2}\right)}`
   * - Multiquadric
     - 'mq'
     - :math:`\sqrt{1+\left(\sigma d\right)^{2}}`
   * - Inverse mMultiquadric
     - 'imq'
     - :math:`1/{\sqrt{1+\left(\sigma d\right)^{2}}}`

Selection of parametric basis functions increase the flexibility of the radial basis function but adds an extra parameter (:math:`\sigma`)to be estimated.

Basic Usage
------------
To generate an RBF model with PySMO, the  *pysmo_surrogate.PysmoRBFTrainer* trainer is instantiated and initialized with the desired configuration arguments, and the training function ``train_surrogate`` is called:

.. code:: python

   # Required imports
   >>> from idaes.core.surrogate.pysmo_surrogate import PysmoRBFTrainer
   >>> import pandas as pd

   # Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   # Define the input and output labels
   input_labels = ['X1', 'X2']
   output_labels = ['Y1','Y2']

   # Create the RBF trainer object
   >>> rbf_trainer = PysmoRBFTrainer(input_labels=input_labels, output_labels=output_labels, training_dataframe = data_training)

   # Set PySMO options
   >>> rbf_trainer.config.basis_function = 'cubic'

   # Train the model
   >>> rbf_train = rbf_trainer.train_surrogate()

Configuration Options
----------------------
``PysmoRBFTrainer`` takes the following optional arguments:

.. list-table:: PySMO RBF options
   :widths: 20 20 60
   :header-rows: 1

   * - **Option**
     - Configuration argument
     - Description
   * - **basis_function**
     - *PysmoRBFTrainer.config.basis_function*
     - Option to specify the type of basis function to be used in the RBF model. Default is 'gaussian'.
   * - **regularization**
     - *PysmoRBFTrainer.config.regularization*
     - | Boolean argument which determines whether model regularization is applied during training. 
       | Default is True.
   * - **solution_method**
     - *PysmoRBFTrainer.config.solution_method*
     - | Method used to solve the parameter estimation problems for the RBF model:
       | BFGS ('BFGS'), maximum likelihood ('algebraic') or Pyomo least squares minimization ('pyomo'). 
       | Default is 'algebraic'.

Output
-------
The result of the ``pysmo_surrogate.PysmoRBFTrainer`` method is a python object containing information about the problem set-up, the optimal radial basis function weights :math:`w_{j}` and different error and quality-of-fit metrics such as the mean-squared-error (MSE) and the :math:`R^{2}` coefficient-of-fit.

Surrogate Visualization
------------------------
For visualizing PySMO-trained surrogates via parity and residual plots, see :ref:`Visualizing Surrogate Model Results<explanations/modeling_extensions/surrogate/plotting/index:Visualizing Surrogate Model Results>`.


Building the IDAES Surrogate Object
------------------------------------
To add the model to an IDAES flowsheet or generate model predictions, the SurrogateTrainer object needs to be transformed into an IDAES SurrogateObject object. This is done by calling ``PySMOSurrogate`` and passing the generated surrogate expressions, along with variable labels and optionally the bounds:

.. code:: python

   >>> surr = PysmoSurrogate(rbf_train, input_labels, output_labels, input_bounds)

The resulting ``PysmoSurrogate`` object may be saved to (and reloaded from) a JSON file; for details, see :ref:`the PySMO main page<explanations/modeling_extensions/surrogate/api/pysmo/index:PySMO: Python-based Surrogate Modeling Objects>`.

Prediction with *PysmoRBFTrainer* models
----------------------------------------------------------
Once the RBF model has been trained and the SurrogateObject object created, predictions for values at previously unsampled points *x_unsampled* can be evaluated by calling SurrogateObject's ``evaluate_surrogate()`` function on the unsampled points:

.. code:: python

   >>> y_unsampled = surr.evaluate_surrogate(x_unsampled)

Flowsheet Integration
----------------------
The result of the RBF training process can be passed directly into a process flowsheet using the IDAES ``SurrogateBlock`` option. The following code snippet demonstrates how a saved RBF model may be integrated directly into an IDAES flowsheet:

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
   >>> surrogates_obj =PysmoSurrogate.load_from_file('rbf_surrogate.json') # rbf_surrogate.json is an existing surrogate JSON file containing the rbf model
   >>> m.fs.surrogate.build_model(surrogates_obj, input_vars=inputs, output_vars=outputs)
   >>> m.fs.surrogate.pprint()

   # Set the variable Y1 as the model objective
   >>> m.fs.obj = Objective(expr=m.fs.Y1, sense=minimize)

   # Solve the model
   >>> solver = SolverFactory('ipopt')
   >>> res = solver.solve(m, tee=True)
   >>> m.fs.display()


For an example of optimizing a flowsheet containing a PySMO-trained RBF surrogate model, see the `Autothermal reformer flowsheet optimization example <https://github.com/IDAES/examples-pse/blob/main/src/Examples/SurrMod/FlowsheetOptimization/PySMO_flowsheet_optimization.ipynb>`_.


References:
----------------
[1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide", https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801

[2] Hongbing Fang & Mark F. Horstemeyer (2006): Global response approximation with radial basis functions, https://www.tandfonline.com/doi/full/10.1080/03052150500422294

[3] Rippa, S. (1999): An algorithm for selecting a good value for the parameter c in radial basis function interpolation, https://doi.org/10.1023/A:1018975909870

[4] Mongillo M.A. (2011) Choosing Basis Functions and Shape Parameters for Radial Basis Function Methods, https://doi.org/10.1137/11S010840