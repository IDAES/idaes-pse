Generating Kriging Models with PySMO
===================================================

The *pysmo.kriging* trains Ordinary Kriging models. Interpolating kriging models assume that the outputs :math:`\hat{y}\in\mathbb{R}^{m\times1}`
are correlated and may be treated as a normally distributed stochastic process. For a set of input measurements
:math:`X=\left\{ x_{1},x_{2},\ldots,x_{m}\right\} ;x_{i}\in\mathbb{R}^{n}`, the output :math:`\hat{y}` is modeled
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
To generate a Kriging model with PySMO, the  *pysmo.kriging* class is first initialized,
and then the function *training* is called on the initialized object:

.. code:: python

   # Required imports
   >>> from idaes.surrogate.pysmo import kriging
   >>> import pandas as pd

   # Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   # Initialize the KrigingModel class, extract the list of features and train the model
   >>> krg_init = kriging.KrigingModel(xy_data, *kwargs)
   >>> features = krg_init.get_feature_vector()
   >>> krg_init.training()

* *xy_data* is a two-dimensional python data structure containing the input and output training data. The output values **MUST** be in the last column.

**Optional Arguments**

* *numerical_gradients*: Whether or not numerical gradients should be used in training. This choice determines the algorithm used to solve the problem.

    - True: The problem is solved with BFGS using central differencing with :math:`\Delta=10^{-6}` to evaluate numerical gradients.
    - False: The problem is solved with Basinhopping, a stochastic optimization algorithm.

* *regularization* - Boolean option which determines whether or not regularization is considered during Kriging training. Default is True.

    - When regularization is turned on, the resulting model is a regressing kriging model.
    - When regularization is turned off, the resulting model is an interpolating kriging model.

*pysmo.kriging* Output
---------------------------------------
The result of *pysmo.kriging* is a python object containing information
about the optimal Kriging hyperparameters :math:`\left(\mu,\sigma^{2},\theta_{1},\ldots,\theta_{n}\right)`
and different error and quality-of-fit metrics such as the mean-squared-error (MSE) and the :math:`R^{2}` coefficient-of-fit.
A Pyomo expression can be generated from the object simply passing a list of variables into the function
*generate_expression*:

.. code:: python

   # Create a python list from the headers of the dataset supplied for training
   >>> list_vars = []
   >>> for i in features.keys():
   >>>     list_vars.append(features[i])

   # Pass list to generate_expression function to obtain a Pyomo expression as output
   >>> print(krg_init.generate_expression(list_vars))

Similar to the *pysmo.polynomial_regression* module, the output of the *generate_expression* function can be passed
into an IDAES or Pyomo module as a constraint, objective or expression.

Prediction with *pysmo.kriging* models
-----------------------------------------------------
Once a Kriging model has been trained, predictions for values at previously unsampled points *x_unsampled* can be evaluated by calling the
*predict_output()* function on the unsampled points:

.. code:: python

   # Create a python list from the headers of the dataset supplied for training
   >>> y_unsampled = kriging_init.predict_output(x_unsampled)

Further details about *pysmo.kriging* module may be found by consulting the examples or reading the paper [...]


.. module:: idaes.surrogate.pysmo.kriging

Available Methods
------------------

.. autoclass:: idaes.surrogate.pysmo.kriging.KrigingModel
    :members: __init__, get_feature_vector, training, predict_output, r2_calculation, generate_expression

References:
----------------
[1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide", https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801

[2] D. R. Jones, A taxonomy of global optimization methods based on response surfaces, Journal of Global Optimization, https://link.springer.com/article/10.1023%2FA%3A1012771025575