Generating Radial Basis Function models with PySMO
===================================================

The *pysmo.radial_basis_function* package has the capability to generate different types of RBF surrogates from data
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

Selection of parametric basis functions increase the flexibility of the radial basis function but adds an extra
parameter (:math:`\sigma`)to be estimated.

Basic Usage
------------
To generate an RBF model with PySMO, the  *pysmo.radial_basis_function* class is first initialized,
and then the function *rbf_training* is called on the initialized object:

.. doctest::

   Required imports
   >>> from idaes.surrogates.pysmo import radial_basis_function
   >>> import pandas as pd

   Load dataset from a csv file
   >>> xy_data = pd.read_csv('data.csv', header=None, index_col=0)

   Initialize the RadialBasisFunctions class, extract the list of features and train the model
   >>> rbf_init = radial_basis_function.RadialBasisFunctions(xy_data, *kwargs)
   >>> features = rbf_init.get_feature_vector()
   >>> rbf_fit = rbf_init.rbf_training()

* *xy_data* is a two-dimensional python data structure containing the input and output training data. The output values **MUST** be in the last column.

**Optional Arguments**

* *basis_function* - option to specify the type of basis function to be used in the RBF model. Default is 'gaussian'.
* *regularization* - boolean which determines whether regularization of the RBF model is considered. Default is True.

    - When regularization is turned on, the resulting model is a regressing RBF model.
    - When regularization is turned off, the resulting model is an interpolating RBF model.


*pysmo.radial_basis_function* Output
---------------------------------------
The result of *pysmo.radial_basis_function* (*rbf_fit* in above example) is a python object containing information
about the optimal radial basis function weights :math:`w_{j}` and different error and quality-of-fit metrics such as
the mean-squared-error (MSE) and the :math:`R^{2}` coefficient-of-fit. A Pyomo expression can be generated from the
object simply passing a list of variables into the function *rbf_generate_expression*:

.. doctest::

   Create a python list from the headers of the dataset supplied for training
   >>> list_vars = []
   >>> for i in features.keys():
   >>>     list_vars.append(features[i])

    Pass list to generate_expression function to obtain a Pyomo expression as output
   >>> print(rbf_fit.rbf_generate_expression(list_vars))

Similar to the *pysmo.polynomial_regression* module, the output of the *rbf_generate_expression* function can be passed
into an IDAES or Pyomo module as a constraint, objective or expression.

Prediction with *pysmo.radial_basis_function* models
-----------------------------------------------------
Once an RBF model has been trained, predictions for values at previously unsampled points *x_unsampled* can be evaluated by calling the
*rbf_predict_output()* function on the resulting Python object and the unsampled points:

.. doctest::

   Create a python list from the headers of the dataset supplied for training
   >>> y_unsampled = rbf_init.rbf_predict_output(rbf_fit, x_unsampled)



Further details about *pysmo.radial_basis_function* module may be found by consulting the examples or reading the paper [...]
