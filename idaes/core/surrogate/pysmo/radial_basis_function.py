#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################


# Imports from the python standard library
from __future__ import division, print_function
from builtins import str
import os.path
import warnings

# Imports from third parties
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pickle
from pyomo.environ import *
import scipy.optimize as opt
from six import string_types

# Imports from IDAES namespace
from idaes.core.surrogate.pysmo.sampling import FeatureScaling as fs

__author__ = "Oluwamayowa Amusat"

"""
The purpose of this file is to perform radial basis functions in Pyomo.
"""


class FeatureScaling:
    """

    A class for scaling and unscaling input and output data. The class contains two main methods: ``data_scaling_minmax`` and ``data_unscaling_minmax``
    """

    def __init__(self):
        pass

    @staticmethod
    def data_scaling_minmax(data):
        """
        ``data_scaling_minmax`` performs column-wise min-max scaling on the input dataset.

        Args:
            data : The input data set to be scaled. Must be a numpy array or dataframe.

        Returns:
            (tuple): tuple containing:
                - **scaled_data**  : A 2-D Numpy Array containing the scaled data. All array values will be between [0, 1].
                - **data_minimum** : A 2-D row vector containing the column-wise minimums of the input data.
                - **data_maximum** : A 2-D row vector containing the column-wise maximums of the input data.

        Raises:
            TypeError: Raised when the input data is not a numpy array or dataframe

        """
        # Confirm that data type is an array or DataFrame
        if isinstance(data, np.ndarray):
            input_data = data
            data_headers = []
        elif isinstance(data, pd.DataFrame):
            input_data = data.values
            data_headers = data.columns.values.tolist()
        else:
            raise TypeError(
                "original_data_input: Pandas dataframe or numpy array required."
            )

        if input_data.ndim == 1:
            input_data = input_data.reshape(len(input_data), 1)
        data_minimum = np.min(input_data, axis=0)
        data_maximum = np.max(input_data, axis=0)
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        scaled_data = (input_data - data_minimum) / scale
        # scaled_data = (input_data - data_minimum)/(data_maximum - data_minimum)
        data_minimum = data_minimum.reshape(1, data_minimum.shape[0])
        data_maximum = data_maximum.reshape(1, data_maximum.shape[0])

        if len(data_headers) > 0:
            scaled_data = pd.DataFrame(scaled_data, columns=data_headers)
        return scaled_data, data_minimum, data_maximum

    @staticmethod
    def data_unscaling_minmax(x_scaled, x_min, x_max):
        """

        ``data_unscaling_minmax`` performs column-wise un-scaling on the a minmax-scaled input dataset.

        Args:
            x_scaled (NumPy Array)  : Data to be un-scaled. Data values should be between 0 and 1.
            x_min (NumPy vector)    : :math:`n \\times 1` vector containing the actual minimum value for each column. Must contain same number of elements as the number of columns in x_scaled.
            x_max (NumPy vector)    : :math:`n \\times 1` vector vector containing the actual minimum value for each column. Must contain same number of elements as the number of columns in x_scaled.

        Returns:
            NumPy Array : A 2-D numpy array containing the scaled data, :math:`x_{min} + x_{scaled} * (x_{max} - x_{min})`

        Raises:
            IndexError: Raised when the dimensions of the arrays are inconsistent.

        """
        # Check if it can be evaluated. Will return index error if dimensions are wrong
        if x_scaled.ndim == 1:  # Check if 1D, and convert to 2D if required.
            x_scaled = x_scaled.reshape(len(x_scaled), 1)
        if (x_scaled.shape[1] != x_min.size) or (x_scaled.shape[1] != x_max.size):
            raise IndexError("Dimensionality problems with data for un-scaling.")
        unscaled_data = x_min + x_scaled * (x_max - x_min)
        return unscaled_data

    # @staticmethod
    # def data_scaling_standardization(data):
    #     # Confirm that data type is an array or DataFrame
    #     if isinstance(data, np.ndarray):
    #         input_data = data
    #     elif isinstance(data, pd.DataFrame):
    #         input_data = data.values
    #     else:
    #         raise TypeError('original_data_input: Pandas dataframe or numpy array required.')
    #
    #     if input_data.ndim == 1:
    #         input_data = input_data.reshape(len(input_data), 1)
    #
    #     data_mean = np.mean(input_data, axis=0)
    #     data_stdev = np.std(input_data, axis=0)
    #     scaled_data = (input_data - data_mean) / data_stdev
    #     data_mean = data_mean.reshape(1, data_mean.shape[0])
    #     data_stdev = data_stdev.reshape(1, data_stdev.shape[0])
    #     return scaled_data, data_mean, data_stdev


class RadialBasisFunctions:
    """
    The RadialBasisFunctions class generates a radial basis function fitting for a training data set.

    The class must first be initialized by calling **RadialBasisFunctions**. Regression is then carried out by calling the method ``training``.

    For a given dataset with n features :math:`x_{1},\ldots,x_{n}`, RadialBasisFunctions is able to consider six types of basis transformations:
        - Linear ('linear')
        - Cubic ('cubic')
        - Gaussian ('gaussian')
        - Multiquadric ('mq')
        - Inverse multiquadric ('imq')
        - Thin-plate spline ('spline')

    ``training`` selects the best hyperparameters (regularization parameter :math:`\lambda` and shape parameter :math:`\sigma`, where necessary) by evaluating the leave-one-out cross-validation error for each (:math:`\lambda,\sigma`) pair.

    It should be noted that the all the training points are treated as centres for the RBF, resulting in a square system.

    **Example:**

    .. code-block:: python

         # Initialize the class
        >>> d = RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True))
        >>> p = d.get_feature_vector()

        # Train RBF model and predict output for an test data x_test
        >>> d.training()
        >>> predictions = d.predict_output(x_test)

    Args:
        XY_data (Numpy Array or Pandas Dataframe)             : The dataset for RBF training. **XY_data** is expected to contain the features (X) and output (Y) data, with the output values (Y) in the last column.

    Further details about the optional inputs may be found under the ``__init__`` method.

    """

    def __init__(
        self,
        XY_data,
        basis_function=None,
        solution_method=None,
        regularization=None,
        fname=None,
        overwrite=False,
    ):
        """

        Initialization of **RadialBasisFunctions** class.

        Args:
            XY_data (Numpy Array or Pandas Dataframe): The dataset for RBF training. **XY_data** is expected to contain feature and output information, with the output values (y) in the last column.

        Keyword Args:
            basis_function(str): The basis function transformation to be applied to the training data. Two classes of basis transformations are available for selection:

                - Fixed basis transformations, which require no shape parameter :math:`\sigma` :

                    (a)  'cubic'      : Cubic basis transformation
                    (b) 'linear'      : Linear basis transformation
                    (c) 'spline'      : Thin-plate spline basis transformation

                - Parametric basis transformations which require a shape parameter:

                    (a) 'gaussian'    : Gaussian basis transformation (Default)
                    (b) 'mq'          : Multiquadric basis transformation
                    (c) 'imq'         : Inverse multiquadric basis transformation


            solution_method(str): The method to be used for solving the RBF least squares optimization problem. Three options are available:

                (a) 'algebraic'  : The explicit algebraic method solves the least squares problem using linear algebra.
                (b) 'BFGS'       : This approach solves the least squares problem using SciPy's BFGS algorithm.
                (c) 'pyomo'      : This option solves the optimization problem in Pyomo with IPOPT as solver. This is the default.

            regularization(bool): This option determines whether or not the regularization parameter :math:`\lambda` is considered during RBF fitting. Default setting is True.


        Returns:
            **self** object with the input information

        Raises:
            ValueError: The input dataset is of the wrong type (not a NumPy array or Pandas Dataframe)

            Exception:
                * **basis_function** entry is not valid.
            Exception:
                * **solution_method** is not 'algebraic', 'pyomo' or 'bfgs'.
            Exception:
                - :math:`\lambda` is not boolean.

        **Example:**

        .. code-block:: python

            # Specify the gaussian basis transformation
            >>> d = RadialBasisFunctions(XY_data, basis_function='gaussian')

        """
        if not isinstance(overwrite, bool):
            raise Exception("overwrite must be boolean.")
        self.overwrite = overwrite
        if fname is None:
            fname = "solution.pickle"
            self.filename = "solution.pickle"
        elif (
            not isinstance(fname, str)
            or os.path.splitext(fname)[-1].lower() != ".pickle"
        ):
            raise Exception(
                'fname must be a string with extension ".pickle". Please correct.'
            )
        if (
            os.path.exists(fname) and overwrite is True
        ):  # Explicit overwrite done by user
            print(
                "Warning:",
                fname,
                "already exists; previous file will be overwritten.\n",
            )
            self.filename = fname
        elif os.path.exists(fname) and overwrite is False:  # User is not overwriting
            self.filename = (
                os.path.splitext(fname)[0]
                + "_v_"
                + pd.Timestamp.today().strftime("%m-%d-%y_%H%M%S")
                + ".pickle"
            )
            print(
                "Warning:",
                fname,
                'already exists; results will be saved to "',
                self.filename,
                '".\n',
            )
            # self.filename = 'solution.pickle'
        elif os.path.exists(fname) is False:
            self.filename = fname

        # Check data types and shapes
        if isinstance(XY_data, pd.DataFrame):
            xy_data = XY_data.values
            self.x_data_columns = list(XY_data.columns)[:-1]
        elif isinstance(XY_data, np.ndarray):
            xy_data = XY_data
            self.x_data_columns = list(range(XY_data.shape[1] - 1))
        else:
            raise ValueError('Pandas dataframe or numpy array required for "XY_data".')

        self.x_data_unscaled = xy_data[:, :-1]
        self.y_data_unscaled = xy_data[:, -1].reshape(xy_data.shape[0], 1)
        xy_data_scaled, self.data_min, self.data_max = fs.data_scaling_minmax(XY_data)
        x_data_scaled = xy_data_scaled[:, :-1]
        y_data_scaled = xy_data_scaled[:, -1]
        self.x_data = x_data_scaled.reshape(self.x_data_unscaled.shape)
        self.y_data = y_data_scaled.reshape(self.y_data_unscaled.shape)
        self.centres = xy_data_scaled[:, :-1]

        if solution_method is None:
            solution_method = "algebraic"
            self.solution_method = solution_method
            print("Default parameter estimation method is used.")
        elif not isinstance(solution_method, string_types):
            raise Exception("Invalid solution method. Must be of type <str>.")
        elif (
            (solution_method.lower() == "algebraic")
            or (solution_method.lower() == "pyomo")
            or (solution_method.lower() == "bfgs")
        ):
            solution_method = solution_method.lower()
            self.solution_method = solution_method
        else:
            raise Exception(
                'Invalid solution method entered. Select one of ALGEBRAIC (solution_method="algebraic") , L-BFGS (solution_method="bfgs") or Pyomo optimization (solution_method="pyomo") methods. '
            )
        print("\nParameter estimation method: ", self.solution_method)

        if basis_function is None:
            basis_function = "gaussian"
            self.basis_function = basis_function
            print("Gaussian basis function is used.")
        elif not isinstance(basis_function, string_types):
            raise Exception("Invalid basis_function. Must be of type <str>.")
        elif (
            (basis_function.lower() == "linear")
            or (basis_function.lower() == "cubic")
            or (basis_function.lower() == "gaussian")
            or (basis_function.lower() == "mq")
            or (basis_function.lower() == "imq")
            or (basis_function.lower() == "spline")
        ):
            basis_function = basis_function.lower()
            self.basis_function = basis_function
        else:
            raise Exception(
                "Invalid basis function entered. See manual for available options. "
            )
        print("Basis function: ", self.basis_function)

        if regularization is None:
            regularization = True
            self.regularization = regularization
        elif not isinstance(regularization, bool):
            raise Exception("Invalid basis_function. Must be boolean")
        elif (regularization is True) or (regularization is False):
            self.regularization = regularization
        print("Regularization done: ", self.regularization)

        # Results
        self.weights = None
        self.sigma = None
        self.regularization_parameter = None
        self.rmse = None
        self.output_predictions = None
        self.condition_number = None
        self.R2 = None
        self.x_data_min = None
        self.x_data_max = None
        self.y_data_min = None
        self.y_data_max = None
        self.solution_status = None

    def r2_distance(self, c):
        """
        The function r2_distance calculates Euclidean distance from the point or array c.

        """
        dist = self.x_data - c
        l2_distance = np.sqrt(np.sum((dist**2), axis=1))
        return l2_distance

    @staticmethod
    def gaussian_basis_transformation(x, shape_parameter):
        """
        The function gaussian_basis_transformation returns the element-by-element Gaussian transformation of the input data x.

        Args:
            x(NumPy Array): Input data to be transformed
            shape parameter(float): Shape parameter of the Gaussian function

        Returns:
            x_mod(NumPy Array): Gaussian transformation of the input data x

        Examples:
            Gaussian of numbers from 0 to 2 for a shape parameter of 2:
                [In]>>  rbf.RadialBasisFunctions.gaussian_basis_transformation(np.arange(3), 2)
                [Out]>> array([1.00000000e+00, 1.83156389e-02, 1.12535175e-07])

        For more information, see Hongbing Fang & Mark F. Horstemeyer (2006): Global response approximation with radial basis functions
        https://www.tandfonline.com/doi/full/10.1080/03052150500422294

        """
        x_mod = np.exp(-1 * ((x * shape_parameter) ** 2))
        return x_mod

    @staticmethod
    def linear_transformation(x):
        """
        The function linear_transformation returns the element-by-element linear transformation of the input data x.

        Args:
            x(NumPy Array): Input data to be transformed

        Returns:
            x_mod(NumPy Array): Linear transformation of the input data x, x_mod = x

        Examples:
            Linear transformation of 0, 1 and 2:
                [In]>>  rbf.RadialBasisFunctions.linear_transformation(np.arange(3))
                [Out]>> array([0, 1, 2], dtype=int32)

        For more information, see Hongbing Fang & Mark F. Horstemeyer (2006): Global response approximation with radial basis functions
        https://www.tandfonline.com/doi/full/10.1080/03052150500422294

        """
        x_mod = x**1
        return x_mod

    @staticmethod
    def cubic_transformation(x):
        """
        The function cubic_transformation returns the element-by-element cubic transformation of the input data x.

        Args:
            x(NumPy Array): Input data to be transformed

        Returns:
            x_mod(NumPy Array): Cubic transformation of the input data x, x_mod = (x ** 3)

        Examples:
            Cubic transformation of 0, 1 and 2:
                [In]>>  rbf.RadialBasisFunctions.cubic_transformation(np.arange(3))
                [Out]>> array([0, 1, 8], dtype=int32)

        For more information, see Hongbing Fang & Mark F. Horstemeyer (2006): Global response approximation with radial basis functions
        https://www.tandfonline.com/doi/full/10.1080/03052150500422294
        """
        x_mod = x**3
        return x_mod

    @staticmethod
    def multiquadric_basis_transformation(x, shape_parameter):
        """
        The function multiquadric_basis_transformation returns the element-by-element Multi-quadric transformation of the input data x.

        Args:
            x(NumPy Array): Input data to be transformed
            shape parameter(float): Shape parameter of the Multiquadric function

        Returns:
            x_mod(NumPy Array): Multiquadric transformation of the input data x; x_mod = sqrt[(1 + (c.x)**2)] where c = shape parameter

        Examples:
            Multiquadric transformation of numbers from 0 to 2 for a shape parameter of 2:
                [In]>>  rbf.RadialBasisFunctions.multiquadric_basis_transformation(np.arange(3), 2)
                [Out]>> array([1.        , 2.23606798, 4.12310563])

        For more information, see
        (1) Hongbing Fang & Mark F. Horstemeyer (2006): Global response approximation with radial basis functions
        https://www.tandfonline.com/doi/full/10.1080/03052150500422294

        (2) Santana-Quintero L.V., Montaño A.A., Coello C.A.C. (2010) A Review of Techniques for Handling Expensive Functions in Evolutionary Multi-Objective Optimization.
        In: Tenne Y., Goh CK. (eds) Computational Intelligence in Expensive Optimization Problems.
        https://link.springer.com/chapter/10.1007/978-3-642-10701-6_2

        """
        x_mod = np.sqrt(((x * shape_parameter) ** 2) + 1)
        # x_mod = np.sqrt(x**2 + shape_parameter**2)  # Alternative implementation
        return x_mod

    @staticmethod
    def inverse_multiquadric_basis_transformation(x, shape_parameter):
        """
        The function inverse_multiquadric_basis_transformation returns the element-by-element  inverse multiquadric transformation of the input data x.
        Direct inverse of the multiquadric basis transformation

        Args:
            x(NumPy Array): Input data to be transformed
            shape_parameter(float): Shape parameter of the inverse multiquadric function

        Returns:
            x_mod(NumPy Array): Inverse multiquadric transformation of the input data x; x_mod =  1 / sqrt[(1 + (c.x)**2)] where c = shape parameter

        Examples:
            Inverse multiquadric transformation of numbers from 0 to 2 for a shape parameter of 2:
                [In]>>  rbf.RadialBasisFunctions.inverse_multiquadric_basis_transformation(np.arange(3), 2)
                [Out]>>array([1.        , 0.4472136 , 0.24253563])

        For more information, see
        (1) Hongbing Fang & Mark F. Horstemeyer (2006): Global response approximation with radial basis functions
        https://www.tandfonline.com/doi/full/10.1080/03052150500422294

        (2) Santana-Quintero L.V., Montaño A.A., Coello C.A.C. (2010) A Review of Techniques for Handling Expensive Functions in Evolutionary Multi-Objective Optimization.
        In: Tenne Y., Goh CK. (eds) Computational Intelligence in Expensive Optimization Problems.
        https://link.springer.com/chapter/10.1007/978-3-642-10701-6_2

        """
        x_mod = 1 / np.sqrt(((x * shape_parameter) ** 2) + 1)
        # x_mod = 1 / (np.sqrt(x ** 2 + shape_parameter ** 2))  # Alternative implementation
        return x_mod

    @staticmethod
    def thin_plate_spline_transformation(x):
        """
        The function thin_plate_spline_transformation returns the element-by-element spline transformation of the input data x.
        Direct inverse of the multiquadric basis transformation

        Args:
            x(NumPy Array): Input data to be transformed

        Returns:
            x_mod(NumPy Array): Spline transformation of the input data x; x_mod =  (x**2).ln(x)

        Except:
            RuntimeWarning: thrown up when ln(x)=0

        Examples:
            Inverse multiquadric transformation of numbers from 1 and 2:
                [In]>>  rbf.RadialBasisFunctions.thin_plate_spline_transformation(np.arange(1,3))
                [Out]>> array([0, 2.77258872])

        For more information, see
        (1) Hongbing Fang & Mark F. Horstemeyer (2006): Global response approximation with radial basis functions
        https://www.tandfonline.com/doi/full/10.1080/03052150500422294

        (2) Santana-Quintero L.V., Montaño A.A., Coello C.A.C. (2010) A Review of Techniques for Handling Expensive Functions in Evolutionary Multi-Objective Optimization.
        In: Tenne Y., Goh CK. (eds) Computational Intelligence in Expensive Optimization Problems.
        https://link.springer.com/chapter/10.1007/978-3-642-10701-6_2

        """
        # x_mod = (x ** 2) * np.log(x)
        # x_mod = np.nan_to_num(x_mod)
        with np.errstate(
            divide="ignore"
        ):  # catch division warnings in log function due to log(0)~=0
            log_x = np.log(x)
        with np.errstate(
            invalid="ignore"
        ):  # catch invalid warnings due to - Inf * 0 evaluations
            x_mod = (x**2) * log_x
        x_mod = np.nan_to_num(x_mod)
        return x_mod

    def basis_generation(self, r):
        """
        The function basis_generation converts the input data to the requisite basis specified by the user.
        This is done in two steps:

        1. The Euclidean distance from each of the points to each of the RBF centres is calculated by calling the r2_distance function.
        2. The distances evaluated in step 1 are transformed to the relevant basis selected by the user.

        Args:
            self(NumPy Array): contains, among other things, the input data
            r(float)        : The shape parameter required for the Gaussian, Multiquadric and Inverse multiquadric transformations.

        Returns:
            x_transformed(NumPy Array): Array of transformed data based on user-defined transformation function

        """

        basis_functions = np.zeros((self.x_data.shape[0], self.centres.shape[0]))
        for i in range(0, self.centres.shape[0]):
            basis_functions[:, i] = self.r2_distance(self.centres[i, :])

        # Initialization of x_transformed
        x_transformed = np.zeros((basis_functions.shape[0], basis_functions.shape[1]))

        if self.basis_function == "gaussian":
            x_transformed = self.gaussian_basis_transformation(basis_functions, r)
        elif self.basis_function == "linear":
            x_transformed = self.linear_transformation(basis_functions)
        elif self.basis_function == "cubic":
            x_transformed = self.cubic_transformation(basis_functions)
        elif self.basis_function == "mq":
            x_transformed = self.multiquadric_basis_transformation(basis_functions, r)
        elif self.basis_function == "imq":
            x_transformed = self.inverse_multiquadric_basis_transformation(
                basis_functions, r
            )
        elif self.basis_function == "spline":
            x_transformed = self.thin_plate_spline_transformation(basis_functions)
        return x_transformed

    @staticmethod
    def cost_function(theta, x, y):
        """
        This function is an implementation of the cost function for linear regression with BFGS:
                cost = [sum of square errors over m samples / (2 * m)]

        This is the objective function for the BFGS optimization problem.

        Args:
            theta        : polynomial coefficients/weights,  (n x 1) in size
            x            : array of features, (m x n) in size
            y            : actual output vector, size (m x 1)

        Returns:
            cost_value   : the cost value for the fit, the objective value of the optimization problem

        """
        y = y.reshape(y.shape[0], 1)
        y_prediction = np.matmul(x, theta)
        y_prediction = y_prediction.reshape(y_prediction.shape[0], 1)
        cost_value = (0.5 / x.shape[0]) * (np.sum((y - y_prediction) ** 2))
        return cost_value

    @staticmethod
    def gradient_function(theta, x, y):
        """
        This function is an implementation of the gradient function for linear regression:
                if
                    cost = [(A.x - y)^2 / 2m]
                then
                    gradient = [((A.x - y)* A) / m]

        This is the gradient function supplied to the BFGS optimization algorithm.

        Args:
            theta        : polynomial coefficients/weights,  (n x 1) in size
            x            : array of features, (m x n) in size
            y            : actual output vector, size (m x 1)

        Returns:
            grad_values   : the cost gradients for the fit, size (n x 1)

        """
        y = y.reshape(y.shape[0], 1)
        y_prediction = np.matmul(x, theta)
        y_prediction = y_prediction.reshape(y_prediction.shape[0], 1)
        t1 = (y_prediction - y) * x
        grad_values = (1 / x.shape[0]) * np.sum(t1, axis=0)
        grad_values = grad_values.reshape(
            theta.size,
        )
        return grad_values

    def bfgs_parameter_optimization(self, x, y):
        """
        This function performs parameter optimization using scipy's BFGS algorithm.
        It takes in the functions pre-defined functions cost_function and gradient_function as the cost and gradient functions.

        Args:
            x            : array of features, (m x n) in size
            y            : actual output vector, size (m x 1)

        Initialization:
            The initial weights are set to zero,
            init_phi = 0

        Returns:
            phi: The optimal linear regression weights found

        """
        init_phi = np.zeros((x.shape[1], 1))
        other_args = (x, y)
        phi = opt.fmin_bfgs(
            self.cost_function,
            init_phi,
            fprime=self.gradient_function,
            args=other_args,
            disp=False,
            gtol=1e-20,
        )
        return phi

    @staticmethod
    def explicit_linear_algebra_solution(x, y):
        """
        The function finds the explicit linear algebra solution to the reqularized problem (X+yI).A = B

        If:
            (x + yI).A = B,

        Then:
            A = inv(X + yI) * B

        where y is the regularization parameter and I is the identity matrix.

        Numpy's inverse and pseudoinverse functions have been used to calculate the inverse here.

        Args:
            x            : regularized array of features (x + yI), (m x n) in size
            y            : actual output vector, size (m x 1)

        Returns:
            phi          : optimal linear regression weights A

        For more details, see to Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide", https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801

        """
        # Find matrix inverse. Use pseudo-inverse if inverse is not available
        try:
            inverse_x = np.linalg.inv(x)
        except np.linalg.LinAlgError as LAE:
            inverse_x = np.linalg.pinv(x)

        phi = np.matmul(inverse_x, y)
        return phi

    @staticmethod
    def pyomo_optimization(x, y):
        """
        Pyomo implementation of least squares optimization problem:

                        Minimize cost = (y' - y) ^ 2
                        subject to: y' = Ax

        The problem is solved within Pyomo's framework using IPOPT as solver.

        Args:
            x            : array of features, (m x n) in size
            y            : actual output vector, size (m x 1)

        Returns:
            phi: The optimal linear regression weights found

        """
        model = ConcreteModel()

        pd.set_option("display.precision", 64)
        x_data = pd.DataFrame(x)
        y_data = pd.DataFrame(y)

        model.M = Set(initialize=x_data.index.values)  # Rows indices passed into set
        model.N = Set(
            initialize=x_data.columns.values
        )  # x column indices passed into set
        model.P = Set(
            initialize=y_data.columns.values
        )  # y column index passed into set

        model.x = Param(model.M, model.N, initialize=x_data.stack().to_dict())
        model.y_real = Param(model.M, model.P, initialize=y_data.stack().to_dict())

        # Define variables
        model.theta = Var(model.N, initialize=0, domain=Reals)
        model.y_predictions = Var(
            model.M, model.P, initialize=y_data.stack().to_dict(), domain=Reals
        )

        # constraint y_p = theta.X
        def xy_product(model, i, k):
            return model.y_predictions[i, k] == sum(
                model.theta[j] * model.x[i, j] for j in model.N for k in model.P
            )

        model.x_theta_product = Constraint(
            model.M, model.P, rule=xy_product, doc="Predicted value calc: y = hx"
        )

        # Cost function - RMSE
        def model_rms_error(model):
            cost_value = (1 / len(model.M)) * sum(
                ((model.y_real[i, k] - model.y_predictions[i, k]) ** 2)
                for i in model.M
                for k in model.P
            )
            return cost_value

        model.prediction_error = Objective(
            rule=model_rms_error, sense=minimize, doc="Minimum RMSE error"
        )

        instance = model
        opt = SolverFactory("ipopt")
        opt.options["max_iter"] = 10000
        opt.options["acceptable_tol"] = 1e-30
        # model.pprint()
        result = opt.solve(instance)  # , tee=True)
        # model.display()

        # Convert theta variable into numpy array
        phi = np.zeros((len(instance.theta), 1))
        iterator = 0
        for s in instance.N:
            phi[iterator, 0] = instance.theta[s].value
            iterator += 1
        return phi

    @staticmethod
    def error_calculation(phi, x, y_data):
        """
        This function calculates the SSE and RMSE errors between the actual and predicted output values,
             ss_error = sum of squared errors / number of samples
             rmse_error = sqrt(sum of squared errors / number of samples)

        Args:
            phi             : weight vector obtained by optimization
            x               : vector of input features
            y               : actual output values

        Returns:
            ss_error        : The average sum of squared errors
            rmse_error      : The root-mean-squared error (RMSE)
            y_prediction    : Predicted values of y, y_prediction = phi.x

        """
        y_prediction = np.matmul(x, phi)
        ss_error = (1 / y_data.shape[0]) * (np.sum((y_data - y_prediction) ** 2))
        rmse_error = np.sqrt(ss_error)
        return ss_error, rmse_error, y_prediction

    @staticmethod
    def r2_calculation(y_true, y_predicted):
        """
        ``r2_calculation`` returns the :math:`R^{2}` as a measure-of-fit between the true and predicted values of the output variable.

        Args:
            y_true(NumPy Array)             : Vector of actual values of the output variable
            y_predicted(NumPy Array)        : Vector of predictions for the output variable based on the surrogate

        Returns:
            float                           : :math:`R^{2}` measure-of-fit between actual and predicted data

        """
        y_true = y_true.reshape(y_true.shape[0], 1)
        y_predicted = y_predicted.reshape(y_predicted.shape[0], 1)
        input_y_mean = np.mean(y_true, axis=0)
        ss_total = np.sum((y_true - input_y_mean) ** 2)
        ss_residual = np.sum((y_predicted - y_true) ** 2)
        r_square = 1 - (ss_residual / ss_total)
        return r_square

    def loo_error_estimation_with_rippa_method(self, sigma, lambda_reg):
        """
        The function loo_error_estimation_with_rippa_method implements the leave-one-out cross-validation (LOOCV) error for square systems

        The LOOCV error is calculated analytically using Rippa's equation:
            Error[k] = alpha[k] / inv(A[kk]),

        where:
            Error[k] = error on leaving out a particular sample k
            alpha[k] = kth radial weight based on data (kth coefficient of full data interpolation)
            A[kk] = kth diagonal element of data matrix.

        Args:
            self                          : contains, among other things, the input data
            sigma(float)                  : shape parameter for the parametric bases (Gaussian, Multiquadric, Inverse multiquadric)
            lambda_reg(float)             : regularization parameter

        Returns:
            condition_number_pure           : condition number of transformed matrix generated from the input data before regularization
            condition_number_regularized    : condition number of transformed matrix generated from the input data after regularization
            loo_error_estimate              : norm of the leave-one-out cross-validation error matrix

        For more information, see
        (1) Rippa, S. (1999) Advances in Computational Mathematics
        https://doi.org/10.1023/A:1018975909870

        (2) Mongillo M.A. (2011) Choosing Basis Functions and Shape Parameters for Radial Basis Function Methods
        https://doi.org/10.1137/11S010840

        """
        x_transformed = self.basis_generation(sigma)
        condition_number_pure = np.linalg.cond(x_transformed)

        x_regularized = x_transformed + (
            lambda_reg * np.eye(x_transformed.shape[0], x_transformed.shape[1])
        )
        condition_number_regularized = np.linalg.cond(x_regularized)

        y_train = self.y_data.reshape(self.y_data.shape[0], 1)

        # SOLVE RADIAL WEIGHTS FOR FULL X DATA
        if self.solution_method == "algebraic":
            radial_weights = self.explicit_linear_algebra_solution(
                x_regularized, y_train
            )
        elif self.solution_method == "pyomo":
            radial_weights = self.pyomo_optimization(x_regularized, y_train)
        elif self.solution_method == "bfgs":
            radial_weights = self.bfgs_parameter_optimization(x_regularized, y_train)
        radial_weights = radial_weights.reshape(radial_weights.shape[0], 1)

        # Evaluate loo-estimate with Rippa formula
        inverse_matrix = np.diag(np.linalg.pinv(x_regularized))
        error_vector = radial_weights.reshape(radial_weights.shape[0], 1) / (
            inverse_matrix.reshape(inverse_matrix.shape[0], 1)
        )
        loo_error_estimate = np.linalg.norm(error_vector)
        return condition_number_pure, condition_number_regularized, loo_error_estimate

    def leave_one_out_crossvalidation(self):
        """
        The function leave_one_out_crossvalidation determines the best hyperparameters (shape and regularization parameters) for a given RBF fitting problem.
        The function cycles through a set of predefined sets to determine the shape parameter and regularization parameter combination which yields the lowest LOOCV error.
        The LOOCV error for each (shape_parameter, regulkarization parameter) pair is evaluated by calling the function loo_error_estimation_with_rippa_method
        The pre-defined shape parameter set considers 24 irregularly spaced values ranging between 0.001 - 1000, while the regularization parameter set considers 21 values ranging between 0.00001 - 1.

        Args:
            self:                           : contains, among other things, the input data

        Returns:
            r_best(float)                 : best found shape parameter
            lambda_best(float)            : best found regularization parameter
            error_best                      : LOOCV error corresponding to be best found hynperparameters

        Note: The optimal shape parameter r_best is only evaluated for parametriuc bases (such as the Gaussian basis). For fixed basis (e.g. linear), the value is returned as zero.

        """
        # Define sigma and lambda ranges
        if (
            (self.basis_function == "gaussian")
            or (self.basis_function == "mq")
            or (self.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]

        if self.regularization is True:
            # reg_parameter = [0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1]
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif self.regularization is False:
            reg_parameter = [0]

        machine_precision = np.finfo(float).eps

        error_vector = np.zeros((len(r_set) * len(reg_parameter), 3))
        counter = 0
        print(
            "==========================================================================================================="
        )
        for i in range(0, len(r_set)):
            sigma = r_set[i]
            for j in range(0, len(reg_parameter)):
                lambda_reg = reg_parameter[j]
                (
                    cond_no_pure,
                    cond_no_reg,
                    cv_error,
                ) = self.loo_error_estimation_with_rippa_method(sigma, lambda_reg)
                error_vector[counter, :] = [sigma, lambda_reg, cv_error]
                counter += 1
                print(
                    sigma,
                    "   |    ",
                    lambda_reg,
                    "   |    ",
                    cv_error,
                    "   |    ",
                    cond_no_pure,
                    "   |    ",
                    cond_no_pure * machine_precision,
                    "   |    ",
                    cond_no_reg,
                    "   |    ",
                    cond_no_reg * machine_precision,
                )
        minimum_value_column = np.argmin(error_vector[:, 2], axis=0)
        r_best = error_vector[minimum_value_column, 0]
        lambda_best = error_vector[minimum_value_column, 1]
        error_best = error_vector[minimum_value_column, 2]
        return r_best, lambda_best, error_best

    def training(self):
        """
        Main function for RBF training.

        To train the RBF:
            (1) The best values of the hyperparameters (:math:`\sigma, \lambda`) are selected via LOOCV.
            (2) The necessary basis transformation at the optimal hyperparameters is generated.
            (3) The condition number for the transformed matrix is calculated.
            (4) The optimal radial weights are evaluated using the selected optimization method.
            (5) The training predictions, prediction errors and r-square coefficient of fit are evaluated by calling the methods ``error_calculation`` and ``r2_calculation``
            (6) A results object is generated by calling the ResultsReport class.

        The LOOCV error for each (:math:`\sigma, \lambda`) pair is evaluated by calling the function ``loo_error_estimation_with_rippa_method``.

        The pre-defined shape parameter set considers 24 irregularly spaced values ranging between 0.001 - 1000, while the regularization parameter set considers 21 values ranging between 0.00001 - 1.

        Returns:
            tuple   : self object (**results**) containing the all information about the best RBF fitting obtained, including:
                - the optimal radial weights  (**results.radial_weights**),
                - when relevant, the optimal shape parameter found :math:`\sigma` (**results.sigma**),
                - when relevant, the optimal regularization parameter found :math:`\lambda` (**results.regularization**),
                - the RBF predictions for the training data (**results.output_predictions**), and
                - the :math:`R^{2}` value on the training data (**results.R2**)

        """

        # Determine best r value
        best_r_value, best_lambda_param, _ = self.leave_one_out_crossvalidation()

        # Generate x matrix
        x_transformed = self.basis_generation(best_r_value)
        x_transformed = x_transformed + (
            best_lambda_param * np.eye(x_transformed.shape[0], x_transformed.shape[1])
        )
        x_condition_number = np.linalg.cond(x_transformed)

        if self.solution_method == "algebraic":
            radial_weights = self.explicit_linear_algebra_solution(
                x_transformed, self.y_data
            )
        elif self.solution_method == "pyomo":
            radial_weights = self.pyomo_optimization(x_transformed, self.y_data)
        elif self.solution_method == "bfgs":
            radial_weights = self.bfgs_parameter_optimization(
                x_transformed, self.y_data
            )
        radial_weights = radial_weights.reshape(radial_weights.shape[0], 1)

        (
            training_ss_error,
            rmse_error,
            y_training_predictions_scaled,
        ) = self.error_calculation(
            radial_weights, self.basis_generation(best_r_value), self.y_data
        )
        r_square = self.r2_calculation(self.y_data, y_training_predictions_scaled)
        y_training_predictions = self.data_min[
            0, -1
        ] + y_training_predictions_scaled * (
            self.data_max[0, -1] - self.data_min[0, -1]
        )

        # Results
        self.weights = radial_weights
        self.sigma = best_r_value
        self.regularization_parameter = best_lambda_param
        self.rmse = rmse_error
        self.output_predictions = y_training_predictions
        self.condition_number = x_condition_number
        self.R2 = r_square
        self.x_data_min = self.data_min[:, :-1]
        self.x_data_max = self.data_max[:, :-1]
        self.y_data_min = self.data_min[:, -1]
        self.y_data_max = self.data_max[:, -1]
        if x_condition_number < (1 / np.finfo(float).eps):
            self.solution_status = "ok"
        else:
            warnings.warn(
                "The parameter matrix A in A.x=B is ill-conditioned (condition number > 1e10). The solution returned may be inaccurate or unstable - inspect rmse error. Regularization (if not already done) may improve solution"
            )
            self.solution_status = "unstable solution"

        self.pickle_save({"model": self})
        return self

    def predict_output(self, x_data):
        """

        The ``predict_output`` method generates output predictions for input data x_data based a previously generated RBF fitting.

        Args:
            x_data(NumPy Array)    : Designs for which the output is to be evaluated/predicted.

        Returns:
             Numpy Array    : Output variable predictions based on the rbf fit.

        """
        radial_weights = self.weights
        centres_matrix = self.centres
        r = self.sigma
        lambda_reg = self.regularization_parameter
        scale = self.x_data_max - self.x_data_min
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_data - self.x_data_min) / scale
        x_data = x_pred_scaled.reshape(x_data.shape)

        basis_vector = np.zeros((x_data.shape[0], centres_matrix.shape[0]))
        # Calculate distances from centres
        for i in range(0, centres_matrix.shape[0]):
            basis_vector[:, i] = np.sqrt(
                np.sum(((x_data - centres_matrix[i, :]) ** 2), axis=1)
            )
        # Initialization of x_transformed
        x_transformed = np.zeros((basis_vector.shape[0], basis_vector.shape[1]))

        # Transform X
        if self.basis_function == "gaussian":
            x_transformed = RadialBasisFunctions.gaussian_basis_transformation(
                basis_vector, r
            )
        elif self.basis_function == "linear":
            x_transformed = RadialBasisFunctions.linear_transformation(basis_vector)
        elif self.basis_function == "cubic":
            x_transformed = RadialBasisFunctions.cubic_transformation(basis_vector)
        elif self.basis_function == "mq":
            x_transformed = RadialBasisFunctions.multiquadric_basis_transformation(
                basis_vector, r
            )
        elif self.basis_function == "imq":
            x_transformed = (
                RadialBasisFunctions.inverse_multiquadric_basis_transformation(
                    basis_vector, r
                )
            )
        elif self.basis_function == "spline":
            x_transformed = RadialBasisFunctions.thin_plate_spline_transformation(
                basis_vector
            )

        # Add regularization shifting?
        x_transformed = x_transformed + (
            0 * np.eye(x_transformed.shape[0], x_transformed.shape[1])
        )
        # x_transformed = x_transformed + (lambda_reg * np.eye(x_transformed.shape[0], x_transformed.shape[1]))
        y_prediction_scaled = np.matmul(x_transformed, radial_weights)
        y_prediction_unscaled = self.y_data_min + y_prediction_scaled * (
            self.y_data_max - self.y_data_min
        )
        return y_prediction_unscaled

    def generate_expression(self, variable_list):
        """
        The ``generate_expression`` method returns the Pyomo expression for the RBF model trained.

        The expression is constructed based on the supplied list of variables **variable_list** and the results of the previous RBF training process.

        Args:
            variable_list(list)           : List of input variables to be used in generating expression. This can be the a list generated from the output of ``get_feature_vector``. The user can also choose to supply a new list of the appropriate length.

        Returns:
            Pyomo Expression              : Pyomo expression of the RBF model based on the variables provided in **variable_list**

        """
        t1 = np.array([variable_list])
        # Reshaping of variable array is necessary when input variables are Pyomo scalar variables
        t1 = t1.reshape(1, len(variable_list)) if t1.ndim > 2 else t1

        basis_vector = []
        # Calculate distances from centres
        for i in range(0, self.centres.shape[0]):
            ans = 0
            for j in range(0, self.centres.shape[1]):
                ans += (
                    (
                        (t1[0, j] - self.x_data_min[0, j])
                        / (self.x_data_max[0, j] - self.x_data_min[0, j])
                    )
                    - self.centres[i, j]
                ) ** 2
            eucl_d = ans**0.5
            basis_vector.append(eucl_d)
        rbf_terms_list = []
        if self.basis_function == "linear":
            for k in range(0, len(basis_vector)):
                rbf_terms_list.append(
                    RadialBasisFunctions.linear_transformation(basis_vector[k])
                )
        elif self.basis_function == "cubic":
            for k in range(0, len(basis_vector)):
                rbf_terms_list.append(
                    RadialBasisFunctions.cubic_transformation(basis_vector[k])
                )
        elif self.basis_function == "gaussian":
            for k in range(0, len(basis_vector)):
                rbf_terms_list.append(exp(-1 * ((self.sigma * basis_vector[k]) ** 2)))
        elif self.basis_function == "mq":
            for k in range(0, len(basis_vector)):
                rbf_terms_list.append(
                    (((basis_vector[k] * self.sigma) ** 2) + 1) ** 0.5
                )
        elif self.basis_function == "imq":
            for k in range(0, len(basis_vector)):
                rbf_terms_list.append(
                    1 / ((((basis_vector[k] * self.sigma) ** 2) + 1) ** 0.5)
                )
        elif self.basis_function == "spline":
            for k in range(0, len(basis_vector)):
                rbf_terms_list.append(((basis_vector[k] ** 2) * log(basis_vector[k])))

        rbf_terms_array = np.asarray(rbf_terms_list)
        rbf_expr = self.y_data_min[0]
        rbf_expr += (self.y_data_max[0] - self.y_data_min[0]) * sum(
            w * t
            for w, t in zip(
                np.nditer(self.weights), np.nditer(rbf_terms_array, flags=["refs_ok"])
            )
        )
        return rbf_expr

    def get_feature_vector(self):
        """

        The ``get_feature_vector`` method generates the list of regression features from the column headers of the input dataset.

        Returns:
            Pyomo IndexedParam  : An indexed parameter list of the variables supplied in the original data

        **Example:**

        .. code-block:: python

            # Create a small dataframe with three columns ('one', 'two', 'three') and two rows (A, B)
            >>> xy_data = pd.DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6])], orient='index', columns=['one', 'two', 'three'])

            # Initialize the **RadialBasisFunctions** class with a linear kernel and print the column headers for the variables
            >>> f = RadialBasisFunctions(xy_data, basis_function='linear')
            >>> p = f.get_feature_vector()
            >>> for i in p.keys():
            >>>     print(i)
            one
            two

        """
        p = Param(self.x_data_columns, mutable=True, initialize=0)
        p.index_set().construct()
        p.construct()
        self.feature_list = p
        return p

    def pickle_save(self, solutions):
        """
        The pickle_save method saves the results of the run in a pickle object.
        """
        try:
            filehandler = open(self.filename, "wb")
            pickle.dump(solutions, filehandler)
            print("\nResults saved in ", str(self.filename))
        except:
            raise IOError("File could not be saved.")

    @staticmethod
    def pickle_load(solution_file):
        """
        pickle_load loads the results of a saved run 'file.obj'.

        Input arguments:
                solution_file            : Pickle object file containing previous solution to be loaded.

        """
        try:
            filehandler = open(solution_file, "rb")
            return pickle.load(filehandler)
        except:
            raise Exception("File could not be loaded.")

    def parity_residual_plots(self):
        """

        inputs:

        Returns:

        """

        fig1 = plt.figure(figsize=(16, 9), tight_layout=True)
        ax = fig1.add_subplot(121)
        ax.plot(self.y_data_unscaled, self.y_data_unscaled, "-")
        ax.plot(self.y_data_unscaled, self.output_predictions, "o")
        ax.set_xlabel(r"True data", fontsize=12)
        ax.set_ylabel(r"Surrogate values", fontsize=12)
        ax.set_title(r"Parity plot", fontsize=12)

        ax2 = fig1.add_subplot(122)
        ax2.plot(
            self.y_data_unscaled,
            self.y_data_unscaled - self.output_predictions,
            "s",
            mfc="w",
            mec="m",
            ms=6,
        )
        ax2.axhline(y=0, xmin=0, xmax=1)
        ax2.set_xlabel(r"True data", fontsize=12)
        ax2.set_ylabel(r"Residuals", fontsize=12)
        ax2.set_title(r"Residual plot", fontsize=12)

        plt.show()

        return

    def _report(self):
        ## Will only work with Python > 3.5
        variable_headers = self.get_feature_vector()
        var_list = []
        for i in variable_headers:
            var_list.append(variable_headers[i])
        eqn = self.generate_expression(var_list)

        double_line = "=" * 120
        s = (
            f"\n{double_line}"
            f"\nResults of radial basis function run:\n"
            f"\nBasis function type               : {self.basis_function}\n"
            f"Shape parameter                    : {self.sigma}\n"
            f"Regularization parameter           : {self.regularization_parameter}\n"
            f"Number of terms in RBF model       : {self.weights.size + 1}\n"  # The additional term is y_min
            f"\nRBF Expression:\n"
            f"--------------------------\n"
            f"\n{eqn}\n"
            f"--------------------------\n"
            f"\nModel training errors:"
            f"\n-----------------------\n"
            f"Mean Squared Error (MSE)         : {self.rmse ** 2}\n"
            f"Root Mean Squared Error (RMSE)   : {self.rmse}\n"
            f"Goodness of fit (R2)             : {self.R2}\n"
            f"\n{double_line}"
        )
        return s

    def print_report(self):
        s = self._report()
        print(s)

    def _repr_pretty_(self, p, cycle=False):
        s = self._report()
        p.text(s)
