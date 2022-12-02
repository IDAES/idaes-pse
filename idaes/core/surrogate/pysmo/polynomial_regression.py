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
from __future__ import division

# from builtins import int, str
import os.path
import warnings

# Imports from third parties
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pickle
from pyomo.environ import *
from pyomo.core.expr.visitor import replace_expressions
import scipy.optimize as opt
from six import string_types

# Imports from IDAES namespace
from idaes.core.surrogate.pysmo.utils import NumpyEvaluator


__author__ = "Oluwamayowa Amusat"

"""
The purpose of this file is to perform polynomial regression in Pyomo.
This will be done in two stages. First, a sampling plan will
be used to select samples for generating a surrogate model.
In the second stage, the surrogate model is constructed by fitting to
different order polynomials. Long term, an iterative adaptive sampling 
approach will be incorporated for model improvement.
Cross-validation is used to select the best model.


FeatureScaling:
Simple script for scaling and un-scaling the input data

Polynomial Regression:
Three approaches are implemented for evaluating polynomial coefficients -
a. Moore-Penrose maximum likelihood method (Forrester et al.)
b. Optimization using the BFGS algorithm.
c. Optimization with Pyomo
The Pyomo optimization approach is enabled as the default at this time.
"""


class FeatureScaling:
    """

    A class for scaling and unscaling input and output data. The class contains two main methods: ``data_scaling`` and ``data_unscaling``
    """

    def __init__(self):
        pass

    @staticmethod
    def data_scaling(data):
        """
        ``data_scaling`` performs column-wise minimax scaling on the input dataset.

        Args:
            data : The input data set to be scaled. Must be a numpy array or dataframe.

        Returns:
            (tuple): tuple containing:
                - **scaled_data**  : A 2-D Numpy Array containing the scaled data. All array values will be between [0, 1].
                - **data_minimum** : A 2-D row vector containing the column-wise minimums of the input data.
                - **data_maximum** : A 2-D row vector containing the column-wise maximums of the input data.

        Raises:
            TypeError:
                Raised when the input data is not a numpy array or dataframe

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
    def data_unscaling(x_scaled, x_min, x_max):
        """

        ``data_unscaling`` performs column-wise un-scaling on the a minmax-scaled input dataset.

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


class PolynomialRegression:
    """

    The PolynomialRegression class performs polynomial regression on a training data set.

    The class must first be initialized by calling PolynomialRegression. Regression is then carried out by calling ``training``.

    For a given dataset with :math:`n` features :math:`x_{1},x_{2},\ldots,x_{n}`, Polyregression is able to consider three types of basis functions:
        (a) Mononomial terms (:math:`x_{i}^{p},p \leq 10`) for all individual features. The maximum degree to be considered can be set by the user (**maximum_polynomial_order**)
        (b) All first order interaction terms :math:`x_{1}x_{2}`, :math:`x_{1}x_{3}` etc. This can be turned on or off by the user (set **multinomials**)
        (c) User defined input features, e.g. :math:`\sin(x_{1})`. These must be Pyomo functions and should be provided as a list by the user calling ``set_additional_terms`` method before the polynomial training is done.

    **Example:**

    .. code-block:: python

        # Initialize the class and set additional terms
        >>> d = PolynomialRegression(full_data, training_data, maximum_polynomial_order=2, max_iter=20, multinomials=1, solution_method='pyomo')
        >>> p = d.get_feature_vector()
        >>> d.set_additional_terms([...extra terms...])

        # Train polynomial model and predict output for an test data x_test
        >>> d.training()
        >>> predictions = d.predict_output(x_test)

    Args:
        regression_data_input(NumPy Array of Pandas Dataframe) : The dataset for regression training. It is expected to contain the features (X) and output (Y) data, with the output values (Y) in the last column.
        original_data_input(NumPy Array of Pandas Dataframe) : If **regression_data_input** was drawn from a larger dataset by some sampling approach, the larger dataset may be provided here.
            When additional data is not available, the same data supplied for training_data can be supplied - this tells the algorithm not to carry out adaptive sampling.
        maximum_polynomial_order(int): The maximum polynomial order to be considered.

    Further details about the optional inputs may be found under the ``__init__`` method.

    """

    def __init__(
        self,
        original_data_input,
        regression_data_input,
        maximum_polynomial_order,
        number_of_crossvalidations=None,
        no_adaptive_samples=None,
        training_split=None,
        max_fraction_training_samples=None,
        max_iter=None,
        solution_method=None,
        multinomials=None,
        fname=None,
        overwrite=False,
    ):
        """
        Initialization of PolynomialRegression class.

        Args:
            regression_data_input(NumPy Array of Pandas Dataframe) : The dataset for regression training. It is expected to contain features and output data, with the output values (Y) in the last column.
            original_data_input(NumPy Array of Pandas Dataframe) : If **regression_data_input** was drawn from a larger dataset by some sampling approach, the larger dataset may be provided here.
            maximum_polynomial_order(int) : The maximum polynomial order to be considered.

        Keyword Args:
            number_of_crossvalidations(int) : The number of polynomial fittings and cross-validations to be carried out for each polynomial function/expression. Must be a positive, non-zero integer. Default=3.

            training_split(float): The training/test split to be used for regression_data_input. Must be between 0 and 1. Default = 0.75

            solution_method(str): The method to be used for solving the least squares optimization problem for polynomial regression. Three options are available:

                (a) "MLE"  : The mle (maximum likelihood estimate) method solves the least squares problem using linear algebra. Details of the method may be found in Forrester et al.
                (b) "BFGS" : This approach solves the least squares problem using scipy's BFGS algorithm.
                (c) "pyomo": This option solves the optimization problem in pyomo with IPOPT as solver. This is the default option.

            multinomials(bool):  This option determines whether or not multinomial terms are considered during polynomial fitting. Takes 0 for No and 1 for Yes. Default = 1.

        Returns:
            **self** object containing all the input information.

        Raises:
            ValueError:
                - The input datasets (**original_data_input** or **regression_data_input**) are of the wrong type (not Numpy arrays or Pandas Dataframes)

            Exception:
                * **maximum_polynomial_order** is not a positive, non-zero integer or **maximum_polynomial_order** is higher than the number of training samples available
            Exception:
                * **solution_method** is not 'mle', 'pyomo' or 'bfgs
            Exception:
                - **multinomials** is not binary (0 or 1)
            Exception:
                - **training_split** is not between 0 and 1
            Exception:
                - **number_of_crossvalidations** is not a positive, non-zero integer
            Exception:
                - **max_fraction_training_samples** is not between 0 and 1
            Exception:
                - **no_adaptive_samples** is not a positive, non-zero integer
            Exception:
                - **max_iter** is not a positive, non-zero integer

            warnings.warn:
                - When the number of cross-validations is too high, i.e. number_of_crossvalidations > 10
        """

        print(
            "\n===========================Polynomial Regression===============================================\n"
        )
        # Checks if fname is provided or exists in the path
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

        if isinstance(original_data_input, pd.DataFrame):
            original_data = original_data_input.values
            # FIXME: if we add an option to specify the response column, this needs to change
            self.regression_data_columns = list(original_data_input.columns)[:-1]
        elif isinstance(original_data_input, np.ndarray):
            original_data = original_data_input
            self.regression_data_columns = list(range(original_data_input.shape[1] - 1))
        else:
            raise ValueError(
                "original_data_input: Pandas dataframe or numpy array required."
            )

        if isinstance(regression_data_input, pd.DataFrame):
            regression_data = regression_data_input.values
        elif isinstance(regression_data_input, np.ndarray):
            regression_data = regression_data_input
        else:
            raise ValueError(
                "regression_data_input: Pandas dataframe or numpy array required."
            )

        # Check for potential dimensionality problems in input data
        if regression_data.shape[0] > original_data.shape[0]:
            raise Exception(
                "The sampled data has more entries than the original dataset."
            )
        elif regression_data.shape[1] != original_data.shape[1]:
            raise Exception(
                "Dimensional discrepancies in the dimensions of the original and regression datasets."
            )
        elif (regression_data.shape[1] == 1) or (original_data.shape[1] == 1):
            raise Exception(
                "Input data requires at least two dimensions (X and Y data)."
            )

        self.original_data = original_data
        self.regression_data = regression_data

        if number_of_crossvalidations is None:
            print("The number of cross-validation cases (3) is used.")
            number_of_crossvalidations = 3
        elif number_of_crossvalidations > 10:
            warnings.warn(
                "The number of cross-validations entered is large. The simulation may take a while to run"
            )
        self.number_of_crossvalidations = number_of_crossvalidations

        if not isinstance(maximum_polynomial_order, int):
            raise Exception("Maximum polynomial order must be an integer")
        elif maximum_polynomial_order > 10:
            warnings.warn(
                "The maximum allowed polynomial order is 10. Value has been adjusted to 10."
            )
            maximum_polynomial_order = 10
        self.max_polynomial_order = maximum_polynomial_order

        self.number_of_x_vars = regression_data.shape[1] - 1

        if training_split is None:
            print("The default training/cross-validation split of 0.75 is used.")
            training_split = 0.75
        elif training_split >= 1 or training_split <= 0:
            raise Exception(
                "Fraction of samples used for training must be between 0 and 1"
            )
        self.fraction_training = training_split

        if no_adaptive_samples is None:
            no_adaptive_samples = 4
        self.no_adaptive_samples = no_adaptive_samples

        self.number_of_samples = regression_data.shape[0]

        if max_fraction_training_samples is None:
            max_fraction_training_samples = 0.5
        elif max_fraction_training_samples > 1 or max_fraction_training_samples < 0:
            raise Exception(
                "The fraction for the maximum number of training samples must be between 0 and 1"
            )
        self.max_fraction_training_samples = max_fraction_training_samples

        if (regression_data.shape[0] < original_data.shape[0]) and max_iter is None:
            max_iter = 10
        if (
            regression_data.shape[0] == original_data.shape[0]
            or no_adaptive_samples == 0
        ):
            print("No iterations will be run.")
            max_iter = 0
        self.max_iter = max_iter

        # Ensure all other key variables are integers
        if not isinstance(self.number_of_crossvalidations, int):
            raise Exception("Number of cross-validations must be an integer")
        elif not isinstance(self.no_adaptive_samples, int):
            raise Exception("Number of adaptive samples must be an integer")
        elif not isinstance(self.max_iter, int):
            raise Exception("Maximum number of iterations must be an integer")
        elif self.max_polynomial_order >= regression_data.shape[0]:
            raise Exception(
                "max_polynomial_order too high for the number of samples supplied"
            )

        if (self.max_polynomial_order <= 0) or (self.number_of_crossvalidations <= 0):
            raise Exception(
                "maximum_polynomial_order and number_of_crossvalidations must be positive, non-zero integers"
            )
        elif (self.no_adaptive_samples < 0) or (self.max_iter < 0):
            raise Exception("no_adaptive_samples and max_iter must be positive")

        if solution_method is None:
            solution_method = "pyomo"
            self.solution_method = solution_method
            print("Default parameter estimation method is used.")
        elif not isinstance(solution_method, string_types):
            raise Exception("Invalid solution method. Must be of type <str>.")
        elif (
            (solution_method.lower() == "mle")
            or (solution_method.lower() == "pyomo")
            or (solution_method.lower() == "bfgs")
        ):
            solution_method = solution_method.lower()
            self.solution_method = solution_method
        else:
            raise Exception(
                'Invalid parameter estimation method entered. Select one of maximum likelihood (solution_method="mle"), Pyomo optimization (solution_method="pyomo") or BFGS (solution_method="bfgs") methods. '
            )
        print("Parameter estimation method: ", self.solution_method, "\n")

        if multinomials is None:
            self.multinomials = 1
        elif multinomials == 1:
            self.multinomials = 1
        elif multinomials == 0:
            self.multinomials = 0
        else:
            raise Exception(
                'Multinomial must be binary: input "1" for "Yes" and "0" for "No". '
            )

        self.feature_list = []
        self.additional_term_expressions = []

        # Results
        self.optimal_weights_array = None
        self.final_polynomial_order = None
        self.errors = None
        self.number_of_iterations = None
        self.iteration_summary = None
        self.additional_features_data = None
        self.final_training_data = None
        self.dataframe_of_optimal_weights_polynomial = None
        self.dataframe_of_optimal_weights_extra_terms = None
        self.extra_terms_feature_vector = None
        self.fit_status = None

    def training_test_data_creation(self, additional_features=None):

        """

        The training_test_data_creation splits data into training and test data sets.

        Given the number of cross-validations and the required training/test split, it:
            - calculates the number of training samples as num_training = int(training_split x total number of samples),
            - shuffles regression_data number_of_crossvalidations times,
            - splits off the top num_training samples in each shuffle as individual training sets, and
            - takes the bottom (total number of samples - num_training) samples in each shuffle to create its corresponding test dataset.

        Args:
            self: containing the number of training samples (self.number_of_samples), training/test split (self.fraction_training) and the required number of cross-validations (self.number_of_crossvalidations).

        Keyword Args:
            additional_features(NumPy Array): A numpy array containing additional features provided by the user. When supplied, additional_features is column-appended to self.regression data before the training and tests sets are created.

        Returns:
            Tuple(training_data, cross_val_data)

            - training_data: Dictionary containing all the training datasets created.

                * When no additional features have been specified, the dictionary has a length of number_of_crossvalidations.
                * When additional features have been specified, the dictionary has a length of 2 * number_of_crossvalidations, with the training data for additional_features stored separately.

            - cross_val_data: Dictionary containing all the test datasets created. Dictionary will have the same length as training_data.

        """
        training_data = {}
        cross_val_data = {}
        num_training = int(np.around(self.number_of_samples * self.fraction_training))
        if num_training == 0:
            raise Exception("The inputted of fraction_training is too low.")
        elif num_training == self.number_of_samples:
            raise Exception("The inputted of fraction_training is too high.")
        for i in range(1, self.number_of_crossvalidations + 1):
            np.random.seed(i)
            if additional_features is None:
                A = np.zeros(
                    (self.regression_data.shape[0], self.regression_data.shape[1])
                )
                A[:, :] = self.regression_data
                np.random.shuffle(
                    A
                )  # Shuffles the rows of the regression data randomly
                training_data["training_set_" + str(i)] = A[0:num_training, :]
                cross_val_data["test_set_" + str(i)] = A[num_training:, :]
            elif additional_features is not None:
                A = np.zeros(
                    (
                        self.regression_data.shape[0],
                        self.regression_data.shape[1] + additional_features.shape[1],
                    )
                )
                A[:, 0 : self.regression_data.shape[1]] = self.regression_data
                A[:, self.regression_data.shape[1] :] = additional_features
                np.random.shuffle(
                    A
                )  # Shuffles the rows of the regression data randomly
                training_data["training_set_" + str(i)] = A[
                    0:num_training, : self.regression_data.shape[1]
                ]
                training_data["training_extras_" + str(i)] = A[
                    0:num_training, self.regression_data.shape[1] :
                ]
                cross_val_data["test_set_" + str(i)] = A[
                    num_training:, : self.regression_data.shape[1]
                ]
                cross_val_data["test_extras_" + str(i)] = A[
                    num_training:, self.regression_data.shape[1] :
                ]
        return training_data, cross_val_data

    @classmethod
    def polygeneration(
        self,
        polynomial_order,
        multinomials,
        x_input_train_data,
        additional_x_training_data=None,
    ):
        """

        This function generates a x-variable vector for the required polynomial order. This is done in four stages:

        - First, generates the pure mononomials are generated by increasing the polynomial degree  by 1 until polynomial_order is reached.
        - Next, the first-order multinomials are generated if self.multinomials = 1. This is implemented in suci a way that each multinomial appears only once, i.e. x_i.x_j = x_j.x_i. The multinomial columns are appended to the enx of the array.
        - Next, a column of ones is inserted in the first column space to represent the constant term.
        - Finally, the columns containing the extra terms supplied by the user are added to the end of the array (when available).

        Thus, the format of the output array is [constant, nmononomials, multinomials, extra terms]

        Args:
            polynomial_order(int):                    The polynomial order currently under consideration
            multinomials(bool):                       Boolean variable that determines whether or not multinomial terms are considered during polynomial fitting.
            x_input_train_data(NumPy Array):           Input data containing features supplied by the user

        Keyword Args:
            additional_x_training_data(NumPy Array):   Array containing additional features supplied by the user

        Returns:
            x_train_data(NumPy Array):                 Array containing all polynomial features to be considered during regression

        Example:
                        if polynomial_order=2, numtinomials=1, x_input_train_data = [x1, x2, x3], additional_x_training_data = [sin(x1), tanh(x3)], then x_train_data will contain the regression features
            x_train_data = [1, x1, x2, x3, x1^2, x2^2, x3^2, x1.x2, x1.x3, x2.x3, sin(x1), tanh(x3)]

        """
        N = x_input_train_data.shape[0]
        x_train_data = x_input_train_data
        # Generate the constant and pure power terms
        for i in range(2, polynomial_order + 1):
            x_train_data = np.concatenate(
                (x_train_data, x_input_train_data**i), axis=1
            )

        if multinomials == 1:
            # Next, generate first order multinomials
            for i in range(0, x_input_train_data.shape[1]):
                for j in range(0, i):
                    x_train_data = np.concatenate(
                        (
                            x_train_data,
                            (
                                x_input_train_data[:, i] * x_input_train_data[:, j]
                            ).reshape(N, 1),
                        ),
                        axis=1,
                    )

        # Concatenate to generate full dataset.
        x_train_data = np.concatenate((np.ones((N, 1)), x_train_data), axis=1)

        # Add additional features if they have been provided:
        if additional_x_training_data is not None:
            x_train_data = np.concatenate(
                (x_train_data, additional_x_training_data), axis=1
            )

        return x_train_data

    @staticmethod
    def cost_function(theta, x, y, reg_parameter):
        """

        This function is an implementation of the cost function for linear regression:
                cost = [sum of square errors over m samples / (2 * m)]  + [reg_parameter * theta*2 / (2 * m)]

        This is the objective function for the BFGS optimization problem.

        Args:
            theta        : polynomial coefficients/weights,  (n x 1) in size
            x            : array of features, (m x n) in size
            y            : actual output vector, size (m x 1)
            reg_parameter: reqularization parameter, set to

        Returns:
            cost_value   : the cost value for the fit, the objective value of the optimization problem

        """

        y = y.reshape(y.shape[0], 1)
        y_prediction = np.matmul(x, theta)
        y_prediction = y_prediction.reshape(y_prediction.shape[0], 1)
        cost_value = (0.5 / x.shape[0]) * (np.sum((y - y_prediction) ** 2))
        cost_penalty = (reg_parameter * 0.5 / x.shape[0]) * (np.sum(theta**2))
        cost_value = cost_value + cost_penalty
        return cost_value

    @staticmethod
    def gradient_function(theta, x, y, reg_parameter):
        """

        This function is an implementation of the gradient function for linear regression:
                if
                    cost = [(A.x - y)^2 / 2m] + [reg_parameter * theta*2/ (2 * m)],
                then
                    gradient = [((A.x - y)* A) / m] + [reg_parameter * theta/ m]

        This is the gradient function supplied to the BFGS optimization algorithm.

        Args:
            theta        : polynomial coefficients/weights,  (n x 1) in size
            x            : array of features, (m x n) in size
            y            : actual output vector, size (m x 1)
            reg_parameter: reqularization parameter

        Returns:
            grad_value   : the cost gradients for the fit, size (n x 1)

        """
        y = y.reshape(y.shape[0], 1)
        y_prediction = np.matmul(x, theta)
        y_prediction = y_prediction.reshape(y_prediction.shape[0], 1)
        t1 = (y_prediction - y) * x
        grad_values = (1 / x.shape[0]) * np.sum(t1, axis=0)
        gradient_penalty = (reg_parameter / x.shape[0]) * theta
        grad_values = grad_values + gradient_penalty
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
            The regularization parameter and initial weights are set to zero,
            reg_parameter = 0
            init_theta = 0

        Returns:
            theta: The optimal linear regression weights found

        """
        init_theta = np.zeros((x.shape[1], 1))
        reg_parameter = 0.0
        other_args = (x, y, reg_parameter)
        theta = opt.fmin_bfgs(
            self.cost_function,
            init_theta,
            fprime=self.gradient_function,
            args=other_args,
        )
        return theta

    @staticmethod
    def MLE_estimate(x, y):
        """

        Maximum likelihood estimate method for solving polynomial regression problems:

            If
                Ax = B,
            then
                x = inv_A * B

            where the inv_A is called the Moore-Penrose inverse.

        Numpy's pseudoinverse function has been used to calculate the inverse here.

        Args:
            x            : array of features, (m x n) in size
            y            : actual output vector, size (m x 1)

        Returns:
            phi: The optimal linear regression weights found

         For more details about the maximum likelihood estimate methos, see to Forrester et al.

        """
        moore_penrose_inverse = np.linalg.pinv(x)  # Moore Penrose inverse of vector x
        phi = np.matmul(moore_penrose_inverse, y)
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
        model.theta = Var(model.N, initialize=0.1, domain=Reals)

        # constraint y_p = theta.X
        def xy_product(model, i, k):
            return sum(model.theta[j] * model.x[i, j] for j in model.N for k in model.P)

        model.y_predictions = Expression(
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
        opt.options["max_iter"] = 10000000
        result = opt.solve(instance)  # , tee=True)

        # Convert theta variable into numpy array
        phi = np.zeros((len(instance.theta), 1))
        iterator = 0
        for s in instance.N:
            phi[iterator, 0] = instance.theta[s].value
            iterator += 1
        return phi

    @staticmethod
    def cross_validation_error_calculation(phi, x_test_data, y_test_data):
        """

        This function calculates the average sum of square errors between the actual and predicted output values,
             ss_error = sum of squared errors / number of samples

        Args:
            phi             : optimal weight vector obtained by optimization
            x_test_data     : vector of features x_test_data
            y_test_data     : actual output values associated with

        Returns:
            ss_error        : The average sum of squared errors

        """
        y_test_prediction = np.matmul(x_test_data, phi)
        ss_error = (1 / y_test_data.shape[0]) * (
            np.sum((y_test_data - y_test_prediction) ** 2)
        )
        return ss_error

    def polyregression(
        self,
        poly_order,
        training_data,
        test_data,
        additional_x_training_data=None,
        additional_x_test_data=None,
    ):
        """

        Function that performs polynomial regression on a given dataset. It returns the estimated parameters and the fitting errors. It

            - calls the method self.polygeneration to generate the required polynomial/feature array based on the current polynomial order poly_order,
            - calls the pre-selected solution algorithm to solve the least squares problem, and
            - calls the cross_validation_error_calculation method to calculate the training and cross-validation errors.


        Args:
            poly_order(int)           : The polynomial order currently being considered - between 1 and max_polynomial_order
            training_data(NumPy Array) : The training data to be regressed
            test_data(NumPy Array)    : The test data to be used to cross-validate the polynomial fit

        Keyword Args:
            additional_x_training_data  : Array containing additional training features based on additional_features list supplied by the user. Will have same number of rows as training_data.
            additional_x_test_data      : Array of additional cross-validation features based on additional_features list supplied by the user. Will have same number of rows as test_data.

        Returns:
            phi_vector                  : the optimal weight vector for the polynomial considered here, returns zeros when problem is underspecified, i.e number of features > number of training samples.
            training_error              : the average SSE estimate in the training dataset, returns Inf when number of features > number of training samples (DoF < 0).
            crossval_error             : the average SSE estimate on the cross-validation dataset, returns Inf when number of features > number of training samples (DoF < 0).

        """
        x_training_data = training_data[:, :-1]
        y_training_data = training_data[:, -1]
        x_test_data = test_data[:, :-1]
        y_test_data = test_data[:, -1]
        x_polynomial_data = self.polygeneration(
            poly_order, self.multinomials, x_training_data, additional_x_training_data
        )

        # Check that the problem has more samples than features - necessary for fitting. If not, return Infinity.
        if x_polynomial_data.shape[0] >= x_polynomial_data.shape[1]:
            if self.solution_method == "mle":
                phi_vector = self.MLE_estimate(
                    x_polynomial_data,
                    y_training_data.reshape(y_training_data.shape[0], 1),
                )
            elif self.solution_method == "bfgs":
                phi_vector = self.bfgs_parameter_optimization(
                    x_polynomial_data, y_training_data
                )
            elif self.solution_method == "pyomo":
                phi_vector = self.pyomo_optimization(x_polynomial_data, y_training_data)
            phi_vector = phi_vector.reshape(
                phi_vector.shape[0], 1
            )  # Pseudo-inverse approach

            x_polynomial_data_test = self.polygeneration(
                poly_order, self.multinomials, x_test_data, additional_x_test_data
            )
            training_error = self.cross_validation_error_calculation(
                phi_vector,
                x_polynomial_data,
                y_training_data.reshape(y_training_data.shape[0], 1),
            )
            crossval_error = self.cross_validation_error_calculation(
                phi_vector,
                x_polynomial_data_test,
                y_test_data.reshape(y_test_data.shape[0], 1),
            )

        else:
            phi_vector = np.zeros((x_polynomial_data.shape[1], 1))
            phi_vector[:, 0] = np.Inf
            training_error = np.Inf
            crossval_error = np.Inf

        # print(poly_order, x_polynomial_data.shape[0], x_polynomial_data.shape[1], training_error, crossval_error)

        return phi_vector, training_error, crossval_error

    def surrogate_performance(
        self, phi_best, order_best, additional_features_array=None
    ):
        """

        This function evaluates the performance of the surrogate model on the entire dataset.
            1. A vector is created to hold the original input data and the predicted y values from the surrogate model is created - comparison_vector
            2. The predicted values from the surrogate model are then evaluated.
            3. The errors on each datapoint(individual error), the mean absolute error and the mean square errors are calculated.
            4. The R-square coefficient is then calculated.
            5. The adjusted R2 is calculated next, taking into account the number of terms in the equation

        The comparison vector is sorted based on the performance of the surrogate model in its prediction - best to worst.
        Note that the error on each data point is based on the error maximization function in ALAMO (Cozad et al., Eq. 7)

        """

        comparison_vector = np.zeros(
            (self.original_data.shape[0], self.original_data.shape[1] + 1)
        )
        comparison_vector[:, : self.original_data.shape[1]] = self.original_data[:, :]

        # Create x terms for the whole input data, and evaluate the predicted y's as phi.X.
        x_evaluation_data = self.polygeneration(
            order_best,
            self.multinomials,
            self.original_data[:, 0 : self.original_data.shape[1] - 1],
            additional_features_array,
        )
        y_prediction = np.matmul(x_evaluation_data, phi_best)
        y_prediction = y_prediction.reshape(y_prediction.shape[0], 1)
        comparison_vector[:, self.original_data.shape[1]] = y_prediction[:, 0]

        # Error calculations:
        den = np.max(comparison_vector[:, -2]) - np.min(comparison_vector[:, -2])
        individual_error = (
            (comparison_vector[:, -1] - comparison_vector[:, -2]) / den
        ) ** 2
        mae_error = (1 / comparison_vector.shape[0]) * np.sum(
            np.abs(comparison_vector[:, -1] - comparison_vector[:, -2])
        )
        mse_error = (1 / comparison_vector.shape[0]) * np.sum(
            (comparison_vector[:, -1] - comparison_vector[:, -2]) ** 2
        )
        # R-squared coefficient calculation:
        input_y_mean = np.mean(comparison_vector[:, -2], axis=0)
        ss_total = np.sum((comparison_vector[:, -2] - input_y_mean) ** 2)
        ss_residual = np.sum((comparison_vector[:, -1] - comparison_vector[:, -2]) ** 2)
        r_square = 1 - (ss_residual / ss_total)

        # Sort comparison vector based on error in predictions
        individual_error = individual_error.reshape(individual_error.shape[0], 1)
        comparison_vector = np.append(comparison_vector, individual_error, 1)
        sorted_comparison_vector = comparison_vector[comparison_vector[:, -1].argsort()]

        # Adjusted R_squared
        samp_size = self.original_data.shape[0]
        no_nonzero_terms = np.count_nonzero(phi_best[1:, 0])
        # Evaluate R2_adjusted only if R2>0: Fit is better than mean value. When fit is worse than hor. line, return 0
        if r_square > 0:
            r2_adjusted = 1 - (
                (1 - r_square) * ((samp_size - 1) / (samp_size - no_nonzero_terms - 1))
            )
        else:
            r2_adjusted = 0

        return sorted_comparison_vector, mae_error, mse_error, r_square, r2_adjusted

    def results_generation(self, beta, order):
        """
        This function prints the results of the fitting to the screen.
        """
        results_df = pd.Series(dtype="float64")
        counter = 1
        print("\n------------------------------------------------------------")
        print("The final coefficients of the regression terms are: \n")
        print("k               |", beta[0, 0])
        results_df = pd.concat([results_df, pd.Series({"k": beta[0, 0]})], axis=0)
        if self.multinomials == 1:
            for i in range(1, order + 1):
                for j in range(1, self.number_of_x_vars + 1):
                    print("(x_", j, ")^", i, "     |", beta[counter, 0])
                    col_name = "(x_" + str(j) + ")^" + str(i)
                    results_df = pd.concat(
                        [results_df, pd.Series({col_name: beta[counter, 0]})], axis=0
                    )
                    counter += 1
            for i in range(1, self.number_of_x_vars + 1):
                for j in range(1, self.number_of_x_vars + 1):
                    if i > j:
                        print("x_", j, ".x_", i, "     |", beta[counter, 0])
                        col_name = "(x_" + str(j) + ")" + ".(x_" + str(i) + ")"
                        results_df = pd.concat(
                            [results_df, pd.Series({col_name: beta[counter, 0]})],
                            axis=0,
                        )
                        counter += 1

        else:
            for i in range(1, order + 1):
                for j in range(1, self.number_of_x_vars + 1):
                    print("(x_", j, ")^", i, "     |", beta[counter, 0])
                    col_name = "(x_" + str(j) + ")^" + str(i)
                    results_df = pd.concat(
                        [results_df, pd.Series({col_name: beta[counter, 0]})], axis=0
                    )
                    counter += 1

        return results_df

    @staticmethod
    def error_plotting(vector_of_results):
        """
        This function generates displays a plot of the different errors
        """
        ax1 = plt.subplot(2, 2, 1)
        ax1.plot(
            vector_of_results[:, 0],
            vector_of_results[:, 2],
            "green",
            vector_of_results[:, 0],
            vector_of_results[:, 3],
            "red",
        )
        ax1.set_title("Training (green) vs Cross-validation error (red)")

        ax2 = plt.subplot(2, 2, 2)
        ax2.plot(vector_of_results[:, 0], vector_of_results[:, 4], "green")
        ax2.set_title("MAE")

        ax3 = plt.subplot(2, 2, 3)
        ax3.plot(vector_of_results[:, 0], vector_of_results[:, 5], "blue")
        ax3.set_title("MSE")

        ax4 = plt.subplot(2, 2, 4)
        ax4.plot(
            vector_of_results[:, 0],
            vector_of_results[:, 6],
            "blue",
            vector_of_results[:, 0],
            vector_of_results[:, 7],
            "red",
        )
        ax4.set_title("R-squared (blue) and Adjusted R-squared (red)")

        plt.show()

        return ax1, ax2, ax3, ax4

    def user_defined_terms(self, additional_regression_features):
        """

        This function generates a 2D array of the additional features from the list supplied by the user.
        Note: It assumes that each list element is 1D

        Args:
            additional_regression_features(list): a list of features to be added to the regression problem. Each element of the list must have the same number of entries as self.number_of_samples

        Returns:
            additional_features_array(NumPy Array): an array of additional training features with len(additional_regression_features) columns to be considered during regression.

        Raises:
            Exception:
                * when additional_regression_features is not a list
            Exception:
                * when the entries in additional_regression_features are not of type 1-D NumPy Array or Pandas Series
            Exception:
                * when the length of the entries in additional_regression_features do not match the number of rows in self.regression_data

        """
        # Check for list
        if not isinstance(additional_regression_features, list):
            raise ValueError("additional_regression_features: list required.")
        # Determine number of additional features, and create 2D array for additional features
        number_additional_features = len(additional_regression_features)
        additional_features_array = np.zeros(
            (self.regression_data.shape[0], number_additional_features)
        )
        for i in range(0, number_additional_features):
            # If entry is an array and is of the right sizes
            if (
                isinstance(additional_regression_features[i], np.ndarray)
                and (
                    len(additional_regression_features[i])
                    == self.regression_data.shape[0]
                )
                and (additional_regression_features[i].ndim == 1)
            ):
                additional_features_array[:, i] = additional_regression_features[i]
            elif (
                isinstance(additional_regression_features[i], pd.DataFrame)
                and (
                    len(additional_regression_features[i])
                    == self.regression_data.shape[0]
                )
                and (additional_regression_features[i].ndim == 1)
            ):
                additional_features_array[:, i] = additional_regression_features[
                    i
                ].values
            elif (
                isinstance(additional_regression_features[i], pd.Series)
                and (
                    len(additional_regression_features[i])
                    == self.regression_data.shape[0]
                )
                and (additional_regression_features[i].ndim == 1)
            ):
                additional_features_array[:, i] = additional_regression_features[
                    i
                ].values
            else:
                raise Exception(
                    "Wrong data dimensions or type - additional_regression_features contain 1-D vectors, have same number of entries as regression_data and be of type pd.Series, pd.Dataframe or np.ndarray."
                )
        return additional_features_array

    def polynomial_regression_fitting(self, additional_regression_features=None):
        """

        polynomial_regression_fitting is the core method which is called in the PolynomialRegression class.
        It ties together all the other functions in the class.

        For each polynomial order, it
                 - calls the function user_defined_terms to generate the array of additional features (when required),
                 - calls the function training_test_data_creation to generate the training and test data sets,
                 - calls the function polyregression to determine the optimal weight vector and the fitting errors,
                 - determines whether the new fit improves is the best so far by the crossvalidation error of the current fit to the previous best,
                 - calls the function surrogate_performance to calculate the errors and R-values of the current fit, and
                 - returns results to user.

        When adaptive sampling is done, the function also
         - selects the adaptive samples to be added to the training data based on the magnitudes of the prediction errors of individual samples in self.original_data, and
         - determines when the the stopping conditions have been satisfied.


        The following stopping conditions are considered when adaptive sampling is in use:
         - The maximum number of training samples allowed by the user has been exceeded
         - Mean absolute error ~= 0
         - Mean squared error ~= 0
         - R^2 = 1
         - The preset iteration number given by the user has been reached
         - All available points in self.original_data have been used for training.

         Keyword Args:
            additional_regression_features(<list>):     Additional features the user wants the algorithm to consider during regression.
                                                        It should be noted that adaptive sampling is not available when additional features have been supplied by the user, i.e. when len(additional_regression_features) > 0.

         Returns:
            results:                                    Python object containing the results of the polynomial regression process including the polynomial order
                                                        (results.polynomial_order), polynomial coefficients (results.optimal_weights_array) and fit errors (results.errors).
                                                        See information on ResultReport class for details on contents.

        """
        # Parameters that represent the best solution found at each iteration based on the cross-validation error
        best_error = 1e20
        train_error_fit = 1e20
        phi_best = 0
        order_best = 0

        if (additional_regression_features is None) or (
            len(additional_regression_features) == 0
        ):
            print(
                "max_fraction_training_samples set at ",
                self.max_fraction_training_samples,
            )
            print(
                "Number of adaptive samples (no_adaptive_samples) set at ",
                self.no_adaptive_samples,
            )
            print("Maximum number of iterations (Max_iter) set at: ", self.max_iter)

            training_data, cross_val_data = self.training_test_data_creation()
            for poly_order in range(1, self.max_polynomial_order + 1):
                for cv_number in range(1, self.number_of_crossvalidations + 1):
                    phi, train_error, cv_error = self.polyregression(
                        poly_order,
                        training_data["training_set_" + str(cv_number)],
                        cross_val_data["test_set_" + str(cv_number)],
                    )
                    if cv_error < best_error:
                        best_error = cv_error
                        phi_best = phi
                        order_best = poly_order
                        train_error_fit = train_error
            print(
                "\nInitial surrogate model is of order",
                order_best,
                " with a cross-val error of %4f" % best_error,
            )
            # Next, Calculate and report errors.
            (
                sorted_comparison_vector,
                mae_error,
                mse_error,
                r_square,
                r_square_adj,
            ) = self.surrogate_performance(phi_best, order_best)
            print(
                "Initial Regression Model Performance:\nOrder: ",
                order_best,
                " / MAE: %4f" % mae_error,
                " / MSE: %4f" % mse_error,
                " / R^2: %4f" % r_square,
                " / Adjusted R^2: %4f" % r_square_adj,
            )

            # Parameters that retain the previous best solutions. They are compared to the best solution at each iteration based on the R-square coefficient.
            (
                order_opt,
                train_error_opt,
                best_error_opt,
                mae_error_opt,
                mse_error_opt,
                r_square_opt,
                r_square_adj_opt,
                phi_opt,
            ) = (
                order_best,
                train_error_fit,
                best_error,
                mae_error,
                mse_error,
                r_square,
                r_square_adj,
                phi_best,
            )

            eps_neg = 1e-6
            eps_pos = 0.999999
            iteration_number = 1
            stopping_criterion = int(
                np.ceil(
                    self.max_fraction_training_samples * self.original_data.shape[0]
                )
            )
            vector_of_results = np.zeros((stopping_criterion, 9))
            while (
                (self.regression_data.shape[0] < stopping_criterion)
                and (mae_error > eps_neg)
                and (mse_error > eps_neg)
                and (r_square < eps_pos)
                and (iteration_number < self.max_iter)
                and (
                    self.regression_data.shape[0] + self.no_adaptive_samples
                    < self.original_data.shape[0]
                )
            ):
                print("\n-------------------------------------------------")
                print("\nIteration ", iteration_number)
                best_error = 1e20

                # Select n_adaptive_samples worst fitting points to be added to the dataset used in the previous evaluation.
                scv_input_data = sorted_comparison_vector[:, :-2]
                sorted_comparison_vector_unique = scv_input_data[
                    np.all(
                        np.any(
                            (scv_input_data - self.regression_data[:, None]), axis=2
                        ),
                        axis=0,
                    )
                ]
                adaptive_samples = sorted_comparison_vector_unique[
                    # PYLINT-WHY: pylint considers self.no_adaptive_samples to be None here
                    # pylint: disable=invalid-unary-operand-type
                    -self.no_adaptive_samples :,
                    :
                    # pylint: enable=invalid-unary-operand-type
                ]
                self.regression_data = np.concatenate(
                    (self.regression_data, adaptive_samples), axis=0
                )
                self.number_of_samples = self.regression_data.shape[
                    0
                ]  # Never forget to update
                print(
                    "\n",
                    self.no_adaptive_samples,
                    " additional points added to training data. New number of training samples: ",
                    self.regression_data.shape[0],
                )

                training_data, cross_val_data = self.training_test_data_creation()

                for poly_order in range(1, self.max_polynomial_order + 1):
                    for cv_number in range(1, self.number_of_crossvalidations + 1):
                        phi, train_error, cv_error = self.polyregression(
                            poly_order,
                            training_data["training_set_" + str(cv_number)],
                            cross_val_data["test_set_" + str(cv_number)],
                        )
                        if cv_error < best_error:
                            best_error = cv_error
                            phi_best = phi
                            order_best = poly_order
                            train_error_fit = train_error
                print(
                    "\nThe best regression model is of order",
                    order_best,
                    " with a cross-val error of %4f" % best_error,
                )

                (
                    sorted_comparison_vector,
                    mae_error,
                    mse_error,
                    r_square,
                    r_square_adj,
                ) = self.surrogate_performance(phi_best, order_best)
                print(
                    "Regression performance on full data in iteration",
                    iteration_number,
                    "\nOrder: ",
                    order_best,
                    " / MAE: %4f" % mae_error,
                    " / MSE: %4f" % mse_error,
                    " / R_sq: %4f" % r_square,
                    " / Adjusted R^2: %4f" % r_square_adj,
                )

                # Determine if solution is improved. If yes, update solution. if no, retain previous best.
                if r_square_adj > r_square_adj_opt:
                    (
                        phi_opt,
                        order_opt,
                        mae_error_opt,
                        mse_error_opt,
                        r_square_opt,
                        r_square_adj_opt,
                        train_error_opt,
                        best_error_opt,
                    ) = (
                        phi_best,
                        order_best,
                        mae_error,
                        mse_error,
                        r_square,
                        r_square_adj,
                        train_error_fit,
                        best_error,
                    )
                    print("New solution found.")
                else:
                    print("Previous solution retained.")
                vector_of_results[iteration_number, :] = [
                    iteration_number,
                    order_opt,
                    train_error_opt,
                    best_error_opt,
                    mae_error_opt,
                    mse_error_opt,
                    r_square_opt,
                    r_square_adj_opt,
                    self.regression_data.shape[0],
                ]
                iteration_number += 1

            # Remove all zero rows in the solution vector
            vector_of_results = vector_of_results[
                ~np.all(vector_of_results == 0, axis=1)
            ]
            # Round phi to 2.d.p and print results to screen
            beta_vector = np.round(phi_opt, 6)
            if r_square_adj_opt < 0.95:
                print("\nPolynomial regression performs poorly for this dataset.")
            else:
                print(
                    "\nPolynomial regression generates a good surrogate model for the input data."
                )
            if iteration_number > 1:
                _, _, _, _ = self.error_plotting(vector_of_results)
            print(
                "\n-------------------------------------------------\n-------------------------------------------------"
            )
            print(
                "Best solution found: ",
                "\nOrder: ",
                order_opt,
                " / MAE: %4f" % mae_error_opt,
                " / MSE: %4f" % mse_error_opt,
                " / R_sq: %4f" % r_square_opt,
                " / Adjusted R^2: %4f" % r_square_adj_opt,
            )
            dataframe_coeffs = self.results_generation(beta_vector, order_opt)

            vector_of_results_df = pd.DataFrame(
                {
                    "Iteration_number": vector_of_results[:, 0],
                    "Polynomial order": vector_of_results[:, 1],
                    "Training error": vector_of_results[:, 2],
                    "Cross-val error": vector_of_results[:, 3],
                    "MAE": vector_of_results[:, 4],
                    "MSE": vector_of_results[:, 5],
                    "R2": vector_of_results[:, 6],
                    "Adjusted R2": vector_of_results[:, 7],
                    "Number of training samples": vector_of_results[:, 8],
                }
            )

            extra_terms_feature_vector = list(
                self.feature_list[i] for i in self.regression_data_columns
            )

            # Results
            self.optimal_weights_array = phi_opt
            self.final_polynomial_order = order_opt
            self.errors = {
                "MAE": mae_error_opt,
                "MSE": mse_error_opt,
                "R2": r_square_opt,
                "Adjusted R2": r_square_adj_opt,
            }
            self.number_of_iterations = iteration_number
            self.iteration_summary = vector_of_results_df
            self.additional_features_data = None
            self.final_training_data = self.regression_data

            self.dataframe_of_optimal_weights_polynomial = dataframe_coeffs
            self.dataframe_of_optimal_weights_extra_terms = []
            self.extra_terms_feature_vector = extra_terms_feature_vector
            if r_square_opt > 0.95:
                self.fit_status = "ok"
            else:
                warnings.warn(
                    "Polynomial regression generates poor fit for the dataset"
                )
                self.fit_status = "poor"

            self.pickle_save({"model": self})
            return self

        else:
            print("No iterations will be run.")
            # Determine number of additional features based on length of list, and convert list into array
            number_additional_features = len(additional_regression_features)
            additional_features_array = self.user_defined_terms(
                additional_regression_features
            )

            training_data, cross_val_data = self.training_test_data_creation(
                additional_features_array
            )
            for poly_order in range(1, self.max_polynomial_order + 1):
                for cv_number in range(1, self.number_of_crossvalidations + 1):
                    phi, train_error, cv_error = self.polyregression(
                        poly_order,
                        training_data["training_set_" + str(cv_number)],
                        cross_val_data["test_set_" + str(cv_number)],
                        training_data["training_extras_" + str(cv_number)],
                        cross_val_data["test_extras_" + str(cv_number)],
                    )
                    if cv_error < best_error:
                        best_error = cv_error
                        phi_best = phi
                        order_best = poly_order
                        train_error_fit = train_error
            print(
                "\nBest surrogate model is of order",
                order_best,
                " with a cross-val S.S. Error  of %4f" % best_error,
            )

            # KEY: Modification of self variable outside initialization. Required to make @surrogate_performance work here.
            self.original_data = self.regression_data
            _, mae_error, mse_error, r_square, _ = self.surrogate_performance(
                phi_best, order_best, additional_features_array
            )

            # Round solution to 6.d.p
            beta_vector = np.round(phi_best, 6)

            # Print results to screen
            dataframe_coeffs = self.results_generation(beta_vector, order_best)

            extra_terms_coeffs = pd.Series(dtype="float64")
            print(
                "\nThe coefficients of the extra terms in additional_regression_features are:\n"
            )
            for af in range(number_additional_features, 0, -1):
                print(
                    "Coeff. additional_regression_features[",
                    number_additional_features - af + 1,
                    "]: ",
                    beta_vector[len(beta_vector) - af, 0],
                )
                col_name = (
                    "Coeff. additional_regression_features["
                    + str(number_additional_features - af + 1)
                    + "]"
                )
                extra_terms_coeffs = pd.concat(
                    [
                        extra_terms_coeffs,
                        pd.Series({col_name: beta_vector[len(beta_vector) - af, 0]}),
                    ],
                    axis=0,
                )

            # Print errors
            print(
                "\nRegression model performance on training data:\nOrder: ",
                order_best,
                " / MAE: %4f" % mae_error,
                " / MSE: %4f" % mse_error,
                " / R^2: %4f" % r_square,
            )

            extra_terms_feature_vector = list(
                self.feature_list[i] for i in self.regression_data_columns
            )

            # Results
            self.optimal_weights_array = phi_best
            self.final_polynomial_order = order_best
            self.errors = {"MAE": mae_error, "MSE": mse_error, "R2": r_square}
            self.number_of_iterations = []
            self.iteration_summary = []
            self.additional_features_data = additional_features_array
            self.final_training_data = self.regression_data
            self.dataframe_of_optimal_weights_polynomial = dataframe_coeffs
            self.dataframe_of_optimal_weights_extra_terms = extra_terms_coeffs
            self.extra_terms_feature_vector = extra_terms_feature_vector
            if r_square > 0.95:
                self.fit_status = "ok"
            else:
                warnings.warn(
                    "Polynomial regression generates poor fit for the dataset"
                )
                self.fit_status = "poor"

            self.pickle_save({"model": self})

            return self

    def get_feature_vector(self):
        """

        The ``get_feature_vector`` method generates the list of regression features from the column headers of the input dataset.

        Returns:
            Pyomo IndexedParam  : An indexed parameter list of the variables supplied in the original data


        **Example:**

        .. code-block:: python

            # Create a small dataframe with three columns ('one', 'two', 'three') and two rows (A, B)
            >>> xy_data = pd.DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6])], orient='index', columns=['one', 'two', 'three'])

            # Initialize the **PolynomialRegression** class and print the column headers for the variables
            >>> f = PolynomialRegression(xy_data, xy_data, maximum_polynomial_order=1, multinomials=True, training_split=0.8)
            >>> p = f.get_feature_vector()
            >>> for i in p.keys():
            >>>     print(i)
            one
            two

        """
        p = Param(self.regression_data_columns, mutable=True, initialize=0)
        p.index_set().construct()
        p.construct()
        self.feature_list = p
        return p

    def set_additional_terms(self, term_list):
        """

        ``set_additional_terms`` accepts additional user-defined features for consideration during regression.

        Args:
            term_list (list) : List of additional terms to be considered as regression features. Each term in the list must be a Pyomo-supported intrinsic function.


        **Example:**

        .. code-block:: python

            # To add the sine and cosine of a variable with header 'X1' in the dataset as additional regression features:
            >>> xy_data = pd.DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6])], orient='index', columns=['X1', 'X2', 'Y'])
            >>> A = PolynomialRegression(xy_data, xy_data, maximum_polynomial_order=5)
            >>> p = A.get_feature_vector()
            >>> A.set_additional_terms([ pyo.sin(p['X1']) , pyo.cos(p['X1']) ])

        """
        self.additional_term_expressions = term_list

    # def fit_surrogate(self):
    def training(self):
        """

        The ``training`` method trains a polynomial model to an input dataset.
        It calls the core method which is called in the PolynomialRegression class (polynomial_regression_fitting).
        It accepts no user input, inheriting the information passed in class initialization.

        Returns:
            tuple   : Python Object (**results**) containing the results of the polynomial regression process including:
                - the polynomial order  (**self.final_polynomial_order**)
                - polynomial coefficients (**self.optimal_weights_array**), and
                - MAE and MSE errors as well as the :math:`R^{2}` (**results.errors**).

        """

        cMap = ComponentMap()
        for i, col in enumerate(self.regression_data_columns):
            cMap[self.feature_list[col]] = self.regression_data[:, i]
        npe = NumpyEvaluator(cMap)
        additional_data = list(
            npe.walk_expression(term) for term in self.additional_term_expressions
        )
        return self.polynomial_regression_fitting(additional_data)

    def generate_expression(self, variable_list):
        """

        The ``generate_expression`` method returns the Pyomo expression for the polynomial model trained.

        The expression is constructed based on a supplied list of variables **variable_list** and the output of ``training``.

        Args:
            variable_list(list)           : List of input variables to be used in generating expression. This can be the a list generated from the results of ``get_feature_vector``. The user can also choose to supply a new list of the appropriate length.

        Returns:
            Pyomo Expression              : Pyomo expression of the polynomial model based on the variables provided in **variable_list**.

        """
        # Reshaping of array necessary when input variables are Pyomo scalar variables
        vl = np.array([variable_list])
        vl = vl.reshape(1, len(variable_list)) if vl.ndim > 2 else vl

        terms = PolynomialRegression.polygeneration(
            self.final_polynomial_order, self.multinomials, vl
        ).transpose()
        n = len(terms)

        ans = sum(
            w * t
            for w, t in zip(
                np.nditer(self.optimal_weights_array),
                np.nditer(terms, flags=["refs_ok"]),
            )
        )

        user_term_map = dict(
            (id(a), b)
            for a, b in zip(
                self.extra_terms_feature_vector,
                variable_list,
            )
        )
        if len(self.additional_term_expressions) > 0:
            for w, expr in zip(
                np.nditer(self.optimal_weights_array[n:]),
                self.additional_term_expressions,
            ):
                ans += float(w) * replace_expressions(expr, user_term_map)
        return ans

    def predict_output(self, x_data):
        """

        The ``predict_output`` method generates output predictions for input data x_data based a previously generated polynomial fitting.

        Args:
            x_data          : Numpy array of designs for which the output is to be evaluated/predicted.

        Returns:
             Numpy Array    : Output variable predictions based on the polynomial fit.

        """
        nf = x_data.shape[1]
        x_list = [i for i in range(0, nf)]
        import pyomo.environ as aml

        m = aml.ConcreteModel()
        i = aml.Set(initialize=x_list)
        m.xx = aml.Var(i)
        m.o2 = aml.Objective(expr=self.generate_expression([m.xx[i] for i in x_list]))
        y_eq = np.zeros((x_data.shape[0], 1))
        for j in range(0, x_data.shape[0]):
            for i in x_list:
                m.xx[i] = x_data[j, i]
            y_eq[j, 0] = aml.value(m.o2([m.xx[i] for i in x_list]))
        return y_eq

    def pickle_save(self, solutions):
        """
        The training method saves the results of the run in a pickle object. It saves an object with two elements: the setup (index[0]) and the results (index[1]).
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
        pickle_load loads the results of a saved run 'file.obj'. It returns an array of two elements: the setup (index[0]) and the results (index[1]).

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
        y_predicted = self.predict_output(self.final_training_data[:, :-1])
        fig1 = plt.figure(figsize=(16, 9), tight_layout=True)
        ax = fig1.add_subplot(121)
        ax.plot(self.final_training_data[:, -1], self.final_training_data[:, -1], "-")
        ax.plot(self.final_training_data[:, -1], y_predicted, "o")
        ax.set_xlabel(r"True data", fontsize=12)
        ax.set_ylabel(r"Surrogate values", fontsize=12)
        ax.set_title(r"Parity plot", fontsize=12)

        ax2 = fig1.add_subplot(122)
        ax2.plot(
            self.final_training_data[:, -1],
            self.final_training_data[:, -1]
            - y_predicted[:,].reshape(
                y_predicted.shape[0],
            ),
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
            f"\nResults of polynomial regression run:\n"
            f"\nPolynomial order                   : {self.final_polynomial_order}\n"
            f"Number of terms in polynomial model: {self.optimal_weights_array.size}\n"
            f"\nPolynomial Expression:\n"
            f"--------------------------\n"
            f"\n{eqn}\n"
            f"--------------------------\n"
            f"\nModel training errors:"
            f"\n-----------------------\n"
            f"Mean Squared Error (MSE)         : {self.errors['MSE']}\n"
            f"Root Mean Squared Error (RMSE)   : {np.sqrt(self.errors['MSE'])}\n"
            f"Mean Absolute error (MSE)        : {self.errors['MAE']}\n"
            f"Goodness of fit (R2)             : {self.errors['R2']}\n"
            f"\n{double_line}"
        )
        return s

    def print_report(self):
        s = self._report()
        print(s)

    def _repr_pretty_(self, p, cycle=False):

        s = self._report()
        p.text(s)

    def confint_regression(self, confidence=0.95):
        """
        The ``confint_regression`` method prints the confidence intervals for the regression patamaters.

        Args:
            confidence      : Required confidence interval level, default = 0.95 (95%)

        """
        from scipy.stats import t

        data = self.final_training_data
        y_pred = self.predict_output(data[:, :-1])
        dof = (
            data.shape[0] - len(self.optimal_weights_array) + 1
        )  # be careful when there are additional features
        ssr = np.sum((data[:, -1] - y_pred[:, 0]) ** 2)
        sig_sq = ssr / dof
        if (self.additional_features_data is None) or (
            len(self.additional_features_data) == 0
        ):
            x_exp = self.polygeneration(
                self.final_polynomial_order, self.multinomials, data[:, :-1]
            )  # will not account for additional features
        else:
            x_exp = self.polygeneration(
                self.final_polynomial_order,
                self.multinomials,
                data[:, :-1],
                additional_x_training_data=self.additional_features_data,
            )

        covar = sig_sq * np.linalg.pinv(x_exp.transpose() @ x_exp)
        ss_reg_params = np.sqrt(
            np.diag(covar)
        )  # standard error for each regression parameter
        t_dist = t.ppf(
            (1 + confidence) / 2, dof
        )  # alternatively, t_dist_data = st.t.interval(0.99, 8)
        # Evaluate confidence intervals, Tabulate and print results
        c_data = np.zeros((self.optimal_weights_array.shape[0], 4))
        c_data[:, 0] = self.optimal_weights_array[:, 0]
        c_data[:, 1] = ss_reg_params[
            :,
        ]
        c_data[:, 2] = (
            self.optimal_weights_array[:, 0]
            - t_dist
            * ss_reg_params[
                :,
            ]
        )
        c_data[:, 3] = (
            self.optimal_weights_array[:, 0]
            + t_dist
            * ss_reg_params[
                :,
            ]
        )

        headers = [
            "Regression coeff.",
            "Std. error",
            "Conf. int. lower",
            "Conf. int. upper",
        ]
        c_data_df = pd.DataFrame(data=c_data, columns=headers)
        print(c_data_df)
        return c_data_df
