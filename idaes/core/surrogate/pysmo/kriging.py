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

__author__ = "Oluwamayowa Amusat"

# Imports from the python standard library
import os.path

# Imports from third parties
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import pickle
from pyomo.core import Param, exp
from scipy.optimize import basinhopping
import scipy.optimize as opt

# Imports from IDAES namespace
from idaes.core.surrogate.pysmo.sampling import FeatureScaling as fs


class MyBounds(object):
    """
    The Class MyBounds tests whether the reguularization parameter value is within the expected range.
     The class is initialized with the preset valies in __init__; the __call__ function returns Booleans indicating whether the regularization parameter value is acceptable.
     The results of the __call__ function is fed into the Basinhopping algorithm using the accept_test parameter.
    """

    def __init__(self, xmax=[1], xmin=[1e-6]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(x[-1] <= self.xmax)
        tmin = bool(x[-1] >= self.xmin)
        return tmax and tmin


class KrigingModel:
    """
    The KrigingModel class trains a Kriging model for a training data set.

    The class must first be initialized by calling **KrigingModel**. Model training is then carried out by calling the ``training`` method.

    **KrigingModel** is able to generate either an interpolating or a regressing Kriging model depending on the settings used during initialization..

    **Example:**

    .. code-block:: python

        # Initialize the class
        >>> d = KrigingModel(training_data, numerical_gradients=True, regularization=True))
        >>> p = d.get_feature_vector()

        # Train Kriging model and predict output for an test data x_test
        >>> d.training()
        >>> predictions = d.predict_output(x_test)

    Args:
        XY_data (NumPy Array or Pandas Dataframe)   : The dataset for Kriging training. **XY_data** is expected to contain both the features (X) and output (Y) information, with the output values (Y) in the last column.

    Further details about the optional inputs may be found under the ``__init__`` method.

    """

    def __init__(
        self,
        XY_data,
        numerical_gradients=True,
        regularization=True,
        fname=None,
        overwrite=False,
    ):
        """
        Initialization of **KrigingModel** class.

        Args:
            XY_data (NumPy Array or Pandas Dataframe)   : The dataset for Kriging training. **XY_data** is expected to contain feature and output data, with the output values (y) in the last column.

        Keyword Args:
            numerical_gradients(bool)               : Whether or not numerical gradients should be used in training. This choice determines the algorithm used to solve the problem.

                                                            - numerical_gradients = True: The problem is solved with BFGS using central differencing with a step size of :math:`10^{-6}` to evaluate numerical gradients.
                                                            - numerical_gradients = False: The problem is solved with Basinhopping, a stochastic optimization algorithm.

            regularization(bool)                    :  This option determines whether or not regularization is considered during Kriging training. Default is True.

                                                            - When regularization is turned off, the model generates an interpolating kriging model.

        Returns:
            self object with the input information and settings.

        Raises:
            ValueError: - The input dataset is of the wrong type (not a NumPy array or Pandas Dataframe)

            Exception:  - numerical_gradients is not boolean

            Exception:  - regularization is not boolean

        **Example:**

        .. code-block:: python

                # Initialize Kriging class with no numerical gradients - solution algorithm will be Basinhopping
                >>> d = KrigingModel(XY_data, numerical_gradients=False)

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

        self.x_data = xy_data[:, :-1].reshape(xy_data.shape[0], xy_data.shape[1] - 1)
        self.y_data = xy_data[:, -1].reshape(xy_data.shape[0], 1)
        self.num_vars = self.x_data.shape[1] + 1  # thetas and reg parameter only
        x_data_scaled, self.x_data_min, self.x_data_max = fs.data_scaling_minmax(
            self.x_data
        )
        self.x_data_scaled = x_data_scaled.reshape(self.x_data.shape)

        if isinstance(numerical_gradients, bool):
            self.num_grads = numerical_gradients
        else:
            raise Exception("numerical_gradients must be boolean.")

        if isinstance(regularization, bool):
            self.regularization = regularization
        else:
            raise Exception("Choice of regularization must be boolean.")

        # Results
        self.optimal_weights = None
        self.optimal_p = None
        self.optimal_mean = None
        self.optimal_variance = None
        self.regularization_parameter = None
        self.optimal_covariance_matrix = None
        self.covariance_matrix_inverse = None
        self.optimal_y_mu = None
        self.output_predictions = None
        self.training_R2 = None
        self.training_rmse = None

    @staticmethod
    def covariance_matrix_generator(x, theta, reg_param, p):
        """
        The covariance_matrix_generator method generates the regularized co-variance matrix for a Kriging model

        Args:
            x                       : scaled features data
            theta                   : Kriging weights
            reg_param               : regularization parameter
            p                       : Kriging exponent, fixed at 2 for smoothness.

        Returns:
            cov_matrix              : Regularized co-variance matrix

        """
        distance_matrix = np.zeros((x.shape[0], x.shape[0]))
        for i in range(0, x.shape[0]):
            distance_matrix[i, :] = (
                np.matmul(((np.abs(x[i, :] - x)) ** p), theta)
            ).transpose()
        cov_matrix = np.exp(-1 * distance_matrix)
        cov_matrix = cov_matrix + reg_param * np.eye(
            cov_matrix.shape[0]
        )  # Regularization parameter addition, see Forrester book
        return cov_matrix

    @staticmethod
    def covariance_inverse_generator(x):
        """
        The covariance_inverse_generator method generates the inverse of the regularized co-variance matrix for a Kriging model

        Args:
            x                       : Regularized co-variance matrix

        Returns:
            inverse_x               : Inverse of regularized co-variance matrix

        """
        try:
            inverse_x = np.linalg.inv(x)
        except np.linalg.LinAlgError as LAE:
            inverse_x = np.linalg.pinv(x)
        return inverse_x

    @staticmethod
    def kriging_mean(cov_inv, y):
        """
        The kriging_mean method calculates the MLE estimate of the mean.

        Args:
            cov_inv (NumPy Array)           : Inverse of the co-variance matrix
            y (NumPy Array)                 : Output values of the training data


        Returns:
            kriging_mean                    : MLE estimate of the Kriging mean

        Reference:
            [1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide",
                https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801

        """
        ones_vec = np.ones((y.shape[0], 1))
        kriging_mean = np.matmul(
            np.matmul(ones_vec.transpose(), cov_inv), y
        ) / np.matmul(np.matmul(ones_vec.transpose(), cov_inv), ones_vec)
        # kriging_mean = np.matmul(ones_vec.transpose(), np.matmul(cov_inv, y)) / np.matmul(ones_vec.transpose(), np.matmul(cov_inv, ones_vec))
        return kriging_mean

    @staticmethod
    def y_mu_calculation(y, mu):
        """
        The y_mu_calculation method calculates the deviation of each output value from the MLE estimate of the mean, mu

        """
        y_mu = y - mu * np.ones((y.shape[0], 1))
        return y_mu

    @staticmethod
    def kriging_sd(cov_inv, y_mu, ns):
        """
        The kriging_sd method calculates the MLE estimate of the Kriging variance.

        Args:
            cov_inv (NumPy Array)           : Inverse of the co-variance matrix
            y_mu (NumPy Array)              : Deviation of y from the Kriging mean estimate (y-mean)


        Returns:
            kriging_sd                      : MLE estimate of the Kriging variance

        Reference:
            [1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide",
                https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801

        """
        sigma_sq = np.matmul(np.matmul(y_mu.transpose(), cov_inv), y_mu) / ns
        # sigma_sq = np.matmul(y_mu.transpose(), np.matmul(cov_inv, y_mu)) / ns
        return sigma_sq

    @staticmethod
    def print_fun(x, f, accepted):
        print("at minimum %.4f accepted %d" % (f, int(accepted)))

    def objective_function(self, var_vector, x, y, p):
        """
        The objective_function method calculates the concentrated likelihood function

        Args:
            var_vector(NumPy Array)        : Numpy array containing the Kriging paramaters (Kriging weights and regularization parameter)
            x(NumPy Array)                 : Scaled version of input features/variables
            y(NumPy Array)                 : Output variable y (unscaled)
            p(float)                      : Kriging model exponent (fixed to 2) to ensure model smoothness

        Returns:
            conc_log_like(float)          : Concentrated likelihood value. Function incurs a large penalty (10000) when co-variance matrix is non-positive definite

        Reference:
            [1] Forrester et al.'s book "Engineering Design via Surrogate Modelling: A Practical Guide",
                https://onlinelibrary.wiley.com/doi/pdf/10.1002/9780470770801

        """
        theta = var_vector[:-1]
        reg_param = var_vector[-1]
        theta = 10**theta  # Assumes log(theta) provided
        ns = y.shape[0]
        cov_mat = self.covariance_matrix_generator(x, theta, reg_param, p)
        try:  # Check Cholesky factorization
            L = np.linalg.cholesky(cov_mat)
            lndetcov = 2 * np.sum(
                np.log(np.abs(np.diag(L)))
            )  # Approximation to 2nd term from Forrester book, making use of the Ch. factorization
            cov_inv = self.covariance_inverse_generator(cov_mat)
            km = self.kriging_mean(cov_inv, y)
            y_mu = self.y_mu_calculation(y, km)
            ssd = self.kriging_sd(cov_inv, y_mu, ns)
            # log_like = (0.5 * ns * np.log(ssd)) + (0.5 * np.log(np.abs(np.linalg.det(cov_mat))))
            log_like = (0.5 * ns * np.log(ssd)) + (0.5 * lndetcov)
            conc_log_like = log_like[0, 0]
        except:  # When Cholesky fails - non-positive definite covariance matrix
            conc_log_like = 1e4
        return conc_log_like

    def numerical_gradient(self, var_vector, x, y, p):
        """
        The numerical_gradient method calculates numerical gradients for the Kriging hyperparameters via central differencing,

        grad(theta) = (f(theta + eps) - f(theta - eps))/(2 * eps)

        Args:
            var_vector(NumPy Array)        : Numpy array containing the Kriging paramaters (Kriging weights and regularization parameter)
            x(NumPy Array)                 : Scaled version of input features/variables
            y(NumPy Array)                 : Output variable y (unscaled)
            p(float)                       : Kriging model exponent (fixed to 2) to ensure model smoothness

        Returns:
            grad_vec(NumPy Array)          : Array of the gradients of the variables in var_vector

        """
        eps = 1e-6
        grad_vec = np.zeros(
            len(var_vector),
        )
        for i in range(0, len(var_vector)):
            var_vector_plus = np.copy(var_vector)
            var_vector_plus[i] = var_vector[i] + eps
            var_vector_minus = np.copy(var_vector)
            var_vector_minus[i] = var_vector[i] - eps
            # of = self.objective_function(var_vector, x, y, p)
            of_plus = self.objective_function(var_vector_plus, x, y, p)
            of_minus = self.objective_function(var_vector_minus, x, y, p)
            grad_current = (of_plus - of_minus) / (2 * eps)
            grad_vec[
                i,
            ] = grad_current
        if self.regularization is False:
            grad_vec[
                -1,
            ] = 0
        return grad_vec

    def parameter_optimization(self, p):
        """
        Parameter (theta) optimization using BFGS or Basinhopping algorithm. This is the core of the Kriging Class.
        Algorithm used will depend on whether the numerical_gradients was set to True or False.
        """
        initial_value_list = np.random.randn(
            self.num_vars - 1,
        )
        initial_value_list = initial_value_list.tolist()
        initial_value_list.append(1e-4)
        initial_value = np.array(initial_value_list)
        initial_value = initial_value.reshape(initial_value.shape[0], 1)
        # Create bounds for variables. All logthetas btw (-4, 4), reg param between (1e-9, 0.1)
        bounds = []
        for i in range(0, len(initial_value_list)):
            if i == len(initial_value_list) - 1:
                if self.regularization is True:
                    bounds.append((1e-6, 0.1))
                else:
                    bounds.append((1e-10, 1e-10))
            else:
                bounds.append((-3, 3))
        bounds = tuple(bounds)

        if self.num_grads:
            print("Optimizing kriging parameters using L-BFGS-B algorithm...")
            other_args = (self.x_data_scaled, self.y_data, p)
            # opt_results = opt.minimize(self.objective_function, initial_value, args=other_args, method='L-BFGS-B', jac=self.numerical_gradient, bounds=bounds, options={'gtol': 1e-7}) #, 'disp': True})
            opt_results1 = opt.minimize(
                self.objective_function,
                initial_value,
                args=other_args,
                method="tnc",
                jac=self.numerical_gradient,
                bounds=bounds,
                options={"gtol": 1e-7},
            )
            opt_results2 = opt.minimize(
                self.objective_function,
                initial_value,
                args=other_args,
                method="L-BFGS-B",
                jac=self.numerical_gradient,
                bounds=bounds,
                options={"gtol": 1e-7},
            )  # , 'disp': True})
            if opt_results1.fun < opt_results2.fun:
                opt_results = opt_results1
            else:
                opt_results = opt_results2
        else:
            print("Optimizing Kriging parameters using Basinhopping algorithm...")
            other_args = {
                "args": (self.x_data_scaled, self.y_data, p),
                "bounds": bounds,
            }
            # other_args = {"args": (self.x_data, self.y_data, p)}
            mybounds = MyBounds()  # Bounds on regularization parameter
            opt_results = basinhopping(
                self.objective_function,
                initial_value_list,
                minimizer_kwargs=other_args,
                niter=250,
                disp=True,
                accept_test=mybounds,
            )  # , interval=5)
        return opt_results

    def optimal_parameter_evaluation(self, var_vector, p):
        """
        The optimal_parameter_evaluation method  evaluates the values of all the parameters of the final Kriging model.
        For an input set of Kriging parameters var_vector and p, it:

            (1) Generates the covariance matrix by calling covariance_matrix_generator
            (2) Finds the co-variance matrix inverse
            (3) Evaluates the Kriging mean and variance
            (4) Evaluates the deviation of each training point from the Kriging mean

        Args:
            var_vector              : Optimal Kriging parameters (weights + regularization parameter)
            p                       : Krigng exponents

        Returns:
            theta                   : Optimal Kriging weights for each variable
            reg_param               : Optimal regularization parameter
            mean                    : Final MLE estimate of the mean
            variance                : Final MLE estimate of the variance
            cov_mat                 : Co-variance matrix of the final model
            cov_inv                 : Inverse of final co-variance matrix
            y_mu                    : Deviation of each output value in the training data from the Kriging mean.

        """
        theta = var_vector[:-1]
        reg_param = var_vector[-1]
        theta = (
            10**theta
        )  # Assumes log(theta) provided. Ensures that theta is always positive
        ns = self.y_data.shape[0]
        cov_mat = self.covariance_matrix_generator(
            self.x_data_scaled, theta, reg_param, p
        )
        cov_inv = self.covariance_inverse_generator(cov_mat)
        mean = self.kriging_mean(cov_inv, self.y_data)
        y_mu = self.y_mu_calculation(self.y_data, mean)
        variance = self.kriging_sd(cov_inv, y_mu, ns)
        print(
            "\nFinal results\n================\nTheta:",
            theta,
            "\nMean:",
            mean,
            "\nRegularization parameter:",
            reg_param,
        )
        return theta, reg_param, mean, variance, cov_mat, cov_inv, y_mu

    @staticmethod
    def error_calculation(theta, p, mean, cov_inv, y_mu, x, y_data):
        """
        This method calculates the SSE and RMSE errors between the actual and predicted output values,
             ss_error = sum of squared errors / number of samples
             rmse_error = sqrt(sum of squared errors / number of samples)

        Args:
            theta           : Kriging weights
            p               : Kriging exponents
            mean            : MLE estimate of the Kriging mean
            cov_inv         : Inverse of the co-variance matrix of the current solution
            y_mu            : Deviation of y valuers from the mean estimate
            x               : Input test data
            y_data          : Actual outputs corresponding to input test data x

        Returns:
            ss_error        : The average sum of squared errors
            rmse_error      : The root-mean-squared error (RMSE)
            y_prediction    : Predicted values of y

        """
        y_prediction = np.zeros((x.shape[0], 1))
        for i in range(0, x.shape[0]):
            cmt = (np.matmul(((np.abs(x[i, :] - x)) ** p), theta)).transpose()
            cov_matrix_tests = np.exp(-1 * cmt)
            y_prediction[i, 0] = mean + np.matmul(
                np.matmul(cov_matrix_tests.transpose(), cov_inv), y_mu
            )
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

    def predict_output(self, x_pred):
        """
        The ``predict_output`` method generates output predictions for input data x_pred based a previously trained Kriging model.

        Args:
            x_pred(NumPy Array)             : Array of designs for which the output is to be evaluated/predicted.

        Returns:
             NumPy Array                    : Output variable predictions based on the Kriging model.

        """
        x_pred_scaled = (x_pred - self.x_data_min) / (self.x_data_max - self.x_data_min)
        x_pred = x_pred_scaled.reshape(x_pred.shape)
        if x_pred.ndim == 1:
            x_pred = x_pred.reshape(1, len(x_pred))
        y_pred = np.zeros((x_pred.shape[0], 1))
        for i in range(0, x_pred.shape[0]):
            cmt = (
                np.matmul(
                    ((np.abs(x_pred[i, :] - self.x_data_scaled)) ** self.optimal_p),
                    self.optimal_weights,
                )
            ).transpose()
            cov_matrix_tests = np.exp(-1 * cmt)
            y_pred[i, 0] = self.optimal_mean + np.matmul(
                np.matmul(cov_matrix_tests.transpose(), self.covariance_matrix_inverse),
                self.optimal_y_mu,
            )
        return y_pred

    def training(self):
        """
        Main function for Kriging training.

        To train the Kriging model:
            (1) The Kriging exponent :math:`\\tau_{i}` is fixed at 2.
            (2) The optimal Kriging hyperparameters :math:`\\left(\mu,\sigma^{2},\\theta_{1},\ldots,\\theta_{n}\\right)` are evaluated by calling the ``optimal_parameter_evaluation`` method using either BFGS or Basinhopping.
            (3) The training predictions, prediction errors and r-square coefficient of fit are evaluated by calling the functions ``error_calculation`` and ``self.r2_calculation``
            (4) A results object is generated by calling the ``ResultsReport`` class.

        Returns:
            tuple   : self object (**results**) containing the all information about the best Kriging model obtained, including:
                - the Kriging model hyperparameters  (**results.optimal_weights**),
                - when relevant, the optimal regularization parameter found :math:`\lambda` (**results.regularization_parameter**),
                - the Kriging mean (**results.optimal_mean**),
                - the Kriging variance (**results.optimal_variance**),
                - the Kriging model regularized co-variance matrix (**results.optimal_covariance_matrix**),
                - the inverse of the co-variance matrix (**results.covariance_matrix_inverse**),
                - the RBF predictions for the training data (**results.output_predictions**),
                - the RMSE of the training output predictions (**results.training_rmse**), and
                - the :math:`R^{2}` value on the training data (**results.R2**)

        """

        # Create p values, for now fixed at p=2. Arraying p makes the code significantly (at least 7x slower)
        p = 2
        # Solve optimization problem
        bh_results = self.parameter_optimization(p)
        # Calculate other variables and parameters
        (
            optimal_theta,
            optimal_reg_param,
            optimal_mean,
            optimal_variance,
            optimal_cov_mat,
            opt_cov_inv,
            optimal_ymu,
        ) = self.optimal_parameter_evaluation(bh_results.x, p)
        # Training performance
        training_ss_error, rmse_error, y_training_predictions = self.error_calculation(
            optimal_theta,
            p,
            optimal_mean,
            opt_cov_inv,
            optimal_ymu,
            self.x_data_scaled,
            self.y_data,
        )
        r2_training = self.r2_calculation(self.y_data, y_training_predictions)

        # Return results
        self.optimal_weights = optimal_theta
        self.optimal_p = p
        self.optimal_mean = optimal_mean
        self.optimal_variance = optimal_variance
        self.regularization_parameter = optimal_reg_param
        self.optimal_covariance_matrix = optimal_cov_mat
        self.covariance_matrix_inverse = opt_cov_inv
        self.optimal_y_mu = optimal_ymu
        self.output_predictions = y_training_predictions
        self.training_R2 = r2_training
        self.training_rmse = rmse_error
        self.pickle_save({"model": self})
        return self

    def generate_expression(self, variable_list):
        """
        The ``generate_expression`` method returns the Pyomo expression for the Kriging model trained.

        The expression is constructed based on the supplied list of variables **variable_list** and the results of the previous Kriging training process.

        Args:
            variable_list(list)           : List of input variables to be used in generating expression. This can be the a list generated from the output of ``get_feature_vector``.  The user can also choose to supply a new list of the appropriate length.

        Returns:
            Pyomo Expression              : Pyomo expression of the Kriging model based on the variables provided in **variable_list**

        """
        t1 = np.array([variable_list])
        # Reshaping of variable array is necessary when input variables are Pyomo scalar variables
        t1 = t1.reshape(1, len(variable_list)) if t1.ndim > 2 else t1

        phi_var = []
        for i in range(0, self.x_data.shape[0]):
            curr_term = sum(
                self.optimal_weights[j]
                * (
                    (
                        (t1[0, j] - self.x_data_min[0, j])
                        / (self.x_data_max[0, j] - self.x_data_min[0, j])
                    )
                    - self.x_data_scaled[i, j]
                )
                ** self.optimal_p
                for j in range(0, self.x_data.shape[1])
            )

            curr_term = exp(-curr_term)
            phi_var.append(curr_term)
        phi_var_array = np.asarray(phi_var)

        phi_inv_times_y_mu = np.matmul(
            self.covariance_matrix_inverse, self.optimal_y_mu
        )
        phi_inv_times_y_mu = phi_inv_times_y_mu.reshape(
            phi_inv_times_y_mu.shape[0],
        )
        kriging_expr = self.optimal_mean[0, 0]
        kriging_expr += sum(
            w * t
            for w, t in zip(
                np.nditer(phi_inv_times_y_mu),
                np.nditer(phi_var_array, flags=["refs_ok"]),
            )
        )
        return kriging_expr

    def get_feature_vector(self):
        """

        The ``get_feature_vector`` method generates the list of regression features from the column headers of the input dataset.

        Returns:
            Pyomo IndexedParam  : An indexed parameter list of the variables supplied in the original data

        """
        p = Param(self.x_data_columns, mutable=True, initialize=0)
        p.index_set().construct()
        p.construct()
        self.feature_list = p
        return p

    def pickle_save(self, solutions):
        """
        The poly_training method saves the results of the run in a pickle object. It saves an object with two elements: the setup (index[0]) and the results (index[1]).
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
        pickle_load loads the results of a saved run 'file.obj'.])

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
        ax.plot(self.y_data, self.y_data, "-")
        ax.plot(self.y_data, self.output_predictions, "o")
        ax.set_xlabel(r"True data", fontsize=12)
        ax.set_ylabel(r"Surrogate values", fontsize=12)
        ax.set_title(r"Parity plot", fontsize=12)

        ax2 = fig1.add_subplot(122)
        ax2.plot(
            self.y_data,
            self.y_data - self.output_predictions,
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
            f"\nResults of Kriging run:\n"
            f"\nKriging mean                     : {self.optimal_mean}\n"
            f"Kriging variance                 : {self.optimal_variance}\n"
            f"Kriging weights                  : {self.optimal_weights}\n"
            f"Regularization parameter         : {self.regularization_parameter}\n"
            f"Number of terms in Kriging model : {self.optimal_y_mu.size + 1}\n"
            f"\nKriging Expression:\n"
            f"--------------------\n"
            f"\n{eqn}\n"
            f"--------------------------\n"
            f"\nModel training errors:"
            f"\n-----------------------\n"
            f"Mean Squared Error (MSE)         : {self.training_rmse ** 2}\n"
            f"Root Mean Squared Error (RMSE)   : {self.training_rmse}\n"
            f"Goodness of fit (R2)             : {self.training_R2}\n"
            f"\n{double_line}"
        )
        return s

    def print_report(self):
        s = self._report()
        print(s)

    def _repr_pretty_(self, p, cycle=False):
        s = self._report()
        p.text(s)
