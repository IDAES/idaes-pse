import numpy as np
from scipy.optimize import basinhopping
import scipy.optimize as opt
import pandas as pd
from pyomo.core import Param, exp
from sampling import FeatureScaling as fs


class ResultReport:
    """
    ===================================================================================================================
    A class for creating an object containing information about the RBF solution to be returned to the user.

        :returns:
        self function containing ten attributes -

            self.weights(<np.ndarray>)              : array containing optimal radial weights (coefficients) for the RBF
            self.sigma(<float>)                     : best shape parameter found for selected parametric basis function. Will return zero for fixed basis functions.
            self.regularization                     : Boolean variable indicating whether regularization was turned on/off for the problem
            self.regularization_parameter           : best regularization parameter found. Will return zero when reqularization is turned off
            self.centres                            : co-ordinates of RBF centres
            self.rmse                               : RMSE error on the training output predictions
            self.output_predictions                 : Predictions from the RBF surrogate for the output variable of the training data
            self.condition_number                   : Condition number of the regularized coefficient matrix used to generate the radial weights in the regression problem
            self.R2                                 : R2 coefficient-of-fit between the actual output and the surrogate predictions
            self.solution_status                    : Judgement on coefficient matrix conditioning. Returns 'ok' when problem is sufficiently well conditioned (condition number < 1 / eps); 'unstable solution' when the problem is ill-conditioned.

    ===================================================================================================================
    """

    def __init__(self, theta, reg_param, mean, variance, cov_mat, cov_inv, ymu, y_training, r2_training, rmse_error, p, x_data, x_data_scaled, x_data_min, x_data_max):
        self.optimal_weights = theta
        self.optimal_p = p
        self.optimal_mean = mean
        self.optimal_variance = variance
        self.regularization_parameter = reg_param
        self.optimal_covariance_matrix = cov_mat
        self.covariance_matrix_inverse = cov_inv
        self.optimal_y_mu = ymu
        self.output_predictions = y_training
        self.training_R2 = r2_training
        self.training_rmse = rmse_error
        self.x_data = x_data
        self.scaled_x = x_data_scaled
        self.x_min = x_data_min
        self.x_max = x_data_max

    def kriging_generate_expression(self, variable_list):
        t1 = np.array([variable_list])
        phi_var = []
        for i in range(0, self.x_data.shape[0]):
            curr_term = sum(self.optimal_weights[j] * ( ((t1[0, j] - self.x_min[0, j])/(self.x_max[0, j] - self.x_min[0, j])) - self.scaled_x[i, j]) ** self.optimal_p for j in range(0, self.x_data.shape[1]))

            curr_term = exp(-curr_term)
            phi_var.append(curr_term)
        phi_var_array = np.asarray(phi_var)

        phi_inv_times_y_mu = np.matmul(self.covariance_matrix_inverse, self.optimal_y_mu)
        phi_inv_times_y_mu = phi_inv_times_y_mu.reshape(phi_inv_times_y_mu.shape[0], )
        kriging_expr = self.optimal_mean[0,0]
        kriging_expr += sum(w * t for w, t in zip(np.nditer(phi_inv_times_y_mu), np.nditer(phi_var_array, flags=['refs_ok']) ))
        return kriging_expr


class MyBounds(object):
    def __init__(self, xmax=[1], xmin=[1e-6]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(x[-1] <= self.xmax)
        tmin = bool(x[-1] >= self.xmin)
        return tmax and tmin


class KrigingModel:

    def __init__(self, XY_data, numerical_gradients=True, regularization=True):

        # Check data types and shapes
        if isinstance(XY_data, pd.DataFrame):
            xy_data = XY_data.values
            self.x_data_columns = list(XY_data.columns)[:-1]
        elif isinstance(XY_data, np.ndarray):
            xy_data = XY_data
            self.x_data_columns = list(range(XY_data.shape[1] - 1))
        else:
            raise ValueError('Pandas dataframe or numpy array required for "XY_data".')

        self.x_data = xy_data[:, :-1].reshape(xy_data.shape[0], xy_data.shape[1]-1)
        self.y_data = xy_data[:, -1].reshape(xy_data.shape[0], 1)
        self.num_vars = self.x_data.shape[1] + 1  # thetas and reg parameter only
        x_data_scaled, self.x_data_min, self.x_data_max = fs.data_scaling_minmax(self.x_data)
        self.x_data_scaled = x_data_scaled.reshape(self.x_data.shape)

        if isinstance(regularization, bool):
            self.num_grads = numerical_gradients
        else:
            raise Exception('numerical_gradients must be boolean.')

        if isinstance(regularization, bool):
            self.regularization = regularization
        else:
            raise Exception('Choice of regularization must be boolean.')


    @staticmethod
    def covariance_matrix_generator(x, theta, reg_param, p):
        distance_matrix = np.zeros((x.shape[0], x.shape[0]))
        for i in range(0, x.shape[0]):
            distance_matrix[i, :] = (np.matmul(((np.abs(x[i, :] - x)) ** p), theta)).transpose()
        cov_matrix = np.exp(-1 * distance_matrix)
        cov_matrix = cov_matrix + reg_param * np.eye(cov_matrix.shape[0])  # Regularization parameter addition, see Forrester book
        return cov_matrix

    @staticmethod
    def covariance_inverse_generator(x):
        try:
            inverse_x = np.linalg.inv(x)
        except np.linalg.LinAlgError as LAE:
            inverse_x = np.linalg.pinv(x)
        return inverse_x

    @staticmethod
    def kriging_mean(cov_inv, y):
        ones_vec = np.ones((y.shape[0], 1))
        kriging_mean = np.matmul(np.matmul(ones_vec.transpose(), cov_inv), y) / np.matmul(np.matmul(ones_vec.transpose(), cov_inv), ones_vec)
        # kriging_mean = np.matmul(ones_vec.transpose(), np.matmul(cov_inv, y)) / np.matmul(ones_vec.transpose(), np.matmul(cov_inv, ones_vec))
        return kriging_mean

    @staticmethod
    def y_mu_calculation(y, mu):
        y_mu = y - mu * np.ones((y.shape[0], 1))
        return y_mu

    @staticmethod
    def kriging_sd(cov_inv, y_mu, ns):
        sigma_sq = np.matmul(np.matmul(y_mu.transpose(), cov_inv), y_mu) / ns
        # sigma_sq = np.matmul(y_mu.transpose(), np.matmul(cov_inv, y_mu)) / ns
        return sigma_sq

    @staticmethod
    def print_fun(x, f, accepted):
        print("at minimum %.4f accepted %d" % (f, int(accepted)))

    def objective_function(self, var_vector, x, y, p):
        theta = var_vector[:-1]
        reg_param = var_vector[-1]
        theta = 10 ** theta  # Assumes log(theta) provided
        ns = y.shape[0]
        cov_mat = self.covariance_matrix_generator(x, theta, reg_param, p)
        try:  # Check Cholesky factorization
            L = np.linalg.cholesky(cov_mat)
            lndetcov = 2 * np.sum(np.log(np.abs(np.diag(L))))  # Approximation to 2nd term from Forrester book, making use of the Ch. factorization
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
        eps = 1e-6
        grad_vec = np.zeros(len(var_vector), )
        for i in range(0, len(var_vector)):
            var_vector_plus = np.copy(var_vector)
            var_vector_plus[i] = var_vector[i] + eps
            var_vector_minus = np.copy(var_vector)
            var_vector_minus[i] = var_vector[i] - eps
            # of = self.objective_function(var_vector, x, y, p)
            of_plus = self.objective_function(var_vector_plus, x, y, p)
            of_minus = self.objective_function(var_vector_minus, x, y, p)
            grad_current = (of_plus - of_minus) / (2 * eps)
            grad_vec[i, ] = grad_current
        if self.regularization is False:
            grad_vec[-1, ] = 0
        return grad_vec

    def parameter_optimization(self, p):
        """
        Parameter (theta) optimization using BFGS algorithm
        """
        initial_value_list = np.random.randn(self.num_vars - 1, )
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
            print('Optimizing kriging parameters using L-BFGS-B algorithm...')
            other_args = (self.x_data_scaled, self.y_data, p)
            opt_results = opt.minimize(self.objective_function, initial_value, args=other_args, method='L-BFGS-B', jac=self.numerical_gradient, bounds=bounds, options={'gtol': 1e-7}) #, 'disp': True})
        else:
            print('Optimizing Kriging parameters using Basinhopping algorithm...')
            other_args = {"args": (self.x_data_scaled, self.y_data, p), 'bounds': bounds}
            # other_args = {"args": (self.x_data, self.y_data, p)}
            mybounds = MyBounds()  # Bounds on regularization parameter
            opt_results = basinhopping(self.objective_function, initial_value_list, minimizer_kwargs=other_args, niter=1000, disp=True, accept_test=mybounds) # , interval=5)
        return opt_results

    def optimal_parameter_evaluation(self, var_vector, p):
        theta = var_vector[:-1]
        reg_param = var_vector[-1]
        theta = 10 ** theta  # Assumes log(theta) provided. Ensures that theta is always positive
        ns = self.y_data.shape[0]
        cov_mat = self.covariance_matrix_generator(self.x_data_scaled, theta, reg_param, p)
        cov_inv = self.covariance_inverse_generator(cov_mat)
        mean = self.kriging_mean(cov_inv, self.y_data)
        y_mu = self.y_mu_calculation(self.y_data, mean)
        variance = self.kriging_sd(cov_inv, y_mu, ns)
        print('\nFinal results\n================\nTheta:', theta, '\nMean:', mean, '\nRegularization parameter:', reg_param)
        return theta, reg_param, mean, variance, cov_mat, cov_inv, y_mu

    @staticmethod
    def error_calculation(theta, p, mean, cov_inv, y_mu, x, y_data):
        """
        """
        y_prediction = np.zeros((x.shape[0], 1))
        for i in range(0, x.shape[0]):
            cmt = (np.matmul(((np.abs(x[i, :] - x)) ** p), theta)).transpose()
            cov_matrix_tests = np.exp(-1 * cmt)
            y_prediction[i, 0] = mean + np.matmul(np.matmul(cov_matrix_tests.transpose(), cov_inv), y_mu)
        ss_error = (1 / y_data.shape[0]) * (np.sum((y_data - y_prediction) ** 2))
        rmse_error = np.sqrt(ss_error)
        return ss_error, rmse_error, y_prediction

    @staticmethod
    def r2_calculation(y_true, y_predicted):
        """

        ===============================================================================================================
        The function r2_calculation returns the R2-coefficient as a measure-of-fit between the true and predicted values of the output parameter.
              R2 = 1 - (ss_residual / ss_total)

        Input arguments:
            y_true(<np.ndarray>)             : Vector of actual values of the output variable
            x_predicted(<np.ndarray>)        : Vector of predictions for the output variable based on the surrogate

        :returns
            r_square(<float>)                : R2 measure-of-fit between actual anf predcited data
         ==============================================================================================================

        """
        y_true = y_true.reshape(y_true.shape[0], 1)
        y_predicted = y_predicted.reshape(y_predicted.shape[0], 1)
        input_y_mean = np.mean(y_true, axis=0)
        ss_total = np.sum((y_true - input_y_mean) ** 2)
        ss_residual = np.sum((y_predicted - y_true) ** 2)
        r_square = 1 - (ss_residual / ss_total)
        return r_square

    def kriging_predict_output(self, kriging_params, x_pred):
        x_pred_scaled = ((x_pred - self.x_data_min) / (self.x_data_max - self.x_data_min))
        x_pred = x_pred_scaled.reshape(x_pred.shape)
        if x_pred.ndim == 1:
            x_pred = x_pred.reshape(1, len(x_pred))
        y_pred = np.zeros((x_pred.shape[0], 1))
        for i in range(0, x_pred.shape[0]):
            cmt = (np.matmul(((np.abs(x_pred[i, :] - self.x_data_scaled)) ** kriging_params.optimal_p), kriging_params.optimal_weights)).transpose()
            cov_matrix_tests = np.exp(-1 * cmt)
            y_pred[i, 0] = kriging_params.optimal_mean + np.matmul(np.matmul(cov_matrix_tests.transpose(), kriging_params.covariance_matrix_inverse), kriging_params.optimal_y_mu)
        return y_pred

    def kriging_training(self):
        # Create p values, for now fixed at p=2. Arraying p makes the code significantly (at least 7x slower)
        p = 2
        # Solve optimization problem
        bh_results = self.parameter_optimization(p)
        # Calculate other variables and parameters
        optimal_theta, optimal_reg_param, optimal_mean, optimal_variance, optimal_cov_mat, opt_cov_inv, optimal_ymu = self.optimal_parameter_evaluation(bh_results.x, p)
        # Training performance
        training_ss_error, rmse_error, y_training_predictions = self.error_calculation(optimal_theta, p, optimal_mean, opt_cov_inv, optimal_ymu, self.x_data_scaled, self.y_data)
        r2_training = self.r2_calculation(self.y_data, y_training_predictions)
        # Return results
        results = ResultReport(optimal_theta, optimal_reg_param, optimal_mean, optimal_variance, optimal_cov_mat, opt_cov_inv, optimal_ymu, y_training_predictions, r2_training, rmse_error, p, self.x_data, self.x_data_scaled, self.x_data_min, self.x_data_max)
        return results

    def get_feature_vector(self):
        p = Param(self.x_data_columns, mutable=True, initialize=0)
        p.index_set().construct()
        p.construct()
        self.feature_list = p
        return p
