#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
import sys
import os
import io
from unittest.mock import patch

sys.path.append(os.path.abspath(".."))  # current folder is ~/tests
from idaes.core.surrogate.pysmo import kriging
from idaes.core.surrogate.pysmo.kriging import KrigingModel
from idaes import logger as idaes_logger
import numpy as np
import pandas as pd
import pytest

# Turn down the logging during the test
kriging.set_log_level(idaes_logger.ERROR)


class TestKrigingModel:
    y = np.array(
        [
            [i, j, ((i + 1) ** 2) + ((j + 1) ** 2)]
            for i in np.linspace(0, 10, 21)
            for j in np.linspace(0, 10, 21)
        ]
    )
    full_data = {"x1": y[:, 0], "x2": y[:, 1], "y": y[:, 2]}
    training_data = [
        [i, j, ((i + 1) ** 2) + ((j + 1) ** 2)]
        for i in np.linspace(0, 10, 5)
        for j in np.linspace(0, 10, 5)
    ]
    test_data = [[i, (i + 1) ** 2] for i in range(10)]
    test_data_large = [[i, (i + 1) ** 2] for i in range(200)]
    test_data_1d = [[(i + 1) ** 2] for i in range(10)]
    test_data_3d = [[i, (i + 1) ** 2, (i + 2) ** 2] for i in range(10)]
    sample_points = [[i, (i + 1) ** 2] for i in range(8)]
    sample_points_large = [[i, (i + 1) ** 2] for i in range(100)]
    sample_points_1d = [[(i + 1) ** 2] for i in range(8)]
    sample_points_3d = [[i, (i + 1) ** 2, (i + 2) ** 2] for i in range(8)]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__01(self, array_type):
        input_array = array_type(self.test_data)
        kriging_class = KrigingModel(input_array)
        assert kriging_class.num_grads is True
        assert kriging_class.regularization is True

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__02(self, array_type):
        input_array = array_type(self.test_data)
        with pytest.raises(ValueError):
            KrigingModel(input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__03(self, array_type):
        input_array = array_type(self.test_data)
        with pytest.raises(Exception):
            KrigingModel(input_array, numerical_gradients=1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__04(self, array_type):
        input_array = array_type(self.test_data)
        with pytest.raises(Exception):
            KrigingModel(input_array, regularization=1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__05(self, array_type):
        input_array = array_type(self.test_data)
        with pytest.raises(Exception):
            KrigingModel(input_array, overwrite=1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__06(self, array_type):
        input_array = array_type(self.test_data)
        with pytest.raises(Exception):
            KrigingModel(input_array, fname="solution.pkl")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__07(self, array_type):
        input_array = array_type(self.test_data)
        with pytest.raises(Exception):
            KrigingModel(input_array, fname=1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__08(self, array_type, tmp_path):
        input_array = array_type(self.test_data)
        file_name = str(tmp_path / "test_filename.pickle")
        kriging_class1 = KrigingModel(input_array, fname=file_name, overwrite=True)
        kriging_class1.training()
        kriging_class2 = KrigingModel(input_array, fname=file_name, overwrite=True)
        assert kriging_class1.filename == kriging_class2.filename

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__09(self, array_type, tmp_path):
        input_array = array_type(self.test_data)
        file_name1 = str(tmp_path / "test_filename1.pickle")
        file_name2 = str(tmp_path / "test_filename2.pickle")
        kriging_class1 = KrigingModel(input_array, fname=file_name1, overwrite=True)
        kriging_class1.training()
        kriging_class2 = KrigingModel(input_array, fname=file_name2, overwrite=True)
        assert kriging_class1.filename == file_name1
        assert kriging_class2.filename == file_name2

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_covariance_matrix_generator(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        p = 2
        theta = np.array([1, 2])
        reg_param = 1.00000000e-06
        cov_matrix = kriging_class.covariance_matrix_generator(
            kriging_class.x_data_scaled, theta, reg_param, p
        )

        cov_matrix_exp = np.array(
            [
                [1.000001, 0.60653066, 0.13533528],
                [0.60653066, 1.000001, 0.60653066],
                [0.13533528, 0.60653066, 1.000001],
            ]
        )
        np.testing.assert_array_equal(
            np.round(cov_matrix, 7), np.round(cov_matrix_exp, 7)
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_covariance_inverse_generator_01(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        cov_matrix = np.array(
            [
                [1.000001, 0.60653066, 0.13533528],
                [0.60653066, 1.000001, 0.60653066],
                [0.13533528, 0.60653066, 1.000001],
            ]
        )
        cov_matrix_inv_exp = np.array(
            [
                [1.82957788, -1.51792604, 0.67306158],
                [-1.51792604, 2.84133453, -1.51792604],
                [0.67306158, -1.51792604, 1.82957788],
            ]
        )

        inverse_x = kriging_class.covariance_inverse_generator(cov_matrix)
        np.testing.assert_array_equal(
            np.round(inverse_x, 7), np.round(cov_matrix_inv_exp, 7)
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_covariance_inverse_generator_02(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        cov_matrix = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        inverse_x = kriging_class.covariance_inverse_generator(cov_matrix)
        np.testing.assert_array_equal(np.round(inverse_x, 7), np.round(cov_matrix, 7))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_kriging_mean(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        cov_matrix_inv = np.array(
            [
                [1.82957788, -1.51792604, 0.67306158],
                [-1.51792604, 2.84133453, -1.51792604],
                [0.67306158, -1.51792604, 1.82957788],
            ]
        )
        kriging_mean = kriging_class.kriging_mean(cov_matrix_inv, kriging_class.y_data)
        kriging_mean_exp = 20.18496
        assert np.round(kriging_mean_exp, 5) == np.round(kriging_mean[0][0], 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_y_mu_calculation(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        kriging_mean = 20.18496
        y_mu = kriging_class.y_mu_calculation(kriging_class.y_data, kriging_mean)
        y_mu_exp = np.array([[-18.18496], [-6.93496], [16.81504]])
        np.testing.assert_array_equal(np.round(y_mu, 5), np.round(y_mu_exp, 5))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_kriging_sd(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        cov_matrix_inv = np.array(
            [
                [1.82957788, -1.51792604, 0.67306158],
                [-1.51792604, 2.84133453, -1.51792604],
                [0.67306158, -1.51792604, 1.82957788],
            ]
        )
        y_mu_exp = np.array([[-18.18496], [-6.93496], [16.81504]])
        sigma_sq = kriging_class.kriging_sd(
            cov_matrix_inv, y_mu_exp, kriging_class.y_data.shape[0]
        )
        sigma_sq_exp = 272.84104637
        assert np.round(sigma_sq_exp, 5) == np.round(sigma_sq[0][0], 5)

    # @pytest.mark.unit
    # @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    # def test_print_fun(self, array_type):
    #     input_array = array_type(self.training_data)
    #     kriging_class = KrigingModel(input_array[0:3], regularization=True)
    #     capturedOutput = io.StringIO()
    #     sys.stdout = capturedOutput
    #     kriging_class.print_fun(1, 2, 3.7)
    #     sys.stdout = sys.__stdout__
    #     assert "at minimum 2.0000 accepted 3\n" == capturedOutput.getvalue()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_objective_function(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        p = 2
        var_vector = np.array([1, 2, 1.00000000e-06])
        conc_log_like = kriging_class.objective_function(
            var_vector, kriging_class.x_data_scaled, kriging_class.y_data, p
        )

        conc_log_like_exp = 8.0408619
        assert np.round(conc_log_like_exp, 5) == np.round(conc_log_like, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_numerical_gradient_01(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=True)
        p = 2
        var_vector = np.array([1, 2, 1.00000000e-06])
        grad_vec = kriging_class.numerical_gradient(
            var_vector, kriging_class.x_data_scaled, kriging_class.y_data, p
        )
        grad_vec_exp = np.array([0, 0, 8.8817842e-10])
        np.testing.assert_array_equal(np.round(grad_vec, 5), np.round(grad_vec_exp, 5))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_numerical_gradient_02(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=False)
        p = 2
        var_vector = np.array([1, 2, 1.00000000e-06])
        grad_vec = kriging_class.numerical_gradient(
            var_vector, kriging_class.x_data_scaled, kriging_class.y_data, p
        )
        grad_vec_exp = np.array([0, 0, 0])
        np.testing.assert_array_equal(np.round(grad_vec, 5), np.round(grad_vec_exp, 5))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_parameter_optimization_01(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3])
        p = 2
        np.random.seed(0)
        opt_results = kriging_class.parameter_optimization(p)
        assert len(opt_results.x) == 3
        assert opt_results.success is True

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_parameter_optimization_02(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], numerical_gradients=False)
        p = 2
        opt_results = kriging_class.parameter_optimization(p)
        assert len(opt_results.x) == 3
        assert not opt_results.minimization_failures

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_optimal_parameter_evaluation(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3])
        p = 2
        var_vector = np.array([1, 2, 1.00000000e-06])
        (
            theta,
            reg_param,
            mean,
            variance,
            cov_mat,
            cov_inv,
            y_mu,
        ) = kriging_class.optimal_parameter_evaluation(var_vector, p)
        np.testing.assert_array_equal(theta, [10**1, 10**2])
        np.testing.assert_array_equal(reg_param, 1.00000000e-06)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_error_calculation(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=False)
        p = 2
        var_vector = np.array([1, 2, 1.00000000e-06])
        (
            theta,
            reg_param,
            mean,
            variance,
            cov_mat,
            cov_inv,
            y_mu,
        ) = kriging_class.optimal_parameter_evaluation(var_vector, p)
        y_prediction_exp = np.zeros((kriging_class.x_data_scaled.shape[0], 1))
        for i in range(0, kriging_class.x_data_scaled.shape[0]):
            cmt = (
                np.matmul(
                    (
                        (
                            np.abs(
                                kriging_class.x_data_scaled[i, :]
                                - kriging_class.x_data_scaled
                            )
                        )
                        ** p
                    ),
                    theta,
                )
            ).transpose()
            cov_matrix_tests = np.exp(-1 * cmt)
            y_prediction_exp[i, 0] = mean + np.matmul(
                np.matmul(cov_matrix_tests.transpose(), cov_inv), y_mu
            )

        ss_error, rmse_error, y_prediction = kriging_class.error_calculation(
            theta,
            p,
            mean,
            cov_inv,
            y_mu,
            kriging_class.x_data_scaled,
            kriging_class.y_data,
        )

        np.testing.assert_array_equal(y_prediction, y_prediction_exp)
        assert (
            np.sum((kriging_class.y_data - y_prediction_exp) ** 2)
            / kriging_class.x_data_scaled.shape[0]
            == ss_error
        )
        assert (
            np.sqrt(
                np.sum((kriging_class.y_data - y_prediction_exp) ** 2)
                / kriging_class.x_data_scaled.shape[0]
            )
            == rmse_error
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_r2_calculation(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=False)
        p = 2
        var_vector = np.array([1, 2, 1.00000000e-06])
        (
            theta,
            reg_param,
            mean,
            variance,
            cov_mat,
            cov_inv,
            y_mu,
        ) = kriging_class.optimal_parameter_evaluation(var_vector, p)

        ss_error, rmse_error, y_prediction = kriging_class.error_calculation(
            theta,
            p,
            mean,
            cov_inv,
            y_mu,
            kriging_class.x_data_scaled,
            kriging_class.y_data,
        )
        r_square = kriging_class.r2_calculation(kriging_class.y_data, y_prediction)
        assert 0.999999999999 == r_square

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_predict_output_01(self, array_type):
        input_array = array_type(self.training_data)
        np.random.seed(0)
        kriging_class = KrigingModel(input_array)
        kriging_class.training()
        y_pred = kriging_class.predict_output(kriging_class.x_data_scaled)
        assert y_pred.shape[0] == kriging_class.x_data_scaled.shape[0]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_predict_output(self, array_type):
        input_array = array_type(self.training_data)
        np.random.seed(0)
        kriging_class = KrigingModel(input_array)
        kriging_class.training()
        y_pred = kriging_class.predict_output(np.array([0.1, 0.2]))
        assert y_pred.shape[0] == 1

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_training(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array[0:3], regularization=False)
        np.random.seed(0)
        p = 2
        bh_results = kriging_class.parameter_optimization(p)
        # Calculate other variables and parameters
        (
            optimal_theta,
            optimal_reg_param,
            optimal_mean,
            optimal_variance,
            optimal_cov_mat,
            opt_cov_inv,
            optimal_ymu,
        ) = kriging_class.optimal_parameter_evaluation(bh_results.x, p)
        # Training performance
        (
            training_ss_error,
            rmse_error,
            y_training_predictions,
        ) = kriging_class.error_calculation(
            optimal_theta,
            p,
            optimal_mean,
            opt_cov_inv,
            optimal_ymu,
            kriging_class.x_data_scaled,
            kriging_class.y_data,
        )
        r2_training = kriging_class.r2_calculation(
            kriging_class.y_data, y_training_predictions
        )

        np.random.seed(0)
        results = kriging_class.training()
        np.testing.assert_array_equal(results.optimal_weights, optimal_theta)
        np.testing.assert_array_equal(
            results.regularization_parameter, optimal_reg_param
        )
        np.testing.assert_array_equal(results.optimal_mean, optimal_mean)
        np.testing.assert_array_equal(results.optimal_variance, optimal_variance)
        np.testing.assert_array_equal(
            results.optimal_covariance_matrix, optimal_cov_mat
        )
        np.testing.assert_array_equal(results.optimal_y_mu, optimal_ymu)
        np.testing.assert_array_equal(
            results.output_predictions, y_training_predictions
        )
        np.testing.assert_array_equal(results.training_R2, r2_training)
        np.testing.assert_array_equal(results.training_rmse, rmse_error)
        np.testing.assert_array_equal(results.optimal_p, p)
        np.testing.assert_array_equal(results.x_data, kriging_class.x_data)
        np.testing.assert_array_equal(
            results.x_data_scaled, kriging_class.x_data_scaled
        )
        np.testing.assert_array_equal(results.x_data_min, kriging_class.x_data_min)
        np.testing.assert_array_equal(results.x_data_max, kriging_class.x_data_max)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_get_feature_vector_01(self, array_type):
        input_array = array_type(self.full_data)
        kriging_class = KrigingModel(input_array, regularization=False)
        p = kriging_class.get_feature_vector()
        expected_dict = {"x1": 0, "x2": 0}
        assert expected_dict == p.extract_values()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_get_feature_vector_02(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array, regularization=False)
        p = kriging_class.get_feature_vector()
        expected_dict = {0: 0, 1: 0}
        assert expected_dict == p.extract_values()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_kriging_generate_expression(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array, regularization=False)
        results = kriging_class.training()
        p = kriging_class.get_feature_vector()
        lv = []
        for i in p.keys():
            lv.append(p[i])
        results.generate_expression(lv)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_pickle_load01(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array, regularization=False)
        kriging_class.training()
        kriging_class.pickle_load(kriging_class.filename)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_pickle_load02(self, array_type):
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array, regularization=False)
        with pytest.raises(Exception):
            kriging_class.pickle_load("file_not_existing.pickle")

    @pytest.mark.unit
    @patch("matplotlib.pyplot.show")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_parity_residual_plots(self, mock_show, array_type):  # noqa
        input_array = array_type(self.training_data)
        kriging_class = KrigingModel(input_array, regularization=False)
        kriging_class.training()
        kriging_class.parity_residual_plots()


if __name__ == "__main__":
    pytest.main()
