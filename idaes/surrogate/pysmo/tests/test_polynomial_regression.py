##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
import sys, os
sys.path.append(os.path.abspath('..'))# current folder is ~/tests

from idaes.surrogate.pysmo.polynomial_regression import PolynomialRegression, FeatureScaling
import numpy as np
import pandas as pd
import pyutilib.th as unittest
from unittest.mock import patch
import pytest

'''
coverage run test_polynomial_regression.py
coverage report -m
coverage html

'''

class FeatureScalingTestCases(unittest.TestCase):
    """
    test_data_scaling_01: Test behaviour when input is a numpy array and with 1D array.
    test_data_scaling_02: Test behaviour when input is a numpy array and with 2D array.
    test_data_scaling_03: Test behaviour when input is a numpy array and with 3D array.
    test_data_scaling_04: Test behaviour when input is a numpy array and with 3D array with a varibale is constant.
    test_data_scaling_05: Test behaviour list input TypeError

    test_data_scaling_06: Test behaviour when input is a Pandas DataFrame and with 1D array.
    test_data_scaling_07: Test behaviour when input is a Pandas DataFrame and with 2D array.
    test_data_scaling_08: Test behaviour when input is a Pandas DataFrame and with 3D array.
    test_data_scaling_09: Test behaviour when input is a Pandas DataFrame and with 3D array with a varibale is constant.

    test_data_unscaling_01: Test behaviour when input is a numpy array and with 1D array.
    test_data_unscaling_02: Test behaviour when input is a numpy array and with 2D array.
    test_data_unscaling_03: Test behaviour when input is a numpy array and with 3D array. 
    test_data_unscaling_04: Test behaviour when input is a numpy array and with 3D array with a varibale is constant.

    test_data_unscaling_05: Test behaviour IndexError when input array size > array size
    test_data_unscaling_06: Test behaviour IndexError when input array size < array size  

    """
    def setUp(self):
        input_array_np_1d =  np.array(range(10))
        input_array_np_2d = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        input_array_np_3d = np.array([[0,10, 11], [1,11, 15], [2,12, 21], [3,13, 29], [4,14, 39], [5,15, 51], [6,16, 65], [7,17, 81], [8,18, 99], [9,19, 119]])
        input_array_np_3d_constant = np.array([[0, 10, 11], [1, 10, 14], [2, 10, 19], [3, 10, 26], [4, 10, 35], [5, 10, 46], [6, 10, 59], [7, 10, 74], [8, 10, 91], [9, 10, 110]])

        input_array_pd_1d = pd.DataFrame( np.array(range(10)))
        input_array_pd_2d = pd.DataFrame([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        input_array_pd_3d = pd.DataFrame([[0,10, 11], [1,11, 15], [2,12, 21], [3,13, 29], [4,14, 39], [5,15, 51], [6,16, 65], [7,17, 81], [8,18, 99], [9,19, 119]])
        input_array_pd_3d_constant = pd.DataFrame([[0, 10, 11], [1, 10, 14], [2, 10, 19], [3, 10, 26], [4, 10, 35], [5, 10, 46], [6, 10, 59], [7, 10, 74], [8, 10, 91], [9, 10, 110]])

        self.test_data_numpy_1d = input_array_np_1d
        self.test_data_numpy_2d = input_array_np_2d
        self.test_data_numpy_3d = input_array_np_3d
        self.test_data_numpy_3d_constant = input_array_np_3d_constant

        self.test_data_pandas_1d = input_array_pd_1d
        self.test_data_pandas_2d = input_array_pd_2d
        self.test_data_pandas_3d = input_array_pd_3d
        self.test_data_pandas_3d_constant = input_array_pd_3d_constant
        
    
    def test_data_scaling_01(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_numpy_1d

        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1.reshape(10, 1))

    def test_data_scaling_02(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9, 100]])
        expected_output_2 = np.array([[0, 1]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_03(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_numpy_3d 
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9, 19, 119]])
        expected_output_2 = np.array([[0, 10, 11]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_04(self):
        # 3D with constant
        # Sample data generated from the expression (x_1 + 1)^2 + 10, x1: between 0 and 9
        input_array = self.test_data_numpy_3d_constant
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9, 10, 110]])
        expected_output_2 = np.array([[0, 10, 11]])
        scale = expected_output_3 - expected_output_2
        scale[scale == 0.0] = 1.0
        expected_output_1 = (input_array - expected_output_2) / scale
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_05(self):
        # TypeError with list
        input_array = self.test_data_numpy_2d.tolist()
        with pytest.raises(TypeError):
            FeatureScaling.data_scaling(input_array)

    def test_data_scaling_06(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_pandas_1d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    def test_data_scaling_07(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array =self.test_data_pandas_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9, 100]])
        expected_output_2 = np.array([[0, 1]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_08(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_pandas_3d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9, 19, 119]])
        expected_output_2 = np.array([[0, 10, 11]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_09(self):
        # 3D with constant
        # Sample data generated from the expression (x_1 + 1)^2 + 10, x1: between 0 and 9
        input_array = self.test_data_pandas_3d_constant
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        expected_output_3 = np.array([[9, 10, 110]])
        expected_output_2 = np.array([[0, 10, 11]])
        scale = expected_output_3 - expected_output_2
        scale[scale == 0.0] = 1.0
        expected_output_1 = (input_array - expected_output_2) / scale
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    

    def test_data_unscaling_01(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_numpy_1d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        output_1 = output_1.reshape(output_1.shape[0], )
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array.reshape(10, 1))

    def test_data_unscaling_02(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)
    
    def test_data_unscaling_03(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_numpy_3d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)
    
    def test_data_unscaling_04(self):
        # 3D with constant
        # Sample data generated from the expression (x_1 + 1)^2 + 10, x1: between 0 and 9
        input_array = self.test_data_numpy_3d_constant
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    def test_data_unscaling_05(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)

        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling(output_1, min_array, max_array)

    def test_data_unscaling_06(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
                
        min_array = np.array([[1,2,3]])
        max_array = np.array([[5,6,7]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling(output_1, min_array, max_array)

class PolynomialRegressionTestCases(unittest.TestCase):
    '''
    test__init__01: Check that all default values are correctly loaded when input is sparse with both input arrays are numpy
    test__init__02: Check that all default values are correctly loaded when input is sparse with both input arrays are pandas
    test__init__03: Check that all default values are correctly loaded when input is sparse with input array is pandas or numpy
    test__init__04: Check that all default values are correctly loaded when input is sparse with input array is pandas or numpy
    test__init__05: Check that all manual values are correctly loaded when input is sparse with both input arrays are numpy
    test__init__06: Test behaviour raise ValueError with original_data_input = list
    test__init__07: Test behaviour raise ValueError with regression_data_input = list
    test__init__08: Test behaviour raise Exception with regression_data_input = list
    test__init__09: Test behaviour raise Exception with the sampled data has more entries than the original dataset.
    test__init__10: Test behaviour raise Exception with Dimensional discrepancies in the dimensions of the original > regression datasets
    test__init__11: Test behaviour raise Exception with Dimensional discrepancies in the dimensions of the original < regression datasets
    test__init__12: Test behaviour raise Exception with Input data requires at least two dimensions (X and Y data)
    test__init__13: Test behaviour raise Warning with number_of_crossvalidations > 10 and correctly set as 11
    test__init__14: Test behaviour raise Exception with non-integer maximum_polynomial_order
    test__init__13: Test behaviour raise Warning with maximum_polynomial_order > 10 and set t0 10
    test__init__16: Test behaviour raise Exception with training_split >=1 
    test__init__17: Test behaviour raise Exception with training_split<=0
    test__init__18: Test behaviour raise Exception with max_fraction_training_samples>1
    test__init__19: Test behaviour raise Exception with max_fraction_training_samples<0
    test__init__20: Check if regression_data.shape[0] == original_data.shape[0], set max_iter =0
    test__init__21: Check if no_adaptive_samples == 0, set max_iter =0
    test__init__22: Test behaviour raise Exception with number_of_crossvalidations is not integer
    test__init__23: Test behaviour raise Exception with no_adaptive_samples is not integer    
    test__init__24: Test behaviour raise Exception with max_iter is not integer
    test__init__25: Test behaviour raise Exception with max_polynomial_order >= num samples    
    test__init__26: Test behaviour raise Exception with invalid solution_method type
    test__init__27: Test behaviour raise Exception with invalid str solution_method 
    test__init__28: Test behaviour raise Exception with multinomials is not 0 or 1
    test__init__29: Test behaviour raise Exception with maximum_polynomial_order< 0
    test__init__30: Test behaviour raise Exception with number_of_crossvalidations< 0
    test__init__31: Test behaviour raise Exception with no_adaptive_samples< 0
    test__init__32: Test behaviour raise Exception with max_iter< 0

    test_training_test_data_creation_01: Test behaviour raise Exception with num_training = 0,training_split=0.01
    test_training_test_data_creation_02: Test behaviour raise Exception with num_training == self.number_of_samples, training_split=0.99
    test_training_test_data_creation_03: Check 1. splited training / test size = cross validation size, 2. size of each train / test, 3. each train / test are correctly splitted with default class values
    test_training_test_data_creation_04: Check 1. splited training / test size = cross validation size, 2. size of each train / test, 3. each train / test are correctly splitted with manual class values
    test_training_test_data_creation_05: Check with additional data, 1. splited training / test size = 2*cross validation size, 2. size of each train, train_extra / test, test_extra, 3. ach train, train_extra / test, test_extra are correctly splitted with default class values

    test_polygeneration_01:
    test_polygeneration_02:
    test_polygeneration_03:
    test_polygeneration_04:
    test_polygeneration_05:
    Tests: : Unit tests for polygeneration function. Five tests covering the range of max_polynomial_order are considered.
        The first test checks matrix/array returned when polynomial order is 1 -  the minimum possible value.
        The second test checks matrix/array returned when polynomial order is 2.
        The third test checks matrix/array returned when polynomial order is 10 (maximum possible value) with multinomials = default (1).
        The fourth test checks matrix/array returned when polynomial order is 10 (maximum possible value) with multinomials = 0.
        The fifth test checks matrix/array returned when polynomial order is 1 and additional terms have been supplied.

    test_cost_function_01:
    test_cost_function_02:
    test_cost_function_03:
        Tests: : Unit tests for cost_function.  The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
            Three demonstration tests are done:
            The first test checks that the correct loss value is returned when theta is initialized as a vector of zeros.
            The second test checks that the correct loss value is returned for a random point within the solution space.
            The third test checks that the loss value is zero when the exact solution is found.

    test_gradient_function_01:
    test_gradient_function_02:
    test_gradient_function_03:
    Tests: : Unit tests for gradient_function. The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
        Three demonstration tests are done:
        The first test checks that the correct gradients returned when theta is initialized as a vector of zeros.
        The second test checks that the correct gradients are returned for a random point within the solution space.
        The third test checks that the gradients are zero when the exact solution is found. 

    test_bfgs_parameter_optimization_01:
    test_bfgs_parameter_optimization_02:
        Tests: : Unit tests for bfgs_parameter_optimization function. Two unit tests are done:
            The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
            The second test is used to check that higher order multi-variable problems are solved correctly.
            - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
            - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.
    
    test_mle_estimate_01:
    test_mle_estimate_02:
    Tests: : Unit tests for MLE_estimate function. The same two unit tests are done:
        The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
        The second test is used to check that higher order multi-variable problems are solved correctly.
        - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
        - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.

    test_pyomo_optimization_01: 
    test_pyomo_optimization_02: 
      Unit tests for pyomo_optimization function. As with the other parameter estimation methods, two unit tests are done:
        The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
        The second test is used to check that higher order multi-variable problems are solved correctly.
        - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
        - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.

    test_cross_validation_error_calculation_01:
    test_cross_validation_error_calculation_02:
    test_cross_validation_error_calculation_03:
      Tests: : Unit tests for cross_validation_error_calculation function.  The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
        The three tests used for the cost_function are repeated for the cross-validation error. 
        The first test checks that the correct cv error is returned when theta is initialized as a vector of zeros.
        The second test checks that the correct cv error is returned for a random point within the solution space.
        The third test checks that the cv error is zero when the exact solution is found.
          In each case, the result is expected to be twice as large as the cost function's expected value: the cost function is half of the total squared error.

    test_polyregression_01:
    test_polyregression_02:
    test_polyregression_03:
    test_polyregression_04:
       Tests: : Unit tests for the polyregression function. We verify that:
           - 1: The correct optimization method is called and the correct optimization solution returned when 'mle' is selected as the solution_method.
           - 2: The correct optimization method is called and the correct optimization solution returned when 'bfgs' is selected as the solution_method.
           - 3: The correct optimization method is called and the correct optimization solution returned when 'pyomo' is selected as the solution_method.
           - 3: No optimization problem is solved when the problem is underspecified - all the outputs default to np.Inf.
           For the tests, the relevant optimization method in each case was replaced by a mock function.

    test_surrogate_performance_01:
    test_surrogate_performance_02:
    test_surrogate_performance_03:
    Tests: : Unit tests for surrogate_performance function.  The single variable problem (x+1) ^ 2 is considered. 
        Two tests considered for the function to validate its MAE, MSE, R2 and adjusted R2 calculations.
        - 01: The errors reported when the initialization values of theta (zeros) are supplied are compared to manually calculated values: see excel sheet.
        - 02: The errors reported when the actual values of theta (zeros) are supplied are checked. The errors are expected to be zero, while R^2 ~= 1.
        - 02: The errors reported when a random value of theta (zeros) is supplied are checked.
    The returned array contains [x,y_real,y_prediction] data and does not need to be checked once the errors are correct.

    test_results_generation_01:
    test_results_generation_02:
    test_results_generation_03:
       Tests: : Unit tests for results_generation. For specified solution arrays and orders, we verify that:
           - 1: The index and correct values are returned for a polynomial order of 1 and no multinomials
           - 2: The index and correct values are returned for a polynomial order of 3 and no multinomials
           - 3: The index and correct values are returned for a polynomial order of 2 and multinomials present

    test_error_plotting:
       Tests: : Unit tests for error_plotting function. For a sample dataset of the expected shape (n x 8), we verify that:
           - 1: The data we supply is the data used in the plot for each of the subplots, and that the right columns are accessed in each case.
           - 2: The right plot titles are displayed for each figure.
           - 3. The vertices of the path (and thus actual line segments) for each plot is the same as the input data; more information on matplotlib's path may be found at https://matplotlib.org/api/path_api.html
    
    test_user_defined_terms_01:
    test_user_defined_terms_02:
    test_user_defined_terms_03:
    test_user_defined_terms_04:
    test_user_defined_terms_05:
       Tests: : Unit tests for user_defined_terms function. For a sample dataset, we verify that:
           - 1&2: The function behaves properly when the right data types and shapes are entered.
           - 3: An exception is raised when the the number of samples in the term_list array is different from the regression_data length.
           - 4: An exception is raised when list element is not 1-dimensional
           - 5. An exception is raised when the term list (additional_terms) is not a list

    test_polynomial_regression_fitting_01: checking the status is 'ok', R2 > 0.95, __init__ in ResultReport class is covered here
    test_polynomial_regression_fitting_02: checking the status is 'poor' with Warning, R2 <= 0.95, __init__ in ResultReport class is covered here
    test_polynomial_regression_fitting_03: checking the status is 'ok', R2 > 0.95, method generate_expression in ResultReport class is covered here
    test_polynomial_regression_fitting_04: checking the status is 'ok', R2 > 0.95 with additional train_data

    test_get_feature_vector_01:
    test_get_feature_vector_02:
       Tests: : Unit tests for get_feature_vector. We verify that:
           - 1: The (key, val) dictionary obtained from the generated IndexParam matches the headers of the data when the input is a dataframe
           - 2: The (key, val) dictionary obtained from the generated IndexParam matches is numerically consistent with that of the input array when a numpy array is supplied. 

    test_set_additional_terms_01: need to check, working with any inputs

    test_poly_training_01: checking the status is 'ok', R2 > 0.95, by running polynomial_regression_fitting, ResultReport class is covered here

    test_generate_expression:   test only while it is running or not (not compared values)
    '''
    def setUp(self):
        # Data generated from the expression (x_1 + 1)^2 + (x_2 + 1) ^ 2 between 0 and 10 for x_1 and x_2
        x1 = np.linspace(0, 10, 21)
        x2 = np.linspace(0, 10, 21)
        y = np.array([ [i,j,((i + 1) ** 2) + ((j + 1) ** 2)] for i in x1 for j in x2])
        self.full_data = pd.DataFrame({'x1': y[:, 0], 'x2': y[:, 1], 'y': y[:, 2]})
        
        # Theoretical sample data from range for(x_1 + 1)^2 + (x_2 + 1) ^ 2 between 0 and 10 for x_1 and x_2
        x1 = np.linspace(0, 10, 5)
        x2 = np.linspace(0, 10, 5)
        self.training_data = np.array([ [i,j,((i + 1) ** 2) + ((j + 1) ** 2)] for i in x1 for j in x2])

        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        self.test_data_numpy = np.array([[i,(i+1)**2] for i in range(10)])
        self.test_data_pandas = pd.DataFrame([[i,(i+1)**2] for i in range(10)])
        self.test_data_large =  np.array([[i,(i+1)**2] for i in range(200)])
        self.test_data_numpy_1d = np.array([[(i+1)**2] for i in range(10)])
        self.test_data_numpy_3d = np.array([[i,(i+1)**2,(i+2)**2] for i in range(10)]) 

        self.sample_points_numpy = np.array([[i,(i+1)**2] for i in range(8)])
        self.sample_points_large =  np.array([[i,(i+1)**2] for i in range(100)])
        self.sample_points_pandas= pd.DataFrame([[i,(i+1)**2] for i in range(8)])
        self.sample_points_numpy_1d = np.array([[(i+1)**2] for i in range(8)])
        self.sample_points_numpy_3d = np.array([[i,(i+1)**2,(i+2)**2] for i in range(8)])



    def test__init__01(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy

        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
        assert PolyClass.max_polynomial_order == 5
        assert PolyClass.number_of_crossvalidations == 3  # Default number of cross-validations
        assert PolyClass.no_adaptive_samples == 4  # Default number of adaptive samples
        assert PolyClass.fraction_training == 0.75  # Default training split
        assert PolyClass.max_fraction_training_samples == 0.5  # Default fraction for the maximum number of training samples
        assert PolyClass.max_iter == 10  # Default maximum number of iterations
        assert PolyClass.solution_method == 'pyomo'  # Default solution_method
        assert PolyClass.multinomials == 1  # Default multinomials
    
    def test__init__02(self):       
        original_data_input= self.test_data_pandas
        regression_data_input = self.sample_points_pandas
        
        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
        assert PolyClass.max_polynomial_order == 5
        assert PolyClass.number_of_crossvalidations == 3  # Default number of cross-validations
        assert PolyClass.no_adaptive_samples == 4  # Default number of adaptive samples
        assert PolyClass.fraction_training == 0.75  # Default training split
        assert PolyClass.max_fraction_training_samples == 0.5  # Default fraction for the maximum number of training samples
        assert PolyClass.max_iter == 10  # Default maximum number of iterations
        assert PolyClass.solution_method == 'pyomo'  # Default solution_method
        assert PolyClass.multinomials == 1  # Default multinomials

    def test__init__03(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_pandas

        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
        assert PolyClass.max_polynomial_order == 5
        assert PolyClass.number_of_crossvalidations == 3  # Default number of cross-validations
        assert PolyClass.no_adaptive_samples == 4  # Default number of adaptive samples
        assert PolyClass.fraction_training == 0.75  # Default training split
        assert PolyClass.max_fraction_training_samples == 0.5  # Default fraction for the maximum number of training samples
        assert PolyClass.max_iter == 10  # Default maximum number of iterations
        assert PolyClass.solution_method == 'pyomo'  # Default solution_method
        assert PolyClass.multinomials == 1  # Default multinomials

    def test__init__04(self):       
        original_data_input= self.test_data_pandas
        regression_data_input = self.sample_points_numpy

        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
        assert PolyClass.max_polynomial_order == 5
        assert PolyClass.number_of_crossvalidations == 3  # Default number of cross-validations
        assert PolyClass.no_adaptive_samples == 4  # Default number of adaptive samples
        assert PolyClass.fraction_training == 0.75  # Default training split
        assert PolyClass.max_fraction_training_samples == 0.5  # Default fraction for the maximum number of training samples
        assert PolyClass.max_iter == 10  # Default maximum number of iterations
        assert PolyClass.solution_method == 'pyomo'  # Default solution_method
        assert PolyClass.multinomials == 1  # Default multinomials
    
    def test__init__05(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy

        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=3, 
                                        number_of_crossvalidations=5,
                                        no_adaptive_samples=6, 
                                        training_split=0.5, 
                                        max_fraction_training_samples=0.4, 
                                        max_iter=20, 
                                        solution_method='MLe', 
                                        multinomials=0)
        assert PolyClass.max_polynomial_order == 3
        assert PolyClass.number_of_crossvalidations == 5  # Default number of cross-validations
        assert PolyClass.no_adaptive_samples == 6  # Default number of adaptive samples
        assert PolyClass.fraction_training == 0.5  # Default training split
        assert PolyClass.max_fraction_training_samples == 0.4  # Default fraction for the maximum number of training samples
        assert PolyClass.max_iter == 20  # Default maximum number of iterations
        assert PolyClass.solution_method == 'mle'  # Default solution_method, doesn't matter lower / upper characters
        assert PolyClass.multinomials == 0  # Default multinomials
    
    def test__init__06(self):       
        original_data_input= list(self.test_data_numpy)
        regression_data_input = self.sample_points_numpy
        with pytest.raises(ValueError):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
    
    def test__init__07(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = list(self.sample_points_numpy)
        with pytest.raises(ValueError):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)

    def test__init__08(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = list(self.sample_points_numpy)
        with pytest.raises(ValueError):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
    
    def test__init__09(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_large
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
    
    def test__init__10(self):       
        original_data_input= self.test_data_numpy_3d
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
    
    def test__init__11(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy_3d
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
    
    def test__init__12(self):       
        original_data_input= self.test_data_numpy_1d
        regression_data_input = self.sample_points_numpy_1d
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)
    
    def test__init__13(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.warns(Warning):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            number_of_crossvalidations=11)
            assert PolyClass.number_of_crossvalidations == 11  # Default number of cross-validations
    
    def test__init__14(self):
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=1.2)
    
    def test__init__15(self):       
        original_data_input= self.test_data_large
        regression_data_input = self.sample_points_large
        with pytest.warns(Warning):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=11)
            assert PolyClass.max_polynomial_order == 10  

    def test__init__16(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            training_split=1)
    
    def test__init__17(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            training_split=-1)

    def test__init__18(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            max_fraction_training_samples=1.2)
    
    def test__init__19(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            max_fraction_training_samples=-1.2)

    def test__init__20(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        PolyClass = PolynomialRegression(regression_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            max_iter=100)
        assert PolyClass.max_iter == 0
    
    def test__init__21(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            no_adaptive_samples=0,
                                            max_iter=100)
        assert PolyClass.max_iter == 0

    def test__init__22(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            number_of_crossvalidations=1.2)
    
    def test__init__23(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            no_adaptive_samples=4.2)
        
    def test__init__24(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            max_iter=4.2)
    
    def test__init__25(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=15)
    
    def test__init__26(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            solution_method=1)
    
    def test__init__27(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            solution_method='idaes')
    
    def test__init__28(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=5,
                                            multinomials=3)

    def test__init__29(self):
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=-2)
    
    def test__init__30(self):
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=3,
                                            number_of_crossvalidations=-3)

    def test__init__31(self):
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=3,
                                            no_adaptive_samples=-3)
    
    def test__init__32(self):
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                            maximum_polynomial_order=3,
                                            max_iter=-3)

        
    def test_training_test_data_creation_01(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        
        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5,
                                        training_split=0.01)
        with pytest.raises(Exception):
            training_data, cross_val_data = PolyClass.training_test_data_creation()
    
    def test_training_test_data_creation_02(self):       
        original_data_input= self.test_data_numpy
        regression_data_input = self.sample_points_numpy
        
        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5,
                                        training_split=0.99)
        with pytest.raises(Exception):
            training_data, cross_val_data = PolyClass.training_test_data_creation()
    
    def test_training_test_data_creation_03(self):       
        original_data_input= self.full_data
        regression_data_input = self.training_data
        
        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)

        training_data, cross_val_data = PolyClass.training_test_data_creation()

        expected_training_size = int(np.around(PolyClass.number_of_samples * PolyClass.fraction_training))
        expected_test_size = PolyClass.regression_data.shape[0] - expected_training_size


        assert len(training_data) == PolyClass.number_of_crossvalidations
        assert len(cross_val_data) == PolyClass.number_of_crossvalidations

        for i in range(1,PolyClass.number_of_crossvalidations+1):
            assert training_data["training_set_"+str(i)].shape[0] == expected_training_size
            assert cross_val_data["test_set_"+str(i)].shape[0] == expected_test_size
        
            concat_01 = (np.concatenate((training_data["training_set_"+str(i)], cross_val_data["test_set_"+str(i)]), axis=0))
            sample_data_sorted = regression_data_input[
                np.lexsort((regression_data_input[:, 2], regression_data_input[:, 1], regression_data_input[:, 0]))]
            concat_01_sorted = concat_01[np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))]
            np.testing.assert_equal(sample_data_sorted, concat_01_sorted)
    
    def test_training_test_data_creation_04(self):       
        original_data_input= self.full_data
        regression_data_input = self.training_data
        
        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=3, 
                                        number_of_crossvalidations=5,
                                        no_adaptive_samples=6, 
                                        training_split=0.5, 
                                        max_fraction_training_samples=0.4, 
                                        max_iter=20, 
                                        solution_method='MLe', 
                                        multinomials=0)

        training_data, cross_val_data = PolyClass.training_test_data_creation()

        expected_training_size = int(np.around(PolyClass.number_of_samples * PolyClass.fraction_training))
        expected_test_size = PolyClass.regression_data.shape[0] - expected_training_size


        assert len(training_data) == PolyClass.number_of_crossvalidations
        assert len(cross_val_data) == PolyClass.number_of_crossvalidations

        for i in range(1,PolyClass.number_of_crossvalidations+1):
            assert training_data["training_set_"+str(i)].shape[0] == expected_training_size
            assert cross_val_data["test_set_"+str(i)].shape[0] == expected_test_size
        
            concat_01 = (np.concatenate((training_data["training_set_"+str(i)], cross_val_data["test_set_"+str(i)]), axis=0))
            sample_data_sorted = regression_data_input[
                np.lexsort((regression_data_input[:, 2], regression_data_input[:, 1], regression_data_input[:, 0]))]
            concat_01_sorted = concat_01[np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))]
            np.testing.assert_equal(sample_data_sorted, concat_01_sorted)    

    def test_training_test_data_creation_05(self):       
        original_data_input= self.full_data
        regression_data_input = self.training_data
        
        PolyClass = PolynomialRegression(original_data_input, regression_data_input,
                                        maximum_polynomial_order=5)

        additional_data_input = np.array([ [i**2,((i + 1) * 2) + ((j + 1) * 2),j**4,((i + 1) * 2) + ((j + 1) ** 2)] for i in range(5) for j in range(5)])        
        training_data, cross_val_data = PolyClass.training_test_data_creation(additional_features=additional_data_input)

        expected_training_size = int(np.around(PolyClass.number_of_samples * PolyClass.fraction_training))
        expected_test_size = PolyClass.regression_data.shape[0] - expected_training_size


        assert len(training_data) == PolyClass.number_of_crossvalidations*2 
        assert len(cross_val_data) == PolyClass.number_of_crossvalidations*2 

        for i in range(1,PolyClass.number_of_crossvalidations+1):
            assert training_data["training_set_"+str(i)].shape[0] == expected_training_size
            assert training_data["training_extras_"+str(i)].shape[0] == expected_training_size
            assert cross_val_data["test_set_"+str(i)].shape[0] == expected_test_size
            assert cross_val_data["test_extras_"+str(i)].shape[0] == expected_test_size
        
            concat_01 = (np.concatenate((training_data["training_set_"+str(i)], cross_val_data["test_set_"+str(i)]), axis=0))
            sample_data_sorted = regression_data_input[
                np.lexsort((regression_data_input[:, 2], regression_data_input[:, 1], regression_data_input[:, 0]))]
            concat_01_sorted = concat_01[np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))]
            np.testing.assert_equal(sample_data_sorted, concat_01_sorted)

            concat_02 = (np.concatenate((training_data["training_extras_"+str(i)], cross_val_data["test_extras_"+str(i)]), axis=0))
            additional_data_sorted = additional_data_input[
                np.lexsort((additional_data_input[:, 3], additional_data_input[:, 2], additional_data_input[:, 1], additional_data_input[:, 0]))]
            concat_02_sorted = concat_02[np.lexsort((concat_02[:, 3],concat_02[:, 2], concat_02[:, 1], concat_02[:, 0]))]
            np.testing.assert_equal(additional_data_sorted, concat_02_sorted)


    def test_polygeneration_01(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=1)
        x_input_train_data = self.training_data[:, :-1]
        poly_degree = 1
        output_1 = data_feed.polygeneration(poly_degree, data_feed.multinomials, x_input_train_data)
        # Create expected solution vector
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = 4  # New number of features should be = 2 * max_polynomial_order + 2 for two input features
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] * x_input_train_data[:, 1]
        # Compare arrays
        np.testing.assert_equal(output_1, expected_output)

    def test_polygeneration_02(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2)
        x_input_train_data = self.training_data[:, :-1]
        poly_degree = 2
        output_1 = data_feed.polygeneration(poly_degree, data_feed.multinomials, x_input_train_data)
        # Create expected solution vector
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = 6  # New number of features should be = 2 * max_polynomial_order + 2 for two input features
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] ** 2
        expected_output[:, 4] = x_input_train_data[:, 1] ** 2
        expected_output[:, 5] = x_input_train_data[:, 0] * x_input_train_data[:, 1]
        # Compare arrays
        np.testing.assert_equal(output_1, expected_output)

    def test_polygeneration_03(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=10)
        x_input_train_data = self.training_data[:, :-1]
        poly_degree = 10
        output_1 = data_feed.polygeneration(poly_degree, data_feed.multinomials, x_input_train_data)
        # Create expected solution vector
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = 22  # New number of features should be = 2 * max_polynomial_order + 2 for two input features
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] ** 2
        expected_output[:, 4] = x_input_train_data[:, 1] ** 2
        expected_output[:, 5] = x_input_train_data[:, 0] ** 3
        expected_output[:, 6] = x_input_train_data[:, 1] ** 3
        expected_output[:, 7] = x_input_train_data[:, 0] ** 4
        expected_output[:, 8] = x_input_train_data[:, 1] ** 4
        expected_output[:, 9] = x_input_train_data[:, 0] ** 5
        expected_output[:, 10] = x_input_train_data[:, 1] ** 5
        expected_output[:, 11] = x_input_train_data[:, 0] ** 6
        expected_output[:, 12] = x_input_train_data[:, 1] ** 6
        expected_output[:, 13] = x_input_train_data[:, 0] ** 7
        expected_output[:, 14] = x_input_train_data[:, 1] ** 7
        expected_output[:, 15] = x_input_train_data[:, 0] ** 8
        expected_output[:, 16] = x_input_train_data[:, 1] ** 8
        expected_output[:, 17] = x_input_train_data[:, 0] ** 9
        expected_output[:, 18] = x_input_train_data[:, 1] ** 9
        expected_output[:, 19] = x_input_train_data[:, 0] ** 10
        expected_output[:, 20] = x_input_train_data[:, 1] ** 10
        expected_output[:, 21] = x_input_train_data[:, 0] * x_input_train_data[:, 1]
        # Compare arrays
        np.testing.assert_equal(output_1, expected_output)

    def test_polygeneration_04(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=10,
                                         multinomials=0)
        x_input_train_data = self.training_data[:, :-1]
        poly_degree = 10
        output_1 = data_feed.polygeneration(poly_degree, data_feed.multinomials, x_input_train_data)
        # Create expected solution vector
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = 21  # New number of features should be = 2 * max_polynomial_order + 2 for two input features
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] ** 2
        expected_output[:, 4] = x_input_train_data[:, 1] ** 2
        expected_output[:, 5] = x_input_train_data[:, 0] ** 3
        expected_output[:, 6] = x_input_train_data[:, 1] ** 3
        expected_output[:, 7] = x_input_train_data[:, 0] ** 4
        expected_output[:, 8] = x_input_train_data[:, 1] ** 4
        expected_output[:, 9] = x_input_train_data[:, 0] ** 5
        expected_output[:, 10] = x_input_train_data[:, 1] ** 5
        expected_output[:, 11] = x_input_train_data[:, 0] ** 6
        expected_output[:, 12] = x_input_train_data[:, 1] ** 6
        expected_output[:, 13] = x_input_train_data[:, 0] ** 7
        expected_output[:, 14] = x_input_train_data[:, 1] ** 7
        expected_output[:, 15] = x_input_train_data[:, 0] ** 8
        expected_output[:, 16] = x_input_train_data[:, 1] ** 8
        expected_output[:, 17] = x_input_train_data[:, 0] ** 9
        expected_output[:, 18] = x_input_train_data[:, 1] ** 9
        expected_output[:, 19] = x_input_train_data[:, 0] ** 10
        expected_output[:, 20] = x_input_train_data[:, 1] ** 10
        # Compare arrays
        np.testing.assert_equal(output_1, expected_output)

    def test_polygeneration_05(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=1)
        x_input_train_data = self.training_data[:, :-1]
        poly_degree = 1
        additional_term = np.sqrt(x_input_train_data)
        output_1 = data_feed.polygeneration(poly_degree, data_feed.multinomials, x_input_train_data, additional_term)
        # Create expected solution vector
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = 6  # New number of features should be = 2 * max_polynomial_order + 4
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] * x_input_train_data[:, 1]
        expected_output[:, 4] = additional_term[:, 0]
        expected_output[:, 5] = additional_term[:, 1]
        # Compare arrays
        np.testing.assert_equal(output_1, expected_output)


    def test_cost_function_01(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.zeros((x_data_nc, 1))
        expected_value = 6613.875  # Calculated externally as sum(y^2) / 2m
        output_1 = PolynomialRegression.cost_function(theta, x_vector, y, reg_parameter=0)
        assert output_1 == expected_value

    def test_cost_function_02(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array([[4.5], [3], [3], [1], [1], [0]])  # coefficients in (x1 + 1.5)^2 + (x2 + 1.5) ^ 2
        expected_value = 90.625  # Calculated externally as sum(dy^2) / 2m
        output_1 = PolynomialRegression.cost_function(theta, x_vector, y, reg_parameter=0)
        assert output_1 == expected_value

    def test_cost_function_03(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array([[2], [2], [2], [1], [1], [0]])  # Actual coefficients in (x1 + 1)^2 + (x2 + 1) ^ 2
        expected_value = 0  # Value should return zero for exact solution
        output_1 = PolynomialRegression.cost_function(theta, x_vector, y, reg_parameter=0)
        assert output_1 == expected_value


    def test_gradient_function_01(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.zeros((x_data_nc,))
        expected_value = np.array(
            [[-97], [-635], [-635], [-5246.875], [-5246.875], [-3925]])  # Calculated externally: see Excel sheet
        expected_value = expected_value.reshape(expected_value.shape[0], )
        output_1 = PolynomialRegression.gradient_function(theta, x_vector, y, reg_parameter=0)
        np.testing.assert_equal(output_1, expected_value)

    def test_gradient_function_02(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array([[4.5], [3], [3], [1], [1], [0]])  # coefficients in (x1 + 1.5)^2 + (x2 + 1.5) ^ 2
        theta = theta.reshape(theta.shape[0], )
        expected_value = np.array(
            [[12.5], [75], [75], [593.75], [593.75], [437.5]])  # Calculated externally: see Excel sheet
        expected_value = expected_value.reshape(expected_value.shape[0], )
        output_1 = PolynomialRegression.gradient_function(theta, x_vector, y, reg_parameter=0)
        np.testing.assert_equal(output_1, expected_value)

    def test_gradient_function_03(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array([[2], [2], [2], [1], [1], [0]])  # Actual coefficients in (x1 + 1)^2 + (x2 + 1) ^ 2
        theta = theta.reshape(theta.shape[0], )
        expected_value = np.array([[0], [0], [0], [0], [0], [0]])  # Calculated externally: see Excel sheet
        expected_value = expected_value.reshape(expected_value.shape[0], )
        output_1 = PolynomialRegression.gradient_function(theta, x_vector, y, reg_parameter=0)
        np.testing.assert_equal(output_1, expected_value)

    
    def test_bfgs_parameter_optimization_01(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = x[:, ] ** 2
        x_vector[:, 1] = x[:, ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.], [2.], [1.]]).reshape(3, )
        data_feed = PolynomialRegression(self.test_data_numpy, input_array, maximum_polynomial_order=5,
                                         solution_method='bfgs')
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == 'bfgs'
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    def test_bfgs_parameter_optimization_02(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        x = self.training_data[:, : -1]
        y = self.training_data[:, -1]
        x_vector = np.zeros((x.shape[0], 6))
        x_vector[:, 0] = x[:, 0] ** 2
        x_vector[:, 1] = x[:, 1] ** 2
        x_vector[:, 2] = x[:, 0]
        x_vector[:, 3] = x[:, 1]
        x_vector[:, 4] = x[:, 1] * x[:, 0]
        x_vector[:, 5] = 1
        expected_value = np.array([[1.], [1.], [2.], [2.], [0.], [2.]]).reshape(6, )
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=4,
                                         solution_method='bfgs')
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == 'bfgs'
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    
    @staticmethod
    def test_mle_estimate_01():
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = x[:, ] ** 2
        x_vector[:, 1] = x[:, ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.], [2.], [1.]]).reshape(3, )
        output_1 = PolynomialRegression.MLE_estimate(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    def test_mle_estimate_02(self):
        x = self.training_data[:, : -1]
        y = self.training_data[:, -1]
        x_vector = np.zeros((x.shape[0], 6))
        x_vector[:, 0] = x[:, 0] ** 2
        x_vector[:, 1] = x[:, 1] ** 2
        x_vector[:, 2] = x[:, 0]
        x_vector[:, 3] = x[:, 1]
        x_vector[:, 4] = x[:, 1] * x[:, 0]
        x_vector[:, 5] = 1
        expected_value = np.array([[1.], [1.], [2.], [2.], [0.], [2.]]).reshape(6, )
        output_1 = PolynomialRegression.MLE_estimate(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    def test_pyomo_optimization_01(self):    
        x_vector = np.array([[i**2, i, 1 ] for i in range(10) ]  )
        y = np.array([[i**2] for i in range(1,11) ]  )
        expected_value = np.array([[1.], [2.], [1.]])
        output_1 = PolynomialRegression.pyomo_optimization(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))
    
    def test_pyomo_optimization_02(self):
        x = self.training_data[:, : -1]
        y = self.training_data[:, -1]
        x_vector = np.zeros((x.shape[0], 6))
        x_vector[:, 0] = x[:, 0] ** 2
        x_vector[:, 1] = x[:, 1] ** 2
        x_vector[:, 2] = x[:, 0]
        x_vector[:, 3] = x[:, 1]
        x_vector[:, 4] = x[:, 1] * x[:, 0]
        x_vector[:, 5] = 1
        expected_value = np.array([[1.], [1.], [2.], [2.], [0.], [2.]])
        output_1 = PolynomialRegression.pyomo_optimization(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))


    def test_cross_validation_error_calculation_01(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1].reshape(self.training_data.shape[0], 1)
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.zeros((x_data_nc, 1))
        expected_value = 2 * 6613.875  # Calculated externally as sum(y^2) / m
        output_1 = PolynomialRegression.cross_validation_error_calculation(theta, x_vector, y)
        assert output_1 == expected_value

    def test_cross_validation_error_calculation_02(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1].reshape(self.training_data.shape[0], 1)
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array([[4.5], [3], [3], [1], [1], [0]])  # coefficients in (x1 + 1.5)^2 + (x2 + 1.5) ^ 2
        expected_value = 2 * 90.625  # Calculated externally as sum(dy^2) / 2m
        output_1 = PolynomialRegression.cross_validation_error_calculation(theta, x_vector, y)
        assert output_1 == expected_value

    def test_cross_validation_error_calculation_03(self):
        x = self.training_data[:, :-1]
        y = self.training_data[:, -1].reshape(self.training_data.shape[0], 1)
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array([[2], [2], [2], [1], [1], [0]])  # Actual coefficients in (x1 + 1)^2 + (x2 + 1) ^ 2
        expected_value = 2 * 0  # Value should return zero for exact solution
        output_1 = PolynomialRegression.cross_validation_error_calculation(theta, x_vector, y)
        assert output_1 == expected_value


    def mock_optimization(self, x, y):
        return 10 * np.ones((x.shape[1], 1))

    @patch.object(PolynomialRegression, 'MLE_estimate', mock_optimization)
    def test_polyregression_01(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=5, solution_method='mle')
        poly_order = 2
        training_data = self.training_data[0:20, :]
        test_data = self.training_data[20:, :]
        expected_output = 10 * np.ones((6, 1))
        output_1, _, _ = data_feed.polyregression(poly_order, training_data, test_data)
        np.testing.assert_array_equal(expected_output, output_1)

    @patch.object(PolynomialRegression, 'bfgs_parameter_optimization', mock_optimization)
    def test_polyregression_02(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=5, solution_method='bfgs')
        poly_order = 2
        training_data = self.training_data[0:20, :]
        test_data = self.training_data[20:, :]
        expected_output = 10 * np.ones((6, 1))
        output_1, _, _ = data_feed.polyregression(poly_order, training_data, test_data)
        np.testing.assert_array_equal(expected_output, output_1)

    @patch.object(PolynomialRegression, 'pyomo_optimization', mock_optimization)
    def test_polyregression_03(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=5)
        poly_order = 2
        training_data = self.training_data[0:20, :]
        test_data = self.training_data[20:, :]
        expected_output = 10 * np.ones((6, 1))
        output_1, _, _ = data_feed.polyregression(poly_order, training_data, test_data)
        np.testing.assert_array_equal(expected_output, output_1)

    def test_polyregression_04(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=5)
        poly_order = 10
        training_data = self.training_data[0:20, :]
        test_data = self.training_data[20:, :]
        expected_output = np.Inf
        output_1, output_2, output_3 = data_feed.polyregression(poly_order, training_data, test_data)
        np.testing.assert_array_equal(expected_output, output_1)
        np.testing.assert_array_equal(expected_output, output_2)
        np.testing.assert_array_equal(expected_output, output_3)


    def test_surrogate_performance_01(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        order_best = 2
        phi_best = np.array([[0.], [0.], [0.]])
        expected_value_1 = 38.5
        expected_value_2 = 2533.3
        expected_value_3 = -1.410256
        expected_value_4 = 0
        data_feed = PolynomialRegression(self.test_data_numpy, input_array, maximum_polynomial_order=5)
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(phi_best, order_best)
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2
        assert np.round(output_3, 4) == np.round(expected_value_3, 4)
        assert np.round(output_4, 4) == np.round(expected_value_4, 4)

    def test_surrogate_performance_02(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        order_best = 2
        phi_best = np.array([[1.], [2.], [1.]])
        expected_value_1 = 0
        expected_value_2 = 0
        expected_value_3 = 1
        expected_value_4 = 1
        data_feed = PolynomialRegression(self.test_data_numpy, input_array, maximum_polynomial_order=5)
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(phi_best, order_best)
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2
        assert np.round(output_3, 4) == np.round(expected_value_3, 4)
        assert np.round(output_4, 4) == np.round(expected_value_4, 4)

    def test_surrogate_performance_03(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        order_best = 2
        phi_best = np.array([[1.], [1.], [1.]])
        expected_value_1 = 4.5
        expected_value_2 = 28.5
        expected_value_3 = 0.972884259
        expected_value_4 = 0.931219147
        data_feed = PolynomialRegression(self.test_data_numpy, input_array, maximum_polynomial_order=5)
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(phi_best, order_best)
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2
        assert np.round(output_3, 4) == np.round(expected_value_3, 4)
        assert np.round(output_4, 4) == np.round(expected_value_4, 4)

  
    def test_results_generation_01(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=3, multinomials=0)
        order = 1
        beta = np.array([ [0], [0], [0] ])
        expected_df = pd.Series()
        row_list = np.array([['k'], ['(x_1)^1'], ['(x_2)^1']])
        expected_df = expected_df.append(pd.Series({row_list[0, 0]: beta[0, 0], row_list[1, 0]: beta[1, 0], row_list[2, 0]: beta[2, 0]}))
        output_df = data_feed.results_generation(beta, order)
        assert output_df.index.to_list() == expected_df.index.to_list()
        assert expected_df.all() == output_df.all()

    def test_results_generation_02(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=3, multinomials=0)
        order = 3
        beta = np.array([ [1], [0.3], [6], [500], [500000], [0.001], [50] ])
        expected_df = pd.Series()
        row_list = np.array([['k'], ['(x_1)^1'], ['(x_2)^1'], ['(x_1)^2'], ['(x_2)^2'], ['(x_1)^3'], ['(x_2)^3']])
        expected_df = expected_df.append(pd.Series({row_list[0, 0]: beta[0, 0], row_list[1, 0]: beta[1, 0], row_list[2, 0]: beta[2, 0], row_list[3, 0]: beta[3, 0], row_list[4, 0]: beta[4, 0], row_list[5, 0]: beta[5, 0], row_list[6, 0]: beta[6, 0]}))
        output_df = data_feed.results_generation(beta, order)
        assert output_df.index.to_list() == expected_df.index.to_list()
        assert expected_df.all() == output_df.all()

    def test_results_generation_03(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=3, multinomials=1)
        order = 2
        beta = np.array([ [1], [0.3], [6], [500], [500000], [0.001]])
        expected_df = pd.Series()
        row_list = np.array([['k'], ['(x_1)^1'], ['(x_2)^1'], ['(x_1)^2'], ['(x_2)^2'], ['(x_1).(x_2)']])
        expected_df = expected_df.append(pd.Series({row_list[0, 0]: beta[0, 0], row_list[1, 0]: beta[1, 0], row_list[2, 0]: beta[2, 0], row_list[3, 0]: beta[3, 0], row_list[4, 0]: beta[4, 0], row_list[5, 0]: beta[5, 0]}))
        output_df = data_feed.results_generation(beta, order)
        assert output_df.index.to_list() == expected_df.index.to_list()
        assert expected_df.all() == output_df.all()

    @patch("matplotlib.pyplot.show")
    def test_error_plotting(self, mock_show):
        mock_show.return_value = None
        # Generate typical data values for eaxch variable in order to test plot function
        plotting_data = np.array(
            [[1, 7, 23.1, 29.2, 0.01, 1300, -1.6, -1.1, 5], [2, 9, 0.0055, 0.0015, 0.006, 159, 0.34, 0.19, 10],
             [3, 6, 0.0007, 0.0009, 0.001, 2.3, 0.998, 0.994, 15], [4, 2, 0.0005, 0.0004, 0.0008, 0.02, 1.0, 1.0, 20]])
        ax1, ax2, ax3, ax4 = PolynomialRegression.error_plotting(plotting_data)
        expected_output_1 = np.array([[1, 23.1], [2, 0.0055], [3, 0.0007], [4, 0.0005]])
        expected_output_2 = np.array([[1, 29.2], [2, 0.0015], [3, 0.0009], [4, 0.0004]])
        expected_output_3 = np.array([[1, 0.01], [2, 0.006], [3, 0.001], [4, 0.0008]])
        expected_output_4 = np.array([[1, 1300], [2, 159], [3, 2.3], [4, 0.02]])
        expected_output_5 = np.array([[1, -1.6], [2, 0.34], [3, 0.998], [4, 1.0]])
        expected_output_6 = np.array([[1, -1.1], [2, 0.19], [3, 0.994], [4, 1.0]])
        # Compare the data input to the actual plotted data
        np.testing.assert_array_equal(ax1.lines[0].get_xydata(), expected_output_1)
        np.testing.assert_array_equal(ax1.lines[1].get_xydata(), expected_output_2)
        np.testing.assert_array_equal(ax2.lines[0].get_xydata(), expected_output_3)
        np.testing.assert_array_equal(ax3.lines[0].get_xydata(), expected_output_4)
        np.testing.assert_array_equal(ax4.lines[0].get_xydata(), expected_output_5)
        np.testing.assert_array_equal(ax4.lines[1].get_xydata(), expected_output_6)
        # Compare the expected subplot titles to the actual display
        assert ax1.title.get_text() == 'Training (green) vs Cross-validation error (red)'
        assert ax2.title.get_text() == 'MAE'
        assert ax3.title.get_text() == 'MSE'
        assert ax4.title.get_text() == 'R-squared (blue) and Adjusted R-squared (red)'
        # Compare the line/curve segments for each line plot to the actual input data
        np.testing.assert_array_equal(ax1.lines[0].get_path()._vertices, expected_output_1)
        np.testing.assert_array_equal(ax1.lines[1].get_path()._vertices, expected_output_2)
        np.testing.assert_array_equal(ax2.lines[0].get_path()._vertices, expected_output_3)
        np.testing.assert_array_equal(ax3.lines[0].get_path()._vertices, expected_output_4)
        np.testing.assert_array_equal(ax4.lines[0].get_path()._vertices, expected_output_5)
        np.testing.assert_array_equal(ax4.lines[1].get_path()._vertices, expected_output_6)


    def test_user_defined_terms_01(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        additional_terms = [np.sin(self.training_data[:, 0]), np.sin(self.training_data[:, 1])]
        returned_features_array = data_feed.user_defined_terms(additional_terms)

        # Create two additional columns
        expected_array = np.zeros((self.training_data.shape[0], 2))
        expected_array[:, 0] = np.sin(self.training_data[:, 0])
        expected_array[:, 1] = np.sin(self.training_data[:, 1])

        np.testing.assert_equal(expected_array, returned_features_array)

    def test_user_defined_terms_02(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.full_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        additional_terms = [np.sin(self.full_data['x1']), np.sin(self.full_data['x2'])]
        returned_features_array = data_feed.user_defined_terms(additional_terms)

        # Create two additional columns
        expected_array = np.zeros((self.full_data.values.shape[0], 2))
        expected_array[:, 0] = np.sin(self.full_data.values[:, 0])
        expected_array[:, 1] = np.sin(self.full_data.values[:, 1])

        np.testing.assert_equal(expected_array, returned_features_array)

    def test_user_defined_terms_03(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        additional_terms = [np.sin(self.full_data['x1']), np.sin(self.full_data['x2'])]
        with pytest.raises(Exception):
            data_feed.user_defined_terms(additional_terms)

    def test_user_defined_terms_04(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        p = np.sin(self.training_data[:, 0]).reshape(self.training_data.shape[0], 1)  # 2-D array, should raise error
        additional_terms = [p]
        with pytest.raises(Exception):
            data_feed.user_defined_terms(additional_terms)

    def test_user_defined_terms_05(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        p = np.sin(self.training_data[:, 0]).reshape(self.training_data.shape[0], 1)
        additional_terms = p  # Additional terms as array, not list
        with pytest.raises(ValueError):
            data_feed.user_defined_terms(additional_terms)


    def test_polynomial_regression_fitting_01(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2)
        data_feed.get_feature_vector()
        results = data_feed.polynomial_regression_fitting()
        assert results.fit_status == 'ok'

    @patch("matplotlib.pyplot.show")
    def test_polynomial_regression_fitting_02(self, mock_show):
        mock_show.return_value = None
        data_feed = PolynomialRegression(self.full_data, self.training_data[:5], maximum_polynomial_order=1)
        data_feed.get_feature_vector()
        with pytest.warns(Warning):
            results = data_feed.polynomial_regression_fitting()
            assert results.fit_status == 'poor'

    def test_polynomial_regression_fitting_03(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2)
        data_feed.get_feature_vector()
        results = data_feed.polynomial_regression_fitting()
        x_input_train_data = self.training_data[:, :-1]
        assert results.fit_status == 'ok'
    
    def test_polynomial_regression_fitting_04(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2)
        data_feed.get_feature_vector()
        additional_regression_features = [np.sin(self.training_data[:, 0]), np.sin(self.training_data[:, 1])]    
        results = data_feed.polynomial_regression_fitting(additional_regression_features)      
        results = data_feed.polynomial_regression_fitting()
        assert results.fit_status == 'ok'


    def test_get_feature_vector_01(self):
        data_feed = PolynomialRegression(self.full_data, self.full_data, maximum_polynomial_order=3, multinomials=0)
        output = data_feed.get_feature_vector()
        expected_dict =  {'x1': 0, 'x2': 0}
        assert expected_dict == output.extract_values()

    def test_get_feature_vector_02(self):
        data_feed = PolynomialRegression(self.training_data, self.training_data, maximum_polynomial_order=3, multinomials=0)
        output = data_feed.get_feature_vector()
        expected_dict =  {0: 0, 1: 0}
        assert expected_dict == output.extract_values()

    
    def test_set_additional_terms_01(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=3)
        data_feed.set_additional_terms('a')
        assert 'a' == data_feed.additional_term_expressions
        data_feed.set_additional_terms(1)
        assert 1 == data_feed.additional_term_expressions
        data_feed.set_additional_terms([1,2])
        assert [1,2] == data_feed.additional_term_expressions
        data_feed.set_additional_terms(np.array([1,2]))
        np.testing.assert_equal(np.array([1,2]), data_feed.additional_term_expressions)
        

    def test_poly_training_01(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2)
        data_feed.get_feature_vector()
        results = data_feed.poly_training()
        assert results.fit_status == 'ok'
        
     
    def test_generate_expression(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2)
        
        p =data_feed.get_feature_vector()
        results = data_feed.poly_training()

        lv =[]
        for i in p.keys():
            lv.append(p[i])
        poly_expr = results.generate_expression((lv))
        
        

if __name__ == '__main__':
    unittest.main()
    



