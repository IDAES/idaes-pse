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
import sys
import os
sys.path.append(os.path.abspath('..'))# current folder is ~/tests

from idaes.surrogate.pysmo.radial_basis_function import RadialBasisFunctions, FeatureScaling
import numpy as np
import pandas as pd
import pyutilib.th as unittest
from unittest.mock import patch
from scipy.spatial import distance
import pytest

'''
coverage run test_radial_basis_function.py
coverage report -m
coverage html

'''

class FeatureScalingTestCases(unittest.TestCase):
    """
    test_data_scaling_minmax_01: Test behaviour when input is a numpy array and with 1D array.
    test_data_scaling_minmax_02: Test behaviour when input is a numpy array and with 2D array.
    test_data_scaling_minmax_03: Test behaviour when input is a numpy array and with 3D array.
    test_data_scaling_minmax_04: Test behaviour when input is a numpy array and with 3D array with a varibale is constant.
    test_data_scaling_minmax_05: Test behaviour list input TypeError

    test_data_scaling_minmax_06: Test behaviour when input is a Pandas DataFrame and with 1D array.
    test_data_scaling_minmax_07: Test behaviour when input is a Pandas DataFrame and with 2D array.
    test_data_scaling_minmax_08: Test behaviour when input is a Pandas DataFrame and with 3D array.
    test_data_scaling_minmax_09: Test behaviour when input is a Pandas DataFrame and with 3D array with a varibale is constant.

    test_data_unscaling_minmax_01: Test behaviour when input is a numpy array and with 1D array.
    test_data_unscaling_minmax_02: Test behaviour when input is a numpy array and with 2D array.
    test_data_unscaling_minmax_03: Test behaviour when input is a numpy array and with 3D array. 
    test_data_unscaling_minmax_04: Test behaviour when input is a numpy array and with 3D array with a varibale is constant.

    test_data_unscaling_minmax_05: Test behaviour IndexError when input array size > array size
    test_data_unscaling_minmax_06: Test behaviour IndexError when input array size < array size  

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
        
    
    def test_data_scaling_minmax_01(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_numpy_1d

        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1.reshape(10, 1))

    def test_data_scaling_minmax_02(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 100]])
        expected_output_2 = np.array([[0, 1]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_minmax_03(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_numpy_3d 
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 19, 119]])
        expected_output_2 = np.array([[0, 10, 11]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_minmax_04(self):
        # 3D with constant
        # Sample data generated from the expression (x_1 + 1)^2 + 10, x1: between 0 and 9
        input_array = self.test_data_numpy_3d_constant
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 10, 110]])
        expected_output_2 = np.array([[0, 10, 11]])
        scale = expected_output_3 - expected_output_2
        scale[scale == 0.0] = 1.0
        expected_output_1 = (input_array - expected_output_2) / scale
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_minmax_05(self):
        # TypeError with list
        input_array = self.test_data_numpy_2d.tolist()
        with pytest.raises(TypeError):
            FeatureScaling.data_scaling_minmax(input_array)

    def test_data_scaling_minmax_06(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_pandas_1d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    def test_data_scaling_minmax_07(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array =self.test_data_pandas_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 100]])
        expected_output_2 = np.array([[0, 1]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_minmax_08(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_pandas_3d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 19, 119]])
        expected_output_2 = np.array([[0, 10, 11]])
        expected_output_1 = (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    
    def test_data_scaling_minmax_09(self):
        # 3D with constant
        # Sample data generated from the expression (x_1 + 1)^2 + 10, x1: between 0 and 9
        input_array = self.test_data_pandas_3d_constant
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 10, 110]])
        expected_output_2 = np.array([[0, 10, 11]])
        scale = expected_output_3 - expected_output_2
        scale[scale == 0.0] = 1.0
        expected_output_1 = (input_array - expected_output_2) / scale
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)
    

    def test_data_unscaling_minmax_01(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_numpy_1d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        output_1 = output_1.reshape(output_1.shape[0], )
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array.reshape(10, 1))

    def test_data_unscaling_minmax_02(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)
    
    def test_data_unscaling_minmax_03(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_numpy_3d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)
    
    def test_data_unscaling_minmax_04(self):
        # 3D with constant
        # Sample data generated from the expression (x_1 + 1)^2 + 10, x1: between 0 and 9
        input_array = self.test_data_numpy_3d_constant
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    def test_data_unscaling_minmax_05(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)

        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)

    def test_data_unscaling_minmax_06(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
                
        min_array = np.array([[1,2,3]])
        max_array = np.array([[5,6,7]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)


class RadialBasisFunctionTestCases(unittest.TestCase):
    '''
    test__init__01: Check that all default values are correctly loaded when input is sparse with both input arrays are numpy
    test__init__02: Check that all default values are correctly loaded when input is sparse with both input arrays are pandas
    test__init__03: Check that all user input values are correctly loaded when input is sparse with both input arrays are numpy
    test__init__04: Test behaviour raise Exception with list inputs
    test__init__05: Test behaviour raise Exception with invalid solution method type
    test__init__06: Test behaviour raise Exception with invalid solution method value
    test__init__07: Test behaviour raise Exception with invalid basis_function type
    test__init__08: Test behaviour raise Exception with invalid basis_function value
    test__init__09: Test behaviour raise Exception with invalid regularization type

    test_r2_distance: 
        Tests: : Unit tests for eucl_distance, a function that evaluates the distance between a point c and a set of points in a numpy array.
            The test checks that the function is able to calculate the distance from a single point (n x 1 row vector) to a set of design points (supplied in an n x m array)

    test_gaussian_basis_transformation: Tests for the basis transformations, test the results returned by each transformation for a wide range of values.
    test_linear_transformation: Tests for the basis transformations, test the results returned by each transformation for a wide range of values.
    test_cubic_transformation: Tests for the basis transformations, test the results returned by each transformation for a wide range of values.
    test_multiquadric_basis_transformation: Tests for the basis transformations, test the results returned by each transformation for a wide range of values.
    test_inverse_multiquadric_basis_transformation: Tests for the basis transformations, test the results returned by each transformation for a wide range of values.
    test_thin_plate_spline_transformation: Tests for the basis transformations, test the results returned by each transformation for a wide range of values.

    test_basis_generation:
        Tests: : Tests for the function basis_generation.
        The test checks that for a small dataset, the expected result is returned by the function for each basis transformation type.

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

    test_explicit_linear_algebra_solution_01:
    test_explicit_linear_algebra_solution_02:
        Tests: : Unit tests for explicit_linear_algebra_solution function. The same two unit tests are done:
            The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
            The second test is used to check that higher order multi-variable problems are solved correctly.
            - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
            - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.

    test_pyomo_optimization_01:
    test_pyomo_optimization_02:
        Tests: : Unit tests for pyomo_optimization function. As with the other parameter estimation methods, two unit tests are done:
            The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
            The second test is used to check that higher order multi-variable problems are solved correctly.
            - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
            - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.

    test_error_calculation_01:
    test_error_calculation_02:
    test_error_calculation_03:
        Tests: : Unit tests for error_calculation function.  The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
            The three tests used for the cost_function are repeated for the RMSE and SSE errors. 
            - The first test checks that the correct cv error is returned when theta is initialized as a vector of zeros.
            - The second test checks that the correct cv error is returned for a random point within the solution space.
            - The third test checks that the cv error is zero when the exact solution is found.
        In each case, the result is expected to be twice as large as the cost function's expected value: the cost function is half of the total squared error.
  
    test_r2_calculation_01:
    test_r2_calculation_02:
        Tests: : Unit tests for r2_calculation function.  The function is tested for two cases where the predicted values 
        deviate from the actual by 5% and 50% respectively. The R2-value in each case is compared against Excel evaluated values.
  
    test_loo_error_estimation_with_rippa_method_01:
    test_loo_error_estimation_with_rippa_method_02:
    test_loo_error_estimation_with_rippa_method_03:
        Tests: : Unit tests for loo_error_estimation_with_rippa_method. 
            - Three tests are run to check that the correct condition number and errors are returned by the method irrespective of solution_method:
                - First test checks the algebraic method
                - Second test checks the pyomo method
                - Third test checks the bfgs solution
            - In each test, the basis_generation method and relevant optimization method are replaced with mocks. The results from the mocks
            are used to compute the condition number and LOOCV errors.

    test_leave_one_out_crossvalidation_01: Test leave_one_out_crossvalidation with all default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_02: Test leave_one_out_crossvalidation with basis_function='cubic', all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_03: Test leave_one_out_crossvalidation with basis_function='linear', all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_04: Test leave_one_out_crossvalidation with basis_function='spline', all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_05: Test leave_one_out_crossvalidation with basis_function='gaussian', all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_06: Test leave_one_out_crossvalidation with basis_function='mq', all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_07: Test leave_one_out_crossvalidation with basis_function='imq', all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_08: Test leave_one_out_crossvalidation with solution_method='algebraic', all other default inputs  - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_09: Test leave_one_out_crossvalidation with solution_method='BFGS', all other default inputs  - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_10: Test leave_one_out_crossvalidation with solution_method='pyomo', all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")
    test_leave_one_out_crossvalidation_11: Test leave_one_out_crossvalidation with regularization=True, all other default inputs - (all calulations test with method "loo_error_estimation_with_rippa_method")

    test_rbf_training_01: Test rbf_training with solution_method='algebraic'. All methods are already tested. Report class is tested here with x_condition_number < (1 / np.finfo(float).eps)
    test_rbf_training_02: Test rbf_training with solution_method='pyomo'. All methods are already tested. Report class is tested here with x_condition_number >= (1 / np.finfo(float).eps), raise warning.
    test_rbf_training_03: Test rbf_training with solution_method='bfgs'. All methods are already tested. Report class is tested here with x_condition_number >= (1 / np.finfo(float).eps), raise warning.

    test_rbf_predict_output_01:
    test_rbf_predict_output_02:
    test_rbf_predict_output_03:
    test_rbf_predict_output_04:
    test_rbf_predict_output_05:
    test_rbf_predict_output_06:
        Tests: : Unit tests for rbf_predict_output, a function that generates predictions for a test dataset. 
             For a small dataset, we verify that the current predictions are produced for each transformation type.

    test_get_feature_vector_01:
    test_get_feature_vector_02:
        Tests: : Unit tests for get_feature_vector. We verify that:
            - 1: The (key, val) dictionary obtained from the generated IndexParam matches the headers of the data when the input is a dataframe
            - 2: The (key, val) dictionary obtained from the generated IndexParam matches is numerically consistent with that of the input array when a numpy array is supplied. 

    test_rbf_generate_expression_01: test only while it is running with self.basis_function == 'linear' or not (not compared values)
    test_rbf_generate_expression_02: test only while it is running with self.basis_function == 'cubic' or not (not compared values)
    test_rbf_generate_expression_03: test only while it is running with self.basis_function == 'gaussian' or not (not compared values)
    test_rbf_generate_expression_04: test only while it is running with self.basis_function == 'mq' or not (not compared values)
    test_rbf_generate_expression_05: test only while it is running with self.basis_function == 'imq' or not (not compared values)
    test_rbf_generate_expression_06: test only while it is running with self.basis_function == 'spline' or not (not compared values)
    ''' 

    def setUp(self):
        # Data generated from the expression (x_1 + 1)^2 + (x_2 + 1) ^ 2 between 0 and 10 for x_1 and x_2
        x1 = np.linspace(1, 10, 21)
        x2 = np.linspace(1, 10, 21)
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
        RbfClass = RadialBasisFunctions(self.test_data_numpy, basis_function=None, solution_method=None, regularization=None)
        assert RbfClass.solution_method == 'algebraic'
        assert RbfClass.basis_function == 'gaussian'
        assert RbfClass.regularization == True
    
    def test__init__02(self):
        RbfClass = RadialBasisFunctions(self.test_data_pandas, basis_function=None, solution_method=None, regularization=None)
        assert RbfClass.solution_method == 'algebraic'
        assert RbfClass.basis_function == 'gaussian'
        assert RbfClass.regularization == True
    
    def test__init__03(self):
        RbfClass = RadialBasisFunctions(self.test_data_numpy, basis_function='LineaR', solution_method='PyoMo', regularization=False)
        assert RbfClass.solution_method == 'pyomo'
        assert RbfClass.basis_function == 'linear'
        assert RbfClass.regularization == False
    
    def test__init__04(self):
        with pytest.raises(Exception):
            RbfClass = RadialBasisFunctions([1,2,3,4], basis_function='LineaR', solution_method='PyoMo', regularization=False)
        
    def test__init__05(self):
        with pytest.raises(Exception):
            RbfClass = RadialBasisFunctions(self.test_data_numpy, basis_function=None, solution_method=1, regularization=None)

    def test__init__06(self):
        with pytest.raises(Exception):
            RbfClass = RadialBasisFunctions(self.test_data_numpy, basis_function=None, solution_method='idaes', regularization=None)
    
    def test__init__07(self):
        with pytest.raises(Exception):
            RbfClass = RadialBasisFunctions(self.test_data_numpy, basis_function=1, solution_method=None, regularization=None)
    
    def test__init__08(self):
        with pytest.raises(Exception):
            RbfClass = RadialBasisFunctions(self.test_data_numpy, basis_function='idaes', solution_method=None, regularization=None)
    
    def test__init__09(self):
        with pytest.raises(Exception):
            RbfClass = RadialBasisFunctions(self.test_data_numpy, basis_function=None, solution_method=None, regularization=1)

    
    def test_r2_distance(self):
        u = np.array([[0.1, 0.9]])
        data_feed = RadialBasisFunctions(self.training_data)
        output = data_feed.r2_distance(u)

        scaled=FeatureScaling.data_scaling_minmax(self.training_data)
        scaled = scaled[0]
        scaled_x = scaled[:,:-1]
        expected_output = np.sqrt(np.sum(np.square(scaled_x-u),axis=1))
        np.testing.assert_almost_equal(expected_output, output, decimal=6)


    def test_gaussian_basis_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        shape_list = [0.001, 1, 1000]
        expected_output_1 = np.exp(-1 * ((d_vec * shape_list[0]) ** 2))
        expected_output_2 = np.exp(-1 * ((d_vec * shape_list[1]) ** 2))
        expected_output_3 = np.exp(-1 * ((d_vec * shape_list[2]) ** 2))
        output_1 = RadialBasisFunctions.gaussian_basis_transformation(d_vec, shape_list[0])
        output_2 = RadialBasisFunctions.gaussian_basis_transformation(d_vec, shape_list[1])
        output_3 = RadialBasisFunctions.gaussian_basis_transformation(d_vec, shape_list[2])
        np.testing.assert_array_equal(expected_output_1, output_1)
        np.testing.assert_array_equal(expected_output_2, output_2)
        np.testing.assert_array_equal(expected_output_3, output_3)

    def test_linear_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        output_1 = RadialBasisFunctions.linear_transformation(d_vec)
        np.testing.assert_array_equal(d_vec, output_1)

    def test_cubic_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        expected_output = d_vec ** 3
        output = RadialBasisFunctions.cubic_transformation(d_vec)
        np.testing.assert_array_equal(expected_output, output)

    def test_multiquadric_basis_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        shape_list = [0.001, 1, 1000]
        expected_output_1 = np.sqrt(((d_vec * shape_list[0]) ** 2) + 1)
        expected_output_2 = np.sqrt(((d_vec * shape_list[1]) ** 2) + 1)
        expected_output_3 = np.sqrt(((d_vec * shape_list[2]) ** 2) + 1)
        output_1 = RadialBasisFunctions.multiquadric_basis_transformation(d_vec, shape_list[0])
        output_2 = RadialBasisFunctions.multiquadric_basis_transformation(d_vec, shape_list[1])
        output_3 = RadialBasisFunctions.multiquadric_basis_transformation(d_vec, shape_list[2])
        np.testing.assert_array_equal(expected_output_1, output_1)
        np.testing.assert_array_equal(expected_output_2, output_2)
        np.testing.assert_array_equal(expected_output_3, output_3)

    def test_inverse_multiquadric_basis_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        shape_list = [0.001, 1, 1000]
        expected_output_1 = 1 / np.sqrt(((d_vec * shape_list[0]) ** 2) + 1)
        expected_output_2 = 1 / np.sqrt(((d_vec * shape_list[1]) ** 2) + 1)
        expected_output_3 = 1 / np.sqrt(((d_vec * shape_list[2]) ** 2) + 1)
        output_1 = RadialBasisFunctions.inverse_multiquadric_basis_transformation(d_vec, shape_list[0])
        output_2 = RadialBasisFunctions.inverse_multiquadric_basis_transformation(d_vec, shape_list[1])
        output_3 = RadialBasisFunctions.inverse_multiquadric_basis_transformation(d_vec, shape_list[2])
        np.testing.assert_array_equal(expected_output_1, output_1)
        np.testing.assert_array_equal(expected_output_2, output_2)
        np.testing.assert_array_equal(expected_output_3, output_3)

    def test_thin_plate_spline_transformation(self):
        d_vec = np.array([[5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70], [50, np.NaN]])
        expected_output = np.nan_to_num(d_vec ** 2 * np.log(d_vec))
        output = RadialBasisFunctions.thin_plate_spline_transformation(d_vec)
        np.testing.assert_array_equal(expected_output, output)
        

    def test_basis_generation(self):

        # Pair-wise distances for the first three samples of training_data
        scaled=FeatureScaling.data_scaling_minmax(self.training_data[0:3, :])
        scaled = scaled[0]
        scaled_x = scaled[:,:-1]
        
        distance_array = distance.cdist(scaled_x, scaled_x, 'euclidean')

        # Linear
        data_feed_01 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='linear')
        expected_output_1 = distance_array
        output_1 = data_feed_01.basis_generation(2)
        np.testing.assert_array_equal(expected_output_1, output_1)

        # Cubic
        data_feed_02 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='cubic')
        expected_output_2 = distance_array ** 3
        output_2 = data_feed_02.basis_generation(2)
        np.testing.assert_array_equal(expected_output_2, output_2)

        # # Spline
        data_feed_03 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='spline')
        expected_output_3 = np.nan_to_num(distance_array ** 2 * np.log(distance_array))
        output_3 = data_feed_03.basis_generation(2)
        np.testing.assert_array_equal(expected_output_3, output_3)

        # # Gaussian
        data_feed_04 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='gaussian')
        shape_value = 2
        expected_output_4 = np.exp(-1 * ((distance_array * shape_value) ** 2))
        output_4 = data_feed_04.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_4, output_4)

        # # Multiquadric
        data_feed_05 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='mq')
        shape_value = 2
        expected_output_5 = np.sqrt(((distance_array * shape_value) ** 2) + 1)
        output_5 = data_feed_05.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_5, output_5)

        # # Inverse multiquadric
        data_feed_06 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='imq')
        shape_value = 2
        expected_output_6 = 1 / np.sqrt(((distance_array * shape_value) ** 2) + 1)
        output_6 = data_feed_06.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_6, output_6)
        

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
        expected_value = 6613.875
        output_1 = RadialBasisFunctions.cost_function(theta, x_vector, y)
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
        theta = np.array([[4.5], [3], [3], [1], [1], [0]])
        expected_value = 90.625  # Calculated externally as sum(dy^2) / 2m
        output_1 = RadialBasisFunctions.cost_function(theta, x_vector, y)
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
        theta = np.array([[2], [2], [2], [1], [1], [0]])
        expected_value = 0
        output_1 = RadialBasisFunctions.cost_function(theta, x_vector, y)
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
            [[-97], [-635], [-635], [-5246.875], [-5246.875], [-3925]])
        expected_value = expected_value.reshape(expected_value.shape[0], )
        output_1 = RadialBasisFunctions.gradient_function(theta, x_vector, y)
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
        output_1 = RadialBasisFunctions.gradient_function(theta, x_vector, y)
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
        output_1 = RadialBasisFunctions.gradient_function(theta, x_vector, y)
        np.testing.assert_equal(output_1, expected_value)


    def test_bfgs_parameter_optimization_01(self):
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = x[:, ] ** 2
        x_vector[:, 1] = x[:, ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.], [2.], [1.]]).reshape(3, )
        data_feed = RadialBasisFunctions(self.test_data_numpy, basis_function='linear', solution_method='bfgs')
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == 'bfgs'
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    def test_bfgs_parameter_optimization_02(self):
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
        data_feed = RadialBasisFunctions(self.full_data, solution_method='bfgs')
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == 'bfgs'
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    
    def test_explicit_linear_algebra_solution_01(self):
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = x[:, ] ** 2
        x_vector[:, 1] = x[:, ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.], [2.], [1.]]).reshape(3, )
        output_1 = RadialBasisFunctions.explicit_linear_algebra_solution(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    def test_explicit_linear_algebra_solution_02(self):
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
        output_1 = RadialBasisFunctions.explicit_linear_algebra_solution(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))


    def test_pyomo_optimization_01(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = x[:, ] ** 2
        x_vector[:, 1] = x[:, ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.], [2.], [1.]])
        output_1 = RadialBasisFunctions.pyomo_optimization(x_vector, y)
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
        output_1 = RadialBasisFunctions.pyomo_optimization(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))


    def test_error_calculation_01(self):
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
        expected_value_1 = 2 * 6613.875  # Calculated externally as sum(y^2) / m
        expected_value_2 = expected_value_1 ** 0.5
        output_1, output_2, _ = RadialBasisFunctions.error_calculation(theta, x_vector, y)
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2

    def test_error_calculation_02(self):
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
        expected_value_1 = 2 * 90.625  # Calculated externally as sum(dy^2) / 2m
        expected_value_2 = expected_value_1 ** 0.5
        output_1, output_2, _ = RadialBasisFunctions.error_calculation(theta, x_vector, y)
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2

    def test_error_calculation_03(self):
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
        expected_value_1 = 2 * 0  # Value should return zero for exact solution
        expected_value_2 = expected_value_1 ** 0.5
        output_1, output_2, _ = RadialBasisFunctions.error_calculation(theta, x_vector, y)
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2


    def test_r2_calculation_01(self):
        y_actual = np.array([[1], [4], [9], [16], [25], [36], [49], [64], [81], [100]])
        y_pred = y_actual * 1.05
        expected_output = 0.993974359  # Evaluated in Excel
        output = RadialBasisFunctions.r2_calculation(y_actual, y_pred)
        assert round(abs(expected_output-output), 7) == 0
        

    def test_r2_calculation_02(self):
        y_actual = np.array([[1], [4], [9], [16], [25], [36], [49], [64], [81], [100]])
        y_pred = y_actual * 1.50
        expected_output = 0.3974358974  # Evaluated in Excel
        output = RadialBasisFunctions.r2_calculation(y_actual, y_pred)
        assert round(abs(expected_output-output), 7) == 0
        
    
    def mock_basis_generation(self, r):
        return np.ones((self.x_data.shape[0], self.x_data.shape[0]))
    def mock_optimization(self, x, y):
        return 500 * np.ones((x.shape[0], 1))
    @patch.object(RadialBasisFunctions, 'basis_generation', mock_basis_generation)
    @patch.object(RadialBasisFunctions, 'explicit_linear_algebra_solution', mock_optimization)
    def test_loo_error_estimation_with_rippa_method_01(self):
        reg_param = 0.1
        shape_factor = 1
        expected_x = np.ones((self.training_data.shape[0], self.training_data.shape[0])) + (reg_param * np.eye(self.training_data.shape[0], self.training_data.shape[0]))
        expected_inverse_x = np.diag(np.linalg.pinv(expected_x))
        expected_radial_weights = 500 * np.ones((self.training_data.shape[0], 1))
        expected_errors = np.linalg.norm(expected_radial_weights / (expected_inverse_x.reshape(expected_inverse_x.shape[0], 1)))

        data_feed = RadialBasisFunctions(self.training_data, solution_method='algebraic')
        _, output_1, output_2 = data_feed.loo_error_estimation_with_rippa_method(shape_factor, reg_param)
        assert output_1 == np.linalg.cond(expected_x)
        np.testing.assert_array_equal(output_2, expected_errors)
        print()

    @patch.object(RadialBasisFunctions, 'basis_generation', mock_basis_generation)
    @patch.object(RadialBasisFunctions, 'pyomo_optimization', mock_optimization)
    def test_loo_error_estimation_with_rippa_method_02(self):
        reg_param = 0.1
        shape_factor = 1
        expected_x = np.ones((self.training_data.shape[0], self.training_data.shape[0])) + (reg_param * np.eye(self.training_data.shape[0], self.training_data.shape[0]))
        expected_inverse_x = np.diag(np.linalg.pinv(expected_x))
        expected_radial_weights = 500 * np.ones((self.training_data.shape[0], 1))
        expected_errors = np.linalg.norm(expected_radial_weights / (expected_inverse_x.reshape(expected_inverse_x.shape[0], 1)))

        data_feed = RadialBasisFunctions(self.training_data, solution_method='pyomo')
        _, output_1, output_2 = data_feed.loo_error_estimation_with_rippa_method(shape_factor, reg_param)
        assert output_1 == np.linalg.cond(expected_x)
        np.testing.assert_array_equal(output_2, expected_errors)
        print()

    @patch.object(RadialBasisFunctions, 'basis_generation', mock_basis_generation)
    @patch.object(RadialBasisFunctions, 'bfgs_parameter_optimization', mock_optimization)
    def test_loo_error_estimation_with_rippa_method_03(self):
        reg_param = 0.1
        shape_factor = 1
        expected_x = np.ones((self.training_data.shape[0], self.training_data.shape[0])) + (reg_param * np.eye(self.training_data.shape[0], self.training_data.shape[0]))
        expected_inverse_x = np.diag(np.linalg.pinv(expected_x))
        expected_radial_weights = 500 * np.ones((self.training_data.shape[0], 1))
        expected_errors = np.linalg.norm(expected_radial_weights / (expected_inverse_x.reshape(expected_inverse_x.shape[0], 1)))

        data_feed = RadialBasisFunctions(self.training_data, solution_method='bfgs')
        _, output_1, output_2 = data_feed.loo_error_estimation_with_rippa_method(shape_factor, reg_param)
        assert output_1 == np.linalg.cond(expected_x)
        np.testing.assert_array_equal(output_2, expected_errors)


    def test_leave_one_out_crossvalidation_01(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function=None,solution_method=None, regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    def test_leave_one_out_crossvalidation_02(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='cubic',solution_method=None, regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
    
        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
    
    def test_leave_one_out_crossvalidation_03(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='linear' ,solution_method=None, regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        
        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
    
    def test_leave_one_out_crossvalidation_04(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='spline' ,solution_method=None, regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
        
    def test_leave_one_out_crossvalidation_05(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='gaussian',solution_method=None, regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
    
    def test_leave_one_out_crossvalidation_06(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='mq',solution_method=None, regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
    
    def test_leave_one_out_crossvalidation_07(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function= 'imq',solution_method=None, regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
    
    def test_leave_one_out_crossvalidation_08(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function=None,solution_method='algebraic', regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
    
    def test_leave_one_out_crossvalidation_09(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function=None,solution_method='BFGS', regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors
    
    def test_leave_one_out_crossvalidation_10(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function=None,solution_method='pyomo', regularization=False)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    def test_leave_one_out_crossvalidation_11(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function=None,solution_method=None, regularization=True)
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()

        if (data_feed.basis_function == 'gaussian') or (data_feed.basis_function == 'mq') or (data_feed.basis_function.lower() == 'imq'):
            r_set = [0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1.0, 2.0, 5.0, 7.5, 10.0, 20.0, 50.0, 75.0, 100.0, 200.0, 500.0, 1000.0]
        else:
            r_set = [0]

        if data_feed.regularization is True:
            reg_parameter = [0.00001, 0.00002, 0.00005, 0.000075, 0.0001, 0.0002, 0.0005, 0.00075, 0.001, 0.002, 0.005, 0.0075, 0.01, 0.02, 0.05, 0.075, 0.1, 0.2, 0.5, 0.75, 1]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _,_,expected_errors = data_feed.loo_error_estimation_with_rippa_method(r_best, lambda_best)

        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors


    def test_rbf_training_01(self):
        data_feed = RadialBasisFunctions(self.test_data_numpy,basis_function=None,solution_method='algebraic', regularization=False)
        results = data_feed.rbf_training()

        best_r_value, best_lambda_param, _ = data_feed.leave_one_out_crossvalidation()

        # Generate x matrix
        x_transformed = data_feed.basis_generation(best_r_value)
        x_transformed = x_transformed + (best_lambda_param * np.eye(x_transformed.shape[0], x_transformed.shape[1]))
        x_condition_number = np.linalg.cond(x_transformed)

        if data_feed.solution_method == 'algebraic':
            radial_weights = data_feed.explicit_linear_algebra_solution(x_transformed, data_feed.y_data)
        elif data_feed.solution_method == 'pyomo':
            radial_weights = data_feed.pyomo_optimization(x_transformed, data_feed.y_data)
        elif data_feed.solution_method == 'bfgs':
            radial_weights = data_feed.bfgs_parameter_optimization(x_transformed, data_feed.y_data)
        radial_weights = radial_weights.reshape(radial_weights.shape[0], 1)

        training_ss_error, rmse_error, y_training_predictions_scaled = data_feed.error_calculation(radial_weights, x_transformed, data_feed.y_data)
        r_square = data_feed.r2_calculation(data_feed.y_data, y_training_predictions_scaled)
        y_training_predictions = data_feed.data_min[0, -1] + y_training_predictions_scaled * (data_feed.data_max[0, -1] - data_feed.data_min[0, -1])
        
        np.testing.assert_array_equal(radial_weights, results.weights)
        np.testing.assert_array_equal(best_r_value, results.sigma)
        np.testing.assert_array_equal(best_lambda_param, results.regularization)
        np.testing.assert_array_equal(data_feed.centres, results.centres )
        np.testing.assert_array_equal(y_training_predictions, results.output_predictions )
        np.testing.assert_array_equal(rmse_error, results.rmse )
        np.testing.assert_array_equal(x_condition_number, results.condition_number )
        np.testing.assert_array_equal(data_feed.regularization, results.regularization )
        np.testing.assert_array_equal(r_square, results.R2 )
        assert data_feed.basis_function == results.basis_function
        np.testing.assert_array_equal(data_feed.data_min[:, :-1], results.x_data_min )
        np.testing.assert_array_equal(data_feed.data_max[:, :-1], results.x_data_max )
        np.testing.assert_array_equal(data_feed.data_min[:, -1], results.y_data_min )
        np.testing.assert_array_equal(data_feed.data_max[:, -1], results.y_data_max )
        assert results.solution_status == 'ok'
    
    def test_rbf_training_02(self):
        data_feed = RadialBasisFunctions(self.test_data_numpy,basis_function=None,solution_method='pyomo', regularization=False)
        results = data_feed.rbf_training()
        with pytest.warns(Warning):
            results = data_feed.rbf_training()          
            assert results.solution_status == 'unstable solution'
    
    def test_rbf_training_03(self):
        data_feed = RadialBasisFunctions(self.test_data_numpy,basis_function=None,solution_method='bfgs', regularization=False)
        results = data_feed.rbf_training()
        with pytest.warns(Warning):
            results = data_feed.rbf_training()          
            assert results.solution_status == 'unstable solution'


    def test_rbf_predict_output_01(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='linear', regularization=False)
        results = data_feed.rbf_training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.rbf_predict_output(results, x_test)
                
        data_minimum =results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum)/scale
        x_test = x_pred_scaled.reshape(x_test.shape)

        distance_vec = distance.cdist(x_test, results.centres, 'euclidean')
        expected_output = np.matmul(distance_vec, results.weights)
        expected_output = results.y_data_min + expected_output * (results.y_data_max - results.y_data_min)
        assert expected_output == output
        
    def test_rbf_predict_output_02(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='cubic', regularization=False)
        results = data_feed.rbf_training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.rbf_predict_output(results, x_test)
                
        data_minimum =results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum)/scale
        x_test = x_pred_scaled.reshape(x_test.shape)

        distance_vec = distance.cdist(x_test, results.centres, 'euclidean')
        expected_output = np.matmul(distance_vec**3, results.weights)
        expected_output = results.y_data_min + expected_output * (results.y_data_max - results.y_data_min)
        assert expected_output == output
        
    def test_rbf_predict_output_03(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='gaussian', regularization=False)
        results = data_feed.rbf_training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.rbf_predict_output(results, x_test)
                
        data_minimum =results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum)/scale
        x_test = x_pred_scaled.reshape(x_test.shape)

        distance_vec = distance.cdist(x_test, results.centres, 'euclidean')
        expected_output = np.matmul(np.exp(-1 * ((distance_vec * results.sigma) ** 2)),  results.weights)
        expected_output = results.y_data_min + expected_output * (results.y_data_max - results.y_data_min)
        assert expected_output == output
    
    def test_rbf_predict_output_04(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='imq', regularization=False)
        results = data_feed.rbf_training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.rbf_predict_output(results, x_test)
                
        data_minimum =results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum)/scale
        x_test = x_pred_scaled.reshape(x_test.shape)

        distance_vec = distance.cdist(x_test, results.centres, 'euclidean')
        expected_output = np.matmul(1 / np.sqrt(((distance_vec * results.sigma) ** 2) + 1),  results.weights)
        expected_output = results.y_data_min + expected_output * (results.y_data_max - results.y_data_min)
        assert expected_output == output
        
    def test_rbf_predict_output_05(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='mq', regularization=False)
        results = data_feed.rbf_training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.rbf_predict_output(results, x_test)
                
        data_minimum =results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum)/scale
        x_test = x_pred_scaled.reshape(x_test.shape)

        distance_vec = distance.cdist(x_test, results.centres, 'euclidean')
        expected_output = np.matmul(np.sqrt(((distance_vec * results.sigma) ** 2) + 1),  results.weights)
        expected_output = results.y_data_min + expected_output * (results.y_data_max - results.y_data_min)
        assert expected_output == output

    def test_rbf_predict_output_06(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='spline', regularization=False)
        results = data_feed.rbf_training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.rbf_predict_output(results, x_test)
                
        data_minimum =results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum)/scale
        x_test = x_pred_scaled.reshape(x_test.shape)

        distance_vec = distance.cdist(x_test, results.centres, 'euclidean')
        expected_output = np.matmul(np.nan_to_num(distance_vec ** 2 * np.log(distance_vec)),  results.weights)
        expected_output = results.y_data_min + expected_output * (results.y_data_max - results.y_data_min)
        assert expected_output == output


    def test_get_feature_vector_01(self):
        data_feed = RadialBasisFunctions(self.full_data, basis_function='linear')
        output = data_feed.get_feature_vector()
        expected_dict = {'x1': 0, 'x2': 0}
        assert expected_dict == output.extract_values()
        
    def test_get_feature_vector_02(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='linear')
        output = data_feed.get_feature_vector()
        expected_dict = {0: 0, 1: 0}
        assert expected_dict == output.extract_values()
        
        
    def test_rbf_generate_expression_01(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='linear',solution_method=None, regularization=False)
        p = data_feed.get_feature_vector()
        results = data_feed.rbf_training()       
        lv =[]
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.rbf_generate_expression((lv))
    
    def test_rbf_generate_expression_02(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='cubic',solution_method=None, regularization=False)
        p = data_feed.get_feature_vector()
        results = data_feed.rbf_training()       
        lv =[]
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.rbf_generate_expression((lv))
    
    def test_rbf_generate_expression_03(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='gaussian',solution_method=None, regularization=False)
        p = data_feed.get_feature_vector()
        results = data_feed.rbf_training()       
        lv =[]
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.rbf_generate_expression((lv))
    
    def test_rbf_generate_expression_04(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='mq',solution_method=None, regularization=False)
        p = data_feed.get_feature_vector()
        results = data_feed.rbf_training()       
        lv =[]
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.rbf_generate_expression((lv))
    
    def test_rbf_generate_expression_05(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='imq',solution_method=None, regularization=False)
        p = data_feed.get_feature_vector()
        results = data_feed.rbf_training()       
        lv =[]
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.rbf_generate_expression((lv))
    
    def test_rbf_generate_expression_06(self):
        data_feed = RadialBasisFunctions(self.training_data,basis_function='spline',solution_method=None, regularization=False)
        p = data_feed.get_feature_vector()
        results = data_feed.rbf_training()       
        lv =[]
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.rbf_generate_expression((lv))
    
    
if __name__ == '__main__':
    unittest.main()
