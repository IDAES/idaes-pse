from pysmo.radial_basis_function import RadialBasisFunctions, FeatureScaling
import numpy as np
import pandas as pd
import pyutilib.th as unittest
from unittest.mock import patch


class FeatureScalingTestCases(unittest.TestCase):
    def setUp(self):
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})
        return self.test_data

    """
    Tests 1 - 7: Tests to check the range of potential inputs for the number of samples:
     - 01: Test behaviour when input is a pandas dataframe
     - 02: Test behaviour when input is a numpy array and with 1D array. Returns a 2D N x 1 array.
     - 03: Test behaviour input is neither of the acceptable types - should throw up a
     - 04: Test behaviour when number_of_samples > acceptable range. Should raise exception.
     - 05: Test behaviour when number_of_samples = None. Should complete successfully with default number_of_samples
     - 06: Test behaviour when number_of_samples = fractional. Should raise exception.
     - 07: Test that data loading behaves as expected when a numpy array is provided.
     - 08: Test that a ValueError is thrown when the input is not a numpy array or Pandas DataFrame
    Tests 3 and 5 also check x_data attribute is returned correctly.
    """

    def test_data_scaling_01(self):
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(self.test_data)
        expected_output_2 = np.array([[0, 1]])
        expected_output_3 = np.array([[9, 100]])
        expected_output_1 = (self.test_data - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    def test_data_scaling_02(self):
        test_data_as_array = self.test_data.values
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(test_data_as_array[:, 0])
        expected_output_2 = np.array([[0]])
        expected_output_3 = np.array([[9]])
        expected_output_1 = (test_data_as_array[:, 0] - expected_output_2) / (expected_output_3 - expected_output_2)
        expected_output_1 = expected_output_1.reshape(expected_output_1.shape[1], 1)  # Because shape is changed backend
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    def test_data_scaling_03(self):
        test_data_as_list = self.test_data.values.tolist()
        with self.assertRaises(TypeError):
            FeatureScaling.data_scaling_minmax(test_data_as_list)

    def test_data_unscaling_01(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1, 6]])
        max_array = np.array([[5, 10]])
        output = FeatureScaling.data_unscaling_minmax(data_array, min_array, max_array)
        expected_output = np.array([[1, 10], [2, 9], [3, 8], [4, 7], [5, 6]])
        np.testing.assert_array_equal(output, expected_output)

    def test_data_unscaling_02(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with self.assertRaises(IndexError):
            FeatureScaling.data_unscaling_minmax(data_array, min_array, max_array)

    def test_data_unscaling_03(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1, 6, 10]])
        max_array = np.array([[5, 10, 15]])
        with self.assertRaises(IndexError):
            FeatureScaling.data_unscaling_minmax(data_array, min_array, max_array)

    def test_data_unscaling_04(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1, 6]])
        max_array = np.array([[5, 10]])
        output = FeatureScaling.data_unscaling_minmax(data_array[:, 0], min_array[:, 0], max_array[:, 0])
        expected_output = np.array([[1], [2], [3], [4], [5]])
        np.testing.assert_array_equal(output, expected_output)


class RadialBasisFunctionTestCases(unittest.TestCase):

    def setUp(self):
        # Data generated from the expression (x_1 + 1)^2 + (x_2 + 1) ^ 2 between 0 and 10 for x_1 and x_2
        x1 = np.linspace(0, 10, 21)
        x2 = np.linspace(0, 10, 21)
        y = np.zeros((len(x1) * len(x2), 3))
        counter = 0
        for i in x1:
            for j in x2:
                y[counter, 0] = i
                y[counter, 1] = j
                y[counter, 2] = ((i + 1) ** 2) + ((j + 1) ** 2)
                counter = counter + 1
        self.full_data = pd.DataFrame({'x1': y[:, 0], 'x2': y[:, 1], 'y': y[:, 2]})

        # Theoretical sample data from range for(x_1 + 1)^2 + (x_2 + 1) ^ 2 between 0 and 10 for x_1 and x_2
        x1 = np.linspace(0, 10, 5)
        x2 = np.linspace(0, 10, 5)
        training_samples = np.zeros((len(x1) * len(x2), 3))
        counter = 0
        for i in x1:
            for j in x2:
                training_samples[counter, 0] = i
                training_samples[counter, 1] = j
                training_samples[counter, 2] = ((i + 1) ** 2) + ((j + 1) ** 2)
                counter = counter + 1
        self.training_data = training_samples

        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = input_array  # pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})

        return self.full_data, self.training_data, self.test_data

    def mock_basis_generation(self, r):
        return np.ones((self.x_data.shape[0], self.x_data.shape[0]))

    def mock_optimization(self, x, y):
        return 500 * np.ones((x.shape[0], 1))

    """
    Tests 01 - 03 : Demonstration of behaviour for different data input types
     - 01: Test behaviour when the input is a dataframe. We check that all attributes related to the data are properly loaded.
     - 02: Test behaviour when the input is a numpy array. We check that all attributes related to the data are properly loaded.
     - 03: Test behaviour when the input is not a numpy array or dataframe. A ValueError is expected to be raised.
    """

    def test_initialization_rbfs_01(self):
        data_feed = RadialBasisFunctions(self.full_data, 'linear')
        expected_output_1 = self.full_data.values[:, :-1]
        expected_output_2 = self.full_data.values[:, -1].reshape(expected_output_1.shape[0], 1)
        expected_output_3 = expected_output_1
        expected_output_4 = list(self.full_data.columns)[:-1]
        np.testing.assert_array_equal(expected_output_1, data_feed.x_data)
        np.testing.assert_array_equal(expected_output_2, data_feed.y_data)
        np.testing.assert_array_equal(expected_output_3, data_feed.centres)
        self.assertListEqual(expected_output_4, data_feed.x_data_columns)

    def test_initialization_rbfs_02(self):
        data_feed = RadialBasisFunctions(self.training_data, 'linear')
        expected_output_1 = self.training_data[:, :-1]
        expected_output_2 = self.training_data[:, -1].reshape(expected_output_1.shape[0], 1)
        expected_output_3 = expected_output_1
        expected_output_4 = list(range(self.training_data.shape[1] - 1))
        np.testing.assert_array_equal(expected_output_1, data_feed.x_data)
        np.testing.assert_array_equal(expected_output_2, data_feed.y_data)
        np.testing.assert_array_equal(expected_output_3, data_feed.centres)
        self.assertListEqual(expected_output_4, data_feed.x_data_columns)

    def test_initialization_rbfs_03(self):
        sample_points = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        with self.assertRaises(ValueError):
            RadialBasisFunctions(self.training_data.tolist(), basis_function='linear')

    """
    Tests 04-08 : Tests to class attribute solution_method.
     - 04: Test behaviour when solution_method input is not a string. Should raise exception.
     - 05: Test behaviour when solution_method = 'alGebraIC'. Attribute should be set to 'algebraic' irrespective of case.
     - 06: Test behaviour when solve_method = 'pYomO'. Attribute should be set to 'pyomo' irrespective of case.
     - 07: Test behaviour when solve_method = 'bFGs'. Attribute should be set to 'bfgs' irrespective of case.
     - 08: Test behaviour when solve_method is an invalid string. Should raise exception.
     - 09: Test behaviour when solve_method is None. Attribute should be set to default ('algebraic').
    """

    def test_initialization_rbfs_04(self):
        solve_method = 1
        with self.assertRaises(Exception):
            RadialBasisFunctions(self.test_data, solution_method=solve_method)

    def test_initialization_rbfs_05(self):
        solve_method = 'aLgeBraiC'
        output_1 = RadialBasisFunctions(self.test_data, solution_method=solve_method)
        self.assertEqual(output_1.solution_method, solve_method.lower())

    def test_initialization_rbfs_06(self):
        solve_method = 'pYomO'
        output_1 = RadialBasisFunctions(self.test_data, solution_method=solve_method)
        self.assertEqual(output_1.solution_method, solve_method.lower())

    def test_initialization_rbfs_07(self):
        solve_method = 'bFGs'
        output_1 = RadialBasisFunctions(self.test_data, solution_method=solve_method)
        self.assertEqual(output_1.solution_method, solve_method.lower())

    def test_initialization_rbfs_08(self):
        solve_method = 'mle'
        with self.assertRaises(Exception):
            RadialBasisFunctions(self.test_data, solution_method=solve_method)

    def test_initialization_rbfs_09(self):
        default_solve = 'algebraic'
        output = RadialBasisFunctions(self.test_data)
        self.assertEqual(output.solution_method, default_solve)

    """
    Tests 10-18 : Tests to class attribute basis_function.
     - 10: Test behaviour when basis_function = 'lineaR'. Attribute should be set to 'linear' irrespective of case.
     - 11: Test behaviour when basis_function = 'splinE'. Attribute should be set to 'spline' irrespective of case.
     - 12: Test behaviour when basis_function = 'cubiC'. Attribute should be set to 'cubic' irrespective of case.
     - 13: Test behaviour when basis_function = 'gaussiaN'. Attribute should be set to 'gaussian' irrespective of case.
     - 14: Test behaviour when basis_function = 'imQ'. Attribute should be set to 'imq' irrespective of case.
     - 15: Test behaviour when basis_function = 'mQ'. Attribute should be set to 'mq' irrespective of case.
     - 16: Test behaviour when basis_function is an invalid string. Should raise exception.
     - 17: Test behaviour when basis_function is None. Attribute should be set to default ('gaussian').
     - 18: Test behaviour when basis_function input is not a string. Should raise exception.
    """

    def test_initialization_rbfs_10(self):
        basis = 'lineaR'
        output = RadialBasisFunctions(self.test_data, basis_function=basis)
        self.assertEqual(output.basis_function, basis.lower())

    def test_initialization_rbfs_11(self):
        basis = 'splinE'
        output = RadialBasisFunctions(self.test_data, basis_function=basis)
        self.assertEqual(output.basis_function, basis.lower())

    def test_initialization_rbfs_12(self):
        basis = 'cubiC'
        output = RadialBasisFunctions(self.test_data, basis_function=basis)
        self.assertEqual(output.basis_function, basis.lower())

    def test_initialization_rbfs_13(self):
        basis = 'gaussiaN'
        output = RadialBasisFunctions(self.test_data, basis_function=basis)
        self.assertEqual(output.basis_function, basis.lower())

    def test_initialization_rbfs_14(self):
        basis = 'imQ'
        output = RadialBasisFunctions(self.test_data, basis_function=basis)
        self.assertEqual(output.basis_function, basis.lower())

    def test_initialization_rbfs_15(self):
        basis = 'mQ'
        output = RadialBasisFunctions(self.test_data, basis_function=basis)
        self.assertEqual(output.basis_function, basis.lower())

    def test_initialization_rbfs_16(self):
        basis = 'pyomo'
        with self.assertRaises(Exception):
            RadialBasisFunctions(self.test_data, basis_function=basis)

    def test_initialization_rbfs_17(self):
        default_basis = 'gaussian'
        output = RadialBasisFunctions(self.test_data)
        self.assertEqual(output.basis_function, default_basis)

    def test_initialization_rbfs_18(self):
        basis = 1
        with self.assertRaises(Exception):
            RadialBasisFunctions(self.test_data, basis_function=basis)

    """
    Tests 10-18 : Tests to class attribute self.regularization.
     - 19: Test behaviour when regularization input = 'True'. self.regularization attribute should be set to 'True'.
     - 20: Test behaviour when regularization input = 'False'. self.regularization attribute should be set to 'False'.
     - 21: Test behaviour when regularization input is None. self.regularization attribute should be set to the default value of 'True'.
     - 22: Test behaviour when a non-boolean string is supplied as regularization input. An exception should be raised.
     - 23: Test behaviour when a non-string is supplied as regularization input. An exception should be raised.
    """

    def test_initialization_rbfs_19(self):
        reg = True
        output = RadialBasisFunctions(self.test_data, regularization=reg)
        self.assertEqual(output.regularization, reg)

    def test_initialization_rbfs_20(self):
        reg = False
        output = RadialBasisFunctions(self.test_data, regularization=reg)
        self.assertEqual(output.regularization, reg)

    def test_initialization_rbfs_21(self):
        default_reg = True
        output = RadialBasisFunctions(self.test_data)
        self.assertEqual(output.regularization, default_reg)

    def test_initialization_rbfs_22(self):
        reg = 'gaussian'
        with self.assertRaises(Exception):
            RadialBasisFunctions(self.test_data, regularization=reg)

    def test_initialization_rbfs_23(self):
        reg = 5
        with self.assertRaises(Exception):
            RadialBasisFunctions(self.test_data, regularization=reg)

    """
    Tests: : Unit tests for eucl_distance, a function that evaluates the distance between a point c and a set of points in a numpy array.
        The test checks that the function is able to calculate the distance from a single point (n x 1 row vector) to a set of design points (supplied in an n x m array)
    """

    def test_r2_distance(self):
        u = np.array([[3.6, 7.2]])
        data_feed = RadialBasisFunctions(self.training_data[0:10, :], basis_function='linear')
        output = data_feed.r2_distance(u)
        expected_output = np.array([8.049844719, 5.9203040462, 4.2190046219, 3.6124783736, 4.5607017004, 7.2835430939, 4.8270073545, 2.4596747752, 1.1401754251, 3.0083217913])
        np.testing.assert_almost_equal(expected_output, output, decimal=6)

    """
    Tests: : The next six tests are for the basis transformations. 
    They test the results returned by each transformation for a wide range of values.
    """

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

    def test_linear_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        output_1 = RadialBasisFunctions.linear_transformation(d_vec)
        np.testing.assert_array_equal(d_vec, output_1)

    def test_cubic_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        expected_output = d_vec ** 3
        output = RadialBasisFunctions.cubic_transformation(d_vec)
        np.testing.assert_array_equal(expected_output, output)

    def test_thin_plate_spline_transformation(self):
        d_vec = np.array([[0, 0], [5e-6, 7e-6], [0.005, 0.007], [0.05, 0.07], [0.5, 0.7], [5, 7], [50, 70]])
        expected_output = np.nan_to_num(d_vec ** 2 * np.log(d_vec))
        output = RadialBasisFunctions.thin_plate_spline_transformation(d_vec)
        np.testing.assert_array_equal(expected_output, output)

    """
    Tests: : Tests for the function basis_generation.
    The test checks that for a small dataset, the expected result is returned by the function for each basis transformation type.
    """

    def test_basis_generation(self):

        # Pair-wise distances for the first three samples of training_data
        distance_array = np.array([[0, 2.5, 5], [2.5, 0, 2.5], [5, 2.5, 0]])

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

        # Spline
        data_feed_03 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='spline')
        expected_output_3 = np.nan_to_num(distance_array ** 2 * np.log(distance_array))
        output_3 = data_feed_03.basis_generation(2)
        np.testing.assert_array_equal(expected_output_3, output_3)

        # Gaussian
        data_feed_04 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='gaussian')
        shape_value = 2
        expected_output_4 = np.exp(-1 * ((distance_array * shape_value) ** 2))
        output_4 = data_feed_04.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_4, output_4)

        # Multiquadric
        data_feed_05 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='mq')
        shape_value = 2
        expected_output_5 = np.sqrt(((distance_array * shape_value) ** 2) + 1)
        output_5 = data_feed_05.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_5, output_5)

        # Inverse multiquadric
        data_feed_06 = RadialBasisFunctions(self.training_data[0:3, :], basis_function='imq')
        shape_value = 2
        expected_output_6 = 1 / np.sqrt(((distance_array * shape_value) ** 2) + 1)
        output_6 = data_feed_06.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_6, output_6)
        print()

    """
    Tests: : Unit tests for cost_function.  The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
        Three demonstration tests are done:
        The first test checks that the correct loss value is returned when theta is initialized as a vector of zeros.
        The second test checks that the correct loss value is returned for a random point within the solution space.
        The third test checks that the loss value is zero when the exact solution is found.
        - 
    """

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
        self.assertEqual(output_1, expected_value)

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
        self.assertEqual(output_1, expected_value)

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
        self.assertEqual(output_1, expected_value)

    """
    Tests: : Unit tests for gradient_function. The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
        Three demonstration tests are done:
        The first test checks that the correct gradients returned when theta is initialized as a vector of zeros.
        The second test checks that the correct gradients are returned for a random point within the solution space.
        The third test checks that the gradients are zero when the exact solution is found. 
    """

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

    """
    Tests: : Unit tests for bfgs_parameter_optimization function. Two unit tests are done:
        The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
        The second test is used to check that higher order multi-variable problems are solved correctly.
        - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
        - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.
    """

    def test_bfgs_parameter_optimization_01(self):
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = x[:, ] ** 2
        x_vector[:, 1] = x[:, ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.], [2.], [1.]]).reshape(3, )
        data_feed = RadialBasisFunctions(self.test_data, basis_function='linear', solution_method='bfgs')
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        self.assertEqual(data_feed.solution_method, 'bfgs')
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
        self.assertEqual(data_feed.solution_method, 'bfgs')
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    """
    Tests: : Unit tests for pyomo_optimization function. As with the other parameter estimation methods, two unit tests are done:
        The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
        The second test is used to check that higher order multi-variable problems are solved correctly.
        - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
        - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.
    """

    @staticmethod
    def test_pyomo_optimization_01():
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

    """
    Tests: : Unit tests for explicit_linear_algebra_solution function. The same two unit tests are done:
        The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
        The second test is used to check that higher order multi-variable problems are solved correctly.
        - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
        - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.
    """

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

    """
    Tests: : Unit tests for error_calculation function.  The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
        The three tests used for the cost_function are repeated for the RMSE and SSE errors. 
         - The first test checks that the correct cv error is returned when theta is initialized as a vector of zeros.
         - The second test checks that the correct cv error is returned for a random point within the solution space.
         - The third test checks that the cv error is zero when the exact solution is found.
    In each case, the result is expected to be twice as large as the cost function's expected value: the cost function is half of the total squared error.
    """

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
        self.assertEqual(output_1, expected_value_1)
        self.assertEqual(output_2, expected_value_2)

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
        self.assertEqual(output_1, expected_value_1)
        self.assertEqual(output_2, expected_value_2)

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
        self.assertEqual(output_1, expected_value_1)
        self.assertEqual(output_2, expected_value_2)

    """
    Tests: : Unit tests for r2_calculation function.  The function is tested for two cases where the predicted values 
    deviate from the actual by 5% and 50% respectively. The R2-value in each case is compared against Excel evaluated values.
    """

    def test_r2_calculation_01(self):
        y_actual = np.array([[1], [4], [9], [16], [25], [36], [49], [64], [81], [100]])
        y_pred = y_actual * 1.05
        expected_output = 0.993974359  # Evaluated in Excel
        output = RadialBasisFunctions.r2_calculation(y_actual, y_pred)
        self.assertAlmostEqual(expected_output, output, places=7)
        print()

    def test_r2_calculation_02(self):
        y_actual = np.array([[1], [4], [9], [16], [25], [36], [49], [64], [81], [100]])
        y_pred = y_actual * 1.50
        expected_output = 0.3974358974  # Evaluated in Excel
        output = RadialBasisFunctions.r2_calculation(y_actual, y_pred)
        self.assertAlmostEqual(expected_output, output, places=7)
        print()

    """
       Tests: : Unit tests for get_feature_vector. We verify that:
           - 1: The (key, val) dictionary obtained from the generated IndexParam matches the headers of the data when the input is a dataframe
           - 2: The (key, val) dictionary obtained from the generated IndexParam matches is numerically consistent with that of the input array when a numpy array is supplied. 
    """

    def test_get_feature_vector_01(self):
        data_feed = RadialBasisFunctions(self.full_data, basis_function='linear')
        output = data_feed.get_feature_vector()
        expected_dict = {'x1': 0, 'x2': 0}
        self.assertDictEqual(expected_dict, output.extract_values())
        print()

    def test_get_feature_vector_02(self):
        data_feed = RadialBasisFunctions(self.training_data, basis_function='linear')
        output = data_feed.get_feature_vector()
        expected_dict = {0: 0, 1: 0}
        self.assertDictEqual(expected_dict, output.extract_values())
        print()

    """
       Tests: : Unit tests for rbf_predict_output, a function that generates predictions for a test dataset. 
       For a small dataset, we verify that the current predictions are produced for each transformation type.
    """

    def test_rbf_predict_output_01(self):
        data_feed = RadialBasisFunctions(self.training_data[0:3, :], basis_function='linear')
        x_test = np.array([[0, 7.5]])
        opt_wts = np.array([[2], [1], [3]])
        lambda_reg = 0
        reg = 1
        distance_vec = np.array([[7.5, 5, 2.5]])
        expected_output = np.matmul(distance_vec, opt_wts)
        output = data_feed.rbf_predict_output(opt_wts, x_test, data_feed.centres, reg, lambda_reg)
        self.assertEqual(expected_output, output)
        print()

    def test_rbf_predict_output_02(self):
        data_feed = RadialBasisFunctions(self.training_data[0:3, :], basis_function='cubic')
        x_test = np.array([[0, 7.5]])
        opt_wts = np.array([[2], [1], [3]])
        lambda_reg = 0
        reg = 1
        distance_vec = np.array([[7.5, 5, 2.5]])
        expected_output = np.matmul(distance_vec**3, opt_wts)
        output = data_feed.rbf_predict_output(opt_wts, x_test, data_feed.centres, reg, lambda_reg)
        self.assertEqual(expected_output, output)
        print()

    def test_rbf_predict_output_03(self):
        data_feed = RadialBasisFunctions(self.training_data[0:3, :], basis_function='gaussian')
        x_test = np.array([[0, 7.5]])
        opt_wts = np.array([[2], [1], [3]])
        lambda_reg = 0
        reg = 1
        distance_vec = np.array([[7.5, 5, 2.5]])
        expected_output = np.matmul(np.exp(-(distance_vec**reg)**2), opt_wts)
        output = data_feed.rbf_predict_output(opt_wts, x_test, data_feed.centres, reg, lambda_reg)
        self.assertEqual(expected_output, output)
        print()

    def test_rbf_predict_output_04(self):
        data_feed = RadialBasisFunctions(self.training_data[0:3, :], basis_function='imq')
        x_test = np.array([[0, 7.5]])
        opt_wts = np.array([[2], [1], [3]])
        lambda_reg = 0
        reg = 1
        distance_vec = np.array([[7.5, 5, 2.5]])
        expected_output = np.matmul(1 / np.sqrt(((distance_vec * reg) ** 2) + 1), opt_wts)
        output = data_feed.rbf_predict_output(opt_wts, x_test, data_feed.centres, reg, lambda_reg)
        self.assertEqual(expected_output, output)
        print()

    def test_rbf_predict_output_05(self):
        data_feed = RadialBasisFunctions(self.training_data[0:3, :], basis_function='mq')
        x_test = np.array([[0, 7.5]])
        opt_wts = np.array([[2], [1], [3]])
        lambda_reg = 0
        reg = 1
        distance_vec = np.array([[7.5, 5, 2.5]])
        expected_output = np.matmul(np.sqrt(((distance_vec * reg) ** 2) + 1), opt_wts)
        output = data_feed.rbf_predict_output(opt_wts, x_test, data_feed.centres, reg, lambda_reg)
        self.assertEqual(expected_output, output)
        print()

    def test_rbf_predict_output_06(self):
        data_feed = RadialBasisFunctions(self.training_data[0:3, :], basis_function='spline')
        x_test = np.array([[0, 7.5]])
        opt_wts = np.array([[2], [1], [3]])
        lambda_reg = 0
        reg = 1
        distance_vec = np.array([[7.5, 5, 2.5]])
        expected_output = np.matmul(np.nan_to_num(distance_vec ** 2 * np.log(distance_vec)), opt_wts)
        output = data_feed.rbf_predict_output(opt_wts, x_test, data_feed.centres, reg, lambda_reg)
        self.assertEqual(expected_output, output)
        print()

    """
    Tests: : Unit tests for loo_error_estimation_with_rippa_method. 
        - Three tests are run to check that the correct condition number and errors are returned by the method irrespective of solution_method:
            - First test checks the algebraic method
            - Second test checks the pyomo method
            - Third test checks the bfgs solution
        - In each test, the basis_generation method and relevant optimization method are replaced with mocks. The results from the mocks
          are used to compute the condition number and LOOCV errors.
    """

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
        self.assertEqual(output_1, np.linalg.cond(expected_x))
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
        self.assertEqual(output_1, np.linalg.cond(expected_x))
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
        self.assertEqual(output_1, np.linalg.cond(expected_x))
        np.testing.assert_array_equal(output_2, expected_errors)
        print()


if __name__ == '__main__':
    unittest.main()
