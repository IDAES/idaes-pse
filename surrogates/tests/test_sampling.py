from pysmo.sampling import LatinHypercubeSampling, UniformSampling, HaltonSampling, HammersleySampling, CVTSampling, SamplingMethods, FeatureScaling
import numpy as np
import pandas as pd
import pyutilib.th as unittest


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


class LatinHypercubeSamplingTestCases(unittest.TestCase):

    def setUp(self):
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})
        return self.test_data

    """
    Tests 1 - 7: Tests to check the range of potential inputs for the number of samples:
     - 01: Test behaviour when number_of_samples = -ve. Should raise exception.
     - 02: Test behaviour when number_of_samples = 0. Should raise exception.
     - 03: Test behaviour when number_of_samples = + integer within acceptable range. Should complete successfully.
     - 04: Test behaviour when number_of_samples > acceptable range. Should raise exception.
     - 05: Test behaviour when number_of_samples = None. Should complete successfully with default number_of_samples
     - 06: Test behaviour when number_of_samples = fractional. Should raise exception.
     - 07: Test that data loading behaves as expected when a numpy array is provided.
     - 08: Test that a ValueError is thrown when the input is not a numpy array or Pandas DataFrame
    Tests 3 and 5 also check x_data attribute is returned correctly.
    """

    def test_initialization_01(self):
        with self.assertRaises(Exception):
            LatinHypercubeSampling(self.test_data, -1)

    def test_initialization_02(self):
        with self.assertRaises(Exception):
            LatinHypercubeSampling(self.test_data, 0)

    def test_initialization_03(self):
        output_1 = LatinHypercubeSampling(self.test_data, 7)
        expected_output_1 = 7
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_1.number_of_samples, expected_output_1)
        self.assertEqual(output_1.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_04(self):
        with self.assertRaises(Exception):
            LatinHypercubeSampling(self.test_data, 11)

    def test_initialization_05(self):
        output_2 = LatinHypercubeSampling(self.test_data)
        expected_output = 5
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_2.number_of_samples, expected_output)
        self.assertEqual(output_2.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_06(self):
        with self.assertRaises(Exception):
            LatinHypercubeSampling(self.test_data, 1.7)

    def test_initialization_07(self):
        test_data_as_array = self.test_data.values
        output_1 = LatinHypercubeSampling(test_data_as_array, 7)
        np.testing.assert_array_equal(output_1.data, test_data_as_array)

    def test_initialization_08(self):
        test_data_as_list = self.test_data.values.tolist()
        with self.assertRaises(ValueError):
            LatinHypercubeSampling(test_data_as_list, 7)

    """
    Test: The random_shuffling function shuffles the columns of arrays. The next test simply checks that the array columns
    have been shuffled by comparing the lists and sets generated from each column before and after shuffling.

    The sets are expected to be the same, but the lists are expected to be different due to the re-arrangement.
    """

    def test_random_shuffling(self):
        x_pure = np.array(
            [[0, 1, 2, 3], [3, 4, 5, 6], [6, 7, 8, 9], [9, 10, 11, 12], [12, 13, 14, 15], [15, 16, 17, 18]])
        x = np.array([[0, 1, 2, 3], [3, 4, 5, 6], [6, 7, 8, 9], [9, 10, 11, 12], [12, 13, 14, 15], [15, 16, 17, 18]])
        output_1 = LatinHypercubeSampling.random_shuffling(x)
        # Test 1: Check that randomization of columns has been done by comparing individual columns before and after
        self.assertNotEqual(output_1[:, 0].tolist(), x_pure[:, 0].tolist())
        self.assertNotEqual(output_1[:, 1].tolist(), x_pure[:, 1].tolist())
        self.assertNotEqual(output_1[:, 2].tolist(), x_pure[:, 2].tolist())
        self.assertNotEqual(output_1[:, 3].tolist(), x_pure[:, 3].tolist())
        # Test 2: Check that the columns still contain exactly the same values as before randomization by comparing the sets
        self.assertEqual(set(output_1[:, 0]), set(x_pure[:, 0]))
        self.assertEqual(set(output_1[:, 1]), set(x_pure[:, 1]))
        self.assertEqual(set(output_1[:, 2]), set(x_pure[:, 2]))
        self.assertEqual(set(output_1[:, 3]), set(x_pure[:, 3]))

    """
    Test: Tests for lhs_points_generation: the function(s) that conduct random number generation for the input data.
     - The first test checks that its performance for the minimum number of sample points (1).
       The shape and returned value are checked against expected values.

     - The second test checks its performance for a typical value of number_of_samples. The first and last values in the
       returned array are verified to be within the expected ranges. The shapes are also checked to be correct.

     - The third test checks its performance for the upper bound - the number of input data points.
    """

    def test_lhs_points_generation_01(self):
        number_of_samples = 1
        data_feed = LatinHypercubeSampling(self.test_data, number_of_samples)  # Three samples requested.
        output_1 = data_feed.lhs_points_generation()
        # Test set 1: Check that the output has the right shape
        self.assertEqual(output_1.shape[0], number_of_samples)
        self.assertEqual(output_1.shape[1], (self.test_data.shape[1] - 1))
        # Test set 2: Check that the values returned are within the expected range
        expected_lower = 0  # Minimum value for x in input_array
        expected_upper = 1  # Maximum value for x in input_array
        self.assertTrue((output_1 > expected_lower) and (output_1 < expected_upper))

    def test_lhs_points_generation_02(self):
        number_of_samples = 3
        data_feed = LatinHypercubeSampling(self.test_data, number_of_samples)  # Three samples requested.
        output_1 = data_feed.lhs_points_generation()
        # Test set 1: Check that the output has the right shape
        self.assertEqual(output_1.shape[0], number_of_samples)
        self.assertEqual(output_1.shape[1], (self.test_data.shape[1] - 1))
        # Test set 2: Check that the values returned are within the expected range.
        expected_lower_firstvalue = 0
        expected_upper_lastvalue = 1
        expected_upper_firstvalue = expected_upper_lastvalue * 1 / number_of_samples
        expected_lower_lastvalue = expected_upper_lastvalue * (number_of_samples - 1) / number_of_samples
        self.assertTrue((output_1[0, 0] > expected_lower_firstvalue) and (output_1[0, 0] < expected_upper_firstvalue))
        self.assertTrue((output_1[number_of_samples - 1, 0] > expected_lower_lastvalue) and (
                output_1[number_of_samples - 1, 0] < expected_upper_lastvalue))

    def test_lhs_points_generation_03(self):
        number_of_samples = 10
        data_feed = LatinHypercubeSampling(self.test_data, number_of_samples)  # Three samples requested.
        output_1 = data_feed.lhs_points_generation()
        # Test set 1: Check that the output has the right shape
        self.assertEqual(output_1.shape[0], number_of_samples)
        self.assertEqual(output_1.shape[1], (self.test_data.shape[1] - 1))
        # Test set 2: Check that the values returned are within the expected range.
        expected_lower_firstvalue = 0
        expected_upper_lastvalue = 1
        expected_upper_firstvalue = expected_upper_lastvalue * 1 / number_of_samples
        expected_lower_lastvalue = expected_upper_lastvalue * (number_of_samples - 1) / number_of_samples
        self.assertTrue((output_1[0, 0] > expected_lower_firstvalue) and (output_1[0, 0] < expected_upper_firstvalue))
        self.assertTrue((output_1[number_of_samples - 1, 0] > expected_lower_lastvalue) and (
                output_1[number_of_samples - 1, 0] < expected_upper_lastvalue))


class UniformSamplingTestCases(unittest.TestCase):

    def setUp(self):
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})
        return self.test_data

    """
    Tests 1 - 12: Tests to check the range of potential inputs for UniformSampling class:
     - 01: Test behaviour when list_of_samples_per_variable is of the wrong type (not a list). Code should raise TypeError.
     - 02: Test behaviour when list_of_samples_per_variable is of the wrong length/dimension. Code should raise ValueError.
     - 03: Test behaviour when any entry in list_of_samples_per_variable < 2. Code should raise ValueError.
     - 04: Test behaviour when list_of_samples_per_variable is not an integer. Code should raise TypeError.
     - 05: Test behaviour when the total number  of samples to be generated(product of entries in list_of_samples_per_variable)
           is greater than the number of samples in the input data. Code should raise Exception.
     - 06: Test that the function behaves as expected and all attributes are properly loaded when all inputs are okay. 
     - 07: Test behaviour when 'edges' entry is non-boolean. Code should raise Exception.
     - 08: Test behaviour when no input is supplied for 'edges'. Default value of 'True' should be loaded.
     - 09: Ensure that 'self.edge' is loaded correctly when Boolean value 'True' is supplied.
     - 10: Ensure that 'self.edge' is loaded correctly when Boolean value 'False' is supplied.
     - 11: Test that the function behaves as expected and all attributes are properly loaded when data_input is an array.
     - 12: Test that a ValueError is thrown when the input is not a numpy array or Pandas DataFrame     
    """

    def test_initialization_01(self):
        with self.assertRaises(TypeError):
            UniformSampling(self.test_data, np.array([2]))

    def test_initialization_02(self):
        with self.assertRaises(ValueError):
            UniformSampling(self.test_data, [2, 2])

    def test_initialization_03(self):
        with self.assertRaises(ValueError):
            UniformSampling(self.test_data, [1])

    def test_initialization_04(self):
        with self.assertRaises(TypeError):
            UniformSampling(self.test_data, [2.3])

    def test_initialization_05(self):
        with self.assertRaises(Exception):
            UniformSampling(self.test_data, [11])

    def test_initialization_06(self):
        output_2 = UniformSampling(self.test_data, [3])
        expected_output = 3
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        expected_output_3 = [3]
        expected_output_4 = self.test_data.columns.values.tolist()
        self.assertEqual(output_2.number_of_samples, expected_output)
        self.assertEqual(output_2.x_data.tolist(), expected_output_2.tolist())
        self.assertEqual(output_2.dim_vector, expected_output_3)
        self.assertEqual(output_2.data_headers, expected_output_4)

    def test_initialization_07(self):
        with self.assertRaises(Exception):
            UniformSampling(self.test_data, [3], edges='x')

    def test_initialization_08(self):
        output = UniformSampling(self.test_data, [3])
        self.assertEqual(output.edge, True)

    def test_initialization_09(self):
        output = UniformSampling(self.test_data, [3], edges=True)
        self.assertEqual(output.edge, True)

    def test_initialization_10(self):
        output = UniformSampling(self.test_data, [3], edges=False)
        self.assertEqual(output.edge, False)

    def test_initialization_11(self):
        test_data_as_array = self.test_data.values
        output = UniformSampling(test_data_as_array, [7])
        np.testing.assert_array_equal(output.data, test_data_as_array)
        self.assertEqual(output.number_of_samples, 7)

    def test_initialization_12(self):
        test_data_as_list = self.test_data.values.tolist()
        with self.assertRaises(ValueError):
            UniformSampling(test_data_as_list,[7])


class HaltonSamplingTestCases(unittest.TestCase):

    def setUp(self):
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})
        return self.test_data

    """
    Tests 1 - 7: Tests to check the range of potential inputs for the number of samples:
     - 01: Test behaviour when number_of_samples = -ve. Should raise exception.
     - 02: Test behaviour when number_of_samples = 0. Should raise exception.
     - 03: Test behaviour when number_of_samples = + integer within acceptable range. Should complete successfully.
     - 04: Test behaviour when number_of_samples > acceptable range. Should raise exception.
     - 05: Test behaviour when number_of_samples = None. Should complete successfully with default number_of_samples
     - 06: Test behaviour when number_of_samples = fractional. Should raise exception.
     - 07: Test that data loading behaves as expected when a numpy array is provided.
     - 08: Test that a ValueError is thrown when the input is not a numpy array or Pandas DataFrame
    Tests 3 and 5 also check x_data attribute is returned correctly.
    """

    def test_initialization_01(self):
        with self.assertRaises(Exception):
            HaltonSampling(self.test_data, -1)

    def test_initialization_02(self):
        with self.assertRaises(Exception):
            HaltonSampling(self.test_data, 0)

    def test_initialization_03(self):
        output_1 = HaltonSampling(self.test_data, 7)
        expected_output_1 = 7
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_1.number_of_samples, expected_output_1)
        self.assertEqual(output_1.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_04(self):
        with self.assertRaises(Exception):
            HaltonSampling(self.test_data, 11)

    def test_initialization_05(self):
        output_2 = HaltonSampling(self.test_data)
        expected_output = 5
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_2.number_of_samples, expected_output)
        self.assertEqual(output_2.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_06(self):
        with self.assertRaises(Exception):
            HaltonSampling(self.test_data, 1.7)

    def test_initialization_07(self):
        test_data_as_array = self.test_data.values
        output_1 = HaltonSampling(test_data_as_array, 7)
        np.testing.assert_array_equal(output_1.data, test_data_as_array)

    def test_initialization_08(self):
        test_data_as_list = self.test_data.values.tolist()
        with self.assertRaises(ValueError):
            HaltonSampling(test_data_as_list, 7)


class HammersleySamplingTestCases(unittest.TestCase):

    def setUp(self):
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})
        return self.test_data

    """
    Tests 1 - 7: Tests to check the range of potential inputs for the number of samples:
     - 01: Test behaviour when number_of_samples = -ve. Should raise exception.
     - 02: Test behaviour when number_of_samples = 0. Should raise exception.
     - 03: Test behaviour when number_of_samples = + integer within acceptable range. Should complete successfully.
     - 04: Test behaviour when number_of_samples > acceptable range. Should raise exception.
     - 05: Test behaviour when number_of_samples = None. Should complete successfully with default number_of_samples
     - 06: Test behaviour when number_of_samples = fractional. Should raise exception.
     - 07: Test that data loading behaves as expected when a numpy array is provided.
     - 08: Test that a ValueError is thrown when the input is not a numpy array or Pandas DataFrame
    Tests 3 and 5 also check x_data attribute is returned correctly.
    """

    def test_initialization_01(self):
        with self.assertRaises(Exception):
            HammersleySampling(self.test_data, -1)

    def test_initialization_02(self):
        with self.assertRaises(Exception):
            HammersleySampling(self.test_data, 0)

    def test_initialization_03(self):
        output_1 = HammersleySampling(self.test_data, 7)
        expected_output_1 = 7
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_1.number_of_samples, expected_output_1)
        self.assertEqual(output_1.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_04(self):
        with self.assertRaises(Exception):
            HammersleySampling(self.test_data, 11)

    def test_initialization_05(self):
        output_2 = HammersleySampling(self.test_data)
        expected_output = 5
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_2.number_of_samples, expected_output)
        self.assertEqual(output_2.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_06(self):
        with self.assertRaises(Exception):
            HammersleySampling(self.test_data, 1.7)

    def test_initialization_07(self):
        test_data_as_array = self.test_data.values
        output_1 = HammersleySampling(test_data_as_array, 7)
        np.testing.assert_array_equal(output_1.data, test_data_as_array)

    def test_initialization_08(self):
        test_data_as_list = self.test_data.values.tolist()
        with self.assertRaises(ValueError):
            HammersleySampling(test_data_as_list, 7)


class CVTSamplingTestCases(unittest.TestCase):

    def setUp(self):
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})
        return self.test_data

    """
    Tests 1 - 7: Tests to check the range of potential inputs for CVTSampling initialization:
     - 01: Test behaviour when number_of_samples = -ve. Should raise exception.
     - 02: Test behaviour when number_of_samples = 0. Should raise exception.
     - 03: Test behaviour when number_of_samples = + integer within acceptable range. Should complete successfully.
     - 04: Test behaviour when number_of_samples > acceptable range. Should raise exception.
     - 05: Test behaviour when number_of_samples = None. Should complete successfully with default number_of_samples
     - 06: Test behaviour when number_of_samples = fractional. Should raise exception.
     - 07: Test that data loading behaves as expected when a numpy array is provided.
     - 08: Test that a ValueError is thrown when the input is not a numpy array or Pandas DataFrame
     - 09&10: Test behaviour when tolerance >= acceptable maximum. Should raise Exception.
     - 11: Test behaviour when tolerance is within acceptable range. Should complete successfully, and the eps attribute should be properly loaded.
     - 12: Test behaviour when tolerance = None. Should complete successfully with the eps attribute set to the value of 1e-7.
     - 13: Test behaviour when tolerance < aavisable value. Should complete loading with the eps attrinute set to the provided value, but should raise warning (not tested).
     - 14: Test behaviour when tolerance input is of wrong type. Should raise Exception.
    """

    def test_initialization_01(self):
        with self.assertRaises(Exception):
            CVTSampling(self.test_data, -1)

    def test_initialization_02(self):
        with self.assertRaises(Exception):
            CVTSampling(self.test_data, 0)

    def test_initialization_03(self):
        output_1 = CVTSampling(self.test_data, 7)
        expected_output_1 = 7
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_1.number_of_centres, expected_output_1)
        self.assertEqual(output_1.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_04(self):
        with self.assertRaises(Exception):
            CVTSampling(self.test_data, 11)

    def test_initialization_05(self):
        output_2 = CVTSampling(self.test_data)
        expected_output = 5
        expected_output_2 = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        self.assertEqual(output_2.number_of_centres, expected_output)
        self.assertEqual(output_2.x_data.tolist(), expected_output_2.tolist())

    def test_initialization_06(self):
        with self.assertRaises(Exception):
            CVTSampling(self.test_data, 1.7)

    def test_initialization_07(self):
        test_data_as_array = self.test_data.values
        output_1 = CVTSampling(test_data_as_array, 7)
        np.testing.assert_array_equal(output_1.data, test_data_as_array)

    def test_initialization_08(self):
        test_data_as_list = self.test_data.values.tolist()
        with self.assertRaises(ValueError):
            CVTSampling(test_data_as_list, 7)

    def test_initialization_09(self):
        with self.assertRaises(Exception):
            CVTSampling(self.test_data, tolerance=1)

    def test_initialization_10(self):
        with self.assertRaises(Exception):
            CVTSampling(self.test_data, tolerance=0.1)

    def test_initialization_11(self):
        tol = 0.01
        output = CVTSampling(self.test_data, tolerance=tol)
        self.assertEqual(output.eps, tol)

    def test_initialization_12(self):
        deft = 1e-07
        output = CVTSampling(self.test_data)
        self.assertEqual(output.eps, deft)

    def test_initialization_13(self):
        tol = 1e-9
        output = CVTSampling(self.test_data, tolerance=tol)
        self.assertEqual(output.eps, tol)

    def test_initialization_14(self):
        with self.assertRaises(Exception):
            CVTSampling(self.test_data, tolerance='0.1')

    """
    Tests: : Unit tests for eucl_distance, a function that evaluates the distance between two points (u, v).
        Four demonstration tests are done:
        The first test checks that the correct result is obtained when both inputs are single value arrays.
        The second test checks that the correct result is returned when both inputs are 2D vectors.
        The third test checks that the correct result is returned when both inputs are arrays of the same size. 
        The fourth test checks that the function is able to calculate the distance from a single point (n x 1 row vector) to a set of design points (supplied in an n x m array)
    """

    def test_eucl_distance_01(self):
        u = np.array( [ [3] ])
        v = np.array( [ [5] ])
        expected_output = 2
        output = CVTSampling.eucl_distance(u, v)
        self.assertEqual(expected_output, output)

    def test_eucl_distance_02(self):
        u = np.array([ [1, 2] ])
        v = np.array([ [3, 4] ])
        expected_output = 8 ** 0.5
        output = CVTSampling.eucl_distance(u, v)
        self.assertEqual(expected_output, output)

    def test_eucl_distance_03(self):
        u = np.array([ [1, 2], [3, 4] ])
        v = np.array([ [5, 6], [7, 8] ])
        expected_output = np.array([32**0.5, 32**0.5])
        output = CVTSampling.eucl_distance(u, v)
        np.testing.assert_array_equal(expected_output, output)

    def test_eucl_distance_04(self):
        u = np.array([ [1, 2]])
        v = np.array([ [5, 6], [7, 8] ])
        expected_output = np.array([32**0.5, 72**0.5])
        output = CVTSampling.eucl_distance(u, v)
        np.testing.assert_array_equal(expected_output, output)

    """
    Tests: : Unit tests for create_centres, a function that generates new mass centroids for the design space based on McQueen's method.
        Four demonstration tests are done:
        The first test checks that the correct result for the new centres is returned when the counter is at its lowest value (1).
        The second test checks that the correct result for the new centres is returned when the counter is at an arbitrary value (10).
        The third test checks that the correct procedure is followed when one of the centres in initial_centres has no close design point to it.
        The fourth test checks that the the approach works as expected for problems with more than two dimensions.
    """

    def test_create_centres_01(self):
        initial_centres = np.array([ [0, 0], [1, 1] ])
        current_random_points = np.array([[0.6, 0.6], [0.3, 0.3]])
        current_centres = np.array([1, 0])
        counter = 1
        expected_output = np.array([ [0.15, 0.15], [0.8, 0.8]])
        output = CVTSampling.create_centres(initial_centres, current_random_points, current_centres, counter)
        np.testing.assert_array_equal(expected_output, output)

    def test_create_centres_02(self):
        initial_centres = np.array([ [0, 0], [1, 1] ])
        current_random_points = np.array([[0.6, 0.6], [0.3, 0.3]])
        current_centres = np.array([1, 0])
        counter = 10
        expected_output = np.array([ [0.3/11, 0.3/11], [10.6/11, 10.6/11]])
        output = CVTSampling.create_centres(initial_centres, current_random_points, current_centres, counter)
        np.testing.assert_array_equal(expected_output, output)

    def test_create_centres_03(self):
        initial_centres = np.array([ [0, 0], [1, 1] ])
        current_random_points = np.array([[0.6, 0.6], [0.8, 0.8]])
        current_centres = np.array([1, 1])
        counter = 5
        expected_output = np.array([ [0.5/6, 0.5/6], [5.7/6, 5.7/6]])
        output = CVTSampling.create_centres(initial_centres, current_random_points, current_centres, counter)
        np.testing.assert_array_equal(expected_output, output)

    def test_create_centres_04(self):
        initial_centres = np.array([ [0, 0, 0], [1, 1, 1] ])
        current_random_points = np.array([[0.1, 0.1, 0.1], [0.3, 0.3, 0.3], [0.5, 0.5, 0.5], [0.7, 0.7, 0.7], [0.9, 0.9, 0.9]])
        current_centres = np.array([0, 0, 0, 1, 1])
        counter = 4
        expected_output = np.array([ [0.3/5, 0.3/5, 0.3/5], [4.8/5, 4.8/5, 4.8/5]])
        output = CVTSampling.create_centres(initial_centres, current_random_points, current_centres, counter)
        np.testing.assert_array_equal(expected_output, output)


class SampliingMethodsTestCases(unittest.TestCase):
    def setUp(self):
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        self.test_data = pd.DataFrame({'x': input_array[:, 0], 'y': input_array[:, 1]})
        return self.test_data

    """
    Test: Tests the function nearest_neighbour for a range of possible points in the input data
    """

    def test_nearest_neighbour(self):
        output_1 = SamplingMethods.nearest_neighbour(self, self.test_data.values, 0.35)
        output_2 = SamplingMethods.nearest_neighbour(self, self.test_data.values, 4.46)
        output_3 = SamplingMethods.nearest_neighbour(self, self.test_data.values, 4.56)
        output_4 = SamplingMethods.nearest_neighbour(self, self.test_data.values, 8.79)
        output_5 = SamplingMethods.nearest_neighbour(self, self.test_data.values, 0)
        expected_output_1 = 0
        expected_output_2 = 4
        expected_output_3 = 5
        expected_output_4 = 9
        expected_output_5 = 0
        expected_output_6 = 1
        expected_output_7 = 25
        expected_output_8 = 36
        expected_output_9 = 100
        expected_output_10 = 1
        self.assertEqual(output_1[0, ], expected_output_1)
        self.assertEqual(output_2[0, ], expected_output_2)
        self.assertEqual(output_3[0, ], expected_output_3)
        self.assertEqual(output_4[0, ], expected_output_4)
        self.assertEqual(output_5[0, ], expected_output_5)
        self.assertEqual(output_1[1, ], expected_output_6)
        self.assertEqual(output_2[1, ], expected_output_7)
        self.assertEqual(output_3[1, ], expected_output_8)
        self.assertEqual(output_4[1, ], expected_output_9)
        self.assertEqual(output_5[1, ], expected_output_10)

    """
    Test: Tests the points_selection function for a range of possible points in the input data.
    Expected to return the closest point in the actual dataset based on L-2 distance.
    Tests 01-05: points_selection is called with every single sampling technique it applies to.
    """

    def test_points_selection_01(self):
        number_of_samples = 5
        x = np.array([[8.79], [0.35], [4.56]])
        data_feed = LatinHypercubeSampling(self.test_data, number_of_samples)
        output = data_feed.points_selection(self.test_data.values, x)
        expected_output = np.array([[9., 100.], [0., 1.], [5., 36.]])
        self.assertEqual(output.tolist(), expected_output.tolist())

    def test_points_selection_02(self):
        number_of_samples = 5
        x = np.array([[8.79], [0.35], [4.56]])
        data_feed = HammersleySampling(self.test_data, number_of_samples)
        output = data_feed.points_selection(self.test_data.values, x)
        expected_output = np.array([[9., 100.], [0., 1.], [5., 36.]])
        self.assertEqual(output.tolist(), expected_output.tolist())

    def test_points_selection_03(self):
        number_of_samples = 5
        x = np.array([[8.79], [0.35], [4.56]])
        data_feed = HaltonSampling(self.test_data, number_of_samples)
        output = data_feed.points_selection(self.test_data.values, x)
        expected_output = np.array([[9., 100.], [0., 1.], [5., 36.]])
        self.assertEqual(output.tolist(), expected_output.tolist())

    def test_points_selection_04(self):
        number_of_samples = 5
        x = np.array([[8.79], [0.35], [4.56]])
        data_feed = CVTSampling(self.test_data, number_of_samples)
        output = data_feed.points_selection(self.test_data.values, x)
        expected_output = np.array([[9., 100.], [0., 1.], [5., 36.]])
        self.assertEqual(output.tolist(), expected_output.tolist())

    def test_points_selection_05(self):
        nc = [5]
        x = np.array([[8.79], [0.35], [4.56]])
        data_feed = UniformSampling(self.test_data, list_of_samples_per_variable=nc)
        output = data_feed.points_selection(self.test_data.values, x)
        expected_output = np.array([[9., 100.], [0., 1.], [5., 36.]])
        self.assertEqual(output.tolist(), expected_output.tolist())

    """
    Tests: : Unit tests for prime_number_generator, a function that generates prime numbers.
        Two demonstration tests are done:
        The first test checks that the function returns the right result for the minimum value of n (n=1).
        The second test checks that the function returns the right result for an arbitrary value of n (n=10).
    """

    def test_prime_number_generator_01(self):
        output = SamplingMethods.prime_number_generator(self, n=1)
        expected_output = [2]
        self.assertEqual(expected_output, output)

    def test_prime_number_generator_02(self):
        output = SamplingMethods.prime_number_generator(self, n=10)
        expected_output = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        self.assertEqual(expected_output, output)

    """
    Tests: : Unit tests for base_conversion, converts numbers from base 10 to any base.
        For a number 'a' to be converted to base 'b', four demonstration tests are done:
        - The first and second tests check that the function returns the right string result for small bases when a > b 
        - The third test checks that the function returns the right string result for arbitrarily large bases when a > b
        - The third test checks that the function returns the right string  result ('a') for a < b
    """

    def test_base_conversion_01(self):
        output = SamplingMethods.base_conversion(self, a=4, b=2)
        expected_output = ['1','0','0']
        self.assertEqual(expected_output, output)

    def test_base_conversion_02(self):
        output = SamplingMethods.base_conversion(self, a=5, b=3)
        expected_output = ['1','2']
        self.assertEqual(expected_output, output)

    def test_base_conversion_03(self):
        output = SamplingMethods.base_conversion(self, a=45, b=33)
        expected_output = ['1','12']
        self.assertEqual(expected_output, output)

    def test_base_conversion_04(self):
        output = SamplingMethods.base_conversion(self, a=3, b=5)
        expected_output = ['3']
        self.assertEqual(expected_output, output)

    """
    Tests: : Unit tests for prime_base_to_decimal, a function that converts decimal numbers from base any base to base 10.
       The test checks that the current values are returned when 0.11 (in string form) is converted from bases 2, 5 and 19 to base 10.
    """

    def test_prime_base_to_decimal(self):
        base_repr = ['0.', '1', '1']
        output_1 = SamplingMethods.prime_base_to_decimal(self, num=base_repr, base = 2)
        output_2 = SamplingMethods.prime_base_to_decimal(self, num=base_repr, base=5)
        output_3 = SamplingMethods.prime_base_to_decimal(self, num=base_repr, base=19)
        expected_output_1 = 0.75
        expected_output_2 = 0.24
        expected_output_3 = (1/19) + (1/361)
        self.assertEqual(expected_output_1, output_1)
        self.assertAlmostEqual(expected_output_2, output_2, places=6)
        self.assertAlmostEqual(expected_output_3, output_3, places=6)

    """
    Tests: : Unit tests for data_sequencing, a function which generates the first no_samples elements of the Halton or Hammersley sequence.
        Three demonstration tests are done:
        - The first test checks that the function returns the first five elements of the van der Corput sequence for base 2 (no_samples > prime_base)
        - The second test checks that the function returns the first element of the van der Corput sequence for base 3 (no_samples < prime_base)
        - The third test checks that the function returns the first five elements of the van der Corput sequence for base 5 (no_samples = prime_base)
    """

    def test_data_sequencing_01(self):
        data_feed = HaltonSampling(self.test_data, number_of_samples=5)
        output = data_feed.data_sequencing(no_samples=5, prime_base=2)
        expected_output = np.array( [0, 0.5, 0.25, 0.75, 0.125] )
        np.testing.assert_array_equal(expected_output, output)

    def test_data_sequencing_02(self):
        data_feed = HaltonSampling(self.test_data, number_of_samples=5)
        output = data_feed.data_sequencing(no_samples=1, prime_base=3)
        expected_output = np.array( [0] )
        np.testing.assert_array_equal(expected_output, output)

    def test_data_sequencing_03(self):
        data_feed = HammersleySampling(self.test_data, number_of_samples=5)
        output = data_feed.data_sequencing(no_samples=5, prime_base=5)
        expected_output = np.array( [0, 0.2, 0.4, 0.6, 0.8] )
        np.testing.assert_array_equal(expected_output, output)


if __name__ == '__main__':
    unittest.main()



