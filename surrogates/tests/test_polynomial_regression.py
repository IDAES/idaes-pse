# from github_files.polynomial_regression.polynomial_regression import PolynomialRegression, FeatureScaling
from pysmo.polynomial_regression import PolynomialRegression, FeatureScaling
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
        output_1, output_2, output_3 = FeatureScaling.data_scaling(self.test_data)
        expected_output_2 = np.array([[0, 1]])
        expected_output_3 = np.array([[9, 100]])
        expected_output_1 = (self.test_data - expected_output_2) / (expected_output_3 - expected_output_2)
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    def test_data_scaling_02(self):
        test_data_as_array = self.test_data.values
        output_1, output_2, output_3 = FeatureScaling.data_scaling(test_data_as_array[:, 0])
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
            FeatureScaling.data_scaling(test_data_as_list)

    def test_data_unscaling_01(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1, 6]])
        max_array = np.array([[5, 10]])
        output = FeatureScaling.data_unscaling(data_array, min_array, max_array)
        expected_output = np.array([[1, 10], [2, 9], [3, 8], [4, 7], [5, 6]])
        np.testing.assert_array_equal(output, expected_output)

    def test_data_unscaling_02(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with self.assertRaises(IndexError):
            FeatureScaling.data_unscaling(data_array, min_array, max_array)

    def test_data_unscaling_03(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1, 6, 10]])
        max_array = np.array([[5, 10, 15]])
        with self.assertRaises(IndexError):
            FeatureScaling.data_unscaling(data_array, min_array, max_array)

    def test_data_unscaling_04(self):
        data_array = np.array([[0, 1], [0.25, 0.75], [0.5, 0.5], [0.75, 0.25], [1, 0]])
        min_array = np.array([[1, 6]])
        max_array = np.array([[5, 10]])
        output = FeatureScaling.data_unscaling(data_array[:, 0], min_array[:, 0], max_array[:, 0])
        expected_output = np.array([[1], [2], [3], [4], [5]])
        np.testing.assert_array_equal(output, expected_output)


# ######################################################################################################################
class PolynomialRegressionTestCases(unittest.TestCase):
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

    def mock_optimization(self, x, y):
        return 10 * np.ones((x.shape[1], 1))

    """
    Test 1: Check that all default values are correctly loaded when input is sparse
    """

    def test_initialization_polyregression_01(self):
        sample_points = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64]])
        maximum_polynomial_order = 5
        output_1 = PolynomialRegression(self.test_data, sample_points,
                                        maximum_polynomial_order=maximum_polynomial_order)
        assert output_1.max_polynomial_order == maximum_polynomial_order
        assert output_1.number_of_crossvalidations == 3  # Default number of cross-validations
        assert output_1.no_adaptive_samples == 4  # Default number of adaptive samples
        assert output_1.fraction_training == 0.75  # Default training split
        assert output_1.max_fraction_training_samples == 0.5  # Default fraction for the maximum number of training samples
        assert output_1.max_iter == 10  # Default maximum number of iterations

    """
    Tests 2 - 6 : Tests to check the behaviour of the class attribute max_polynomial_order to different input values
     - 02: Test behaviour when maximum_polynomial_order = -ve. Should raise exception.
     - 03: Test behaviour when maximum_polynomial_order = 0. Should raise exception.
     - 04: Test behaviour when maximum_polynomial_order = fractional. Should raise exception.
     - 05: Test behaviour when maximum_polynomial_order > number of samples. Should raise exception.
     - 06: Test behaviour when maximum_polynomial_order = 5 (acceptable value). Assert that attribute is correctly loaded.
    - 06: Test behaviour when maximum_polynomial_order input > 10. Assert that attribute is set to maximum allowable (10) and warning displayed.
    """

    def test_initialization_polyregression_02(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = -1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=maximum_polynomial_order)

    def test_initialization_polyregression_03(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 0
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=maximum_polynomial_order)

    def test_initialization_polyregression_04(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 1.5
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=maximum_polynomial_order)

    def test_initialization_polyregression_05(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 4
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=maximum_polynomial_order)

    def test_initialization_polyregression_06(self):
        sample_points = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        maximum_polynomial_order = 5
        output_1 = PolynomialRegression(self.test_data, sample_points,
                                        maximum_polynomial_order=maximum_polynomial_order)
        self.assertEqual(output_1.max_polynomial_order, maximum_polynomial_order)

    def test_initialization_polyregression_38(self):
        maximum_polynomial_order = 15
        output_1 = PolynomialRegression(self.full_data, self.training_data,
                                        maximum_polynomial_order=maximum_polynomial_order)
        self.assertEqual(output_1.max_polynomial_order, 10)  # Reset to 10

    """
    Tests 7-12 : Tests to check the behaviour of the class attribute max_iter to different input values
     - 07: Test behaviour when max_iter input = -ve. Should raise exception.
     - 08: Test behaviour when max_iter input = 0. Attribute should be set to zero.
     - 09: Test behaviour when max_iter input = fractional. Should raise exception.
     - 10: Test behaviour when max_iter input  = + integer. Attribute should be set to value.
     - 11: Test behaviour when max_iter input  = None. Attribute should be set to default.
     - 12: Test behaviour when sample data and original data are the same size. Max_iter attribute should be set to zero irrespective of supplied value.
    """

    def test_initialization_polyregression_07(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 2
        max_iter = -1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=maximum_polynomial_order,
                                 max_iter=max_iter)

    def test_initialization_polyregression_08(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 2
        max_iter = 0
        output_1 = PolynomialRegression(self.test_data, sample_points,
                                        maximum_polynomial_order=maximum_polynomial_order, max_iter=max_iter)
        self.assertEqual(output_1.max_iter, max_iter)

    def test_initialization_polyregression_09(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 2
        max_iter = 1.75
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=maximum_polynomial_order,
                                 max_iter=max_iter)

    def test_initialization_polyregression_10(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 2
        max_iter = 20
        output_1 = PolynomialRegression(self.test_data, sample_points,
                                        maximum_polynomial_order=maximum_polynomial_order, max_iter=max_iter)
        self.assertEqual(output_1.max_iter, max_iter)

    def test_initialization_polyregression_11(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        maximum_polynomial_order = 2
        output_1 = PolynomialRegression(self.test_data, sample_points,
                                        maximum_polynomial_order=maximum_polynomial_order)
        self.assertEqual(output_1.max_iter, 10)  # Should return default value

    def test_initialization_polyregression_12(self):
        sample_points = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        maximum_polynomial_order = 2
        max_iter = 5
        output_1 = PolynomialRegression(self.test_data, sample_points,
                                        maximum_polynomial_order=maximum_polynomial_order, max_iter=max_iter)
        self.assertEqual(output_1.max_iter,
                         0)  # Should be zero when the original data and sample data are the same size

    """
    Tests 13-18 : Tests to check the behaviour of the class attribute number_of_crossvalidations to different input values
     - 13: Test behaviour when number_of_crossvalidations input = -ve. Should raise exception.
     - 14: Test behaviour when number_of_crossvalidations input = 0. Should raise exception since minimum allowable value is 1.
     - 15: Test behaviour when number_of_crossvalidations input  = None. Attribute should be set to default.
     - 16: Test behaviour when number_of_crossvalidations input  = + integer. Attribute should be set to value.
     - 17: Test behaviour when number_of_crossvalidations input = fractional. Should raise exception.
     - 18: Test behaviour when number_of_crossvalidations input > 10. Attribute should be set to value but warning given.
    """

    def test_initialization_polyregression_13(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        num_cv = -1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 number_of_crossvalidations=num_cv)

    def test_initialization_polyregression_14(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        num_cv = 0
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 number_of_crossvalidations=num_cv)

    def test_initialization_polyregression_15(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)
        self.assertEqual(output_1.number_of_crossvalidations, 3)  # Should return default value

    def test_initialization_polyregression_16(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        num_cv = 5
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        number_of_crossvalidations=num_cv)
        self.assertEqual(output_1.number_of_crossvalidations, num_cv)

    def test_initialization_polyregression_17(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        num_cv = 1.5
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 number_of_crossvalidations=num_cv)

    def test_initialization_polyregression_18(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        num_cv = 15
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        number_of_crossvalidations=num_cv)
        self.assertEqual(output_1.number_of_crossvalidations,
                         num_cv)  # Should return default value, check that warning is returned?

    """
    Tests 19-23 : Tests to check the behaviour of the class attribute fraction_training (training_split) to different input values
     - 19: Test behaviour when fraction_training input = -ve. Should raise exception.
     - 20: Test behaviour when fraction_training input = 0. Should raise exception; the value should be greater than zero.
     - 21: Test behaviour when fraction_training input  = None. Attribute should be set to default value of 0.75.
     - 22: Test behaviour when fraction_training input  = valid fraction close to UB of 1. Attribute should be set to value.
     - 23: Test behaviour when fraction_training input = 1. Should raise exception; the value should be less than one.
    """

    def test_initialization_polyregression_19(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        t_s = -0.1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 training_split=t_s)

    def test_initialization_polyregression_20(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        t_s = 0
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 training_split=t_s)

    def test_initialization_polyregression_21(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)
        self.assertEqual(output_1.fraction_training, 0.75)  # Default value

    def test_initialization_polyregression_22(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        t_s = 0.99
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        training_split=t_s)
        self.assertEqual(output_1.fraction_training, t_s)

    def test_initialization_polyregression_23(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        t_s = 1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 training_split=t_s)

    """
    Tests 24-28 : Tests to check the behaviour of the class attribute max_fraction_training_samples to different input values
     - 24: Test behaviour when max_fraction_training_samples input = -ve. Should raise exception.
     - 25: Test behaviour when max_fraction_training_samples input = 0. Attribute should be set to zero.
     - 26: Test behaviour when max_fraction_training_samples input  = None. Attribute should be set to default value of 0.5.
     - 27: Test behaviour when max_fraction_training_samples input = 1 (Upper bound). Attribute should be set to one.
     - 28: Test behaviour when max_fraction_training_samples input > 1. Should raise exception.
    """

    def test_initialization_polyregression_24(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        max_tf = -0.1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 max_fraction_training_samples=max_tf)

    def test_initialization_polyregression_25(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        max_tf = 0
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        max_fraction_training_samples=max_tf)
        self.assertEqual(output_1.max_fraction_training_samples, max_tf)

    def test_initialization_polyregression_26(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)
        self.assertEqual(output_1.max_fraction_training_samples, 0.5)  # Default value

    def test_initialization_polyregression_27(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        max_tf = 1
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        max_fraction_training_samples=max_tf)
        self.assertEqual(output_1.max_fraction_training_samples, max_tf)

    def test_initialization_polyregression_28(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        max_tf = 1.01
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 max_fraction_training_samples=max_tf)

    """
    Tests 29-33 : Tests to check the behaviour of the class attribute no_adaptive_samples to different input values
     - 29: Test behaviour when no_adaptive_samples input = -ve. Should raise exception.
     - 30: Test behaviour when no_adaptive_samples input = 0. Attribute should be set to zero. Max_iter must also forced to zero irrespective of entered value.
     - 31: Test behaviour when no_adaptive_samples input  = None. Attribute should be set to default value of 4.
     - 32: Test behaviour when no_adaptive_samples input = +ve integer. Attribute should be set to value.
     - 33: Test behaviour when no_adaptive_samples input = Fractional. Should raise exception.
    """

    def test_initialization_polyregression_29(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        no_as = -1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 no_adaptive_samples=no_as)

    def test_initialization_polyregression_30(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        no_as = 0
        max_iter = 5
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        no_adaptive_samples=no_as, max_iter=max_iter)
        self.assertEqual(output_1.no_adaptive_samples, no_as)
        self.assertEqual(output_1.max_iter, 0)

    def test_initialization_polyregression_31(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)
        self.assertEqual(output_1.no_adaptive_samples, 4)  # Default value

    def test_initialization_polyregression_32(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        no_as = 5
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        no_adaptive_samples=no_as)
        self.assertEqual(output_1.no_adaptive_samples, no_as)

    def test_initialization_polyregression_33(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        no_as = 5.01
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 no_adaptive_samples=no_as)

    """
    Tests 34-37 : Demonstration of behaviour when potential dimensionality challenges may arise
     - 34: Test behaviour when there are more samples for training than input data points. Code should raise exception.
     - 35: Test behaviour when the number of samples for training is the same as the number of input data points. In this case, the number of iterations must be set to zero.
     - 36: Test behaviour when the x + y variables in the sample and original datasets is different. Code should raise exception. 
     - 37: Test behaviour when when there is only one column in the input or sample data. Code should raise exception.
     - 45: Ensures that a ValueError is thrown when the input training data is not a a numpy array or dataframe.
     - 46: Ensures that a ValueError is thrown when the original (input) dataset is not a numpy array or dataframe.
    """

    def test_initialization_polyregression_34(self):
        sample_points = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100], [10, 121]])
        max_poly_order = 2
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)

    def test_initialization_polyregression_35(self):
        sample_points = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        sample_points_as_df = pd.DataFrame({'x': sample_points[:, 0], 'y': sample_points[:, 1]})
        max_poly_order = 2
        output_1 = PolynomialRegression(self.test_data, sample_points_as_df, maximum_polynomial_order=max_poly_order)
        self.assertEqual(output_1.max_iter, 0)

    def test_initialization_polyregression_36(self):
        sample_points = np.array(
            [[0, 1, 1], [1, 4, 5], [2, 9, 11], [3, 16, 19], [4, 25, 29]])  # More columns than original data
        max_poly_order = 2
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)

    def test_initialization_polyregression_37(self):
        sample_points = np.array([[0], [1], [2], [3], [4]])  # More columns than original data
        max_poly_order = 2
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data[:, 0].reshape(self.test_data.shape[0], 1), sample_points,
                                 maximum_polynomial_order=max_poly_order)

    def test_initialization_polyregression_45(self):
        sample_points = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        sample_points_as_list = sample_points.tolist()
        max_poly_order = 2
        with self.assertRaises(ValueError):
            PolynomialRegression(self.test_data, sample_points_as_list, maximum_polynomial_order=max_poly_order)

    def test_initialization_polyregression_46(self):
        sample_points = np.array(
            [[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        test_data_as_list = self.test_data.tolist()
        max_poly_order = 2
        with self.assertRaises(ValueError):
            PolynomialRegression(test_data_as_list, sample_points, maximum_polynomial_order=max_poly_order)

    """
    Tests 39-44 : Tests to class attribute solve_method.
     - 39: Test behaviour when solve_method input is not a string. Should raise exception.
     - 40: Test behaviour when solve_method = 'mLe'. Attribute should be set to 'mle' irrespective of case.
     - 41: Test behaviour when solve_method = 'pYomO'. Attribute should be set to 'pyomo' irrespective of case.
     - 42: Test behaviour when solve_method = 'bFGs'. Attribute should be set to 'bfgs' irrespective of case.
     - 43: Test behaviour when solve_method is an invalid string. Should raise exception.
     - 44: Test behaviour when solve_method is None. Attribute should be set to default ('pyomo').
    """

    def test_initialization_polyregression_39(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        solve_method = 1
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 solution_method=solve_method)

    def test_initialization_polyregression_40(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        solve_method = 'mLe'
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        solution_method=solve_method)
        self.assertEqual(output_1.solution_method, solve_method.lower())

    def test_initialization_polyregression_41(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        solve_method = 'pYomO'
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        solution_method=solve_method)
        self.assertEqual(output_1.solution_method, solve_method.lower())

    def test_initialization_polyregression_42(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        solve_method = 'bFGs'
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        solution_method=solve_method)
        self.assertEqual(output_1.solution_method, solve_method.lower())

    def test_initialization_polyregression_43(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        solve_method = 'xyz'
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 solution_method=solve_method)

    def test_initialization_polyregression_44(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        default_solve = 'pyomo'
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)
        self.assertEqual(output_1.solution_method, default_solve)

    """
    Tests 47-51 : Tests to class attribute mumtinomials.
     - 47: Test behaviour when multinomials input is a string. Should raise exception.
     - 48: Test behaviour when multinomials input = 0. Attribute should be set to value.
     - 49: Test behaviour when multinomials input = 1. Attribute should be set to value.
     - 50: Test behaviour when multinomials input = None. Attribute should be set to default (1).
     - 51: Test behaviour when multinomials is out of range. Exception should be raised.
    """

    def test_initialization_polyregression_47(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        mn = 'a'
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 multinomials=mn)

    def test_initialization_polyregression_48(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        mn = 0
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        multinomials=mn)
        self.assertEqual(output_1.multinomials, mn)

    def test_initialization_polyregression_49(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        mn = 1
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                        multinomials=mn)
        self.assertEqual(output_1.multinomials, mn)

    def test_initialization_polyregression_50(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        output_1 = PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order)
        self.assertEqual(output_1.multinomials, 1)

    def test_initialization_polyregression_51(self):
        sample_points = np.array([[8, 81], [2, 9], [5, 36]])
        max_poly_order = 2
        mn = 2
        with self.assertRaises(Exception):
            PolynomialRegression(self.test_data, sample_points, maximum_polynomial_order=max_poly_order,
                                 multinomials=mn)

    """
    Tests: : Unit tests for polygeneration function. Five tests covering the range of max_polynomial_order are considered.
        The first test checks matrix/array returned when polynomial order is 1 -  the minimum possible value.
        The second test checks matrix/array returned when polynomial order is 2.
        The third test checks matrix/array returned when polynomial order is 10 (maximum possible value) with multinomials = default (1).
        The fourth test checks matrix/array returned when polynomial order is 10 (maximum possible value) with multinomials = 0.
        The fifth test checks matrix/array returned when polynomial order is 1 and additional terms have been supplied.
        - 
    """

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
        expected_value = 6613.875  # Calculated externally as sum(y^2) / 2m
        output_1 = PolynomialRegression.cost_function(theta, x_vector, y, reg_parameter=0)
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
        theta = np.array([[4.5], [3], [3], [1], [1], [0]])  # coefficients in (x1 + 1.5)^2 + (x2 + 1.5) ^ 2
        expected_value = 90.625  # Calculated externally as sum(dy^2) / 2m
        output_1 = PolynomialRegression.cost_function(theta, x_vector, y, reg_parameter=0)
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
        theta = np.array([[2], [2], [2], [1], [1], [0]])  # Actual coefficients in (x1 + 1)^2 + (x2 + 1) ^ 2
        expected_value = 0  # Value should return zero for exact solution
        output_1 = PolynomialRegression.cost_function(theta, x_vector, y, reg_parameter=0)
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

    """
    Tests: : Unit tests for bfgs_parameter_optimization function. Two unit tests are done:
        The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
        The second test is used to check that higher order multi-variable problems are solved correctly.
        - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
        - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.
    """

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
        data_feed = PolynomialRegression(self.test_data, input_array, maximum_polynomial_order=5,
                                         solution_method='bfgs')
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        self.assertEqual(data_feed.solution_method, 'bfgs')
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
        self.assertEqual(data_feed.solution_method, 'bfgs')
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    """
    Tests: : Unit tests for MLE_estimate function. The same two unit tests are done:
        The first test is used to check that a single variable problem y = (ax + b)^2 is solved correctly.
        The second test is used to check that higher order multi-variable problems are solved correctly.
        - 01: y = a.x ^ 2 + b.x + c, find values a-c given data to 4.d.p.
        - 02: y = a.x_1 ^ 2 + b.x_2 ^ 2 + c.x_1 + e.x_2 + f.x_1 * x_2 + g, find values of a-g given data to 4.d.p.
    """

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

    """
    Tests: : Unit tests for cross_validation_error_calculation function.  The second order problem (ax1 + b)^2 + (cx2 + d)^2 is considered.
        The three tests used for the cost_function are repeated for the cross-validation error. 
        The first test checks that the correct cv error is returned when theta is initialized as a vector of zeros.
        The second test checks that the correct cv error is returned for a random point within the solution space.
        The third test checks that the cv error is zero when the exact solution is found.
    In each case, the result is expected to be twice as large as the cost function's expected value: the cost function is half of the total squared error.
    """

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
        self.assertEqual(output_1, expected_value)

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
        self.assertEqual(output_1, expected_value)

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
        self.assertEqual(output_1, expected_value)

    """
    Tests: : Unit tests for surrogate_performance function.  The single variable problem (x+1) ^ 2 is considered. 
        Two tests considered for the function to validate its MAE, MSE, R2 and adjusted R2 calculations.
        - 01: The errors reported when the initialization values of theta (zeros) are supplied are compared to manually calculated values: see excel sheet.
        - 02: The errors reported when the actual values of theta (zeros) are supplied are checked. The errors are expected to be zero, while R^2 ~= 1.
        - 02: The errors reported when a random value of theta (zeros) is supplied are checked.
    The returned array contains [x,y_real,y_prediction] data and does not need to be checked once the errors are correct.

    """

    def test_surrogate_performance_01(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        order_best = 2
        phi_best = np.array([[0.], [0.], [0.]])
        expected_value_1 = 38.5
        expected_value_2 = 2533.3
        expected_value_3 = -1.410256
        expected_value_4 = 0
        data_feed = PolynomialRegression(self.test_data, input_array, maximum_polynomial_order=5)
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(phi_best, order_best)
        self.assertEqual(output_1, expected_value_1)
        self.assertEqual(output_2, expected_value_2)
        self.assertEqual(np.round(output_3, 4), np.round(expected_value_3, 4))
        self.assertEqual(np.round(output_4, 4), np.round(expected_value_4, 4))

    def test_surrogate_performance_02(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        order_best = 2
        phi_best = np.array([[1.], [2.], [1.]])
        expected_value_1 = 0
        expected_value_2 = 0
        expected_value_3 = 1
        expected_value_4 = 1
        data_feed = PolynomialRegression(self.test_data, input_array, maximum_polynomial_order=5)
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(phi_best, order_best)
        self.assertEqual(output_1, expected_value_1)
        self.assertEqual(output_2, expected_value_2)
        self.assertEqual(np.round(output_3, 4), np.round(expected_value_3, 4))
        self.assertEqual(np.round(output_4, 4), np.round(expected_value_4, 4))

    def test_surrogate_performance_03(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        input_array = np.array([[0, 1], [1, 4], [2, 9], [3, 16], [4, 25], [5, 36], [6, 49], [7, 64], [8, 81], [9, 100]])
        order_best = 2
        phi_best = np.array([[1.], [1.], [1.]])
        expected_value_1 = 4.5
        expected_value_2 = 28.5
        expected_value_3 = 0.972884259
        expected_value_4 = 0.931219147
        data_feed = PolynomialRegression(self.test_data, input_array, maximum_polynomial_order=5)
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(phi_best, order_best)
        self.assertEqual(output_1, expected_value_1)
        self.assertEqual(output_2, expected_value_2)
        self.assertEqual(np.round(output_3, 4), np.round(expected_value_3, 4))
        self.assertEqual(np.round(output_4, 4), np.round(expected_value_4, 4))

    """
       Tests: : Unit tests for error_plotting function. For a sample dataset of the expected shape (n x 8), we verify that:
           - 1: The data we supply is the data used in the plot for each of the subplots, and that the right columns are accessed in each case.
           - 2: The right plot titles are displayed for each figure.
           - 3. The vertices of the path (and thus actual line segments) for each plot is the same as the input data; more information on matplotlib's path may be found at https://matplotlib.org/api/path_api.html

       """

    def test_error_plotting(self):
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
        self.assertTrue(ax1.title.get_text() == 'Training (green) vs Cross-validation error (red)')
        self.assertTrue(ax2.title.get_text() == 'MAE')
        self.assertTrue(ax3.title.get_text() == 'MSE')
        self.assertTrue(ax4.title.get_text() == 'R-squared (blue) and Adjusted R-squared (red)')
        # Compare the line/curve segments for each line plot to the actual input data
        np.testing.assert_array_equal(ax1.lines[0].get_path()._vertices, expected_output_1)
        np.testing.assert_array_equal(ax1.lines[1].get_path()._vertices, expected_output_2)
        np.testing.assert_array_equal(ax2.lines[0].get_path()._vertices, expected_output_3)
        np.testing.assert_array_equal(ax3.lines[0].get_path()._vertices, expected_output_4)
        np.testing.assert_array_equal(ax4.lines[0].get_path()._vertices, expected_output_5)
        np.testing.assert_array_equal(ax4.lines[1].get_path()._vertices, expected_output_6)

    """
    Tests: : Unit tests for training_test_data_creation, the function that splits the sample data into training and cross-validation sets.
    Five tests are performed:
     - 01: The function is tested for the minimum valid number of cross-validations (1). We check:
        - The length of the returned dictionaries for the training and cross-validation sets - must be 1.
        - The the number of entries in the created training and cross-validation sets - must be fn. of training_split entered.
        - When the training and test sets are recombined and sorted, they give the sample data set inputed.

     - 02:  The function is tested for a larger number of cross-validations. In addition to the requirements in test 1, 
         we also check that the cross-validation sets created are different from each other.

     - 03/04: The function is tested for cases in which the training_split is too low or high. Exceptions need to be raised in those cases.

     - 05: Checks that the expected arrays are obtained when additional_features is specified

    """

    def test_training_test_data_creation_01(self):
        num_cv = 1
        split = 0.4
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        training_data, cross_val_data = data_feed.training_test_data_creation()
        expected_training_size = int(np.around(self.training_data.shape[0] * split))
        expected_test_size = self.training_data.shape[0] - expected_training_size

        # Confirm number of created training and test data sets in dictionaries
        self.assertEqual(len(training_data), num_cv)
        self.assertEqual(len(cross_val_data), num_cv)

        # Then confirm shapes
        self.assertEqual(training_data["training_set_1"].shape[0], expected_training_size)
        self.assertEqual(cross_val_data["test_set_1"].shape[0], expected_test_size)

        # Finally, check that concatenation generates the original sample data.
        # ## First, re-combine the training and cross-validation datasets for each cross-validation case
        concat_01 = (np.concatenate((training_data["training_set_1"], cross_val_data["test_set_1"]), axis=0))
        # ## Next, sort the sample and concatenated data in the same order. np.lexsort allows the specification of column order during sorting
        sample_data_sorted = self.training_data[
            np.lexsort((self.training_data[:, 2], self.training_data[:, 1], self.training_data[:, 0]))]
        concat_01_sorted = concat_01[np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))]
        # ## Finally, assert that the re-formed arrays are the same as the sorted input data array
        np.testing.assert_equal(sample_data_sorted, concat_01_sorted)

    def test_training_test_data_creation_02(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        training_data, cross_val_data = data_feed.training_test_data_creation()
        expected_training_size = int(np.around(self.training_data.shape[0] * split))
        expected_test_size = self.training_data.shape[0] - expected_training_size

        # Confirm number of created training and test data sets in dictionaries
        self.assertEqual(len(training_data), num_cv)
        self.assertEqual(len(cross_val_data), num_cv)

        # Then confirm shapes
        self.assertEqual(training_data["training_set_1"].shape[0], expected_training_size)
        self.assertEqual(training_data["training_set_2"].shape[0], expected_training_size)
        self.assertEqual(training_data["training_set_3"].shape[0], expected_training_size)
        self.assertEqual(cross_val_data["test_set_1"].shape[0], expected_test_size)
        self.assertEqual(cross_val_data["test_set_2"].shape[0], expected_test_size)
        self.assertEqual(cross_val_data["test_set_3"].shape[0], expected_test_size)

        # Check that the training and cross-validation sets are different in each case
        with self.assertRaises(AssertionError):
            np.testing.assert_equal(training_data["training_set_1"], training_data["training_set_2"])
        with self.assertRaises(AssertionError):
            np.testing.assert_equal(training_data["training_set_2"], training_data["training_set_3"])
        with self.assertRaises(AssertionError):
            np.testing.assert_equal(training_data["training_set_1"], training_data["training_set_3"])
        with self.assertRaises(AssertionError):
            np.testing.assert_equal(cross_val_data["test_set_1"], cross_val_data["test_set_2"])
        with self.assertRaises(AssertionError):
            np.testing.assert_equal(cross_val_data["test_set_1"], cross_val_data["test_set_3"])
        with self.assertRaises(AssertionError):
            np.testing.assert_equal(cross_val_data["test_set_2"], cross_val_data["test_set_3"])

        # Finally, check that concatenation generates the original sample data.
        # ## First, re-combine the training and cross-validation datasets for each cross-validation case
        concat_01 = (np.concatenate((training_data["training_set_1"], cross_val_data["test_set_1"]), axis=0))
        concat_02 = (np.concatenate((training_data["training_set_2"], cross_val_data["test_set_2"]), axis=0))
        concat_03 = (np.concatenate((training_data["training_set_3"], cross_val_data["test_set_3"]), axis=0))
        # ## Next, sort the sample and concatenated data in the same order. np.lexsort allows the specification of column order during sorting
        sample_data_sorted = self.training_data[
            np.lexsort((self.training_data[:, 2], self.training_data[:, 1], self.training_data[:, 0]))]
        concat_01_sorted = concat_01[np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))]
        concat_02_sorted = concat_02[np.lexsort((concat_02[:, 2], concat_02[:, 1], concat_02[:, 0]))]
        concat_03_sorted = concat_03[np.lexsort((concat_03[:, 2], concat_03[:, 1], concat_03[:, 0]))]
        # ## Finally, assert that the re-formed arrays are the same as the sorted input data array
        np.testing.assert_equal(sample_data_sorted, concat_01_sorted)
        np.testing.assert_equal(sample_data_sorted, concat_02_sorted)
        np.testing.assert_equal(sample_data_sorted, concat_03_sorted)

    def test_training_test_data_creation_03(self):
        num_cv = 3
        split = 0.02
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        with self.assertRaises(Exception):
            data_feed.training_test_data_creation()

    def test_training_test_data_creation_04(self):
        num_cv = 3
        split = 0.9999
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        with self.assertRaises(Exception):
            data_feed.training_test_data_creation()

    def test_training_test_data_creation_05(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)

        # Create two additional columns
        additional_data = np.zeros((self.training_data.shape[0], 2))
        additional_data[:, 0] = np.sin(self.training_data[:, 0])
        additional_data[:, 1] = np.sin(self.training_data[:, 1])
        training_data, cross_val_data = data_feed.training_test_data_creation(additional_features=additional_data)

        expected_training_size = int(np.around(self.training_data.shape[0] * split))
        expected_test_size = self.training_data.shape[0] - expected_training_size

        # Confirm number of created training and test data sets in dictionaries = 2 * no cross-validations
        self.assertEqual(len(training_data), num_cv * 2)  # Training and additional data for each CV
        self.assertEqual(len(cross_val_data), num_cv * 2)

        # Then confirm shapes
        self.assertEqual(training_data["training_set_1"].shape, (expected_training_size, self.training_data.shape[1]))
        self.assertEqual(training_data["training_set_2"].shape, (expected_training_size, self.training_data.shape[1]))
        self.assertEqual(training_data["training_set_3"].shape, (expected_training_size, self.training_data.shape[1]))
        self.assertEqual(cross_val_data["test_set_1"].shape, (expected_test_size, self.training_data.shape[1]))
        self.assertEqual(cross_val_data["test_set_2"].shape, (expected_test_size, self.training_data.shape[1]))
        self.assertEqual(cross_val_data["test_set_3"].shape, (expected_test_size, self.training_data.shape[1]))
        self.assertEqual(training_data["training_extras_1"].shape, (expected_training_size, additional_data.shape[1]))
        self.assertEqual(training_data["training_extras_2"].shape, (expected_training_size, additional_data.shape[1]))
        self.assertEqual(training_data["training_extras_3"].shape, (expected_training_size, additional_data.shape[1]))
        self.assertEqual(cross_val_data["test_extras_1"].shape, (expected_test_size, additional_data.shape[1]))
        self.assertEqual(cross_val_data["test_extras_2"].shape, (expected_test_size, additional_data.shape[1]))
        self.assertEqual(cross_val_data["test_extras_3"].shape, (expected_test_size, additional_data.shape[1]))

        # Finally, check that concatenation generates the original sample data and additional_features_array.
        # ## First, re-combine the training and cross-validation datasets for each cross-validation case
        concat_01 = (np.concatenate((training_data["training_set_1"], cross_val_data["test_set_1"]), axis=0))
        concat_02 = (np.concatenate((training_data["training_set_2"], cross_val_data["test_set_2"]), axis=0))
        concat_03 = (np.concatenate((training_data["training_set_3"], cross_val_data["test_set_3"]), axis=0))
        # ## Next, sort the sample and concatenated data in the same order. np.lexsort allows the specification of column order during sorting
        sample_data_sorted = self.training_data[
            np.lexsort((self.training_data[:, 2], self.training_data[:, 1], self.training_data[:, 0]))]
        concat_01_sorted = concat_01[np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))]
        concat_02_sorted = concat_02[np.lexsort((concat_02[:, 2], concat_02[:, 1], concat_02[:, 0]))]
        concat_03_sorted = concat_03[np.lexsort((concat_03[:, 2], concat_03[:, 1], concat_03[:, 0]))]
        # ## Finally, assert that the re-formed arrays are the same as the sorted input data array
        np.testing.assert_equal(sample_data_sorted, concat_01_sorted)
        np.testing.assert_equal(sample_data_sorted, concat_02_sorted)
        np.testing.assert_equal(sample_data_sorted, concat_03_sorted)

        # ## First, re-combine the training and cross-validation datasets for each cross-validation case
        concat_04 = (np.concatenate((training_data["training_extras_1"], cross_val_data["test_extras_1"]), axis=0))
        concat_05 = (np.concatenate((training_data["training_extras_2"], cross_val_data["test_extras_2"]), axis=0))
        concat_06 = (np.concatenate((training_data["training_extras_3"], cross_val_data["test_extras_3"]), axis=0))
        # ## Next, sort the sample and concatenated data in the same order. np.lexsort allows the specification of column order during sorting
        extra_data_sorted = additional_data[np.lexsort((additional_data[:, 1], additional_data[:, 0]))]
        concat_04_sorted = concat_04[np.lexsort((concat_04[:, 1], concat_04[:, 0]))]
        concat_05_sorted = concat_05[np.lexsort((concat_05[:, 1], concat_05[:, 0]))]
        concat_06_sorted = concat_06[np.lexsort((concat_06[:, 1], concat_06[:, 0]))]
        # ## Finally, assert that the re-formed arrays are the same as the sorted input data array
        np.testing.assert_equal(extra_data_sorted, concat_04_sorted)
        np.testing.assert_equal(extra_data_sorted, concat_05_sorted)
        np.testing.assert_equal(extra_data_sorted, concat_06_sorted)

    """
       Tests: : Unit tests for user_defined_terms function. For a sample dataset, we verify that:
           - 1&2: The function behaves properly when the right data types and shapes are entered.
           - 3: An exception is raised when the the number of samples in the term_list array is different from the regression_data length.
           - 4: An exception is raised when list element is not 1-dimensional
           - 5. An exception is raised when the term list (additional_terms) is not a list
    """

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
        with self.assertRaises(Exception):
            data_feed.user_defined_terms(additional_terms)

    def test_user_defined_terms_04(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        p = np.sin(self.training_data[:, 0]).reshape(self.training_data.shape[0], 1)  # 2-D array, should raise error
        additional_terms = [p]
        with self.assertRaises(Exception):
            data_feed.user_defined_terms(additional_terms)

    def test_user_defined_terms_05(self):
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=2,
                                         training_split=split, number_of_crossvalidations=num_cv)
        p = np.sin(self.training_data[:, 0]).reshape(self.training_data.shape[0], 1)
        additional_terms = p  # Additional terms as array, not list
        with self.assertRaises(ValueError):
            data_feed.user_defined_terms(additional_terms)

    """
       Tests: : Unit tests for results_generation. For specified solution arrays and orders, we verify that:
           - 1: The index and correct values are returned for a polynomial order of 1 and no multinomials
           - 2: The index and correct values are returned for a polynomial order of 3 and no multinomials
           - 3: The index and correct values are returned for a polynomial order of 2 and multinomials present
    """

    def test_results_generation_01(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=3, multinomials=0)
        order = 1
        beta = np.array([ [0], [0], [0] ])
        expected_df = pd.Series()
        row_list = np.array([['k'], ['(x_1)^1'], ['(x_2)^1']])
        expected_df = expected_df.append(pd.Series({row_list[0, 0]: beta[0, 0], row_list[1, 0]: beta[1, 0], row_list[2, 0]: beta[2, 0]}))
        output_df = data_feed.results_generation(beta, order)
        self.assertEqual(output_df.index.to_list(), expected_df.index.to_list())
        self.assertEqual(expected_df.all(), output_df.all())

    def test_results_generation_02(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=3, multinomials=0)
        order = 3
        beta = np.array([ [1], [0.3], [6], [500], [500000], [0.001], [50] ])
        expected_df = pd.Series()
        row_list = np.array([['k'], ['(x_1)^1'], ['(x_2)^1'], ['(x_1)^2'], ['(x_2)^2'], ['(x_1)^3'], ['(x_2)^3']])
        expected_df = expected_df.append(pd.Series({row_list[0, 0]: beta[0, 0], row_list[1, 0]: beta[1, 0], row_list[2, 0]: beta[2, 0], row_list[3, 0]: beta[3, 0], row_list[4, 0]: beta[4, 0], row_list[5, 0]: beta[5, 0], row_list[6, 0]: beta[6, 0]}))
        output_df = data_feed.results_generation(beta, order)
        self.assertEqual(output_df.index.to_list(), expected_df.index.to_list())
        self.assertEqual(expected_df.all(), output_df.all())

    def test_results_generation_03(self):
        data_feed = PolynomialRegression(self.full_data, self.training_data, maximum_polynomial_order=3, multinomials=1)
        order = 2
        beta = np.array([ [1], [0.3], [6], [500], [500000], [0.001]])
        expected_df = pd.Series()
        row_list = np.array([['k'], ['(x_1)^1'], ['(x_2)^1'], ['(x_1)^2'], ['(x_2)^2'], ['(x_1).(x_2)']])
        expected_df = expected_df.append(pd.Series({row_list[0, 0]: beta[0, 0], row_list[1, 0]: beta[1, 0], row_list[2, 0]: beta[2, 0], row_list[3, 0]: beta[3, 0], row_list[4, 0]: beta[4, 0], row_list[5, 0]: beta[5, 0]}))
        output_df = data_feed.results_generation(beta, order)
        self.assertEqual(output_df.index.to_list(), expected_df.index.to_list())
        self.assertEqual(expected_df.all(), output_df.all())

    """
       Tests: : Unit tests for get_feature_vector. We verify that:
           - 1: The (key, val) dictionary obtained from the generated IndexParam matches the headers of the data when the input is a dataframe
           - 2: The (key, val) dictionary obtained from the generated IndexParam matches is numerically consistent with that of the input array when a numpy array is supplied. 
    """

    def test_get_feature_vector_01(self):
        data_feed = PolynomialRegression(self.full_data, self.full_data, maximum_polynomial_order=3, multinomials=0)
        output = data_feed.get_feature_vector()
        expected_dict =  {'x1': 0, 'x2': 0}
        self.assertDictEqual(expected_dict, output.extract_values())
        print()

    def test_get_feature_vector_02(self):
        data_feed = PolynomialRegression(self.training_data, self.training_data, maximum_polynomial_order=3, multinomials=0)
        output = data_feed.get_feature_vector()
        expected_dict =  {0: 0, 1: 0}
        self.assertDictEqual(expected_dict, output.extract_values())
        print()

    """   
       Tests: : Unit tests for the polyregression function. We verify that:
           - 1: The correct optimization method is called and the correct optimization solution returned when 'mle' is selected as the solution_method.
           - 2: The correct optimization method is called and the correct optimization solution returned when 'bfgs' is selected as the solution_method.
           - 3: The correct optimization method is called and the correct optimization solution returned when 'pyomo' is selected as the solution_method.
           - 3: No optimization problem is solved when the problem is underspecified - all the outputs default to np.Inf.
           For the tests, the relevant optimization method in each case was replaced by a mock function.
    """

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

if __name__ == '__main__':
    unittest.main()



