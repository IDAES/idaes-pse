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

from idaes.core.surrogate.pysmo.sampling import (
    LatinHypercubeSampling,
    UniformSampling,
    HammersleySampling,
    CVTSampling,
    SamplingMethods,
    FeatureScaling,
)
import numpy as np
import pandas as pd
import pyomo.common.unittest as unittest
import pytest


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
        input_array_np_1d = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        input_array_np_2d = np.array(
            [
                [0, 1],
                [1, 4],
                [2, 9],
                [3, 16],
                [4, 25],
                [5, 36],
                [6, 49],
                [7, 64],
                [8, 81],
                [9, 100],
            ]
        )
        input_array_np_3d = np.array(
            [
                [0, 10, 11],
                [1, 11, 15],
                [2, 12, 21],
                [3, 13, 29],
                [4, 14, 39],
                [5, 15, 51],
                [6, 16, 65],
                [7, 17, 81],
                [8, 18, 99],
                [9, 19, 119],
            ]
        )
        input_array_np_3d_constant = np.array(
            [
                [0, 10, 11],
                [1, 10, 14],
                [2, 10, 19],
                [3, 10, 26],
                [4, 10, 35],
                [5, 10, 46],
                [6, 10, 59],
                [7, 10, 74],
                [8, 10, 91],
                [9, 10, 110],
            ]
        )

        input_array_pd_1d = pd.DataFrame(
            [[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]]
        )
        input_array_pd_2d = pd.DataFrame(
            [
                [0, 1],
                [1, 4],
                [2, 9],
                [3, 16],
                [4, 25],
                [5, 36],
                [6, 49],
                [7, 64],
                [8, 81],
                [9, 100],
            ]
        )
        input_array_pd_3d = pd.DataFrame(
            [
                [0, 10, 11],
                [1, 11, 15],
                [2, 12, 21],
                [3, 13, 29],
                [4, 14, 39],
                [5, 15, 51],
                [6, 16, 65],
                [7, 17, 81],
                [8, 18, 99],
                [9, 19, 119],
            ]
        )
        input_array_pd_3d_constant = pd.DataFrame(
            [
                [0, 10, 11],
                [1, 10, 14],
                [2, 10, 19],
                [3, 10, 26],
                [4, 10, 35],
                [5, 10, 46],
                [6, 10, 59],
                [7, 10, 74],
                [8, 10, 91],
                [9, 10, 110],
            ]
        )

        self.test_data_numpy_1d = input_array_np_1d
        self.test_data_numpy_2d = input_array_np_2d
        self.test_data_numpy_3d = input_array_np_3d
        self.test_data_numpy_3d_constant = input_array_np_3d_constant

        self.test_data_pandas_1d = input_array_pd_1d
        self.test_data_pandas_2d = input_array_pd_2d
        self.test_data_pandas_3d = input_array_pd_3d
        self.test_data_pandas_3d_constant = input_array_pd_3d_constant

    @pytest.mark.unit
    def test_data_scaling_minmax_01(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_numpy_1d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = (input_array - expected_output_2) / (
            expected_output_3 - expected_output_2
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
    def test_data_scaling_minmax_02(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 100]])
        expected_output_2 = np.array([[0, 1]])
        expected_output_1 = (input_array - expected_output_2) / (
            expected_output_3 - expected_output_2
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
    def test_data_scaling_minmax_03(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_numpy_3d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 19, 119]])
        expected_output_2 = np.array([[0, 10, 11]])
        expected_output_1 = (input_array - expected_output_2) / (
            expected_output_3 - expected_output_2
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
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

    @pytest.mark.unit
    def test_data_scaling_minmax_05(self):
        # TypeError with list
        input_array = self.test_data_numpy_2d.tolist()
        with pytest.raises(TypeError):
            FeatureScaling.data_scaling_minmax(input_array)

    @pytest.mark.unit
    def test_data_scaling_minmax_06(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_pandas_1d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = (input_array - expected_output_2) / (
            expected_output_3 - expected_output_2
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
    def test_data_scaling_minmax_07(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_pandas_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 100]])
        expected_output_2 = np.array([[0, 1]])
        expected_output_1 = (input_array - expected_output_2) / (
            expected_output_3 - expected_output_2
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
    def test_data_scaling_minmax_08(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_pandas_3d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 19, 119]])
        expected_output_2 = np.array([[0, 10, 11]])
        expected_output_1 = (input_array - expected_output_2) / (
            expected_output_3 - expected_output_2
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
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

    @pytest.mark.unit
    def test_data_unscaling_minmax_01(self):
        # 1D
        # Sample data generated from between 0 and 9
        input_array = self.test_data_numpy_1d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_02(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_03(self):
        # 3D
        # Sample data generated from the expression (x_1 + 1)^2 + x_2, x1: between 0 and 9, x2: between 10 and 19
        input_array = self.test_data_numpy_3d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_04(self):
        # 3D with constant
        # Sample data generated from the expression (x_1 + 1)^2 + 10, x1: between 0 and 9
        input_array = self.test_data_numpy_3d_constant
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_05(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)

        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_06(self):
        # 2D
        # Sample data generated from the expression (x_1 + 1)^2 between 0 and 9
        input_array = self.test_data_numpy_2d
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)

        min_array = np.array([[1, 2, 3]])
        max_array = np.array([[5, 6, 7]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)


class SamplingMethodsTestCases(unittest.TestCase):
    def setUp(self):
        input_array_np_1d = np.array([[0], [1], [2], [3], [4], [5], [6], [7], [8], [9]])
        input_array_np_2d = np.array(
            [
                [0, 1],
                [1, 4],
                [2, 9],
                [3, 16],
                [4, 25],
                [5, 36],
                [6, 49],
                [7, 64],
                [8, 81],
                [9, 100],
            ]
        )
        input_array_np_3d = np.array(
            [
                [0, 10, 11],
                [1, 11, 15],
                [2, 12, 21],
                [3, 13, 29],
                [4, 14, 39],
                [5, 15, 51],
                [6, 16, 65],
                [7, 17, 81],
                [8, 18, 99],
                [9, 19, 119],
            ]
        )
        input_array_pd = pd.DataFrame(
            {
                "x1": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "x2": [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
                "y1": [11, 14, 19, 26, 35, 46, 59, 74, 91, 110],
                "y2": [21, 24, 29, 36, 45, 56, 69, 84, 101, 120],
            }
        )

        self.test_data_numpy_1d = input_array_np_1d
        self.test_data_numpy_2d = input_array_np_2d
        self.test_data_numpy_3d = input_array_np_3d
        self.test_data_pandas = input_array_pd

    """
   
    test_nearest_neighbour_01: Test behaviour with array N x d = (10,3) data and a=[-0.5,1]
    test_nearest_neighbour_02: Test behaviour with array N x d = (10,2) data and a=[-0.5]
    test_nearest_neighbour_03: Test behaviour with array N x d = (10,1) data and a=[]
    test_nearest_neighbour_04: working even with diffent input size
    test_nearest_neighbour_05: Test behaviour raise ValueError if dimension of point is not matching with array N x d = (10,3) 
    
    test_points_selection_01: Test behaviour with array N x d = (10,3) data and a=[[-0.5,10],[10,100]]
    test_points_selection_02: Test behaviour with array N x d = (10,2) data and a=[[-0.5],[10]]
    test_points_selection_03: Test behaviour with array N x d = (10,1) data and a=[[],[]] <---both return the first element
    test_points_selection_04: Test behaviour raise ValueError if dimension of point is not matching with array N x d = (10,3); small
    test_points_selection_05: Test behaviour raise ValueError if dimension of point is not matching with array N x d = (10,3); large
    
    test_sample_point_selection_01: selection - Test behaviour with array N x d = (10,3) data and a=[[0,0],[10,19]]
    test_sample_point_selection_02: selection - Test behaviour with array N x d = (10,2) data and a=[[0],[7]]
    test_sample_point_selection_03: selection - Test behaviour with array N x d = (10,1) data and a=[[],[]] <---both return the first element, so return only 1 []
    test_sample_point_selection_04: selection - Test behaviour raise ValueError if dimension of point is not matching with array N x d = (10,3); large
    test_sample_point_selection_05: selection - Test behaviour raise ValueError if dimension of point is not matching with array N x d = (10,3); large
    
    test_points_selection_06: creation - Test behaviour with points dimension should be d with array N x d = (10,3) 
    test_points_selection_07: creation - Test behaviour with points dimension should be d with array N x d = (10,2) 
    test_points_selection_08: creation - Test behaviour with points dimension should be d with array N x d = (10,1) 
    test_points_selection_09: creation - raise IndexError if dimension of point is not matching with array N x d = (10,3); small
    test_points_selection_10: creation - raise IndexError if dimension of point is not matching with array N x d = (10,3); large
    
    test_points_selection_01: Test behaviour with n = 3
    test_points_selection_02: Test behaviour with n = 1
    test_points_selection_03: Test behaviour with n = 0
    test_points_selection_04: Test behaviour with n = -1
    test_points_selection_04: Test behaviour with n = 2.9
    
    test_base_conversion_01: Test behaviour with 5 to base 2
    test_base_conversion_02: Test behaviour with 57 to base 47
    test_base_conversion_03: Test behaviour with negative base - works, returns always 0
    test_base_conversion_04: Test behaviour raise ZeroDivisionError with 0 base
    test_base_conversion_05: Test behaviour with 1 base, infinity loop
    
    test_prime_base_to_decimal_01: Test behaviour with 0.01 in base 2 to base 10
    test_prime_base_to_decimal_02: Test behaviour with 0.01 in base 20 to base 10
    test_prime_base_to_decimal_03: working with base 1
    test_prime_base_to_decimal_04: working with base -1

    test_selection_columns_preprocessing_01 - 11: Test behaviour of function for a dataframe input under a variety of xlabels and ylabels conditions. 
        See test docstrings for more specific details.
    test_selection_columns_preprocessing_14 - 26: Repeats above tests 1-13 for case where input is a numpy array rather than a dataframe
    """

    def _create_sampling(self, input_array, sample_points):
        sampling = SamplingMethods()
        sampling.data_headers = [i for i in range(0, sample_points.shape[1] + 1)]
        sampling.x_data = np.zeros((2, input_array.shape[1] - 1))
        return sampling

    @pytest.mark.unit
    def test_nearest_neighbour_01(self):
        input_array = self.test_data_numpy_3d
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [-0.5, 1])
        np.testing.assert_array_equal(closest_point, input_array[0, :])

    @pytest.mark.unit
    def test_nearest_neighbour_02(self):
        input_array = self.test_data_numpy_2d
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [-0.5])
        np.testing.assert_array_equal(closest_point, input_array[0, :])

    @pytest.mark.unit
    def test_nearest_neighbour_03(self):
        input_array = self.test_data_numpy_1d
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [])
        np.testing.assert_array_equal(closest_point, input_array[0, :])

    @pytest.mark.unit
    def test_nearest_neighbour_04(self):
        input_array = self.test_data_numpy_3d
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [0.5])

    @pytest.mark.unit
    def test_nearest_neighbour_05(self):
        input_array = self.test_data_numpy_3d
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        with pytest.raises(ValueError):
            closest_point = SamplingMethods.nearest_neighbour(
                self, input_array, [0.5, 0.9, 10]
            )

    @pytest.mark.unit
    def test_points_selection_01(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[-0.5, 10], [10, 100]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        equivalent_points = sampling_methods.points_selection(
            input_array, generated_sample_points
        )
        np.testing.assert_array_equal(equivalent_points[0], input_array[0, :])
        np.testing.assert_array_equal(equivalent_points[1], input_array[-1, :])

    @pytest.mark.unit
    def test_points_selection_02(self):
        input_array = self.test_data_numpy_2d
        generated_sample_points = np.array([[-0.5], [10]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        equivalent_points = sampling_methods.points_selection(
            input_array, generated_sample_points
        )
        np.testing.assert_array_equal(equivalent_points[0], input_array[0, :])
        np.testing.assert_array_equal(equivalent_points[1], input_array[-1, :])

    @pytest.mark.unit
    def test_points_selection_03(self):
        input_array = self.test_data_numpy_1d
        generated_sample_points = np.array([[], []])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        equivalent_points = sampling_methods.points_selection(
            input_array, generated_sample_points
        )
        np.testing.assert_array_equal(equivalent_points[0], input_array[0, :])
        np.testing.assert_array_equal(equivalent_points[1], input_array[0, :])

    @pytest.mark.unit
    def test_points_selection_04(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[0.5], [10]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            equivalent_points = sampling_methods.points_selection(
                input_array, generated_sample_points
            )

    @pytest.mark.unit
    def test_points_selection_05(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[0.5, 0.7, 10], [10, 0.9, 20]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            equivalent_points = sampling_methods.points_selection(
                input_array, generated_sample_points
            )

    @pytest.mark.unit
    def test_sample_point_selection_01(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[0, 0], [10, 19]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="selection"
        )
        np.testing.assert_array_equal(unique_sample_points[0], input_array[0, :])
        np.testing.assert_array_equal(unique_sample_points[1], input_array[-1, :])

    @pytest.mark.unit
    def test_sample_point_selection_02(self):
        input_array = self.test_data_numpy_2d
        generated_sample_points = np.array([[0], [7]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="selection"
        )
        np.testing.assert_array_equal(unique_sample_points[0], input_array[0, :])
        np.testing.assert_array_equal(unique_sample_points[1], input_array[-1, :])

    @pytest.mark.unit
    def test_sample_point_selection_03(self):
        input_array = self.test_data_numpy_1d
        generated_sample_points = np.array([[], []])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="selection"
        )
        np.testing.assert_array_equal(unique_sample_points[0], input_array[0, :])

    @pytest.mark.unit
    def test_sample_point_selection_04(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[0.5], [7]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            unique_sample_points = sampling_methods.sample_point_selection(
                input_array, generated_sample_points, sampling_type="selection"
            )

    @pytest.mark.unit
    def test_sample_point_selection_05(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[0.5, 1, 10], [7, 19, 20]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            unique_sample_points = sampling_methods.sample_point_selection(
                input_array, generated_sample_points, sampling_type="selection"
            )

    @pytest.mark.unit
    def test_sample_point_selection_06(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[0.5, 11, 3], [7, 19, 4]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="creation"
        )
        min_, max_ = input_array[0, :], input_array[1, :]
        testing = min_ + generated_sample_points * (max_ - min_)
        np.testing.assert_array_equal(testing, unique_sample_points)

    @pytest.mark.unit
    def test_sample_point_selection_07(self):
        input_array = self.test_data_numpy_2d
        generated_sample_points = np.array([[0.5, 1], [7, 19]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="creation"
        )
        min_, max_ = input_array[0, :], input_array[1, :]
        testing = min_ + generated_sample_points * (max_ - min_)
        np.testing.assert_array_equal(testing, unique_sample_points)

    @pytest.mark.unit
    def test_sample_point_selection_08(self):
        input_array = self.test_data_numpy_1d
        generated_sample_points = np.array([[0.5], [7]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        sampling_methods.x_data = np.zeros((2, input_array.shape[1] - 1))
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="creation"
        )
        min_, max_ = input_array[0, :], input_array[1, :]
        testing = min_ + generated_sample_points * (max_ - min_)
        np.testing.assert_array_equal(testing, unique_sample_points)

    @pytest.mark.unit
    def test_sample_point_selection_09(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[], []])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(IndexError):
            unique_sample_points = sampling_methods.sample_point_selection(
                input_array, generated_sample_points, sampling_type="creation"
            )

    @pytest.mark.unit
    def test_sample_point_selection_10(self):
        input_array = self.test_data_numpy_3d
        generated_sample_points = np.array([[0.5, 1, 10, 11], [7, 19, 10, 12]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(IndexError):
            unique_sample_points = sampling_methods.sample_point_selection(
                input_array, generated_sample_points, sampling_type="creation"
            )

    @pytest.mark.unit
    def test_prime_number_generator_01(self):
        prime_list = SamplingMethods.prime_number_generator(self, 3)
        np.testing.assert_array_equal(prime_list, [2, 3, 5])

    @pytest.mark.unit
    def test_prime_number_generator_02(self):
        prime_list = SamplingMethods.prime_number_generator(self, 1)
        np.testing.assert_array_equal(prime_list, [2])

    @pytest.mark.unit
    def test_prime_number_generator_03(self):
        prime_list = SamplingMethods.prime_number_generator(self, 0)
        np.testing.assert_array_equal(prime_list, [])

    @pytest.mark.unit
    def test_prime_number_generator_04(self):
        prime_list = SamplingMethods.prime_number_generator(self, -1)
        np.testing.assert_array_equal(prime_list, [])

    @pytest.mark.unit
    def test_prime_number_generator_05(self):
        prime_list = SamplingMethods.prime_number_generator(self, 2.9)
        np.testing.assert_array_equal(prime_list, [2, 3, 5])

    @pytest.mark.unit
    def test_base_conversion_01(self):
        string_representation = SamplingMethods.base_conversion(self, 5, 2)
        assert string_representation == ["1", "0", "1"]

    @pytest.mark.unit
    def test_base_conversion_02(self):
        string_representation = SamplingMethods.base_conversion(self, 57, 47)
        assert string_representation == ["1", "10"]

    @pytest.mark.unit
    def test_base_conversion_03(self):
        string_representation = SamplingMethods.base_conversion(self, 10, -1)
        assert string_representation == ["0"]

    @pytest.mark.unit
    def test_base_conversion_04(self):
        with pytest.raises(ZeroDivisionError):
            string_representation = SamplingMethods.base_conversion(self, 10, 0)

    # @pytest.mark.unit
    # def test_base_conversion_05(self):
    #     string_representation = SamplingMethods.base_conversion(self, 10, 1)

    @pytest.mark.unit
    def test_prime_base_to_decimal_01(self):
        string_representation = SamplingMethods.prime_base_to_decimal(
            self, ["0", "0", "1"], 2
        )
        assert 0.25 == string_representation

    @pytest.mark.unit
    def test_prime_base_to_decimal_02(self):
        string_representation = SamplingMethods.prime_base_to_decimal(
            self, ["0", "0", "1"], 20
        )
        assert 0.0025 == string_representation

    @pytest.mark.unit
    def test_prime_base_to_decimal_03(self):
        string_representation = SamplingMethods.prime_base_to_decimal(
            self, ["0", "0", "1"], 1
        )

    @pytest.mark.unit
    def test_prime_base_to_decimal_04(self):
        string_representation = SamplingMethods.prime_base_to_decimal(
            self, ["0", "0", "1"], -1
        )

    @pytest.mark.unit
    def test_data_sequencing(self):
        sampling_methods = SamplingMethods()
        sequence_decimal = sampling_methods.data_sequencing(3, 2)

    @pytest.mark.unit
    def test_selection_columns_preprocessing_01(self):
        """
        Test behaviour with no labels specified
        """
        sampling_methods = SamplingMethods()
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, None, None
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - 1
        assert (
            sampling_methods.data_headers
            == self.test_data_pandas.columns.values.tolist()
        )
        assert (
            sampling_methods.data_headers_xvars
            == self.test_data_pandas.columns.values.tolist()[:-1]
        )
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, :-1]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_02(self):
        """
        Test behaviour with only xlabels specified with multiple elements, same order as in dataframe
        """
        sampling_methods = SamplingMethods()
        x_lab = ["x1", "x2"]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, None
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert (
            sampling_methods.data_headers
            == self.test_data_pandas.columns.values.tolist()
        )
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, :-2]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_03(self):
        """
        Test behaviour with only xlabels specified with multiple elements, different order from dataframe
        """
        sampling_methods = SamplingMethods()
        x_lab = ["y1", "x1"]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, None
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == ["y1", "x1", "x2", "y2"]
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data,
            self.test_data_pandas[["y1", "x1", "x2", "y2"]].values,
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_04(self):
        """
        Test behaviour with only xlabels specified with only 1 element
        """
        sampling_methods = SamplingMethods()
        x_lab = ["y2"]
        expected_order = ["y2", "x1", "x2", "y1"]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, None
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_05(self):
        """
        Test behaviour with only ylabels specified with multiple elements, same order as in dataframe
        """
        sampling_methods = SamplingMethods()
        y_lab = ["y1", "y2"]
        expected_order = ["x1", "x2", "y1", "y2"]
        expected_x_lab = ["x1", "x2"]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, None, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - len(
            y_lab
        )
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == expected_x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[expected_x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_06(self):
        """
        Test behaviour with only ylabels specified with multiple elements, different order from dataframe
        """
        sampling_methods = SamplingMethods()
        y_lab = ["y1", "x1"]
        expected_order = ["x2", "y2", "y1", "x1"]
        expected_x_lab = ["x2", "y2"]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, None, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - len(
            y_lab
        )
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == expected_x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[expected_x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_07(self):
        """
        Test behaviour with only ylabels specified with only 1 element
        """
        sampling_methods = SamplingMethods()
        y_lab = ["x2"]
        expected_order = ["x1", "y1", "y2", "x2"]
        expected_x_lab = ["x1", "y1", "y2"]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, None, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - len(
            y_lab
        )
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == expected_x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[expected_x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_11(self):
        """
        Test behaviour with both xlabels and ylabels specified; multiple x and single y
        """
        sampling_methods = SamplingMethods()
        x_lab = ["x1", "x2"]
        y_lab = ["y2"]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_08(self):
        """
        Test behaviour with both xlabels and ylabels specified; single x and single y; ordered
        """
        sampling_methods = SamplingMethods()
        x_lab = ["x1"]
        y_lab = ["y1"]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == len(x_lab) + len(y_lab)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_09(self):
        """
        Test behaviour with both xlabels and ylabels specified; single x and multiple y
        """
        sampling_methods = SamplingMethods()
        x_lab = ["x1"]
        y_lab = ["y1", "y2"]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_10(self):
        """
        Test behaviour with both xlabels and ylabels specified; single x and single y; unordered
        """
        sampling_methods = SamplingMethods()
        x_lab = ["y1"]
        y_lab = ["x2"]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_12(self):
        """
        Test behaviour with both xlabels and ylabels specified; multiple x and y, all initial columns present
        """
        sampling_methods = SamplingMethods()
        x_lab = ["y1", "x2"]
        y_lab = ["x1", "y2"]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas, x_lab, y_lab
        )
        assert sampling_methods.df_flag == True
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas[expected_order].values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas[x_lab].values
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_13(self):
        """
        Test behaviour when non-existent column name is supplied in label - should raise ValueError
        """
        sampling_methods = SamplingMethods()
        x_lab = ["x1"]
        y_lab = ["y3"]
        expected_order = x_lab + y_lab

        with pytest.raises(IndexError):
            sampling_methods.selection_columns_preprocessing(
                self.test_data_pandas, x_lab, y_lab
            )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_14(self):
        sampling_methods = SamplingMethods()
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, None, None
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - 1
        assert sampling_methods.data_headers == [
            i for i in range(0, self.test_data_pandas.shape[1])
        ]
        assert sampling_methods.data_headers_xvars == [
            i for i in range(0, self.test_data_pandas.shape[1] - 1)
        ]
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, :-1]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_15(self):
        sampling_methods = SamplingMethods()
        x_lab = [0, 1]
        expected_order = [0, 1, 2, 3]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, None
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_16(self):
        sampling_methods = SamplingMethods()
        x_lab = [2, 0]
        expected_order = [2, 0, 1, 3]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, None
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_17(self):
        sampling_methods = SamplingMethods()
        x_lab = [3]
        expected_order = [3, 0, 1, 2]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, None
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_18(self):
        sampling_methods = SamplingMethods()
        y_lab = [2, 3]
        expected_order = [0, 1, 2, 3]
        expected_x_lab = [0, 1]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, None, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - len(
            y_lab
        )
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == expected_x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, expected_x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_19(self):
        sampling_methods = SamplingMethods()
        y_lab = [2, 0]
        expected_order = [1, 3, 2, 0]
        expected_x_lab = [1, 3]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, None, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - len(
            y_lab
        )
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == expected_x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, expected_x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_20(self):
        sampling_methods = SamplingMethods()
        y_lab = [1]
        expected_order = [0, 2, 3, 1]
        expected_x_lab = [0, 2, 3]
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, None, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == self.test_data_pandas.shape[1]
        assert sampling_methods.x_data.shape[1] == self.test_data_pandas.shape[1] - len(
            y_lab
        )
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == expected_x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, expected_x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_21(self):
        sampling_methods = SamplingMethods()
        x_lab = [0, 1]
        y_lab = [3]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_22(self):
        sampling_methods = SamplingMethods()
        x_lab = [0]
        y_lab = [2]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == len(x_lab) + len(y_lab)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_23(self):
        sampling_methods = SamplingMethods()
        x_lab = [0]
        y_lab = [2, 3]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_24(self):
        sampling_methods = SamplingMethods()
        x_lab = [2]
        y_lab = [1]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_25(self):
        sampling_methods = SamplingMethods()
        x_lab = [2, 1]
        y_lab = [0, 3]
        expected_order = x_lab + y_lab
        sampling_methods.selection_columns_preprocessing(
            self.test_data_pandas.values, x_lab, y_lab
        )
        assert sampling_methods.df_flag == False
        assert sampling_methods.data.shape[1] == len(expected_order)
        assert sampling_methods.x_data.shape[1] == len(x_lab)
        assert sampling_methods.data_headers == expected_order
        assert sampling_methods.data_headers_xvars == x_lab
        np.testing.assert_array_equal(
            sampling_methods.data, self.test_data_pandas.values[:, expected_order]
        )
        np.testing.assert_array_equal(
            sampling_methods.x_data, self.test_data_pandas.values[:, x_lab]
        )

    @pytest.mark.unit
    def test_selection_columns_preprocessing_26(self):
        sampling_methods = SamplingMethods()
        x_lab = [0]
        y_lab = [4]
        expected_order = x_lab + y_lab

        with pytest.raises(IndexError):
            sampling_methods.selection_columns_preprocessing(
                self.test_data_pandas.values, x_lab, y_lab
            )


class LatinHypercubeSamplingTestCases(unittest.TestCase):
    """
    test__init__selection_01: input numpy array - Test behaviour generate LatinHypercubeSampling object with selection, default number of sample = 5
    test__init__selection_02: input Pandas DataFrame - Test behaviour generate LatinHypercubeSampling object with selection, default number of sample
    test__init__selection_03: input numpy array - Test behaviour generate LatinHypercubeSampling object with selection, with a selected number of sample
    test__init__selection_04: input Pandas DataFrame - Test behaviour generate LatinHypercubeSampling object with selection, with a selected number of sample
    test__init__selection_05: input numpy array - Test behaviour raise exception with a selected number of sample = 0
    test__init__selection_06: input numpy array - Test behaviour raise exception with a selected number of sample = -1
    test__init__selection_07: input numpy array - Test behaviour raise exception with a selected number of sample > input size
    test__init__selection_08: input numpy array - Test behaviour raise exception with a selected number of sample = 1.1 (non-integer)
    test__init__selection_09: input list - Test behaviour raise ValueError with list input

    test__init__creation_01: input list - Test behaviour generate LatinHypercubeSampling object with default sampling_type, default number of sample = 5
    test__init__creation_02: input list - Test behaviour generate LatinHypercubeSampling object with sampling_type = creation, default number of sample = 5
    test__init__creation_03: input list - Test behaviour generate LatinHypercubeSampling object with creation with a selected number of sample
    test__init__creation_04: input list - Test behaviour raise exception with a selected number of sample = 0
    test__init__creation_05: input list - Test behaviour raise exception with a selected number of sample = -1
    test__init__creation_06: input list - Test behaviour raise exception with a selected number of sample = 1.1 (non-integer)
    test__init__creation_07: input numpy - Test behaviour raise ValueError with numpy input
    test__init__creation_08: input numpy - Test behaviour raise ValueError with pandas input

    test__init__creation_selection_01 - Test behaviour raise Exception with sampling_type = non string
    test__init__creation_selection_02 - Test behaviour raise Exception with sampling_type = incorrect string

    test_variable_sample_creation: Test behaviour, sampled values are in the range (min, max), number of samples = 5, 10, 1
    test_lhs_points_generation: Test behaviour, sampled values are in the range (0, 1) , number of samples = 5, 10, 1, 2 d
    test_random_shuffling: Test behaviour, random_shuffling = sampled values after soring, they are in the range (0, 1) , number of samples = 5, 10, 1, 2 d
    test_sample_points_01: Test behaviour with selection, sample points are unique, all in the input array , number of samples = 5, 10, 1, 2 d
    test_sample_points_02: Test behaviour with creation, sample points are unique, all in the input range (min,max) , number of samples = 5, 10, 1, 2 d

    """

    def setUp(self):
        input_array_np = np.array(
            [
                [0, 10, 11],
                [1, 10, 14],
                [2, 10, 19],
                [3, 10, 26],
                [4, 10, 35],
                [5, 10, 46],
                [6, 10, 59],
                [7, 10, 74],
                [8, 10, 91],
                [9, 10, 110],
            ]
        )
        input_array_pd = pd.DataFrame(
            {
                "x1": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "x2": [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
                "y": [11, 14, 19, 26, 35, 46, 59, 74, 91, 110],
            }
        )
        input_array_list = [[1, 10, 3], [2, 11, 4.5]]
        self.test_data_numpy = input_array_np
        self.test_data_pandas = input_array_pd
        self.test_data_list = input_array_list

    @pytest.mark.unit
    def test__init__selection_01(self):
        input_array = self.test_data_numpy
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)
        np.testing.assert_array_equal(LHSClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_02(self):
        input_array = self.test_data_pandas
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)
        input_array = input_array.to_numpy()
        np.testing.assert_array_equal(LHSClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_03(self):
        input_array = self.test_data_numpy
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=6, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 6)
        np.testing.assert_array_equal(LHSClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_04(self):
        input_array = self.test_data_pandas
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=6, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 6)
        input_array = input_array.to_numpy()
        np.testing.assert_array_equal(LHSClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_05(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=0, sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__selection_06(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=-1, sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__selection_07(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=101, sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__selection_08(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=1.1, sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__selection_09(self):
        input_array = self.test_data_list
        with pytest.raises(ValueError):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__creation_01(self):
        input_array = self.test_data_list
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type=None
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)

    @pytest.mark.unit
    def test__init__creation_02(self):
        input_array = self.test_data_list
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)

    @pytest.mark.unit
    def test__init__creation_03(self):
        input_array = self.test_data_list
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=100, sampling_type="creation"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 100)

    @pytest.mark.unit
    def test__init__creation_04(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=0, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_05(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=-1, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_06(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=1.1, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_07(self):
        input_array = self.test_data_numpy
        with pytest.raises(ValueError):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_08(self):
        input_array = self.test_data_pandas
        with pytest.raises(ValueError):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_selection_01(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type=1
            )

    @pytest.mark.unit
    def test__init__creation_selection_02(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="jp"
            )

    @pytest.mark.unit
    def test_variable_sample_creation(self):
        input_array = self.test_data_numpy
        for num_samples in [None, 10, 1]:
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            minimum, maximum = 10, 100
            out_var_samples = LHSClass.variable_sample_creation(minimum, maximum)
            assert (out_var_samples >= minimum).all() and (
                out_var_samples <= maximum
            ).all()
            np.testing.assert_array_equal(
                LHSClass.number_of_samples, out_var_samples.shape[0]
            )

    @pytest.mark.unit
    def test_lhs_points_generation(self):
        input_array = self.test_data_numpy
        for num_samples in [None, 10, 1]:
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            out_sample_points_vector = LHSClass.lhs_points_generation()
            assert (out_sample_points_vector >= 0).all() and (
                out_sample_points_vector <= 1
            ).all()
            np.testing.assert_array_equal(
                LHSClass.number_of_samples, out_sample_points_vector.shape[0]
            )
            np.testing.assert_array_equal(
                input_array.shape[1] - 1, out_sample_points_vector.shape[1]
            )

    @pytest.mark.unit
    def test_random_shuffling(self):
        input_array = self.test_data_numpy
        for num_samples in [None, 10, 1]:
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            out_sample_points_vector = LHSClass.lhs_points_generation()
            vector_of_points = LHSClass.random_shuffling(out_sample_points_vector)

            sidx1 = out_sample_points_vector.argsort(axis=0)
            out1 = out_sample_points_vector[sidx1, np.arange(sidx1.shape[1])]
            sidx2 = vector_of_points.argsort(axis=0)
            out2 = vector_of_points[sidx2, np.arange(sidx2.shape[1])]

            assert (out_sample_points_vector >= 0).all() and (
                out_sample_points_vector <= 1
            ).all()
            np.testing.assert_array_equal(out1, out2)
            np.testing.assert_array_equal(
                LHSClass.number_of_samples, out_sample_points_vector.shape[0]
            )
            np.testing.assert_array_equal(
                input_array.shape[1] - 1, out_sample_points_vector.shape[1]
            )

    @pytest.mark.unit
    def test_sample_points_01(self):
        for num_samples in [None, 10, 1]:
            input_array = self.test_data_numpy
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            unique_sample_points = LHSClass.sample_points()
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            out_testing = [
                unique_sample_points[i, :] in input_array
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)

    @pytest.mark.unit
    def test_sample_points_02(self):
        for num_samples in [None, 10, 1]:
            input_array = self.test_data_list
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=num_samples, sampling_type="creation"
            )
            unique_sample_points = LHSClass.sample_points()
            input_array = np.array(input_array)
            for i in range(input_array.shape[1]):
                var_range = input_array[:, i]
                assert (unique_sample_points[:, i] >= var_range[0]).all() and (
                    unique_sample_points[:, i] <= var_range[1]
                ).all()
            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0).shape,
                unique_sample_points.shape,
            )


class UniformSamplingTestCases(unittest.TestCase):
    """
    test__init__selection_01: input numpy array - Test behaviour generate UniformSampling object with selection, default edge, list_of_samples_per_variable = [2,5]
    test__init__selection_02: input Pandas DataFrame - Test behaviour generate UniformSampling object with selection, default edge, list_of_samples_per_variable = [2,5]
    test__init__selection_03: input numpy array - Test behaviour raise TypeError with a list_of_samples_per_variable is numpy, default edge
    test__init__selection_04: input numpy array - Test behaviour raise TypeError with a list_of_samples_per_variable is pandas, default edge
    test__init__selection_05: input numpy array - Test behaviour raise ValueError with a list_of_samples_per_variable < number of variables, default edge
    test__init__selection_06: input numpy array - Test behaviour raise ValueError with a list_of_samples_per_variable > number of variables, default edge
    test__init__selection_07: input numpy array - Test behaviour raise ValueError with a min(list_of_samples_per_variable) < 2, default edge
    test__init__selection_08: input numpy array - Test behaviour raise TypeError with a list_of_samples_per_variable is non integer, default edge
    test__init__selection_09: input numpy array - Test behaviour raise Exception with a list_of_samples_per_variable = [2,50], 2*50 > number of input data, default edge
    test__init__selection_10: input numpy array - Test behaviour generate UniformSampling object with selection, edge = True, list_of_samples_per_variable = [2,5]
    test__init__selection_11: input numpy array - Test behaviour generate UniformSampling object with selection, edge = True, list_of_samples_per_variable = [2,5]
    test__init__selection_12: input numpy array - Test behaviour raise Exception with edge = 1
    test__init__selection_13: input numpy array - Test behaviour raise Exception with edge = 'str'

    test__init__creation_01: input list - Test behaviour generate UniformSampling object with default sampling_type, default edge, list_of_samples_per_variable = [2,7,5]
    test__init__creation_02: input list - Test behaviour generate UniformSampling object with sampling_type = creation, default edge, list_of_samples_per_variable = [2,7,5]
    test__init__creation_03: input list - Test behaviour raise exception with a list_of_samples_per_variable = [1,7,5]
    test__init__creation_04: input list - Test behaviour raise exception with a list_of_samples_per_variable = [-1,7,5]
    test__init__creation_06: input list - Test behaviour raise exception with a list_of_samples_per_variable = [1.1,7,5] (non-integer)
    test__init__creation_07: input numpy - Test behaviour raise ValueError with numpy input
    test__init__creation_08: input Pandas - Test behaviour raise ValueError with Pandas input

    test__init__creation_selection_01 - Test behaviour raise Exception with sampling_type = non string
    test__init__creation_selection_02 - Test behaviour raise Exception with sampling_type = incorrect string

    test_sample_points_01: Test behaviour with selection, sample points are unique, all in the input array , number of samples = 5, 10, 1
    test_sample_points_02: Test behaviour with creation, sample points are unique, all in the input range (min,max) , number of samples = [2,5,9],[3,2,10],[4,2,28]
    """

    def setUp(self):
        input_array_np = np.array(
            [
                [0, 10, 11],
                [1, 10, 14],
                [2, 10, 19],
                [3, 10, 26],
                [4, 10, 35],
                [5, 10, 46],
                [6, 10, 59],
                [7, 10, 74],
                [8, 10, 91],
                [9, 10, 110],
            ]
        )
        input_array_pd = pd.DataFrame(
            {
                "x1": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "x2": [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
                "y": [11, 14, 19, 26, 35, 46, 59, 74, 91, 110],
            }
        )
        input_array_list = [[1, 10, 3], [2, 11, 4.5]]
        self.test_data_numpy = input_array_np
        self.test_data_pandas = input_array_pd
        self.test_data_list = input_array_list

    @pytest.mark.unit
    def test__init__selection_01(self):
        input_array = self.test_data_numpy
        UniClass = UniformSampling(input_array, [2, 5], sampling_type="selection")
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_02(self):
        input_array = self.test_data_pandas
        UniClass = UniformSampling(input_array, [2, 5], sampling_type="selection")
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        input_array = input_array.to_numpy()
        np.testing.assert_array_equal(UniClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_03(self):
        input_array = self.test_data_numpy
        with pytest.raises(TypeError):
            UniClass = UniformSampling(
                input_array, np.array([2, 5]), sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__selection_04(self):
        input_array = self.test_data_numpy
        with pytest.raises(TypeError):
            UniClass = UniformSampling(
                input_array, pd.DataFrame([2, 5]), sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__selection_05(self):
        input_array = self.test_data_numpy
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [2], sampling_type="selection")

    @pytest.mark.unit
    def test__init__selection_06(self):
        input_array = self.test_data_numpy
        with pytest.raises(ValueError):
            UniClass = UniformSampling(
                input_array, [2, 5, 5], sampling_type="selection"
            )

    @pytest.mark.unit
    def test__init__selection_07(self):
        input_array = self.test_data_numpy
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [-2, 5], sampling_type="selection")

    @pytest.mark.unit
    def test__init__selection_08(self):
        input_array = self.test_data_numpy
        with pytest.raises(TypeError):
            UniClass = UniformSampling(input_array, [2.1, 5], sampling_type="selection")

    @pytest.mark.unit
    def test__init__selection_09(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 50], sampling_type="selection")

    @pytest.mark.unit
    def test__init__selection_10(self):
        input_array = self.test_data_numpy
        UniClass = UniformSampling(
            input_array, [2, 5], sampling_type="selection", edges=True
        )
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_11(self):
        input_array = self.test_data_numpy
        UniClass = UniformSampling(
            input_array, [2, 5], sampling_type="selection", edges=False
        )
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, input_array[:, :-1])

    @pytest.mark.unit
    def test__init__selection_12(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [2, 5], sampling_type="selection", edges=1
            )

    @pytest.mark.unit
    def test__init__selection_13(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [2, 5], sampling_type="selection", edges="x"
            )

    @pytest.mark.unit
    def test__init__creation_01(self):
        input_array = self.test_data_list
        UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 2 * 7 * 5)

    @pytest.mark.unit
    def test__init__creation_02(self):
        input_array = self.test_data_list
        UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type="creation")
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 2 * 7 * 5)

    @pytest.mark.unit
    def test__init__creation_03(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [1, 7, 5], sampling_type="creation")

    @pytest.mark.unit
    def test__init__creation_04(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [-1, 7, 5], sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_05(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [1.1, 7, 5], sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_06(self):
        input_array = self.test_data_numpy
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="creation")

    @pytest.mark.unit
    def test__init__creation_07(self):
        input_array = self.test_data_pandas
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="creation")

    @pytest.mark.unit
    def test__init__creation_selection_01(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type=1)

    @pytest.mark.unit
    def test__init__creation_selection_02(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="jp")

    @pytest.mark.unit
    def test_sample_points_01(self):
        for num_samples in [[2, 5], [3, 2], [4, 2]]:
            input_array = self.test_data_numpy
            UniClass = UniformSampling(
                input_array, num_samples, sampling_type="selection"
            )
            unique_sample_points = UniClass.sample_points()
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            out_testing = [
                unique_sample_points[i, :] in input_array
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)

    @pytest.mark.unit
    def test_sample_points_02(self):
        input_array = self.test_data_list
        for num_samples in [[2, 5, 9], [3, 2, 10], [4, 2, 28]]:
            input_array = self.test_data_list
            UniClass = UniformSampling(
                input_array, num_samples, sampling_type="creation"
            )
            unique_sample_points = UniClass.sample_points()
            input_array = np.array(input_array)
            for i in range(input_array.shape[1]):
                var_range = input_array[:, i]
                assert (unique_sample_points[:, i] >= var_range[0]).all() and (
                    unique_sample_points[:, i] <= var_range[1]
                ).all()
            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0).shape,
                unique_sample_points.shape,
            )


class HammersleySamplingTestCases(unittest.TestCase):
    """
    __init__ = __init__ in LatinHypercubeSampling except Dimensionality problem:
        test__init__selection: Test behaviour with  dimensionality > 10

    test_sample_points_01: Test behaviour with selection, sample points are unique, all in the input array , number of samples = 5, 10, 1
    test_sample_points_02: Test behaviour with creation, sample points are unique, all in the input range (min,max) , number of samples = 5, 10, 1

    """

    def setUp(self):
        input_array_np = np.array(
            [
                [0, 10, 11],
                [1, 10, 14],
                [2, 10, 19],
                [3, 10, 26],
                [4, 10, 35],
                [5, 10, 46],
                [6, 10, 59],
                [7, 10, 74],
                [8, 10, 91],
                [9, 10, 110],
            ]
        )
        input_array_np_large = np.array(
            [
                [0, 10, 11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [1, 10, 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [2, 10, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [3, 10, 26, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [4, 10, 35, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [5, 10, 46, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [6, 10, 59, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [7, 10, 74, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [8, 10, 91, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                [9, 10, 110, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
            ]
        )
        input_array_pd = pd.DataFrame(
            {
                "x1": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "x2": [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
                "y": [11, 14, 19, 26, 35, 46, 59, 74, 91, 110],
            }
        )
        input_array_list = [[1, 10, 3], [2, 11, 4.5]]
        self.test_data_numpy = input_array_np
        self.test_data_numpy_large = input_array_np_large
        self.test_data_pandas = input_array_pd
        self.test_data_list = input_array_list

    @pytest.mark.unit
    def test__init__selection(self):
        input_array = self.test_data_numpy_large
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    def test_sample_points_01(self):
        for num_samples in [None, 10, 1]:
            input_array = self.test_data_numpy
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            unique_sample_points = HammersleyClass.sample_points()
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            out_testing = [
                unique_sample_points[i, :] in input_array
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)

    @pytest.mark.unit
    def test_sample_points_02(self):
        for num_samples in [None, 10, 1]:
            input_array = self.test_data_list
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=num_samples, sampling_type="creation"
            )
            unique_sample_points = HammersleyClass.sample_points()
            input_array = np.array(input_array)
            for i in range(input_array.shape[1]):
                var_range = input_array[:, i]
                assert (unique_sample_points[:, i] >= var_range[0]).all() and (
                    unique_sample_points[:, i] <= var_range[1]
                ).all()
            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0).shape,
                unique_sample_points.shape,
            )


class CVTSamplingTestCases(unittest.TestCase):
    """
    test__init__selection_01: input numpy array - Test behaviour generate object with selection, default number of sample = 5, default tolerance
    test__init__selection_02: input Pandas DataFrame - Test behaviour generate object with selection, default number of sample, default tolerance
    test__init__selection_03: input numpy array - Test behaviour generate object with selection, with a selected number of sample, default tolerance
    test__init__selection_04: input Pandas DataFrame - Test behaviour generate object with selection, with a selected number of sample, default tolerance
    test__init__selection_05: input numpy array - Test behaviour raise exception with a selected number of sample = 0, default tolerance
    test__init__selection_06: input numpy array - Test behaviour raise exception with a selected number of sample = -1, default tolerance
    test__init__selection_07: input numpy array - Test behaviour raise exception with a selected number of sample > input size, default tolerance
    test__init__selection_08: input numpy array - Test behaviour raise exception with a selected number of sample = 1.1 (non-integer), default tolerance
    test__init__selection_09: input list - Test behaviour raise ValueError with list input
    test__init__selection_10: input numpy array - Test behaviour raise exception with tolerance > 0.1
    test__init__selection_11: input numpy array - Test behaviour raise Warning with tolerance < 0.9
    test__init__selection_12: input numpy array - Test behaviour tolerance   0.09
    test__init__selection_13: input numpy array - Test behaviour tolerance  -0.09

    test__init__creation_01: input list - Test behaviour generate object with default sampling_type, default number of sample = 5, default tolerance
    test__init__creation_02: input list - Test behaviour generate object with sampling_type = creation, default number of sample = 5, default tolerance
    test__init__creation_03: input list - Test behaviour generate object with creation with a selected number of sample, default tolerance
    test__init__creation_04: input list - Test behaviour raise exception with a selected number of sample = 0
    test__init__creation_05: input list - Test behaviour raise exception with a selected number of sample = -1
    test__init__creation_06: input list - Test behaviour raise exception with a selected number of sample = 1.1 (non-integer)
    test__init__creation_07: input numpy - Test behaviour raise ValueError with numpy input
    test__init__creation_08: input numpy - Test behaviour raise ValueError with pandas input
    test__init__creation_09: input list - Test behaviour raise exception with tolerance > 0.1
    test__init__creation_10: input list - Test behaviour raise Warning with tolerance < 0.9
    test__init__creation_11: input list - Test behaviour tolerance   0.09
    test__init__creation_12: input list - Test behaviour tolerance  -0.09

    test__init__creation_selection_01 - Test behaviour raise Exception with sampling_type = non string
    test__init__creation_selection_02 - Test behaviour raise Exception with sampling_type = incorrect string


    test_random_sample_selection_01 - Test random_sample_selection with size = (5,2)
    test_random_sample_selection_02 - Test random_sample_selection with size = (0,2)
    test_random_sample_selection_03 - Test random_sample_selection with size = (2,0)
    test_random_sample_selection_04 - Test behaviour raise ValueError with  size = (5,-1)
    test_random_sample_selection_05 - Test behaviour raise TypeError with size =(5,1.1)

    Tests: : Unit tests for eucl_distance, a function that evaluates the distance between two points (u, v).
        Four demonstration tests are done:
    test_eucl_distance_01 - The first test checks that the correct result is obtained when both inputs are single value arrays.
    test_eucl_distance_02 - The second test checks that the correct result is returned when both inputs are 2D vectors.
    test_eucl_distance_03 - The third test checks that the correct result is returned when both inputs are arrays of the same size.
    test_eucl_distance_04 - The fourth test checks that the function is able to calculate the distance from a single point (n x 1 row vector) to a set of design points (supplied in an n x m array)

    Unit tests for create_centres, a function that generates new mass centroids for the design space based on McQueen's method.
        Four demonstration tests are done:
    test_create_centres_01 - The first test checks that the correct result for the new centres is returned when the counter is at its lowest value (1).
    test_create_centres_02 - The second test checks that the correct result for the new centres is returned when the counter is at an arbitrary value (10).
    test_create_centres_03 - The third test checks that the correct procedure is followed when one of the centres in initial_centres has no close design point to it.
    test_create_centres_04 - The fourth test checks that the the approach works as expected for problems with more than two dimensions.

    test_sample_points_01: Test behaviour with selection, sample points are unique, all in the input array , number of samples = 5, 10, 1
    test_sample_points_02: Test behaviour with creation, sample points are unique, all in the input range (min,max) , number of samples = 5, 10, 1

    """

    def setUp(self):
        input_array_np = np.array(
            [
                [0, 10, 11],
                [1, 10, 14],
                [2, 10, 19],
                [3, 10, 26],
                [4, 10, 35],
                [5, 10, 46],
                [6, 10, 59],
                [7, 10, 74],
                [8, 10, 91],
                [9, 10, 110],
            ]
        )
        input_array_pd = pd.DataFrame(
            {
                "x1": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                "x2": [10, 10, 10, 10, 10, 10, 10, 10, 10, 10],
                "y": [11, 14, 19, 26, 35, 46, 59, 74, 91, 110],
            }
        )
        input_array_list = [[1, 10, 3], [2, 11, 4.5]]
        self.test_data_numpy = input_array_np
        self.test_data_pandas = input_array_pd
        self.test_data_list = input_array_list

    @pytest.mark.unit
    def test__init__selection_01(self):
        input_array = self.test_data_numpy
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=None,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)
        np.testing.assert_array_equal(CVTClass.x_data, input_array[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)

    @pytest.mark.unit
    def test__init__selection_02(self):
        input_array = self.test_data_pandas
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=None,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)
        input_array = input_array.to_numpy()
        np.testing.assert_array_equal(CVTClass.x_data, input_array[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)

    @pytest.mark.unit
    def test__init__selection_03(self):
        input_array = self.test_data_numpy
        CVTClass = CVTSampling(
            input_array, number_of_samples=6, tolerance=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 6)
        np.testing.assert_array_equal(CVTClass.x_data, input_array[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)

    @pytest.mark.unit
    def test__init__selection_04(self):
        input_array = self.test_data_pandas
        CVTClass = CVTSampling(
            input_array, number_of_samples=6, tolerance=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 6)
        input_array = input_array.to_numpy()
        np.testing.assert_array_equal(CVTClass.x_data, input_array[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)

    @pytest.mark.unit
    def test__init__selection_05(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=0,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    def test__init__selection_06(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=-1,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    def test__init__selection_07(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=101,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    def test__init__selection_08(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=1.1,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    def test__init__selection_09(self):
        input_array = self.test_data_list
        with pytest.raises(ValueError):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    def test__init__selection_10(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=0.11,
                sampling_type="selection",
            )

    @pytest.mark.unit
    def test__init__selection_11(self):
        input_array = self.test_data_numpy
        with pytest.warns(Warning):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=1e-10,
                sampling_type="selection",
            )

    @pytest.mark.unit
    def test__init__selection_12(self):
        input_array = self.test_data_numpy
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=0.09,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.eps, 0.09)

    @pytest.mark.unit
    def test__init__selection_13(self):
        input_array = self.test_data_numpy
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=-0.09,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.eps, -0.09)

    @pytest.mark.unit
    def test__init__creation_01(self):
        input_array = self.test_data_list
        CVTClass = CVTSampling(
            input_array, number_of_samples=None, tolerance=None, sampling_type=None
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)

    @pytest.mark.unit
    def test__init__creation_02(self):
        input_array = self.test_data_list
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=None,
            sampling_type="creation",
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)

    @pytest.mark.unit
    def test__init__creation_03(self):
        input_array = self.test_data_list
        CVTClass = CVTSampling(
            input_array, number_of_samples=100, tolerance=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 100)

    @pytest.mark.unit
    def test__init__creation_04(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=0,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    def test__init__creation_05(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=-1,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    def test__init__creation_06(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=1.1,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    def test__init__creation_07(self):
        input_array = self.test_data_numpy
        with pytest.raises(ValueError):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    def test__init__creation_08(self):
        input_array = self.test_data_pandas
        with pytest.raises(ValueError):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    def test__init__creation_09(self):
        input_array = self.test_data_list
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=0.11,
                sampling_type="creation",
            )

    @pytest.mark.unit
    def test__init__creation_10(self):
        input_array = self.test_data_list
        with pytest.warns(Warning):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=1e-10,
                sampling_type="creation",
            )

    @pytest.mark.unit
    def test__init__creation_11(self):
        input_array = self.test_data_list
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=0.09,
            sampling_type="creation",
        )
        np.testing.assert_array_equal(CVTClass.eps, 0.09)

    @pytest.mark.unit
    def test__init__creation_11(self):
        input_array = self.test_data_list
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=-0.09,
            sampling_type="creation",
        )
        np.testing.assert_array_equal(CVTClass.eps, -0.09)

    @pytest.mark.unit
    def test__init__creation_selection_01(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, tolerance=None, sampling_type=1
            )

    @pytest.mark.unit
    def test__init__creation_selection_02(self):
        input_array = self.test_data_numpy
        with pytest.raises(Exception):
            CVTClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, tolerance=None, sampling_type="jp"
            )

    @pytest.mark.unit
    def test_random_sample_selection_01(self):
        size = (5, 2)
        out_random_points = CVTSampling.random_sample_selection(size[0], size[1])
        assert (out_random_points >= 0).all() and (out_random_points <= 1).all()
        assert out_random_points.shape == size

    @pytest.mark.unit
    def test_random_sample_selection_02(self):
        size = (0, 2)
        out_random_points = CVTSampling.random_sample_selection(size[0], size[1])
        assert (out_random_points >= 0).all() and (out_random_points <= 1).all()
        assert out_random_points.shape == size

    @pytest.mark.unit
    def test_random_sample_selection_03(self):
        size = (2, 0)
        out_random_points = CVTSampling.random_sample_selection(size[0], size[1])
        assert (out_random_points >= 0).all() and (out_random_points <= 1).all()
        assert out_random_points.shape == size

    @pytest.mark.unit
    def test_random_sample_selection_04(self):
        size = (5, -1)
        with pytest.raises(ValueError):
            out_random_points = CVTSampling.random_sample_selection(size[0], size[1])

    @pytest.mark.unit
    def test_random_sample_selection_05(self):
        size = (5, 1.1)
        with pytest.raises(TypeError):
            out_random_points = CVTSampling.random_sample_selection(size[0], size[1])

    @pytest.mark.unit
    def test_eucl_distance_01(self):
        u = np.array([[3]])
        v = np.array([[5]])
        expected_output = 2
        output = CVTSampling.eucl_distance(u, v)
        assert expected_output == output

    @pytest.mark.unit
    def test_eucl_distance_02(self):
        u = np.array([[1, 2]])
        v = np.array([[3, 4]])
        expected_output = 8**0.5
        output = CVTSampling.eucl_distance(u, v)
        assert expected_output == output

    @pytest.mark.unit
    def test_eucl_distance_03(self):
        u = np.array([[1, 2], [3, 4]])
        v = np.array([[5, 6], [7, 8]])
        expected_output = np.array([32**0.5, 32**0.5])
        output = CVTSampling.eucl_distance(u, v)
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    def test_eucl_distance_04(self):
        u = np.array([[1, 2]])
        v = np.array([[5, 6], [7, 8]])
        expected_output = np.array([32**0.5, 72**0.5])
        output = CVTSampling.eucl_distance(u, v)
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    def test_create_centres_01(self):
        initial_centres = np.array([[0, 0], [1, 1]])
        current_random_points = np.array([[0.6, 0.6], [0.3, 0.3]])
        current_centres = np.array([1, 0])
        counter = 1
        expected_output = np.array([[0.15, 0.15], [0.8, 0.8]])
        output = CVTSampling.create_centres(
            initial_centres, current_random_points, current_centres, counter
        )
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    def test_create_centres_02(self):
        initial_centres = np.array([[0, 0], [1, 1]])
        current_random_points = np.array([[0.6, 0.6], [0.3, 0.3]])
        current_centres = np.array([1, 0])
        counter = 10
        expected_output = np.array([[0.3 / 11, 0.3 / 11], [10.6 / 11, 10.6 / 11]])
        output = CVTSampling.create_centres(
            initial_centres, current_random_points, current_centres, counter
        )
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    def test_create_centres_03(self):
        initial_centres = np.array([[0, 0], [1, 1]])
        current_random_points = np.array([[0.6, 0.6], [0.8, 0.8]])
        current_centres = np.array([1, 1])
        counter = 5
        expected_output = np.array([[0.5 / 6, 0.5 / 6], [5.7 / 6, 5.7 / 6]])
        output = CVTSampling.create_centres(
            initial_centres, current_random_points, current_centres, counter
        )
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    def test_create_centres_04(self):
        initial_centres = np.array([[0, 0, 0], [1, 1, 1]])
        current_random_points = np.array(
            [
                [0.1, 0.1, 0.1],
                [0.3, 0.3, 0.3],
                [0.5, 0.5, 0.5],
                [0.7, 0.7, 0.7],
                [0.9, 0.9, 0.9],
            ]
        )
        current_centres = np.array([0, 0, 0, 1, 1])
        counter = 4
        expected_output = np.array(
            [[0.3 / 5, 0.3 / 5, 0.3 / 5], [4.8 / 5, 4.8 / 5, 4.8 / 5]]
        )
        output = CVTSampling.create_centres(
            initial_centres, current_random_points, current_centres, counter
        )
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.component
    def test_sample_points_01(self):
        for num_samples in [None, 10, 1]:
            input_array = self.test_data_numpy
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=num_samples,
                tolerance=None,
                sampling_type="selection",
            )
            unique_sample_points = CVTClass.sample_points()
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            out_testing = [
                unique_sample_points[i, :] in input_array
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)

    @pytest.mark.component
    def test_sample_points_02(self):
        for num_samples in [None, 10, 1]:
            input_array = self.test_data_list
            CVTClass = CVTSampling(
                input_array, number_of_samples=num_samples, sampling_type="creation"
            )
            unique_sample_points = CVTClass.sample_points()
            input_array = np.array(input_array)
            for i in range(input_array.shape[1]):
                var_range = input_array[:, i]
                assert (unique_sample_points[:, i] >= var_range[0]).all() and (
                    unique_sample_points[:, i] <= var_range[1]
                ).all()
            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0).shape,
                unique_sample_points.shape,
            )


if __name__ == "__main__":
    unittest.main()
