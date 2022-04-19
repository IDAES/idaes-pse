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
# third-party
import numpy as np
import pandas as pd
import pytest
import sys, os

sys.path.append(os.path.abspath(".."))  # current folder is ~/tests
# this
from idaes.core.surrogate.pysmo.sampling import (
    LatinHypercubeSampling,
    UniformSampling,
    HaltonSampling,
    HammersleySampling,
    CVTSampling,
    SamplingMethods,
    FeatureScaling,
)


class TestFeatureScaling:
    test_data_1d = [[x] for x in range(10)]
    test_data_2d = [[x, (x + 1) ** 2] for x in range(10)]
    test_data_3d = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    test_data_3d_constant = [[x, 10, (x + 1) ** 2 + 10] for x in range(10)]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_1d(self, array_type):
        input_array = array_type(self.test_data_1d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = (input_array - expected_output_2) / (
            expected_output_3 - expected_output_2
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(
            output_1, np.array(expected_output_1).reshape(10, 1)
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_2d(self, array_type):
        input_array = array_type(self.test_data_2d)
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
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_3d(self, array_type):
        input_array = array_type(self.test_data_3d)
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
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_3d_constant(self, array_type):
        input_array = array_type(self.test_data_3d_constant)
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
    def test_data_scaling_minmax_2d_list(self):
        input_array = self.test_data_2d
        with pytest.raises(TypeError):
            FeatureScaling.data_scaling_minmax(input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_01(self):
        input_array = np.array(self.test_data_1d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        output_1 = output_1.reshape(
            output_1.shape[0],
        )
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array.reshape(10, 1))

    @pytest.mark.unit
    def test_data_unscaling_minmax_02(self):
        input_array = np.array(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_03(self):
        input_array = np.array(self.test_data_3d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_04(self):
        input_array = np.array(self.test_data_3d_constant)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_05(self):
        input_array = np.array(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)

    @pytest.mark.unit
    def test_data_unscaling_minmax_06(self):
        input_array = np.array(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        min_array = np.array([[1, 2, 3]])
        max_array = np.array([[5, 6, 7]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)


class TestSamplingMethods:
    test_data_1d = [[x] for x in range(10)]
    test_data_2d = [[x, x + 10] for x in range(10)]
    test_data_3d = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]

    def _create_sampling(self, input_array, sample_points):
        sampling = SamplingMethods()
        sampling.data_headers = [i for i in range(0, sample_points.shape[1] + 1)]
        sampling.x_data = np.zeros((2, input_array.shape[1] - 1))
        return sampling

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_nearest_neighbour_01(self, array_type):
        input_array = array_type(self.test_data_3d)
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [-0.5, 1])
        np.testing.assert_array_equal(closest_point, input_array[0, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_nearest_neighbour_02(self, array_type):
        input_array = array_type(self.test_data_2d)
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [-0.5])
        np.testing.assert_array_equal(closest_point, input_array[0, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_nearest_neighbour_03(self, array_type):
        input_array = array_type(self.test_data_1d)
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [])
        np.testing.assert_array_equal(closest_point, input_array[0, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_nearest_neighbour_04(self, array_type):
        input_array = array_type(self.test_data_3d)
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        closest_point = SamplingMethods.nearest_neighbour(self, input_array, [0.5])
        np.testing.assert_array_equal(closest_point, input_array[0, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_nearest_neighbour_05(self, array_type):
        input_array = array_type(self.test_data_3d)
        self.x_data = np.zeros((2, input_array.shape[1] - 1))
        with pytest.raises(ValueError):
            closest_point = SamplingMethods.nearest_neighbour(
                self, input_array, [0.5, 0.9, 10]
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_points_selection_01(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[-0.5, 10], [10, 100]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        equivalent_points = sampling_methods.points_selection(
            input_array, generated_sample_points
        )
        np.testing.assert_array_equal(equivalent_points[0], input_array[0, :])
        np.testing.assert_array_equal(equivalent_points[1], input_array[-1, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_points_selection_02(self, array_type):
        input_array = array_type(self.test_data_2d)
        generated_sample_points = np.array([[-0.5], [10]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        equivalent_points = sampling_methods.points_selection(
            input_array, generated_sample_points
        )
        np.testing.assert_array_equal(equivalent_points[0], input_array[0, :])
        np.testing.assert_array_equal(equivalent_points[1], input_array[-1, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_points_selection_03(self, array_type):
        input_array = array_type(self.test_data_1d)
        generated_sample_points = np.array([[], []])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        equivalent_points = sampling_methods.points_selection(
            input_array, generated_sample_points
        )
        np.testing.assert_array_equal(equivalent_points[0], input_array[0, :])
        np.testing.assert_array_equal(equivalent_points[1], input_array[0, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_points_selection_04(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[0.5], [10]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            equivalent_points = sampling_methods.points_selection(
                input_array, generated_sample_points
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_points_selection_05(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[0.5, 0.7, 10], [10, 0.9, 20]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            equivalent_points = sampling_methods.points_selection(
                input_array, generated_sample_points
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_01(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[0, 0], [10, 19]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="selection"
        )
        np.testing.assert_array_equal(unique_sample_points[0], input_array[0, :])
        np.testing.assert_array_equal(unique_sample_points[1], input_array[-1, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_02(self, array_type):
        input_array = array_type(self.test_data_2d)
        generated_sample_points = np.array([[0], [7]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="selection"
        )
        np.testing.assert_array_equal(unique_sample_points[0], input_array[0, :])
        np.testing.assert_array_equal(unique_sample_points[1], input_array[-1, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_03(self, array_type):
        input_array = array_type(self.test_data_1d)
        generated_sample_points = np.array([[], []])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="selection"
        )
        np.testing.assert_array_equal(unique_sample_points[0], input_array[0, :])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_04(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[0.5], [7]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            unique_sample_points = sampling_methods.sample_point_selection(
                input_array, generated_sample_points, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_05(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[0.5, 1, 10], [7, 19, 20]])
        sampling_methods = self._create_sampling(input_array, generated_sample_points)
        with pytest.raises(ValueError):
            unique_sample_points = sampling_methods.sample_point_selection(
                input_array, generated_sample_points, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_06(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[0.5, 11, 3], [7, 19, 4]])
        sampling_methods = SamplingMethods()
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="creation"
        )
        min_, max_ = input_array[0, :], input_array[1, :]
        testing = min_ + generated_sample_points * (max_ - min_)
        np.testing.assert_array_equal(testing, unique_sample_points)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_07(self, array_type):
        input_array = array_type(self.test_data_2d)
        generated_sample_points = np.array([[0.5, 1], [7, 19]])
        sampling_methods = SamplingMethods()
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="creation"
        )
        min_, max_ = input_array[0, :], input_array[1, :]
        testing = min_ + generated_sample_points * (max_ - min_)
        np.testing.assert_array_equal(testing, unique_sample_points)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_08(self, array_type):
        input_array = array_type(self.test_data_1d)
        generated_sample_points = np.array([[0.5], [7]])
        sampling_methods = SamplingMethods()
        unique_sample_points = sampling_methods.sample_point_selection(
            input_array, generated_sample_points, sampling_type="creation"
        )
        min_, max_ = input_array[0, :], input_array[1, :]
        testing = min_ + generated_sample_points * (max_ - min_)
        np.testing.assert_array_equal(testing, unique_sample_points)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_09(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[], []])
        sampling_methods = SamplingMethods()
        with pytest.raises(IndexError):
            unique_sample_points = sampling_methods.sample_point_selection(
                input_array, generated_sample_points, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_point_selection_10(self, array_type):
        input_array = array_type(self.test_data_3d)
        generated_sample_points = np.array([[0.5, 1, 10, 11], [7, 19, 10, 12]])
        sampling_methods = SamplingMethods()
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


class TestLatinHypercubeSampling:
    input_array = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    input_array_list = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(2)]
    y = np.array(
        [
            [i, j, ((i + 1) ** 2) + ((j + 1) ** 2)]
            for i in np.linspace(0, 10, 21)
            for j in np.linspace(0, 10, 21)
        ]
    )
    full_data = {"x1": y[:, 0], "x2": y[:, 1], "y": y[:, 2]}

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_01(self, array_type):
        input_array = array_type(self.input_array)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)
        np.testing.assert_array_equal(LHSClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_02(self, array_type):
        input_array = array_type(self.input_array)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=6, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 6)
        np.testing.assert_array_equal(LHSClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_03(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=0, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_04(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=-1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_05(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=101, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_06(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=1.1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_07(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(ValueError):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_01(self, array_type):
        input_array = array_type(self.input_array_list)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type=None
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_02(self, array_type):
        input_array = array_type(self.input_array_list)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_03(self, array_type):
        input_array = array_type(self.input_array_list)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=100, sampling_type="creation"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_04(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=0, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_05(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=-1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_06(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=1.1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_07(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test__init__creation_08(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_09(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_10(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_11(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_12(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_13(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_selection_01(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type=1
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_selection_02(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="jp"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_variable_sample_creation(self, array_type):
        input_array = array_type(self.input_array)
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
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_lhs_points_generation(self, array_type):
        input_array = array_type(self.input_array)
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
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_random_shuffling(self, array_type):
        input_array = array_type(self.input_array)
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
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_01(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
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
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_02(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array_list)
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_sample_points_03(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.full_data)
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            unique_sample_points = LHSClass.sample_points()
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            unique_sample_points = np.array(unique_sample_points)
            out_testing = [
                unique_sample_points[i, :] in np.array(input_array)
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)


class TestUniformSampling:
    input_array = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    input_array_list = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(2)]
    y = np.array(
        [
            [i, j, ((i + 1) ** 2) + ((j + 1) ** 2)]
            for i in np.linspace(0, 10, 21)
            for j in np.linspace(0, 10, 21)
        ]
    )
    full_data = {"x1": y[:, 0], "x2": y[:, 1], "y": y[:, 2]}

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_01(self, array_type):
        input_array = array_type(self.input_array)
        UniClass = UniformSampling(input_array, [2, 5], sampling_type="selection")
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_02(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError):
            UniClass = UniformSampling(
                input_array, np.array([2, 5]), sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_03(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError):
            UniClass = UniformSampling(
                input_array, pd.DataFrame([2, 5]), sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_04(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [2], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_05(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            UniClass = UniformSampling(
                input_array, [2, 5, 5], sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_06(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [-2, 5], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_07(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError):
            UniClass = UniformSampling(input_array, [2.1, 5], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_08(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 50], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_09(self, array_type):
        input_array = array_type(self.input_array)
        UniClass = UniformSampling(
            input_array, [2, 5], sampling_type="selection", edges=True
        )
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_10(self, array_type):
        input_array = array_type(self.input_array)
        UniClass = UniformSampling(
            input_array, [2, 5], sampling_type="selection", edges=False
        )
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_11(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [2, 5], sampling_type="selection", edges=1
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_12(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [2, 5], sampling_type="selection", edges="x"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_13(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_01(self, array_type):
        input_array = array_type(self.input_array_list)
        UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 2 * 7 * 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_02(self, array_type):
        input_array = array_type(self.input_array_list)
        UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type="creation")
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 2 * 7 * 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_03(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [1, 7, 5], sampling_type="creation")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_04(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [-1, 7, 5], sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_05(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            UniClass = UniformSampling(
                input_array, [1.1, 7, 5], sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_06(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="creation")

    @pytest.mark.unit
    def test__init__creation_08(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def test__init__creation_09(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def test__init__creation_10(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def test__init__creation_11(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def test__init__creation_12(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_selection_01(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type=1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_selection_02(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="jp")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_01(self, array_type):
        for num_samples in [[2, 5], [3, 2], [4, 2]]:
            input_array = array_type(self.input_array)
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
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_02(self, array_type):
        for num_samples in [[2, 5, 9], [3, 2, 10], [4, 2, 28]]:
            input_array = array_type(self.input_array_list)
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_03(self, array_type):
        for num_samples in [[2, 5], [3, 2], [4, 2]]:
            input_array = array_type(self.input_array)
            UniClass = UniformSampling(
                input_array, num_samples, sampling_type="selection", edges=False
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
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_sample_points_04(self, array_type):
        for num_samples in [[2, 5], [3, 2], [4, 2]]:
            input_array = array_type(self.full_data)
            UniClass = UniformSampling(
                input_array, num_samples, sampling_type="selection"
            )
            unique_sample_points = UniClass.sample_points()
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            unique_sample_points = np.array(unique_sample_points)
            out_testing = [
                unique_sample_points[i, :] in np.array(input_array)
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)


class TestHaltonSampling:
    input_array = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    input_array_large = [
        [
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
        ]
        for x in range(10)
    ]
    input_array_high = [
        [
            x,
            x * 2,
            x**2,
            x**3,
            x**4,
            x**5,
            x**6,
            x**7,
            x**8,
            x**9,
            x**10,
            x**11,
        ]
        for x in range(10)
    ]
    input_array_list = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(2)]
    y = np.array(
        [
            [i, j, ((i + 1) ** 2) + ((j + 1) ** 2)]
            for i in np.linspace(0, 10, 21)
            for j in np.linspace(0, 10, 21)
        ]
    )
    full_data = {"x1": y[:, 0], "x2": y[:, 1], "y": y[:, 2]}

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_01(self, array_type):
        input_array = array_type(self.input_array)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 5)
        np.testing.assert_array_equal(HaltonClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_02(self, array_type):
        input_array = array_type(self.input_array)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=6, sampling_type="selection"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 6)
        np.testing.assert_array_equal(HaltonClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_03(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=0, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_04(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=-1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_05(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=101, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_06(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=1.1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_07(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(ValueError):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_08(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type=[1, 2, 3]
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_09(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="choose"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_10(self, array_type):
        input_array = array_type(self.input_array_high)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_01(self, array_type):
        input_array = array_type(self.input_array_list)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=None, sampling_type=None
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_02(self, array_type):
        input_array = array_type(self.input_array_list)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_03(self, array_type):
        input_array = array_type(self.input_array_list)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=100, sampling_type="creation"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_04(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=0, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_05(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=-1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_06(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=1.1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_07(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_09(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_10(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_11(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_12(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_13(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection(self, array_type):
        input_array = array_type(self.input_array_high)
        with pytest.raises(Exception):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_01(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            unique_sample_points = HaltonClass.sample_points()
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
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_02(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array_list)
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=num_samples, sampling_type="creation"
            )
            unique_sample_points = HaltonClass.sample_points()
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_sample_points_03(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.full_data)
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            unique_sample_points = HaltonClass.sample_points()
            unique_sample_points = np.array(unique_sample_points)
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            out_testing = [
                unique_sample_points[i, :] in np.array(input_array)
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)


class TestHammersleySampling:
    input_array = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    input_array_large = [
        [
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
            x,
            x + 10,
            (x + 1) ** 2 + x + 10,
        ]
        for x in range(10)
    ]
    input_array_high = [
        [
            x,
            x * 2,
            x**2,
            x**3,
            x**4,
            x**5,
            x**6,
            x**7,
            x**8,
            x**9,
            x**10,
            x**11,
        ]
        for x in range(10)
    ]
    input_array_list = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(2)]
    y = np.array(
        [
            [i, j, ((i + 1) ** 2) + ((j + 1) ** 2)]
            for i in np.linspace(0, 10, 21)
            for j in np.linspace(0, 10, 21)
        ]
    )
    full_data = {"x1": y[:, 0], "x2": y[:, 1], "y": y[:, 2]}

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_01(self, array_type):
        input_array = array_type(self.input_array)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 5)
        np.testing.assert_array_equal(
            HammersleyClass.x_data, np.array(input_array)[:, :-1]
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_02(self, array_type):
        input_array = array_type(self.input_array)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=6, sampling_type="selection"
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 6)
        np.testing.assert_array_equal(
            HammersleyClass.x_data, np.array(input_array)[:, :-1]
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_03(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=0, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_04(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=-1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_05(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=101, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_06(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=1.1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_07(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(ValueError):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_08(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type=[1, 2, 3]
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_09(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="choose"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_10(self, array_type):
        input_array = array_type(self.input_array_high)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_01(self, array_type):
        input_array = array_type(self.input_array_list)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=None, sampling_type=None
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_02(self, array_type):
        input_array = array_type(self.input_array_list)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_03(self, array_type):
        input_array = array_type(self.input_array_list)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=100, sampling_type="creation"
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_04(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=0, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_05(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=-1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_06(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=1.1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_06(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_09(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_10(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_11(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_12(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_13(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection(self, array_type):
        input_array = array_type(self.input_array_large)
        with pytest.raises(Exception):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_01(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
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
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_02(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array_list)
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_sample_points_03(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.full_data)
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=num_samples, sampling_type="selection"
            )
            unique_sample_points = HammersleyClass.sample_points()
            unique_sample_points = np.array(unique_sample_points)
            expected_testing = np.array(
                [True] * unique_sample_points.shape[0], dtype=bool
            )
            out_testing = [
                unique_sample_points[i, :] in np.array(input_array)
                for i in range(unique_sample_points.shape[0])
            ]

            np.testing.assert_array_equal(
                np.unique(unique_sample_points, axis=0), unique_sample_points
            )
            np.testing.assert_array_equal(expected_testing, out_testing)


class TestCVTSampling:
    input_array = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    input_array_list = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(2)]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_01(self, array_type):
        input_array = array_type(self.input_array)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=None,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)
        np.testing.assert_array_equal(CVTClass.x_data, np.array(input_array)[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_02(self, array_type):
        input_array = array_type(self.input_array)
        CVTClass = CVTSampling(
            input_array, number_of_samples=6, tolerance=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 6)
        np.testing.assert_array_equal(CVTClass.x_data, np.array(input_array)[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_03(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=0,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_04(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=-1,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_05(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=101,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_06(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=1.1,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_07(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(ValueError):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_08(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=0.11,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_09(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.warns(Warning):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=1e-10,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_10(self, array_type):
        input_array = array_type(self.input_array)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=0.09,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.eps, 0.09)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_11(self, array_type):
        input_array = array_type(self.input_array)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=-0.09,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.eps, -0.09)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_01(self, array_type):
        input_array = array_type(self.input_array_list)
        CVTClass = CVTSampling(
            input_array, number_of_samples=None, tolerance=None, sampling_type=None
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_02(self, array_type):
        input_array = array_type(self.input_array_list)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=None,
            sampling_type="creation",
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_03(self, array_type):
        input_array = array_type(self.input_array_list)
        CVTClass = CVTSampling(
            input_array, number_of_samples=100, tolerance=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_04(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=0,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_05(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=-1,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_06(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=1.1,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_07(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_08(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=0.11,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_09(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.warns(Warning):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=1e-10,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_10(self, array_type):
        input_array = array_type(self.input_array_list)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=0.09,
            sampling_type="creation",
        )
        np.testing.assert_array_equal(CVTClass.eps, 0.09)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_11(self, array_type):
        input_array = array_type(self.input_array_list)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=-0.09,
            sampling_type="creation",
        )
        np.testing.assert_array_equal(CVTClass.eps, -0.09)

    @pytest.mark.unit
    def test__init__creation_13(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(Exception):
            LHSClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_14(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_15(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_16(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_17(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_selection_01(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, tolerance=None, sampling_type=1
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_selection_02(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(Exception):
            CVTClass = CVTSampling(
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_01(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_02(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array_list)
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
    pytest.main()
