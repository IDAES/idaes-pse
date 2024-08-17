#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
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
    CustomSampling,
    SamplingMethods,
    FeatureScaling,
)
import idaes.logger as idaeslog


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
    def test__init__selection_right_behaviour_with_none_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)
        np.testing.assert_array_equal(LHSClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=6, sampling_type="selection"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 6)
        np.testing.assert_array_equal(LHSClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_random_seed(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        rand_seed = 100
        LHSClass = LatinHypercubeSampling(
            input_array,
            number_of_samples=6,
            sampling_type="selection",
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 6)
        np.testing.assert_array_equal(LHSClass.x_data, np.array(input_array)[:, :-1])
        assert LHSClass.seed_value == rand_seed

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_float_random_seed(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        rand_seed = 15.1
        LHSClass = LatinHypercubeSampling(
            input_array,
            number_of_samples=6,
            sampling_type="selection",
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 6)
        np.testing.assert_array_equal(LHSClass.x_data, np.array(input_array)[:, :-1])
        assert LHSClass.seed_value == int(rand_seed)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_zero_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=0, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=-1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_excess_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="LHS sample size cannot be greater than number of samples in the input data set",
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=101, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=1.1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match='Pandas dataframe or numpy array required for sampling_type "selection."',
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_random_seed(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError, match="Random seed must be an integer."):
            LHSClass = LatinHypercubeSampling(
                input_array,
                number_of_samples=5,
                sampling_type="selection",
                rand_seed="1.2",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_none_samplingtype(self, array_type):
        input_array = array_type(self.input_array_list)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type=None
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_none_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        LHSClass = LatinHypercubeSampling(
            input_array, number_of_samples=100, sampling_type="creation"
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_specified_seed(self, array_type):
        input_array = array_type(self.input_array_list)
        rand_seed = 50
        LHSClass = LatinHypercubeSampling(
            input_array,
            number_of_samples=100,
            sampling_type="creation",
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(LHSClass.data, input_array)
        np.testing.assert_array_equal(LHSClass.number_of_samples, 100)
        assert LHSClass.seed_value == rand_seed

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_zero_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=0, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=-1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=1.1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='List entry of two elements expected for sampling_type "creation."',
        ):
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
    def test__init__creation_missing_bounds(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_lb(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_ub(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_unequal_length_list_bounds(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_equal_input_output_bounds_all(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(ValueError, match="Invalid entry: both lists are equal."):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_nonstring(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="Invalid sampling type entry. Must be of type <str>."
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type=1
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_undefined_string(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.',
        ):
            LHSClass = LatinHypercubeSampling(
                input_array, number_of_samples=None, sampling_type="jp"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__test_single_equal_ub_lb(self, array_type):
        input_array = array_type([[0, 0, 0], [0, 1, 1]])
        with pytest.raises(
            ValueError,
            match="Invalid entry: at least one variable contains the same value for the lower and upper bounds.",
        ):
            LHSClass = LatinHypercubeSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
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

    @pytest.mark.integration
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_equality_fixed_seed(self, array_type):
        rand_seed = 1000
        for num_samples in [None, 1, 10, 100]:  # Test for different number of samples
            input_array = array_type(self.input_array_list)

            LHSClass_A = LatinHypercubeSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                rand_seed=rand_seed,
            )
            unique_sample_points_A = LHSClass_A.sample_points()

            LHSClass_B = LatinHypercubeSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                rand_seed=rand_seed,
            )
            unique_sample_points_B = LHSClass_B.sample_points()

            np.testing.assert_array_equal(
                unique_sample_points_A, unique_sample_points_B
            )


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
    def test__init__selection_right_behaviour(self, array_type):
        input_array = array_type(self.input_array)
        UniClass = UniformSampling(input_array, [2, 5], sampling_type="selection")
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_wrong_type_for_list_of_samples_per_variable_01(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="list_of_samples_per_variable: list required."
        ):
            UniClass = UniformSampling(
                input_array, np.array([2, 5]), sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_wrong_type_for_list_of_samples_per_variable_02(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="list_of_samples_per_variable: list required."
        ):
            UniClass = UniformSampling(
                input_array, pd.DataFrame([2, 5]), sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_wrong_length_for_list_of_samples_per_variable_01(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="Length of list_of_samples_per_variable must equal the number of variables.",
        ):
            UniClass = UniformSampling(input_array, [2], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_wrong_length_for_list_of_samples_per_variable_02(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="Length of list_of_samples_per_variable must equal the number of variables.",
        ):
            UniClass = UniformSampling(
                input_array, [2, 5, 5], sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_negative_entry_in_list_of_samples_per_variable(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="All variables must have at least two points per dimension",
        ):
            UniClass = UniformSampling(input_array, [-2, 5], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_fractional_entry_in_list_of_samples_per_variable(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match="All values in list must be integers"):
            UniClass = UniformSampling(input_array, [2.1, 5], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_excess_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="Sample size cannot be greater than number of samples in the input data set",
        ):
            UniClass = UniformSampling(input_array, [2, 50], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_assert_correct_behaviour_edge_true(self, array_type):
        input_array = array_type(self.input_array)
        UniClass = UniformSampling(
            input_array, [2, 5], sampling_type="selection", edges=True
        )
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, np.array(input_array)[:, :-1])
        assert UniClass.edge == True

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_assert_correct_behaviour_edge_false(self, array_type):
        input_array = array_type(self.input_array)
        UniClass = UniformSampling(
            input_array, [2, 5], sampling_type="selection", edges=False
        )
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 10)
        np.testing.assert_array_equal(UniClass.x_data, np.array(input_array)[:, :-1])
        assert UniClass.edge == False

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_nonboolean_edge_entry_01(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match='Invalid "edges" entry. Must be boolean'):
            UniClass = UniformSampling(
                input_array, [2, 5], sampling_type="selection", edges=1
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_nonboolean_edge_entry_02(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match='Invalid "edges" entry. Must be boolean'):
            UniClass = UniformSampling(
                input_array, [2, 5], sampling_type="selection", edges="x"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match='Pandas dataframe or numpy array required for sampling_type "selection."',
        ):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="selection")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_none_samplingtype(self, array_type):
        input_array = array_type(self.input_array_list)
        UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 2 * 7 * 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_specified_sampling(self, array_type):
        input_array = array_type(self.input_array_list)
        UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type="creation")
        np.testing.assert_array_equal(UniClass.data, input_array)
        np.testing.assert_array_equal(UniClass.number_of_samples, 2 * 7 * 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_wrong_entry_in_list_of_samples_per_variable(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="All variables must have at least two points per dimension",
        ):
            UniClass = UniformSampling(input_array, [1, 7, 5], sampling_type="creation")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation__negative_entry_in_list_of_samples_per_variable(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="All variables must have at least two points per dimension",
        ):
            UniClass = UniformSampling(
                input_array, [-1, 7, 5], sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_invalid_entry_in_list_of_samples_per_variable(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="All variables must have at least two points per dimension",
        ):
            UniClass = UniformSampling(
                input_array, [1.1, 7, 5], sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='List entry of two elements expected for sampling_type "creation."',
        ):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="creation")

    @pytest.mark.unit
    def test__init__creation_missing_bounds(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_lb(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_ub(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def test__init__creation_unequal_length_list_bounds(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    def est__init__creation_equal_input_output_bounds_all(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(ValueError, match="Invalid entry: both lists are equal."):
            UniClass = UniformSampling(input_array, [2, 7, 5], sampling_type=None)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_nonstring(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="Invalid sampling type entry. Must be of type <str>."
        ):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type=1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_undefined_string(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.',
        ):
            UniClass = UniformSampling(input_array, [2, 5], sampling_type="jp")

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_equal_input_output_bounds_one(self, array_type):
        input_array = array_type([[0, 0, 0], [0, 1, 1]])
        with pytest.raises(
            ValueError,
            match="Invalid entry: at least one variable contains the same value for the lower and upper bounds.",
        ):
            UniClass = UniformSampling(
                input_array,
                [2, 7, 5],
                sampling_type="creation",
            )

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
    def test__init__selection_right_behaviour_with_none_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=None, sampling_type="selection"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 5)
        np.testing.assert_array_equal(HaltonClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=6, sampling_type="selection"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 6)
        np.testing.assert_array_equal(HaltonClass.x_data, np.array(input_array)[:, :-1])

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_zero_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=0, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=-1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_excess_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="Sample size cannot be greater than number of samples in the input data set",
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=101, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=1.1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match='Pandas dataframe or numpy array required for sampling_type "selection."',
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_nonstring(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="Invalid sampling type entry. Must be of type <str>."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type=[1, 2, 3]
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_undefined_string(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.',
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="choose"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_method_dimensionality_exceeded(self, array_type):
        input_array = array_type(self.input_array_high)
        with pytest.raises(
            Exception,
            match="Dimensionality problem: This method is not available for problems with dimensionality > 10: the performance of the method degrades substantially at higher dimensions",
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_none_samplingtype(self, array_type):
        input_array = array_type(self.input_array_list)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=None, sampling_type=None
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_hahaviour_with_none_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_hahaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        HaltonClass = HaltonSampling(
            input_array, number_of_samples=100, sampling_type="creation"
        )
        np.testing.assert_array_equal(HaltonClass.data, input_array)
        np.testing.assert_array_equal(HaltonClass.number_of_samples, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_zero_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=0, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=-1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=1.1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError,
            match='List entry of two elements expected for sampling_type "creation."',
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_missing_bounds(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_lb(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_ub(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_unequal_length_list_bounds(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_equal_input_output_bounds_all(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(ValueError, match="Invalid entry: both lists are equal."):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_dimensionality_exceeded(self, array_type):
        input_array = array_type(self.input_array_high)
        with pytest.raises(
            Exception,
            match="Dimensionality problem: This method is not available for problems with dimensionality > 10: the performance of the method degrades substantially at higher dimensions",
        ):
            HaltonClass = HaltonSampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__test_single_equal_ub_lb(self, array_type):
        input_array = array_type([[0, 0, 0], [0, 1, 1]])
        with pytest.raises(
            ValueError,
            match="Invalid entry: at least one variable contains the same value for the lower and upper bounds.",
        ):
            HaltonClass = HaltonSampling(
                input_array,
                number_of_samples=5,
                sampling_type="creation",
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
    def test__init__selection_right_behaviour_with_none_no_samples(self, array_type):
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
    def test__init__selection_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
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
    def test__init__selection_zero_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=0, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=-1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_excess_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="Sample size cannot be greater than number of samples in the input data set",
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=101, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=1.1, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match='Pandas dataframe or numpy array required for sampling_type "selection."',
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_nonstring(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="Invalid sampling type entry. Must be of type <str>."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type=[1, 2, 3]
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_undefined_string(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.',
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="choose"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_method_dimensionality_exceeded(self, array_type):
        input_array = array_type(self.input_array_high)
        with pytest.raises(
            Exception,
            match="Dimensionality problem: This method is not available for problems with dimensionality > 10: the performance of the method degrades substantially at higher dimensions",
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_hahaviour_with_none_samplingtype(self, array_type):
        input_array = array_type(self.input_array_list)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=None, sampling_type=None
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_none_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 5)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        HammersleyClass = HammersleySampling(
            input_array, number_of_samples=100, sampling_type="creation"
        )
        np.testing.assert_array_equal(HammersleyClass.data, input_array)
        np.testing.assert_array_equal(HammersleyClass.number_of_samples, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_zero_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=0, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=-1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=1.1, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='List entry of two elements expected for sampling_type "creation."',
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_missing_bounds(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_lb(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_ub(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_unequal_length_list_bounds(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_equal_input_output_bounds_all(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(ValueError, match="Invalid entry: both lists are equal."):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_method_dimensionality_exceeded(self, array_type):
        input_array = array_type(self.input_array_large)
        with pytest.raises(
            Exception,
            match="Dimensionality problem: This method is not available for problems with dimensionality > 10: the performance of the method degrades substantially at higher dimensions",
        ):
            HammersleyClass = HammersleySampling(
                input_array, number_of_samples=None, sampling_type="selection"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_equal_input_output_bounds_one(self, array_type):
        input_array = array_type([[0, 0, 0], [0, 1, 1]])
        with pytest.raises(
            ValueError,
            match="Invalid entry: at least one variable contains the same value for the lower and upper bounds.",
        ):
            HammersleyClass = HammersleySampling(
                input_array,
                number_of_samples=5,
                sampling_type=None,
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
    def test__init__selection_right_behaviour_with_none_no_samples(self, array_type):
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
    def test__init__selection_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
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
    def test__init__selection_right_behaviour_with_specified_random_seed(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        rand_seed = 100
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=6,
            tolerance=None,
            sampling_type="selection",
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 6)
        np.testing.assert_array_equal(CVTClass.x_data, np.array(input_array)[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)
        assert CVTClass.seed_value == rand_seed

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_float_random_seed(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        rand_seed = 2.2
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=6,
            tolerance=None,
            sampling_type="selection",
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 6)
        np.testing.assert_array_equal(CVTClass.x_data, np.array(input_array)[:, :-1])
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)
        assert CVTClass.seed_value == int(rand_seed)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_zero_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=0,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=-1,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_excess_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="CVT sample size cannot be greater than number of samples in the input data set",
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=101,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=1.1,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match='Pandas dataframe or numpy array required for sampling_type "selection."',
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=None,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_tolerance_too_loose(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="Tolerance must be less than 0.1 to achieve good results"
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=0.11,
                sampling_type="selection",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_tolerance_too_tight(self, array_type, caplog):
        caplog.set_level(idaeslog.WARNING)
        warning_msg = (
            "Tolerance too tight. CVT algorithm may take long time to converge."
        )
        input_array = array_type(self.input_array)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=1e-10,
            sampling_type="selection",
        )
        assert warning_msg in caplog.text
        for record in caplog.records:
            assert record.levelno == idaeslog.WARNING

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_random_seed(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError, match="Random seed must be an integer."):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=5,
                sampling_type="selection",
                rand_seed="1.2",
                tolerance=None,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_valid_tolerance(self, array_type):
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
    def test__init__selection_none_tolerance(self, array_type):
        input_array = array_type(self.input_array)
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=None,
            sampling_type="selection",
        )
        np.testing.assert_array_equal(CVTClass.eps, 1e-7)

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
    def test__init__creation_right_behaviour_with_none_samplingtype(self, array_type):
        input_array = array_type(self.input_array_list)
        CVTClass = CVTSampling(
            input_array, number_of_samples=None, tolerance=None, sampling_type=None
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 5)
        assert CVTClass.sampling_type == "creation"

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_none_no_samples(self, array_type):
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
    def test__init__creation_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        CVTClass = CVTSampling(
            input_array, number_of_samples=100, tolerance=None, sampling_type="creation"
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 100)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_specified_seed(self, array_type):
        input_array = array_type(self.input_array_list)
        rand_seed = 50
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=100,
            tolerance=None,
            sampling_type="creation",
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(CVTClass.data, input_array)
        np.testing.assert_array_equal(CVTClass.number_of_centres, 100)
        assert CVTClass.seed_value == rand_seed

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_zero_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=0,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=-1,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=1.1,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='List entry of two elements expected for sampling_type "creation."',
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=None,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_test_tolerance_too_loose(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="Tolerance must be less than 0.1 to achieve good results"
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=None,
                tolerance=0.11,
                sampling_type="creation",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_tolerance_too_tight(self, array_type, caplog):
        caplog.set_level(idaeslog.WARNING)
        input_array = array_type(self.input_array_list)
        warning_msg = (
            "Tolerance too tight. CVT algorithm may take long time to converge."
        )
        CVTClass = CVTSampling(
            input_array,
            number_of_samples=None,
            tolerance=1e-10,
            sampling_type="creation",
        )
        assert warning_msg in caplog.text
        for record in caplog.records:
            assert record.levelno == idaeslog.WARNING

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_valid_tolerance(self, array_type):
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
    def test__init__creation_missing_bounds(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            LHSClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_lb(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_ub(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_unequal_length_list_bounds(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    def test__init__creation_equal_input_output_bounds_all(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(ValueError, match="Invalid entry: both lists are equal."):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, sampling_type="creation"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_nonstring(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="Invalid sampling type entry. Must be of type <str>."
        ):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, tolerance=None, sampling_type=1
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_undefined_string(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.',
        ):
            CVTClass = CVTSampling(
                input_array, number_of_samples=None, tolerance=None, sampling_type="jp"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_equal_input_output_bounds_one(self, array_type):
        input_array = array_type([[0, 0, 0], [0, 1, 1]])
        with pytest.raises(
            ValueError,
            match="Invalid entry: at least one variable contains the same value for the lower and upper bounds.",
        ):
            CVTClass = CVTSampling(
                input_array,
                number_of_samples=5,
                tolerance=None,
                sampling_type=None,
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
    def test_random_sample_selection_test_negative(self):
        size = (5, -1)
        with pytest.raises(ValueError, match="negative dimensions are not allowed"):
            out_random_points = CVTSampling.random_sample_selection(size[0], size[1])

    @pytest.mark.unit
    def test_random_sample_selection_test_float(self):
        size = (5, 1.1)
        with pytest.raises(
            TypeError, match="'float' object cannot be interpreted as an integer"
        ):
            out_random_points = CVTSampling.random_sample_selection(size[0], size[1])

    @pytest.mark.unit
    def test_eucl_distance_single_values(self):
        u = np.array([[3]])
        v = np.array([[5]])
        expected_output = 2
        output = CVTSampling.eucl_distance(u, v)
        assert expected_output == output

    @pytest.mark.unit
    def test_eucl_distance_1d_arrays(self):
        u = np.array([[1, 2]])
        v = np.array([[3, 4]])
        expected_output = 8**0.5
        output = CVTSampling.eucl_distance(u, v)
        assert expected_output == output

    @pytest.mark.unit
    def test_eucl_distance_2d_arrays(self):
        u = np.array([[1, 2], [3, 4]])
        v = np.array([[5, 6], [7, 8]])
        expected_output = np.array([32**0.5, 32**0.5])
        output = CVTSampling.eucl_distance(u, v)
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    def test_eucl_distance_1d_2d_arrays(self):
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

    @pytest.mark.integration
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_equality_fixed_seed(self, array_type):
        rand_seed = 1000
        for num_samples in [None, 1, 10, 100]:  # Test for different number of samples
            input_array = array_type(self.input_array_list)

            CVTClass_A = CVTSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                rand_seed=rand_seed,
            )
            unique_sample_points_A = CVTClass_A.sample_points()

            CVTClass_B = CVTSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                rand_seed=rand_seed,
            )
            unique_sample_points_B = CVTClass_B.sample_points()

            np.testing.assert_array_equal(
                unique_sample_points_A, unique_sample_points_B
            )


class TestCustomSampling:
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
    def test__init__selection_right_behaviour_with_none_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        CSClass = CustomSampling(
            input_array,
            number_of_samples=None,
            sampling_type="selection",
            list_of_distributions=["uniform", "normal"],
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 5)
        np.testing.assert_array_equal(CSClass.x_data, np.array(input_array)[:, :-1])
        assert CSClass.dist_vector == ["uniform", "normal"]
        assert CSClass.normal_bounds_enforced == False

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        CSClass = CustomSampling(
            input_array,
            number_of_samples=6,
            sampling_type="selection",
            list_of_distributions=["uniform", "normal"],
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 6)
        np.testing.assert_array_equal(CSClass.x_data, np.array(input_array)[:, :-1])
        assert CSClass.dist_vector == ["uniform", "normal"]
        assert CSClass.normal_bounds_enforced == False

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_bounds_option_true(self, array_type):
        input_array = array_type(self.input_array)
        CSClass = CustomSampling(
            input_array,
            number_of_samples=6,
            sampling_type="selection",
            list_of_distributions=["uniform", "normal"],
            strictly_enforce_gaussian_bounds=True,
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 6)
        np.testing.assert_array_equal(CSClass.x_data, np.array(input_array)[:, :-1])
        assert CSClass.dist_vector == ["uniform", "normal"]
        assert CSClass.normal_bounds_enforced == True

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_bounds_option_false(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        CSClass = CustomSampling(
            input_array,
            number_of_samples=6,
            sampling_type="selection",
            list_of_distributions=["uniform", "normal"],
            strictly_enforce_gaussian_bounds=False,
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 6)
        np.testing.assert_array_equal(CSClass.x_data, np.array(input_array)[:, :-1])
        assert CSClass.dist_vector == ["uniform", "normal"]
        assert CSClass.normal_bounds_enforced == False

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_random_seed(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        rand_seed = 1000
        CSClass = CustomSampling(
            input_array,
            number_of_samples=6,
            sampling_type="selection",
            list_of_distributions=["uniform", "normal"],
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 6)
        np.testing.assert_array_equal(CSClass.x_data, np.array(input_array)[:, :-1])
        assert CSClass.dist_vector == ["uniform", "normal"]
        assert CSClass.seed_value == rand_seed

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_right_behaviour_with_specified_float_random_seed(
        self, array_type
    ):
        input_array = array_type(self.input_array)
        rand_seed = 1.2
        CSClass = CustomSampling(
            input_array,
            number_of_samples=6,
            sampling_type="selection",
            list_of_distributions=["uniform", "normal"],
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 6)
        np.testing.assert_array_equal(CSClass.x_data, np.array(input_array)[:, :-1])
        assert CSClass.dist_vector == ["uniform", "normal"]
        assert CSClass.seed_value == int(rand_seed)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_zero_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=0,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=-1,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_excess_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match="Sample size cannot be greater than number of samples in the input data set",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=101,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=1.1,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match='Pandas dataframe or numpy array required for sampling_type "selection."',
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_list_distributions_length_exceeds_inputs(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="Length of list_of_distributions must equal the number of variables.",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_list_distributions_length_less_than_inputs(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="Length of list_of_distributions must equal the number of variables.",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="selection",
                list_of_distributions=["uniform"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_empty_distributions(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(ValueError, match="list_of_distributions cannot be empty."):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="selection",
                list_of_distributions=None,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_distribution_not_list(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            TypeError, match="Error with list_of_distributions: list required."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="selection",
                list_of_distributions=("uniform", "normal"),
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_distribution_entry_not_string(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(TypeError, match="All values in list must be strings"):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="selection",
                list_of_distributions=["uniform", 1.0],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_distribution_not_available(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="list_of_distributions only supports 'random', 'normal' and 'uniform' sampling options.",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="selection",
                list_of_distributions=["uniform", "binomial"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__selection_non_integer_random_seed(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(ValueError, match="Random seed must be an integer."):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=5,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
                rand_seed="1.2",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_hahaviour_with_none_samplingtype(self, array_type):
        input_array = array_type(self.input_array_list)
        CSClass = CustomSampling(
            input_array,
            number_of_samples=None,
            sampling_type=None,
            list_of_distributions=["uniform", "normal", "random"],
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 5)
        assert CSClass.dist_vector == ["uniform", "normal", "random"]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_none_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        CSClass = CustomSampling(
            input_array,
            number_of_samples=None,
            sampling_type="creation",
            list_of_distributions=["uniform", "normal", "random"],
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 5)
        assert CSClass.dist_vector == ["uniform", "normal", "random"]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_specified_no_samples(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        CSClass = CustomSampling(
            input_array,
            number_of_samples=100,
            sampling_type="creation",
            list_of_distributions=["uniform", "normal", "random"],
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 100)
        assert CSClass.dist_vector == ["uniform", "normal", "random"]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_right_behaviour_with_specified_seed(self, array_type):
        input_array = array_type(self.input_array_list)
        rand_seed = 50
        CSClass = CustomSampling(
            input_array,
            number_of_samples=100,
            sampling_type="creation",
            list_of_distributions=["uniform", "normal", "random"],
            rand_seed=rand_seed,
        )
        np.testing.assert_array_equal(CSClass.data, input_array)
        np.testing.assert_array_equal(CSClass.number_of_samples, 100)
        assert CSClass.dist_vector == ["uniform", "normal", "random"]
        assert CSClass.seed_value == rand_seed

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_zero_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=0,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_negative_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError, match="number_of_samples must a positive, non-zero integer."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=-1,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_non_integer_no_samples(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(TypeError, match="number_of_samples must be an integer."):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=1.1,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__creation_wrong_input_data_type(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='List entry of two elements expected for sampling_type "creation."',
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    def test__init__creation_missing_bounds(self):
        input_array = [[2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_lb(self):
        input_array = [np.array([1, 10, 3]), [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    def test__init__creation_wrong_data_input_format_ub(self):
        input_array = [[1, 10, 3], np.array([2, 11, 4.5])]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    def test__init__creation_unequal_length_list_bounds(self):
        input_array = [[1, 10], [2, 11, 4.5]]
        with pytest.raises(
            ValueError, match="data_input must contain two lists of equal lengths."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    def test__init__creation_equal_input_output_bounds_all(self):
        input_array = [[2, 11, 4.5], [2, 11, 4.5]]
        with pytest.raises(ValueError, match="Invalid entry: both lists are equal."):
            csClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_list_distributions_length_less_than_inputs(
        self, array_type
    ):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="Length of list_of_distributions must equal the number of variables.",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_list_distributions_length_exceeds_inputs(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="Length of list_of_distributions must equal the number of variables.",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "uniform", "normal", "random"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_empty_distributions(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(ValueError, match="list_of_distributions cannot be empty."):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=None,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_distribution_not_list(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            TypeError, match="Error with list_of_distributions: list required."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=("uniform", "normal", "random"),
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_distribution_entry_not_string(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(TypeError, match="All values in list must be strings"):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", 1.0, "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__selection_distribution_not_available(self, array_type):
        input_array = array_type(self.input_array_list)
        with pytest.raises(
            ValueError,
            match="list_of_distributions only supports 'random', 'normal' and 'uniform' sampling options.",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="creation",
                list_of_distributions=["uniform", "gaussian", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_equal_input_output_bounds_one(self, array_type):
        input_array = array_type([[0, 0, 0], [0, 1, 1]])
        with pytest.raises(Exception):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type=None,
                list_of_distributions=["uniform", "normal", "random"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test__init__creation_nonboolean_bounds_option(self, array_type):
        input_array = array_type([[0, 0, 0], [1, 1, 1]])
        with pytest.raises(
            TypeError,
            match="Invalid 'strictly_enforce_gaussian_bounds' entry. Must be boolean.",
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type=None,
                list_of_distributions=["uniform", "normal", "random"],
                strictly_enforce_gaussian_bounds="False",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_nonstring(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            TypeError, match="Invalid sampling type entry. Must be of type <str>."
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type=1,
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__samplingtype_undefined_string(self, array_type):
        input_array = array_type(self.input_array)
        with pytest.raises(
            ValueError,
            match='Invalid sampling type requirement entered. Enter "creation" for sampling from a range or "selection" for selecting samples from a dataset.',
        ):
            CSClass = CustomSampling(
                input_array,
                number_of_samples=None,
                sampling_type="jp",
                list_of_distributions=["uniform", "normal"],
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_generate_from_dist_uniform(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )
            dist_type = "uniform"
            dist_res, scaled_samples = CSClass.generate_from_dist(dist_type)
            assert type(scaled_samples) == np.ndarray
            assert scaled_samples.shape == (CSClass.number_of_samples,)
            assert dist_res.__name__ == dist_type

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_generate_from_dist_normal(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )
            dist_type = "normal"
            dist_res, scaled_samples = CSClass.generate_from_dist(dist_type)
            assert type(scaled_samples) == np.ndarray
            assert scaled_samples.shape == (CSClass.number_of_samples,)
            assert dist_res.__name__ == dist_type

    @pytest.mark.unit
    def test_generate_from_dist_normal_unenforced_gaussian_bounds(self):
        CSClass = CustomSampling(
            [[0], [1]],
            number_of_samples=10000,
            sampling_type="creation",
            list_of_distributions=["normal"],
        )
        dist_type = "normal"
        dist_res, scaled_samples = CSClass.generate_from_dist(dist_type)
        assert dist_res.__name__ == dist_type
        assert scaled_samples.min() < 0
        assert scaled_samples.max() > 1

    @pytest.mark.unit
    def test_generate_from_dist_normal_enforced_gaussian_bounds(self):
        CSClass = CustomSampling(
            [[0], [1]],
            number_of_samples=10000,
            sampling_type="creation",
            list_of_distributions=["normal"],
            strictly_enforce_gaussian_bounds=True,
        )
        dist_type = "normal"
        dist_res, scaled_samples = CSClass.generate_from_dist(dist_type)
        assert dist_res.__name__ == dist_type
        assert scaled_samples.min() >= 0
        assert scaled_samples.max() <= 1

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_generate_from_dist_random(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )
            dist_type = "random"
            dist_res, scaled_samples = CSClass.generate_from_dist(dist_type)
            assert type(scaled_samples) == np.ndarray
            assert scaled_samples.shape == (CSClass.number_of_samples,)
            assert dist_res.__name__ == dist_type

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_generate_from_dist_all_types(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["uniform", "normal"],
            )
            for dist_type in ["uniform", "random", "uniform"]:
                dist_res, scaled_samples = CSClass.generate_from_dist(dist_type)
                assert type(scaled_samples) == np.ndarray
                assert scaled_samples.shape == (CSClass.number_of_samples,)
                assert dist_res.__name__ == dist_type

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_01(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["random", "normal"],
            )
            unique_sample_points = CSClass.sample_points()
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
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                list_of_distributions=["random", "normal", "uniform"],
                strictly_enforce_gaussian_bounds=True,
            )
            unique_sample_points = CSClass.sample_points()
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
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["random", "normal"],
            )
            unique_sample_points = CSClass.sample_points()
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_with_list_input_creation_mode(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array_list)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                list_of_distributions=["random", "normal", "uniform"],
            )
            unique_sample_points = CSClass.sample_points()
            assert len(CSClass.dist_vector) == len(input_array[0])
            assert unique_sample_points.shape[0] == CSClass.number_of_samples
            assert unique_sample_points.shape[1] == len(input_array[0])
            assert type(unique_sample_points) == np.ndarray

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_sample_points_with_pandas_dataframe_input_selection_mode(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.full_data)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["random", "normal"],
            )
            unique_sample_points = CSClass.sample_points()
            assert len(CSClass.dist_vector) == input_array.shape[1] - 1
            assert unique_sample_points.shape[1] == input_array.shape[1]
            assert type(unique_sample_points) == pd.DataFrame

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_sample_points_with_numpy_array_input_selection_mode(self, array_type):
        for num_samples in [None, 10, 1]:
            input_array = array_type(self.input_array)
            CSClass = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="selection",
                list_of_distributions=["random", "uniform"],
            )
            unique_sample_points = CSClass.sample_points()
            assert len(CSClass.dist_vector) == input_array.shape[1] - 1
            assert unique_sample_points.shape[1] == input_array.shape[1]
            assert type(unique_sample_points) == np.ndarray

    @pytest.mark.integration
    @pytest.mark.parametrize("array_type", [list])
    def test_sample_points_equality_fixed_seed(self, array_type):
        rand_seed = 1000
        for num_samples in [None, 1, 10, 100]:  # Test for different number of samples
            input_array = array_type(self.input_array_list)
            CSClass_A = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                list_of_distributions=["random", "normal", "uniform"],
                rand_seed=rand_seed,
            )
            unique_sample_points_A = CSClass_A.sample_points()

            CSClass_B = CustomSampling(
                input_array,
                number_of_samples=num_samples,
                sampling_type="creation",
                list_of_distributions=["random", "normal", "uniform"],
                rand_seed=rand_seed,
            )
            unique_sample_points_B = CSClass_B.sample_points()

            np.testing.assert_array_equal(
                unique_sample_points_A, unique_sample_points_B
            )


if __name__ == "__main__":
    pytest.main()
