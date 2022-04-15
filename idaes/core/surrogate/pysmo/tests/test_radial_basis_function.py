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
import sys
import os
from unittest.mock import patch

sys.path.append(os.path.abspath(".."))  # current folder is ~/tests\
from idaes.core.surrogate.pysmo.radial_basis_function import (
    RadialBasisFunctions,
    FeatureScaling,
)
import numpy as np
import pandas as pd
from scipy.spatial import distance
import pytest


class TestFeatureScaling:
    test_data_1d = [[x] for x in range(10)]
    test_data_2d = [[x, (x + 1) ** 2] for x in range(10)]
    test_data_3d = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    test_data_3d_constant = [[x, 10, (x + 1) ** 2 + 10] for x in range(10)]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_01(self, array_type):
        input_array = array_type(self.test_data_1d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9]])
        expected_output_2 = np.array([[0]])
        expected_output_1 = np.array(
            (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1.reshape(10, 1))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_02(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 100]])
        expected_output_2 = np.array([[0, 1]])
        expected_output_1 = np.array(
            (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_03(self, array_type):
        input_array = array_type(self.test_data_3d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        expected_output_3 = np.array([[9, 19, 119]])
        expected_output_2 = np.array([[0, 10, 11]])
        expected_output_1 = np.array(
            (input_array - expected_output_2) / (expected_output_3 - expected_output_2)
        )
        np.testing.assert_array_equal(output_3, expected_output_3)
        np.testing.assert_array_equal(output_2, expected_output_2)
        np.testing.assert_array_equal(output_1, expected_output_1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_minmax_04(self, array_type):
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
    @pytest.mark.parametrize("array_type", [list])
    def test_data_scaling_minmax_05(self, array_type):
        input_array = array_type(self.test_data_2d)
        with pytest.raises(TypeError):
            FeatureScaling.data_scaling_minmax(input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_unscaling_minmax_01(self, array_type):
        input_array = array_type(self.test_data_1d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        output_1 = np.array(output_1).reshape(
            output_1.shape[0],
        )
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, np.array(input_array).reshape(10, 1))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_unscaling_minmax_02(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_unscaling_minmax_03(self, array_type):
        input_array = array_type(self.test_data_3d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_unscaling_minmax_04(self, array_type):
        input_array = array_type(self.test_data_3d_constant)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)
        un_output_1 = FeatureScaling.data_unscaling_minmax(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_unscaling_minmax_05(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)

        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_unscaling_minmax_06(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling_minmax(input_array)

        min_array = np.array([[1, 2, 3]])
        max_array = np.array([[5, 6, 7]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling_minmax(output_1, min_array, max_array)


class TestRadialBasisFunction:
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
        RbfClass = RadialBasisFunctions(
            input_array, basis_function=None, solution_method=None, regularization=None
        )
        assert RbfClass.solution_method == "algebraic"
        assert RbfClass.basis_function == "gaussian"
        assert RbfClass.regularization == True

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__02(self, array_type):
        input_array = array_type(self.test_data)
        RbfClass = RadialBasisFunctions(
            input_array,
            basis_function="LineaR",
            solution_method="PyoMo",
            regularization=False,
        )
        assert RbfClass.solution_method == "pyomo"
        assert RbfClass.basis_function == "linear"
        assert RbfClass.regularization == False

    @pytest.mark.unit
    def test__init__03(self):
        with pytest.raises(Exception):
            RbfClass = RadialBasisFunctions(
                [1, 2, 3, 4],
                basis_function="LineaR",
                solution_method="PyoMo",
                regularization=False,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__04(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array, basis_function=None, solution_method=1, regularization=None
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__05(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array,
                basis_function=None,
                solution_method="idaes",
                regularization=None,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__06(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array, basis_function=1, solution_method=None, regularization=None
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__07(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array,
                basis_function="idaes",
                solution_method=None,
                regularization=None,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__08(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array, basis_function=None, solution_method=None, regularization=1
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__09(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array,
                basis_function="LineaR",
                solution_method="PyoMo",
                regularization=False,
                overwrite=1,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__10(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array,
                basis_function="LineaR",
                solution_method="PyoMo",
                regularization=False,
                fname="solution.pkl",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__11(self, array_type):
        with pytest.raises(Exception):
            input_array = array_type(self.test_data)
            RbfClass = RadialBasisFunctions(
                input_array,
                basis_function="LineaR",
                solution_method="PyoMo",
                regularization=False,
                fname=1,
            )

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__12(self, array_type):
        file_name = "test_filename.pickle"
        input_array = array_type(self.test_data)
        RbfClass1 = RadialBasisFunctions(
            input_array,
            basis_function="LineaR",
            solution_method="PyoMo",
            regularization=False,
            fname=file_name,
            overwrite=True,
        )
        p = RbfClass1.get_feature_vector()
        results = RbfClass1.rbf_training()
        RbfClass2 = RadialBasisFunctions(
            input_array,
            basis_function="LineaR",
            solution_method="PyoMo",
            regularization=False,
            fname=file_name,
            overwrite=True,
        )
        assert RbfClass1.filename == RbfClass2.filename

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test__init__14(self, array_type):
        input_array = array_type(self.test_data)
        file_name1 = "test_filename1.pickle"
        file_name2 = "test_filename2.pickle"
        RbfClass1 = RadialBasisFunctions(
            input_array,
            basis_function="LineaR",
            solution_method="PyoMo",
            regularization=False,
            fname=file_name1,
            overwrite=True,
        )
        p = RbfClass1.get_feature_vector()
        RbfClass1.training()
        RbfClass2 = RadialBasisFunctions(
            input_array,
            basis_function="LineaR",
            solution_method="PyoMo",
            regularization=False,
            fname=file_name2,
            overwrite=True,
        )
        assert RbfClass1.filename == file_name1
        assert RbfClass2.filename == file_name2

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_r2_distance(self, array_type):
        input_array = array_type(self.training_data)
        u = np.array([[0.1, 0.9]])
        data_feed = RadialBasisFunctions(input_array)
        output = data_feed.r2_distance(u)
        scaled = FeatureScaling.data_scaling_minmax(input_array)
        scaled = scaled[0]
        scaled_x = np.array(scaled)[:, :-1]
        expected_output = np.sqrt(np.sum(np.square(scaled_x - u), axis=1))
        np.testing.assert_almost_equal(expected_output, output, decimal=6)

    @pytest.mark.unit
    def test_gaussian_basis_transformation(self):
        d_vec = np.array(
            [
                [0, 0],
                [5e-6, 7e-6],
                [0.005, 0.007],
                [0.05, 0.07],
                [0.5, 0.7],
                [5, 7],
                [50, 70],
            ]
        )
        shape_list = [0.001, 1, 1000]
        expected_output_1 = np.exp(-1 * ((d_vec * shape_list[0]) ** 2))
        expected_output_2 = np.exp(-1 * ((d_vec * shape_list[1]) ** 2))
        expected_output_3 = np.exp(-1 * ((d_vec * shape_list[2]) ** 2))
        output_1 = RadialBasisFunctions.gaussian_basis_transformation(
            d_vec, shape_list[0]
        )
        output_2 = RadialBasisFunctions.gaussian_basis_transformation(
            d_vec, shape_list[1]
        )
        output_3 = RadialBasisFunctions.gaussian_basis_transformation(
            d_vec, shape_list[2]
        )
        np.testing.assert_array_equal(expected_output_1, output_1)
        np.testing.assert_array_equal(expected_output_2, output_2)
        np.testing.assert_array_equal(expected_output_3, output_3)

    @pytest.mark.unit
    def test_linear_transformation(self):
        d_vec = np.array(
            [
                [0, 0],
                [5e-6, 7e-6],
                [0.005, 0.007],
                [0.05, 0.07],
                [0.5, 0.7],
                [5, 7],
                [50, 70],
            ]
        )
        output_1 = RadialBasisFunctions.linear_transformation(d_vec)
        np.testing.assert_array_equal(d_vec, output_1)

    @pytest.mark.unit
    def test_cubic_transformation(self):
        d_vec = np.array(
            [
                [0, 0],
                [5e-6, 7e-6],
                [0.005, 0.007],
                [0.05, 0.07],
                [0.5, 0.7],
                [5, 7],
                [50, 70],
            ]
        )
        expected_output = d_vec**3
        output = RadialBasisFunctions.cubic_transformation(d_vec)
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    def test_multiquadric_basis_transformation(self):
        d_vec = np.array(
            [
                [0, 0],
                [5e-6, 7e-6],
                [0.005, 0.007],
                [0.05, 0.07],
                [0.5, 0.7],
                [5, 7],
                [50, 70],
            ]
        )
        shape_list = [0.001, 1, 1000]
        expected_output_1 = np.sqrt(((d_vec * shape_list[0]) ** 2) + 1)
        expected_output_2 = np.sqrt(((d_vec * shape_list[1]) ** 2) + 1)
        expected_output_3 = np.sqrt(((d_vec * shape_list[2]) ** 2) + 1)
        output_1 = RadialBasisFunctions.multiquadric_basis_transformation(
            d_vec, shape_list[0]
        )
        output_2 = RadialBasisFunctions.multiquadric_basis_transformation(
            d_vec, shape_list[1]
        )
        output_3 = RadialBasisFunctions.multiquadric_basis_transformation(
            d_vec, shape_list[2]
        )
        np.testing.assert_array_equal(expected_output_1, output_1)
        np.testing.assert_array_equal(expected_output_2, output_2)
        np.testing.assert_array_equal(expected_output_3, output_3)

    @pytest.mark.unit
    def test_inverse_multiquadric_basis_transformation(self):
        d_vec = np.array(
            [
                [0, 0],
                [5e-6, 7e-6],
                [0.005, 0.007],
                [0.05, 0.07],
                [0.5, 0.7],
                [5, 7],
                [50, 70],
            ]
        )
        shape_list = [0.001, 1, 1000]
        expected_output_1 = 1 / np.sqrt(((d_vec * shape_list[0]) ** 2) + 1)
        expected_output_2 = 1 / np.sqrt(((d_vec * shape_list[1]) ** 2) + 1)
        expected_output_3 = 1 / np.sqrt(((d_vec * shape_list[2]) ** 2) + 1)
        output_1 = RadialBasisFunctions.inverse_multiquadric_basis_transformation(
            d_vec, shape_list[0]
        )
        output_2 = RadialBasisFunctions.inverse_multiquadric_basis_transformation(
            d_vec, shape_list[1]
        )
        output_3 = RadialBasisFunctions.inverse_multiquadric_basis_transformation(
            d_vec, shape_list[2]
        )
        np.testing.assert_array_equal(expected_output_1, output_1)
        np.testing.assert_array_equal(expected_output_2, output_2)
        np.testing.assert_array_equal(expected_output_3, output_3)

    @pytest.mark.unit
    def test_thin_plate_spline_transformation(self):
        d_vec = np.array(
            [
                [5e-6, 7e-6],
                [0.005, 0.007],
                [0.05, 0.07],
                [0.5, 0.7],
                [5, 7],
                [50, 70],
                [50, np.NaN],
            ]
        )
        expected_output = np.nan_to_num(d_vec**2 * np.log(d_vec))
        output = RadialBasisFunctions.thin_plate_spline_transformation(d_vec)
        np.testing.assert_array_equal(expected_output, output)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_basis_generation(self, array_type):
        input_array = array_type(self.training_data)
        scaled = FeatureScaling.data_scaling_minmax(input_array[0:3])
        scaled = scaled[0]
        scaled_x = np.array(scaled)[:, :-1]
        distance_array = distance.cdist(scaled_x, scaled_x, "euclidean")

        # Linear
        data_feed_01 = RadialBasisFunctions(input_array[0:3], basis_function="linear")
        expected_output_1 = distance_array
        output_1 = data_feed_01.basis_generation(2)
        np.testing.assert_array_equal(expected_output_1, output_1)

        # Cubic
        data_feed_02 = RadialBasisFunctions(input_array[0:3], basis_function="cubic")
        expected_output_2 = distance_array**3
        output_2 = data_feed_02.basis_generation(2)
        np.testing.assert_array_equal(expected_output_2, output_2)

        # # Spline
        data_feed_03 = RadialBasisFunctions(input_array[0:3], basis_function="spline")
        expected_output_3 = np.nan_to_num(distance_array**2 * np.log(distance_array))
        output_3 = data_feed_03.basis_generation(2)
        np.testing.assert_array_equal(expected_output_3, output_3)

        # # Gaussian
        data_feed_04 = RadialBasisFunctions(input_array[0:3], basis_function="gaussian")
        shape_value = 2
        expected_output_4 = np.exp(-1 * ((distance_array * shape_value) ** 2))
        output_4 = data_feed_04.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_4, output_4)

        # # Multiquadric
        data_feed_05 = RadialBasisFunctions(input_array[0:3], basis_function="mq")
        shape_value = 2
        expected_output_5 = np.sqrt(((distance_array * shape_value) ** 2) + 1)
        output_5 = data_feed_05.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_5, output_5)

        # # Inverse multiquadric
        data_feed_06 = RadialBasisFunctions(input_array[0:3], basis_function="imq")
        shape_value = 2
        expected_output_6 = 1 / np.sqrt(((distance_array * shape_value) ** 2) + 1)
        output_6 = data_feed_06.basis_generation(shape_value)
        np.testing.assert_array_equal(expected_output_6, output_6)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cost_function_01(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cost_function_02(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cost_function_03(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
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

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_gradient_function_01(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
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
            [[-97], [-635], [-635], [-5246.875], [-5246.875], [-3925]]
        )
        expected_value = expected_value.reshape(
            expected_value.shape[0],
        )
        output_1 = RadialBasisFunctions.gradient_function(theta, x_vector, y)
        np.testing.assert_equal(output_1, expected_value)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_gradient_function_02(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array(
            [[4.5], [3], [3], [1], [1], [0]]
        )  # coefficients in (x1 + 1.5)^2 + (x2 + 1.5) ^ 2
        theta = theta.reshape(
            theta.shape[0],
        )
        expected_value = np.array(
            [[12.5], [75], [75], [593.75], [593.75], [437.5]]
        )  # Calculated externally: see Excel sheet
        expected_value = expected_value.reshape(
            expected_value.shape[0],
        )
        output_1 = RadialBasisFunctions.gradient_function(theta, x_vector, y)
        np.testing.assert_equal(output_1, expected_value)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_gradient_function_03(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array(
            [[2], [2], [2], [1], [1], [0]]
        )  # Actual coefficients in (x1 + 1)^2 + (x2 + 1) ^ 2
        theta = theta.reshape(
            theta.shape[0],
        )
        expected_value = np.array(
            [[0], [0], [0], [0], [0], [0]]
        )  # Calculated externally: see Excel sheet
        expected_value = expected_value.reshape(
            expected_value.shape[0],
        )
        output_1 = RadialBasisFunctions.gradient_function(theta, x_vector, y)
        np.testing.assert_equal(output_1, expected_value)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_bfgs_parameter_optimization_01(self, array_type):
        input_array = np.array(
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
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = (
            x[
                :,
            ]
            ** 2
        )
        x_vector[:, 1] = x[
            :,
        ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.0], [2.0], [1.0]]).reshape(
            3,
        )
        data_feed = RadialBasisFunctions(
            array_type(self.test_data), basis_function="linear", solution_method="bfgs"
        )
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == "bfgs"
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array])
    @pytest.mark.parametrize("array_type2", [pd.DataFrame])
    def test_bfgs_parameter_optimization_02(self, array_type1, array_type2):
        input_array = array_type1(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
        x_vector = np.zeros((x.shape[0], 6))
        x_vector[:, 0] = x[:, 0] ** 2
        x_vector[:, 1] = x[:, 1] ** 2
        x_vector[:, 2] = x[:, 0]
        x_vector[:, 3] = x[:, 1]
        x_vector[:, 4] = x[:, 1] * x[:, 0]
        x_vector[:, 5] = 1
        expected_value = np.array([[1.0], [1.0], [2.0], [2.0], [0.0], [2.0]]).reshape(
            6,
        )
        data_feed = RadialBasisFunctions(
            array_type2(self.full_data), solution_method="bfgs"
        )
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == "bfgs"
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    def test_explicit_linear_algebra_solution_01(self):
        input_array = np.array(
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
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = (
            x[
                :,
            ]
            ** 2
        )
        x_vector[:, 1] = x[
            :,
        ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.0], [2.0], [1.0]]).reshape(
            3,
        )
        output_1 = RadialBasisFunctions.explicit_linear_algebra_solution(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_explicit_linear_algebra_solution_02(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
        x_vector = np.zeros((x.shape[0], 6))
        x_vector[:, 0] = x[:, 0] ** 2
        x_vector[:, 1] = x[:, 1] ** 2
        x_vector[:, 2] = x[:, 0]
        x_vector[:, 3] = x[:, 1]
        x_vector[:, 4] = x[:, 1] * x[:, 0]
        x_vector[:, 5] = 1
        expected_value = np.array([[1.0], [1.0], [2.0], [2.0], [0.0], [2.0]]).reshape(
            6,
        )
        output_1 = RadialBasisFunctions.explicit_linear_algebra_solution(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    def test_pyomo_optimization_01(self):
        input_array = np.array(
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
        x = input_array[:, 0]
        y = input_array[:, 1]
        x_vector = np.zeros((x.shape[0], 3))
        x_vector[:, 0] = (
            x[
                :,
            ]
            ** 2
        )
        x_vector[:, 1] = x[
            :,
        ]
        x_vector[:, 2] = 1
        expected_value = np.array([[1.0], [2.0], [1.0]])
        output_1 = RadialBasisFunctions.pyomo_optimization(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_pyomo_optimization_02(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1]
        x_vector = np.zeros((x.shape[0], 6))
        x_vector[:, 0] = x[:, 0] ** 2
        x_vector[:, 1] = x[:, 1] ** 2
        x_vector[:, 2] = x[:, 0]
        x_vector[:, 3] = x[:, 1]
        x_vector[:, 4] = x[:, 1] * x[:, 0]
        x_vector[:, 5] = 1
        expected_value = np.array([[1.0], [1.0], [2.0], [2.0], [0.0], [2.0]])
        output_1 = RadialBasisFunctions.pyomo_optimization(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_error_calculation_01(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1].reshape(input_array.shape[0], 1)
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
        expected_value_2 = expected_value_1**0.5
        output_1, output_2, _ = RadialBasisFunctions.error_calculation(
            theta, x_vector, y
        )
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_error_calculation_02(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1].reshape(input_array.shape[0], 1)
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array(
            [[4.5], [3], [3], [1], [1], [0]]
        )  # coefficients in (x1 + 1.5)^2 + (x2 + 1.5) ^ 2
        expected_value_1 = 2 * 90.625  # Calculated externally as sum(dy^2) / 2m
        expected_value_2 = expected_value_1**0.5
        output_1, output_2, _ = RadialBasisFunctions.error_calculation(
            theta, x_vector, y
        )
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_error_calculation_03(self, array_type):
        input_array = array_type(self.training_data)
        x = input_array[:, :-1]
        y = input_array[:, -1].reshape(input_array.shape[0], 1)
        x_data_nr = x.shape[0]
        x_data_nc = 6
        x_vector = np.zeros((x_data_nr, x_data_nc))
        x_vector[:, 0] = 1
        x_vector[:, 1] = x[:, 0]
        x_vector[:, 2] = x[:, 1]
        x_vector[:, 3] = x[:, 0] ** 2
        x_vector[:, 4] = x[:, 1] ** 2
        x_vector[:, 5] = x[:, 0] * x[:, 1]
        theta = np.array(
            [[2], [2], [2], [1], [1], [0]]
        )  # Actual coefficients in (x1 + 1)^2 + (x2 + 1) ^ 2
        expected_value_1 = 2 * 0  # Value should return zero for exact solution
        expected_value_2 = expected_value_1**0.5
        output_1, output_2, _ = RadialBasisFunctions.error_calculation(
            theta, x_vector, y
        )
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2

    @pytest.mark.unit
    def test_r2_calculation_01(self):
        y_actual = np.array([[1], [4], [9], [16], [25], [36], [49], [64], [81], [100]])
        y_pred = y_actual * 1.05
        expected_output = 0.993974359  # Evaluated in Excel
        output = RadialBasisFunctions.r2_calculation(y_actual, y_pred)
        assert round(abs(expected_output - output), 7) == 0

    @pytest.mark.unit
    def test_r2_calculation_02(self):
        y_actual = np.array([[1], [4], [9], [16], [25], [36], [49], [64], [81], [100]])
        y_pred = y_actual * 1.50
        expected_output = 0.3974358974  # Evaluated in Excel
        output = RadialBasisFunctions.r2_calculation(y_actual, y_pred)
        assert round(abs(expected_output - output), 7) == 0

    def mock_basis_generation(self, r):
        return np.ones((self.x_data.shape[0], self.x_data.shape[0]))

    def mock_optimization(self, x, y):
        return 500 * np.ones((x.shape[0], 1))

    @patch.object(RadialBasisFunctions, "basis_generation", mock_basis_generation)
    @patch.object(
        RadialBasisFunctions, "explicit_linear_algebra_solution", mock_optimization
    )
    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_loo_error_estimation_with_rippa_method_01(self, array_type):
        input_array = array_type(self.training_data)
        reg_param = 0.1
        shape_factor = 1
        expected_x = np.ones((input_array.shape[0], input_array.shape[0])) + (
            reg_param * np.eye(input_array.shape[0], input_array.shape[0])
        )
        expected_inverse_x = np.diag(np.linalg.pinv(expected_x))
        expected_radial_weights = 500 * np.ones((input_array.shape[0], 1))
        expected_errors = np.linalg.norm(
            expected_radial_weights
            / (expected_inverse_x.reshape(expected_inverse_x.shape[0], 1))
        )

        data_feed = RadialBasisFunctions(input_array, solution_method="algebraic")
        _, output_1, output_2 = data_feed.loo_error_estimation_with_rippa_method(
            shape_factor, reg_param
        )
        assert output_1 == np.linalg.cond(expected_x)
        np.testing.assert_array_equal(output_2, expected_errors)

    @patch.object(RadialBasisFunctions, "basis_generation", mock_basis_generation)
    @patch.object(RadialBasisFunctions, "pyomo_optimization", mock_optimization)
    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_loo_error_estimation_with_rippa_method_02(self, array_type):
        input_array = array_type(self.training_data)
        reg_param = 0.1
        shape_factor = 1
        expected_x = np.ones((input_array.shape[0], input_array.shape[0])) + (
            reg_param * np.eye(input_array.shape[0], input_array.shape[0])
        )
        expected_inverse_x = np.diag(np.linalg.pinv(expected_x))
        expected_radial_weights = 500 * np.ones((input_array.shape[0], 1))
        expected_errors = np.linalg.norm(
            expected_radial_weights
            / (expected_inverse_x.reshape(expected_inverse_x.shape[0], 1))
        )

        data_feed = RadialBasisFunctions(input_array, solution_method="pyomo")
        _, output_1, output_2 = data_feed.loo_error_estimation_with_rippa_method(
            shape_factor, reg_param
        )
        assert output_1 == np.linalg.cond(expected_x)
        np.testing.assert_array_equal(output_2, expected_errors)

    @patch.object(RadialBasisFunctions, "basis_generation", mock_basis_generation)
    @patch.object(
        RadialBasisFunctions, "bfgs_parameter_optimization", mock_optimization
    )
    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_loo_error_estimation_with_rippa_method_03(self, array_type):
        input_array = array_type(self.training_data)
        reg_param = 0.1
        shape_factor = 1
        expected_x = np.ones((input_array.shape[0], input_array.shape[0])) + (
            reg_param * np.eye(input_array.shape[0], input_array.shape[0])
        )
        expected_inverse_x = np.diag(np.linalg.pinv(expected_x))
        expected_radial_weights = 500 * np.ones((input_array.shape[0], 1))
        expected_errors = np.linalg.norm(
            expected_radial_weights
            / (expected_inverse_x.reshape(expected_inverse_x.shape[0], 1))
        )

        data_feed = RadialBasisFunctions(input_array, solution_method="bfgs")
        _, output_1, output_2 = data_feed.loo_error_estimation_with_rippa_method(
            shape_factor, reg_param
        )
        assert output_1 == np.linalg.cond(expected_x)
        np.testing.assert_array_equal(output_2, expected_errors)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_01(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function=None, solution_method=None, regularization=False
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_02(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="cubic",
            solution_method=None,
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_03(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="linear",
            solution_method=None,
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_04(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="spline",
            solution_method=None,
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_05(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="gaussian",
            solution_method=None,
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_06(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="mq", solution_method=None, regularization=False
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_07(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="imq",
            solution_method=None,
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_08(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function=None,
            solution_method="algebraic",
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_09(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function=None,
            solution_method="BFGS",
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_10(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function=None,
            solution_method="pyomo",
            regularization=False,
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_leave_one_out_crossvalidation_11(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function=None, solution_method=None, regularization=True
        )
        r_best, lambda_best, error_best = data_feed.leave_one_out_crossvalidation()
        if (
            (data_feed.basis_function == "gaussian")
            or (data_feed.basis_function == "mq")
            or (data_feed.basis_function.lower() == "imq")
        ):
            r_set = [
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1.0,
                2.0,
                5.0,
                7.5,
                10.0,
                20.0,
                50.0,
                75.0,
                100.0,
                200.0,
                500.0,
                1000.0,
            ]
        else:
            r_set = [0]
        if data_feed.regularization is True:
            reg_parameter = [
                0.00001,
                0.00002,
                0.00005,
                0.000075,
                0.0001,
                0.0002,
                0.0005,
                0.00075,
                0.001,
                0.002,
                0.005,
                0.0075,
                0.01,
                0.02,
                0.05,
                0.075,
                0.1,
                0.2,
                0.5,
                0.75,
                1,
            ]
        elif data_feed.regularization is False:
            reg_parameter = [0]
        _, _, expected_errors = data_feed.loo_error_estimation_with_rippa_method(
            r_best, lambda_best
        )
        assert (r_best in r_set) == True
        assert (lambda_best in reg_parameter) == True
        assert error_best == expected_errors

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_training_01(self, array_type):
        input_array = array_type(self.test_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function=None,
            solution_method="algebraic",
            regularization=False,
        )
        results = data_feed.training()
        best_r_value, best_lambda_param, _ = data_feed.leave_one_out_crossvalidation()
        x_transformed = data_feed.basis_generation(best_r_value)
        x_transformed = x_transformed + (
            best_lambda_param * np.eye(x_transformed.shape[0], x_transformed.shape[1])
        )
        x_condition_number = np.linalg.cond(x_transformed)
        if data_feed.solution_method == "algebraic":
            radial_weights = data_feed.explicit_linear_algebra_solution(
                x_transformed, data_feed.y_data
            )
        elif data_feed.solution_method == "pyomo":
            radial_weights = data_feed.pyomo_optimization(
                x_transformed, data_feed.y_data
            )
        elif data_feed.solution_method == "bfgs":
            radial_weights = data_feed.bfgs_parameter_optimization(
                x_transformed, data_feed.y_data
            )
        radial_weights = radial_weights.reshape(radial_weights.shape[0], 1)
        (
            training_ss_error,
            rmse_error,
            y_training_predictions_scaled,
        ) = data_feed.error_calculation(radial_weights, x_transformed, data_feed.y_data)
        r_square = data_feed.r2_calculation(
            data_feed.y_data, y_training_predictions_scaled
        )
        y_training_predictions = data_feed.data_min[
            0, -1
        ] + y_training_predictions_scaled * (
            data_feed.data_max[0, -1] - data_feed.data_min[0, -1]
        )
        np.testing.assert_array_equal(radial_weights, results.weights)
        np.testing.assert_array_equal(best_r_value, results.sigma)
        np.testing.assert_array_equal(best_lambda_param, results.regularization)
        np.testing.assert_array_equal(data_feed.centres, results.centres)
        np.testing.assert_array_equal(
            y_training_predictions, results.output_predictions
        )
        np.testing.assert_array_equal(rmse_error, results.rmse)
        np.testing.assert_array_equal(x_condition_number, results.condition_number)
        np.testing.assert_array_equal(data_feed.regularization, results.regularization)
        np.testing.assert_array_equal(r_square, results.R2)
        assert data_feed.basis_function == results.basis_function
        np.testing.assert_array_equal(data_feed.data_min[:, :-1], results.x_data_min)
        np.testing.assert_array_equal(data_feed.data_max[:, :-1], results.x_data_max)
        np.testing.assert_array_equal(data_feed.data_min[:, -1], results.y_data_min)
        np.testing.assert_array_equal(data_feed.data_max[:, -1], results.y_data_max)
        assert results.solution_status == "ok"

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_training_02(self, array_type):
        input_array = array_type(self.test_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function=None,
            solution_method="pyomo",
            regularization=False,
        )
        data_feed.training()
        with pytest.warns(Warning):
            results = data_feed.training()
            assert data_feed.solution_status == "unstable solution"

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_training_03(self, array_type):
        input_array = array_type(self.test_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function=None,
            solution_method="bfgs",
            regularization=False,
        )
        data_feed.training()
        with pytest.warns(Warning):
            data_feed.training()
            assert data_feed.solution_status == "unstable solution"

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_predict_output_01(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="linear", regularization=False
        )
        results = data_feed.training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.predict_output(x_test)
        data_minimum = results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum) / scale
        x_test = x_pred_scaled.reshape(x_test.shape)
        distance_vec = distance.cdist(x_test, results.centres, "euclidean")
        expected_output = np.matmul(distance_vec, results.weights)
        expected_output = results.y_data_min + expected_output * (
            results.y_data_max - results.y_data_min
        )
        assert expected_output == output

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_predict_output_02(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="cubic", regularization=False
        )
        results = data_feed.training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.predict_output(x_test)
        data_minimum = data_feed.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum) / scale
        x_test = x_pred_scaled.reshape(x_test.shape)
        distance_vec = distance.cdist(x_test, results.centres, "euclidean")
        expected_output = np.matmul(distance_vec**3, results.weights)
        expected_output = results.y_data_min + expected_output * (
            results.y_data_max - results.y_data_min
        )
        assert expected_output == output

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_predict_output_03(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="gaussian", regularization=False
        )
        results = data_feed.training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.predict_output(x_test)
        data_minimum = results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum) / scale
        x_test = x_pred_scaled.reshape(x_test.shape)
        distance_vec = distance.cdist(x_test, results.centres, "euclidean")
        expected_output = np.matmul(
            np.exp(-1 * ((distance_vec * results.sigma) ** 2)), results.weights
        )
        expected_output = results.y_data_min + expected_output * (
            results.y_data_max - results.y_data_min
        )
        assert expected_output == output

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_predict_output_04(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="imq", regularization=False
        )
        results = data_feed.training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.predict_output(x_test)
        data_minimum = data_feed.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum) / scale
        x_test = x_pred_scaled.reshape(x_test.shape)
        distance_vec = distance.cdist(x_test, results.centres, "euclidean")
        expected_output = np.matmul(
            1 / np.sqrt(((distance_vec * results.sigma) ** 2) + 1), results.weights
        )
        expected_output = results.y_data_min + expected_output * (
            results.y_data_max - results.y_data_min
        )
        assert expected_output == output

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_predict_output_05(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="mq", regularization=False
        )
        results = data_feed.training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.predict_output(x_test)
        data_minimum = results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum) / scale
        x_test = x_pred_scaled.reshape(x_test.shape)
        distance_vec = distance.cdist(x_test, results.centres, "euclidean")
        expected_output = np.matmul(
            np.sqrt(((distance_vec * results.sigma) ** 2) + 1), results.weights
        )
        expected_output = results.y_data_min + expected_output * (
            results.y_data_max - results.y_data_min
        )
        assert expected_output == output

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_predict_output_06(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="spline", regularization=False
        )
        results = data_feed.training()
        x_test = np.array([[0, 7.5]])
        output = data_feed.predict_output(x_test)
        data_minimum = results.x_data_min
        data_maximum = results.x_data_max
        scale = data_maximum - data_minimum
        scale[scale == 0.0] = 1.0
        x_pred_scaled = (x_test - data_minimum) / scale
        x_test = x_pred_scaled.reshape(x_test.shape)
        distance_vec = distance.cdist(x_test, results.centres, "euclidean")
        expected_output = np.matmul(
            np.nan_to_num(distance_vec**2 * np.log(distance_vec)), results.weights
        )
        expected_output = results.y_data_min + expected_output * (
            results.y_data_max - results.y_data_min
        )
        assert expected_output == output

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_get_feature_vector_01(self, array_type):
        input_array = array_type(self.full_data)
        data_feed = RadialBasisFunctions(input_array, basis_function="linear")
        output = data_feed.get_feature_vector()
        expected_dict = {"x1": 0, "x2": 0}
        assert expected_dict == output.extract_values()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_get_feature_vector_02(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(input_array, basis_function="linear")
        output = data_feed.get_feature_vector()
        expected_dict = {0: 0, 1: 0}
        assert expected_dict == output.extract_values()

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_generate_expression_01(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="linear",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        data_feed.training()
        lv = []
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = data_feed.generate_expression((lv))

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_generate_expression_02(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="cubic",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        data_feed.training()
        lv = []
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = data_feed.generate_expression((lv))

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_generate_expression_03(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="gaussian",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        results = data_feed.training()
        lv = []
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.generate_expression((lv))

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_generate_expression_04(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array, basis_function="mq", solution_method=None, regularization=False
        )
        p = data_feed.get_feature_vector()
        results = data_feed.training()
        lv = []
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.generate_expression((lv))

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_generate_expression_05(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="imq",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        results = data_feed.training()
        lv = []
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.generate_expression((lv))

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_rbf_generate_expression_06(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="spline",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        results = data_feed.training()
        lv = []
        for i in p.keys():
            lv.append(p[i])
        rbf_expr = results.generate_expression((lv))

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_pickle_load01(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="spline",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        data_feed.training()
        data_feed.pickle_load(data_feed.filename)

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_pickle_load02(self, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="spline",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        data_feed.training()
        with pytest.raises(Exception):
            data_feed.pickle_load("file_not_existing.pickle")

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @patch("matplotlib.pyplot.show")
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_parity_residual_plots(self, mock_show, array_type):
        input_array = array_type(self.training_data)
        data_feed = RadialBasisFunctions(
            input_array,
            basis_function="spline",
            solution_method=None,
            regularization=False,
        )
        p = data_feed.get_feature_vector()
        data_feed.training()
        data_feed.parity_residual_plots()


if __name__ == "__main__":
    pytest.main()
