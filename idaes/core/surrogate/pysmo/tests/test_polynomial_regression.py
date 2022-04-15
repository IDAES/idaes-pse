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
import sys, os
from unittest.mock import patch

sys.path.append(os.path.abspath(".."))  # current folder is ~/tests
from idaes.core.surrogate.pysmo.polynomial_regression import (
    PolynomialRegression,
    FeatureScaling,
)
import numpy as np
import pandas as pd
import pytest


class TestFeatureScaling:
    test_data_1d = [[x] for x in range(10)]
    test_data_2d = [[x, (x + 1) ** 2] for x in range(10)]
    test_data_3d = [[x, x + 10, (x + 1) ** 2 + x + 10] for x in range(10)]
    test_data_3d_constant = [[x, 10, (x + 1) ** 2 + 10] for x in range(10)]

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_data_scaling_01(self, array_type):
        input_array = array_type(self.test_data_1d)

        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
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
    def test_data_scaling_02(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
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
    def test_data_scaling_03(self, array_type):
        input_array = array_type(self.test_data_3d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
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
    def test_data_scaling_04(self, array_type):
        input_array = array_type(self.test_data_3d_constant)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
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
    def test_data_scaling_05(self, array_type):
        input_array = array_type(self.test_data_2d)
        with pytest.raises(TypeError):
            FeatureScaling.data_scaling(input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_data_unscaling_01(self, array_type):
        input_array = array_type(self.test_data_1d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        output_1 = output_1.reshape(
            output_1.shape[0],
        )
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array.reshape(10, 1))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_data_unscaling_02(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_data_unscaling_03(self, array_type):
        input_array = array_type(self.test_data_3d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_data_unscaling_04(self, array_type):
        input_array = array_type(self.test_data_3d_constant)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)
        un_output_1 = FeatureScaling.data_unscaling(output_1, output_2, output_3)
        np.testing.assert_array_equal(un_output_1, input_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_data_unscaling_05(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)

        min_array = np.array([[1]])
        max_array = np.array([[5]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling(output_1, min_array, max_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_data_unscaling_06(self, array_type):
        input_array = array_type(self.test_data_2d)
        output_1, output_2, output_3 = FeatureScaling.data_scaling(input_array)

        min_array = np.array([[1, 2, 3]])
        max_array = np.array([[5, 6, 7]])
        with pytest.raises(IndexError):
            FeatureScaling.data_unscaling(output_1, min_array, max_array)


class TestPolynomialRegression:
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
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__01(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)

        PolyClass = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=5
        )
        assert PolyClass.max_polynomial_order == 5
        assert (
            PolyClass.number_of_crossvalidations == 3
        )  # Default number of cross-validations
        assert PolyClass.no_adaptive_samples == 4  # Default number of adaptive samples
        assert PolyClass.fraction_training == 0.75  # Default training split
        assert (
            PolyClass.max_fraction_training_samples == 0.5
        )  # Default fraction for the maximum number of training samples
        assert PolyClass.max_iter == 10  # Default maximum number of iterations
        assert PolyClass.solution_method == "pyomo"  # Default solution_method
        assert PolyClass.multinomials == 1  # Default multinomials

    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    @pytest.mark.unit
    def test__init__02(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        PolyClass = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            number_of_crossvalidations=5,
            no_adaptive_samples=6,
            training_split=0.5,
            max_fraction_training_samples=0.4,
            max_iter=20,
            solution_method="MLe",
            multinomials=0,
        )
        assert PolyClass.max_polynomial_order == 3
        assert (
            PolyClass.number_of_crossvalidations == 5
        )  # Default number of cross-validations
        assert PolyClass.no_adaptive_samples == 6  # Default number of adaptive samples
        assert PolyClass.fraction_training == 0.5  # Default training split
        assert (
            PolyClass.max_fraction_training_samples == 0.4
        )  # Default fraction for the maximum number of training samples
        assert PolyClass.max_iter == 20  # Default maximum number of iterations
        assert (
            PolyClass.solution_method == "mle"
        )  # Default solution_method, doesn't matter lower / upper characters
        assert PolyClass.multinomials == 0  # Default multinomials

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [list])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__03(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(ValueError):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=5
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [list])
    def test__init__04(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(ValueError):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=5
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__05(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points_large)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=5
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__06(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data_3d)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=5
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__07(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points_3d)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=5
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__08(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data_1d)
        regression_data_input = array_type2(self.sample_points_1d)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=5
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__09(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.warns(Warning):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                number_of_crossvalidations=11,
            )
            assert (
                PolyClass.number_of_crossvalidations == 11
            )  # Default number of cross-validations

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__10(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=1.2
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__11(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data_large)
        regression_data_input = array_type2(self.sample_points_large)
        with pytest.warns(Warning):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=11
            )
            assert PolyClass.max_polynomial_order == 10

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__12(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                training_split=1,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__13(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                training_split=-1,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__14(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                max_fraction_training_samples=1.2,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__15(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                max_fraction_training_samples=-1.2,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__16(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        PolyClass = PolynomialRegression(
            regression_data_input,
            regression_data_input,
            maximum_polynomial_order=5,
            max_iter=100,
        )
        assert PolyClass.max_iter == 0

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__17(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        PolyClass = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=5,
            no_adaptive_samples=0,
            max_iter=100,
        )
        assert PolyClass.max_iter == 0

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__18(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                number_of_crossvalidations=1.2,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__19(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                no_adaptive_samples=4.2,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__20(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                max_iter=4.2,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__21(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=15
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__22(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                solution_method=1,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__23(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                solution_method="idaes",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__24(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=5,
                multinomials=3,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__25(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input, regression_data_input, maximum_polynomial_order=-2
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__26(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=3,
                number_of_crossvalidations=-3,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__27(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=3,
                no_adaptive_samples=-3,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__28(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=3,
                max_iter=-3,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__29(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=3,
                overwrite=1,
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__30(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=3,
                fname="solution.pkl",
            )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__31(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        with pytest.raises(Exception):
            PolyClass = PolynomialRegression(
                original_data_input,
                regression_data_input,
                maximum_polynomial_order=3,
                fname=1,
            )

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__32(self, array_type1, array_type2):
        file_name = "sol_check.pickle"
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)

        PolyClass1 = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            fname=file_name,
            overwrite=True,
        )
        PolyClass1.get_feature_vector()
        results = PolyClass1.polynomial_regression_fitting()
        PolyClass2 = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            fname=file_name,
            overwrite=True,
        )
        assert PolyClass1.filename == PolyClass2.filename

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__33(self, array_type1, array_type2):
        file_name1 = "sol_check1.pickle"
        file_name2 = "sol_check2.pickle"
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)

        PolyClass1 = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            fname=file_name1,
            overwrite=True,
        )
        PolyClass1.get_feature_vector()
        results = PolyClass1.polynomial_regression_fitting()
        PolyClass2 = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            fname=file_name2,
            overwrite=True,
        )
        assert PolyClass1.filename == file_name1
        assert PolyClass2.filename == file_name2

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test__init__34(self, array_type1, array_type2):
        file_name = "sol_check.pickle"
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        PolyClass1 = PolynomialRegression(
            original_data_input,
            regression_data_input,
            fname=file_name,
            maximum_polynomial_order=3,
            overwrite=True,
        )
        PolyClass1.get_feature_vector()
        results = PolyClass1.polynomial_regression_fitting()
        PolygClass2 = PolynomialRegression(
            original_data_input,
            regression_data_input,
            fname=file_name,
            maximum_polynomial_order=3,
            overwrite=True,
        )
        assert PolyClass1.filename == PolygClass2.filename

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test_training_test_data_creation_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)

        PolyClass = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=5,
            training_split=0.01,
        )
        with pytest.raises(Exception):
            training_data, cross_val_data = PolyClass.training_test_data_creation()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [np.array, pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array, pd.DataFrame])
    def test_training_test_data_creation_02(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)

        PolyClass = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=5,
            training_split=0.99,
        )
        with pytest.raises(Exception):
            training_data, cross_val_data = PolyClass.training_test_data_creation()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_training_test_data_creation_03(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)

        PolyClass = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=5
        )

        training_data, cross_val_data = PolyClass.training_test_data_creation()

        expected_training_size = int(
            np.around(PolyClass.number_of_samples * PolyClass.fraction_training)
        )
        expected_test_size = PolyClass.regression_data.shape[0] - expected_training_size

        assert len(training_data) == PolyClass.number_of_crossvalidations
        assert len(cross_val_data) == PolyClass.number_of_crossvalidations

        for i in range(1, PolyClass.number_of_crossvalidations + 1):
            assert (
                training_data["training_set_" + str(i)].shape[0]
                == expected_training_size
            )
            assert cross_val_data["test_set_" + str(i)].shape[0] == expected_test_size

            concat_01 = np.concatenate(
                (
                    training_data["training_set_" + str(i)],
                    cross_val_data["test_set_" + str(i)],
                ),
                axis=0,
            )
            sample_data_sorted = regression_data_input[
                np.lexsort(
                    (
                        regression_data_input[:, 2],
                        regression_data_input[:, 1],
                        regression_data_input[:, 0],
                    )
                )
            ]
            concat_01_sorted = concat_01[
                np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))
            ]
            np.testing.assert_equal(sample_data_sorted, concat_01_sorted)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_training_test_data_creation_04(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        PolyClass = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            number_of_crossvalidations=5,
            no_adaptive_samples=6,
            training_split=0.5,
            max_fraction_training_samples=0.4,
            max_iter=20,
            solution_method="MLe",
            multinomials=0,
        )
        training_data, cross_val_data = PolyClass.training_test_data_creation()
        expected_training_size = int(
            np.around(PolyClass.number_of_samples * PolyClass.fraction_training)
        )
        expected_test_size = PolyClass.regression_data.shape[0] - expected_training_size
        assert len(training_data) == PolyClass.number_of_crossvalidations
        assert len(cross_val_data) == PolyClass.number_of_crossvalidations
        for i in range(1, PolyClass.number_of_crossvalidations + 1):
            assert (
                training_data["training_set_" + str(i)].shape[0]
                == expected_training_size
            )
            assert cross_val_data["test_set_" + str(i)].shape[0] == expected_test_size
            concat_01 = np.concatenate(
                (
                    training_data["training_set_" + str(i)],
                    cross_val_data["test_set_" + str(i)],
                ),
                axis=0,
            )
            sample_data_sorted = regression_data_input[
                np.lexsort(
                    (
                        regression_data_input[:, 2],
                        regression_data_input[:, 1],
                        regression_data_input[:, 0],
                    )
                )
            ]
            concat_01_sorted = concat_01[
                np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))
            ]
            np.testing.assert_equal(sample_data_sorted, concat_01_sorted)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_training_test_data_creation_05(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        PolyClass = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=5
        )
        additional_data_input = np.array(
            [
                [
                    i**2,
                    ((i + 1) * 2) + ((j + 1) * 2),
                    j**4,
                    ((i + 1) * 2) + ((j + 1) ** 2),
                ]
                for i in range(5)
                for j in range(5)
            ]
        )
        training_data, cross_val_data = PolyClass.training_test_data_creation(
            additional_features=additional_data_input
        )
        expected_training_size = int(
            np.around(PolyClass.number_of_samples * PolyClass.fraction_training)
        )
        expected_test_size = PolyClass.regression_data.shape[0] - expected_training_size
        assert len(training_data) == PolyClass.number_of_crossvalidations * 2
        assert len(cross_val_data) == PolyClass.number_of_crossvalidations * 2
        for i in range(1, PolyClass.number_of_crossvalidations + 1):
            assert (
                training_data["training_set_" + str(i)].shape[0]
                == expected_training_size
            )
            assert (
                training_data["training_extras_" + str(i)].shape[0]
                == expected_training_size
            )
            assert cross_val_data["test_set_" + str(i)].shape[0] == expected_test_size
            assert (
                cross_val_data["test_extras_" + str(i)].shape[0] == expected_test_size
            )
            concat_01 = np.concatenate(
                (
                    training_data["training_set_" + str(i)],
                    cross_val_data["test_set_" + str(i)],
                ),
                axis=0,
            )
            sample_data_sorted = regression_data_input[
                np.lexsort(
                    (
                        regression_data_input[:, 2],
                        regression_data_input[:, 1],
                        regression_data_input[:, 0],
                    )
                )
            ]
            concat_01_sorted = concat_01[
                np.lexsort((concat_01[:, 2], concat_01[:, 1], concat_01[:, 0]))
            ]
            np.testing.assert_equal(sample_data_sorted, concat_01_sorted)
            concat_02 = np.concatenate(
                (
                    training_data["training_extras_" + str(i)],
                    cross_val_data["test_extras_" + str(i)],
                ),
                axis=0,
            )
            additional_data_sorted = additional_data_input[
                np.lexsort(
                    (
                        additional_data_input[:, 3],
                        additional_data_input[:, 2],
                        additional_data_input[:, 1],
                        additional_data_input[:, 0],
                    )
                )
            ]
            concat_02_sorted = concat_02[
                np.lexsort(
                    (concat_02[:, 3], concat_02[:, 2], concat_02[:, 1], concat_02[:, 0])
                )
            ]
            np.testing.assert_equal(additional_data_sorted, concat_02_sorted)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polygeneration_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        x_input_train_data = regression_data_input[:, :-1]
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=1
        )
        poly_degree = 1
        output_1 = data_feed.polygeneration(
            poly_degree, data_feed.multinomials, x_input_train_data
        )
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = 4  # New number of features should be = 2 * max_polynomial_order + 2 for two input features
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] * x_input_train_data[:, 1]
        np.testing.assert_equal(output_1, expected_output)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polygeneration_02(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        x_input_train_data = regression_data_input[:, :-1]
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=2
        )
        poly_degree = 2
        output_1 = data_feed.polygeneration(
            poly_degree, data_feed.multinomials, x_input_train_data
        )
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = 6  # New number of features should be = 2 * max_polynomial_order + 2 for two input features
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] ** 2
        expected_output[:, 4] = x_input_train_data[:, 1] ** 2
        expected_output[:, 5] = x_input_train_data[:, 0] * x_input_train_data[:, 1]
        np.testing.assert_equal(output_1, expected_output)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polygeneration_03(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        x_input_train_data = regression_data_input[:, :-1]
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=10
        )
        poly_degree = 10
        output_1 = data_feed.polygeneration(
            poly_degree, data_feed.multinomials, x_input_train_data
        )
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
        np.testing.assert_equal(output_1, expected_output)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polygeneration_04(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        x_input_train_data = regression_data_input[:, :-1]
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=10,
            multinomials=0,
        )
        poly_degree = 10
        output_1 = data_feed.polygeneration(
            poly_degree, data_feed.multinomials, x_input_train_data
        )
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
        np.testing.assert_equal(output_1, expected_output)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polygeneration_05(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        x_input_train_data = regression_data_input[:, :-1]
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=1
        )
        poly_degree = 1
        additional_term = np.sqrt(x_input_train_data)
        output_1 = data_feed.polygeneration(
            poly_degree, data_feed.multinomials, x_input_train_data, additional_term
        )
        expected_output_nr = x_input_train_data.shape[0]
        expected_output_nc = (
            6  # New number of features should be = 2 * max_polynomial_order + 4
        )
        expected_output = np.zeros((expected_output_nr, expected_output_nc))
        expected_output[:, 0] = 1
        expected_output[:, 1] = x_input_train_data[:, 0]
        expected_output[:, 2] = x_input_train_data[:, 1]
        expected_output[:, 3] = x_input_train_data[:, 0] * x_input_train_data[:, 1]
        expected_output[:, 4] = additional_term[:, 0]
        expected_output[:, 5] = additional_term[:, 1]
        np.testing.assert_equal(output_1, expected_output)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cost_function_01(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        output_1 = PolynomialRegression.cost_function(
            theta, x_vector, y, reg_parameter=0
        )
        assert output_1 == expected_value

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cost_function_02(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        expected_value = 90.625  # Calculated externally as sum(dy^2) / 2m
        output_1 = PolynomialRegression.cost_function(
            theta, x_vector, y, reg_parameter=0
        )
        assert output_1 == expected_value

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cost_function_03(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        expected_value = 0  # Value should return zero for exact solution
        output_1 = PolynomialRegression.cost_function(
            theta, x_vector, y, reg_parameter=0
        )
        assert output_1 == expected_value

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_gradient_function_01(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        )  # Calculated externally: see Excel sheet
        expected_value = expected_value.reshape(
            expected_value.shape[0],
        )
        output_1 = PolynomialRegression.gradient_function(
            theta, x_vector, y, reg_parameter=0
        )
        np.testing.assert_equal(output_1, expected_value)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_gradient_function_02(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        output_1 = PolynomialRegression.gradient_function(
            theta, x_vector, y, reg_parameter=0
        )
        np.testing.assert_equal(output_1, expected_value)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_gradient_function_03(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        output_1 = PolynomialRegression.gradient_function(
            theta, x_vector, y, reg_parameter=0
        )
        np.testing.assert_equal(output_1, expected_value)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_bfgs_parameter_optimization_01(self, array_type):
        original_data_input = array_type(self.test_data)
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
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
        data_feed = PolynomialRegression(
            original_data_input,
            input_array,
            maximum_polynomial_order=5,
            solution_method="bfgs",
        )
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == "bfgs"
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_bfgs_parameter_optimization_02(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=4,
            solution_method="bfgs",
        )
        output_1 = data_feed.bfgs_parameter_optimization(x_vector, y)
        assert data_feed.solution_method == "bfgs"
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    def test_mle_estimate_01(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
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
        output_1 = PolynomialRegression.MLE_estimate(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_mle_estimate_02(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
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
        output_1 = PolynomialRegression.MLE_estimate(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    def test_pyomo_optimization_01(self):
        x_vector = np.array([[i**2, i, 1] for i in range(10)])
        y = np.array([[i**2] for i in range(1, 11)])
        expected_value = np.array([[1.0], [2.0], [1.0]])
        output_1 = PolynomialRegression.pyomo_optimization(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_pyomo_optimization_02(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1]
        x_vector = np.zeros((x.shape[0], 6))
        x_vector[:, 0] = x[:, 0] ** 2
        x_vector[:, 1] = x[:, 1] ** 2
        x_vector[:, 2] = x[:, 0]
        x_vector[:, 3] = x[:, 1]
        x_vector[:, 4] = x[:, 1] * x[:, 0]
        x_vector[:, 5] = 1
        expected_value = np.array([[1.0], [1.0], [2.0], [2.0], [0.0], [2.0]])
        output_1 = PolynomialRegression.pyomo_optimization(x_vector, y)
        np.testing.assert_array_equal(expected_value, np.round(output_1, 4))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cross_validation_error_calculation_01(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1].reshape(regression_data_input.shape[0], 1)
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
        output_1 = PolynomialRegression.cross_validation_error_calculation(
            theta, x_vector, y
        )
        assert output_1 == expected_value

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cross_validation_error_calculation_02(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1].reshape(regression_data_input.shape[0], 1)
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
        expected_value = 2 * 90.625  # Calculated externally as sum(dy^2) / 2m
        output_1 = PolynomialRegression.cross_validation_error_calculation(
            theta, x_vector, y
        )
        assert output_1 == expected_value

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_cross_validation_error_calculation_02(self, array_type):
        regression_data_input = array_type(self.training_data)
        x = regression_data_input[:, :-1]
        y = regression_data_input[:, -1].reshape(regression_data_input.shape[0], 1)
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
        expected_value = 2 * 0  # Value should return zero for exact solution
        output_1 = PolynomialRegression.cross_validation_error_calculation(
            theta, x_vector, y
        )
        assert output_1 == expected_value

    def mock_optimization(self, x, y):
        return 10 * np.ones((x.shape[1], 1))

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    @patch.object(PolynomialRegression, "MLE_estimate", mock_optimization)
    def test_polyregression_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=5,
            solution_method="mle",
        )
        poly_order = 2
        training_data = regression_data_input[0:20, :]
        test_data = regression_data_input[20:, :]
        expected_output = 10 * np.ones((6, 1))
        output_1, _, _ = data_feed.polyregression(poly_order, training_data, test_data)
        np.testing.assert_array_equal(expected_output, output_1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    @patch.object(
        PolynomialRegression, "bfgs_parameter_optimization", mock_optimization
    )
    def test_polyregression_02(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=5,
            solution_method="bfgs",
        )
        poly_order = 2
        training_data = regression_data_input[0:20, :]
        test_data = regression_data_input[20:, :]
        expected_output = 10 * np.ones((6, 1))
        output_1, _, _ = data_feed.polyregression(poly_order, training_data, test_data)
        np.testing.assert_array_equal(expected_output, output_1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    @patch.object(PolynomialRegression, "pyomo_optimization", mock_optimization)
    def test_polyregression_03(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=5
        )
        poly_order = 2
        training_data = regression_data_input[0:20, :]
        test_data = regression_data_input[20:, :]
        expected_output = 10 * np.ones((6, 1))
        output_1, _, _ = data_feed.polyregression(poly_order, training_data, test_data)
        np.testing.assert_array_equal(expected_output, output_1)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polyregression_04(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=5
        )
        poly_order = 10
        training_data = regression_data_input[0:20, :]
        test_data = regression_data_input[20:, :]
        expected_output = np.Inf
        output_1, output_2, output_3 = data_feed.polyregression(
            poly_order, training_data, test_data
        )
        np.testing.assert_array_equal(expected_output, output_1)
        np.testing.assert_array_equal(expected_output, output_2)
        np.testing.assert_array_equal(expected_output, output_3)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_surrogate_performance_01(self, array_type):
        original_data_input = array_type(self.test_data)
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
        order_best = 2
        phi_best = np.array([[0.0], [0.0], [0.0]])
        expected_value_1 = 38.5
        expected_value_2 = 2533.3
        expected_value_3 = -1.410256
        expected_value_4 = 0
        data_feed = PolynomialRegression(
            original_data_input, input_array, maximum_polynomial_order=5
        )
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(
            phi_best, order_best
        )
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2
        assert np.round(output_3, 4) == np.round(expected_value_3, 4)
        assert np.round(output_4, 4) == np.round(expected_value_4, 4)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_surrogate_performance_02(self, array_type):
        original_data_input = array_type(self.test_data)
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
        order_best = 2
        phi_best = np.array([[1.0], [2.0], [1.0]])
        expected_value_1 = 0
        expected_value_2 = 0
        expected_value_3 = 1
        expected_value_4 = 1
        data_feed = PolynomialRegression(
            original_data_input, input_array, maximum_polynomial_order=5
        )
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(
            phi_best, order_best
        )
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2
        assert np.round(output_3, 4) == np.round(expected_value_3, 4)
        assert np.round(output_4, 4) == np.round(expected_value_4, 4)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array, pd.DataFrame])
    def test_surrogate_performance_03(self, array_type):
        original_data_input = array_type(self.test_data)
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
        order_best = 2
        phi_best = np.array([[1.0], [1.0], [1.0]])
        expected_value_1 = 4.5
        expected_value_2 = 28.5
        expected_value_3 = 0.972884259
        expected_value_4 = 0.9651369
        data_feed = PolynomialRegression(
            original_data_input, input_array, maximum_polynomial_order=5
        )
        _, output_1, output_2, output_3, output_4 = data_feed.surrogate_performance(
            phi_best, order_best
        )
        assert output_1 == expected_value_1
        assert output_2 == expected_value_2
        assert np.round(output_3, 4) == np.round(expected_value_3, 4)
        assert np.round(output_4, 4) == np.round(expected_value_4, 4)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_results_generation_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            multinomials=0,
        )
        order = 1
        beta = np.array([[0], [0], [0]])
        expected_df = pd.Series()
        row_list = np.array([["k"], ["(x_1)^1"], ["(x_2)^1"]])
        expected_df = expected_df.append(
            pd.Series(
                {
                    row_list[0, 0]: beta[0, 0],
                    row_list[1, 0]: beta[1, 0],
                    row_list[2, 0]: beta[2, 0],
                }
            )
        )
        output_df = data_feed.results_generation(beta, order)
        assert output_df.index.to_list() == expected_df.index.to_list()
        assert expected_df.all() == output_df.all()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_results_generation_02(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            multinomials=0,
        )
        order = 3
        beta = np.array([[1], [0.3], [6], [500], [500000], [0.001], [50]])
        expected_df = pd.Series()
        row_list = np.array(
            [
                ["k"],
                ["(x_1)^1"],
                ["(x_2)^1"],
                ["(x_1)^2"],
                ["(x_2)^2"],
                ["(x_1)^3"],
                ["(x_2)^3"],
            ]
        )
        expected_df = expected_df.append(
            pd.Series(
                {
                    row_list[0, 0]: beta[0, 0],
                    row_list[1, 0]: beta[1, 0],
                    row_list[2, 0]: beta[2, 0],
                    row_list[3, 0]: beta[3, 0],
                    row_list[4, 0]: beta[4, 0],
                    row_list[5, 0]: beta[5, 0],
                    row_list[6, 0]: beta[6, 0],
                }
            )
        )
        output_df = data_feed.results_generation(beta, order)
        assert output_df.index.to_list() == expected_df.index.to_list()
        assert expected_df.all() == output_df.all()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_results_generation_03(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            multinomials=1,
        )
        order = 2
        beta = np.array([[1], [0.3], [6], [500], [500000], [0.001]])
        expected_df = pd.Series()
        row_list = np.array(
            [["k"], ["(x_1)^1"], ["(x_2)^1"], ["(x_1)^2"], ["(x_2)^2"], ["(x_1).(x_2)"]]
        )
        expected_df = expected_df.append(
            pd.Series(
                {
                    row_list[0, 0]: beta[0, 0],
                    row_list[1, 0]: beta[1, 0],
                    row_list[2, 0]: beta[2, 0],
                    row_list[3, 0]: beta[3, 0],
                    row_list[4, 0]: beta[4, 0],
                    row_list[5, 0]: beta[5, 0],
                }
            )
        )
        output_df = data_feed.results_generation(beta, order)
        assert output_df.index.to_list() == expected_df.index.to_list()
        assert expected_df.all() == output_df.all()

    @pytest.mark.unit
    @patch("matplotlib.pyplot.show")
    def test_error_plotting(self, mock_show):
        mock_show.return_value = None
        # Generate typical data values for eaxch variable in order to test plot function
        plotting_data = np.array(
            [
                [1, 7, 23.1, 29.2, 0.01, 1300, -1.6, -1.1, 5],
                [2, 9, 0.0055, 0.0015, 0.006, 159, 0.34, 0.19, 10],
                [3, 6, 0.0007, 0.0009, 0.001, 2.3, 0.998, 0.994, 15],
                [4, 2, 0.0005, 0.0004, 0.0008, 0.02, 1.0, 1.0, 20],
            ]
        )
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
        assert (
            ax1.title.get_text() == "Training (green) vs Cross-validation error (red)"
        )
        assert ax2.title.get_text() == "MAE"
        assert ax3.title.get_text() == "MSE"
        assert ax4.title.get_text() == "R-squared (blue) and Adjusted R-squared (red)"
        # Compare the line/curve segments for each line plot to the actual input data
        np.testing.assert_array_equal(
            ax1.lines[0].get_path()._vertices, expected_output_1
        )
        np.testing.assert_array_equal(
            ax1.lines[1].get_path()._vertices, expected_output_2
        )
        np.testing.assert_array_equal(
            ax2.lines[0].get_path()._vertices, expected_output_3
        )
        np.testing.assert_array_equal(
            ax3.lines[0].get_path()._vertices, expected_output_4
        )
        np.testing.assert_array_equal(
            ax4.lines[0].get_path()._vertices, expected_output_5
        )
        np.testing.assert_array_equal(
            ax4.lines[1].get_path()._vertices, expected_output_6
        )

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_user_defined_terms_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=2,
            training_split=split,
            number_of_crossvalidations=num_cv,
        )
        additional_terms = [
            np.sin(regression_data_input[:, 0]),
            np.sin(regression_data_input[:, 1]),
        ]
        returned_features_array = data_feed.user_defined_terms(additional_terms)

        # Create two additional columns
        expected_array = np.zeros((regression_data_input.shape[0], 2))
        expected_array[:, 0] = np.sin(regression_data_input[:, 0])
        expected_array[:, 1] = np.sin(regression_data_input[:, 1])

        np.testing.assert_equal(expected_array, returned_features_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [pd.DataFrame])
    def test_user_defined_terms_02(self, array_type):
        original_data_input = array_type(self.full_data)
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(
            original_data_input,
            original_data_input,
            maximum_polynomial_order=2,
            training_split=split,
            number_of_crossvalidations=num_cv,
        )
        additional_terms = [
            np.sin(original_data_input["x1"]),
            np.sin(original_data_input["x2"]),
        ]
        returned_features_array = data_feed.user_defined_terms(additional_terms)

        # Create two additional columns
        expected_array = np.zeros((original_data_input.values.shape[0], 2))
        expected_array[:, 0] = np.sin(original_data_input.values[:, 0])
        expected_array[:, 1] = np.sin(original_data_input.values[:, 1])

        np.testing.assert_equal(expected_array, returned_features_array)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_user_defined_terms_03(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=2,
            training_split=split,
            number_of_crossvalidations=num_cv,
        )
        additional_terms = [
            np.sin(original_data_input["x1"]),
            np.sin(original_data_input["x2"]),
        ]
        with pytest.raises(Exception):
            data_feed.user_defined_terms(additional_terms)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_user_defined_terms_04(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=2,
            training_split=split,
            number_of_crossvalidations=num_cv,
        )
        p = np.sin(regression_data_input[:, 0]).reshape(
            regression_data_input.shape[0], 1
        )  # 2-D array, should raise error
        additional_terms = [p]
        with pytest.raises(Exception):
            data_feed.user_defined_terms(additional_terms)

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_user_defined_terms_05(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        num_cv = 3
        split = 0.75
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=2,
            training_split=split,
            number_of_crossvalidations=num_cv,
        )
        p = np.sin(regression_data_input[:, 0]).reshape(
            regression_data_input.shape[0], 1
        )
        additional_terms = p  # Additional terms as array, not list
        with pytest.raises(ValueError):
            data_feed.user_defined_terms(additional_terms)

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polynomial_regression_fitting_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=2
        )
        data_feed.get_feature_vector()
        results = data_feed.polynomial_regression_fitting()
        assert results.fit_status == "ok"

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    @patch("matplotlib.pyplot.show")
    def test_polynomial_regression_fitting_02(
        self, mock_show, array_type1, array_type2
    ):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        mock_show.return_value = None
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input[:5], maximum_polynomial_order=1
        )
        data_feed.get_feature_vector()
        with pytest.warns(Warning):
            results = data_feed.polynomial_regression_fitting()
            assert results.fit_status == "poor"

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polynomial_regression_fitting_03(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=2
        )
        data_feed.get_feature_vector()
        results = data_feed.polynomial_regression_fitting()
        x_input_train_data = regression_data_input[:, :-1]
        assert results.fit_status == "ok"

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_polynomial_regression_fitting_04(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=2
        )
        data_feed.get_feature_vector()
        additional_regression_features = [
            np.sin(regression_data_input[:, 0]),
            np.sin(regression_data_input[:, 1]),
        ]
        results = data_feed.polynomial_regression_fitting(
            additional_regression_features
        )
        results = data_feed.polynomial_regression_fitting()
        assert results.fit_status == "ok"

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_get_feature_vector_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            multinomials=0,
        )
        output = data_feed.get_feature_vector()
        expected_dict = {"x1": 0, "x2": 0}
        assert expected_dict == output.extract_values()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type", [np.array])
    def test_get_feature_vector_02(self, array_type):
        regression_data_input = array_type(self.training_data)
        data_feed = PolynomialRegression(
            regression_data_input,
            regression_data_input,
            maximum_polynomial_order=3,
            multinomials=0,
        )
        output = data_feed.get_feature_vector()
        expected_dict = {0: 0, 1: 0}
        assert expected_dict == output.extract_values()

    @pytest.mark.unit
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_set_additional_terms_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=3
        )
        data_feed.set_additional_terms("a")
        assert "a" == data_feed.additional_term_expressions
        data_feed.set_additional_terms(1)
        assert 1 == data_feed.additional_term_expressions
        data_feed.set_additional_terms([1, 2])
        assert [1, 2] == data_feed.additional_term_expressions
        data_feed.set_additional_terms(np.array([1, 2]))
        np.testing.assert_equal(np.array([1, 2]), data_feed.additional_term_expressions)

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_poly_training_01(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=2
        )
        data_feed.get_feature_vector()
        data_feed.training()
        assert data_feed.fit_status == "ok"

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [pd.DataFrame])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_generate_expression(self, array_type1, array_type2):
        original_data_input = array_type1(self.full_data)
        regression_data_input = array_type2(self.training_data)
        data_feed = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=2
        )

        p = data_feed.get_feature_vector()
        data_feed.training()

        lv = []
        for i in p.keys():
            lv.append(p[i])
        poly_expr = data_feed.generate_expression((lv))

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [np.array])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_pickle_load01(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        PolyClass = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=3
        )
        PolyClass.get_feature_vector()
        PolyClass.training()
        PolyClass.pickle_load(PolyClass.filename)

    @pytest.mark.unit
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [np.array])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_pickle_load02(self, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        PolyClass = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=3
        )
        with pytest.raises(Exception):
            PolyClass.pickle_load("file_not_existing.pickle")

    @pytest.mark.unit
    @patch("matplotlib.pyplot.show")
    @pytest.fixture(scope="module")
    @pytest.mark.parametrize("array_type1", [np.array])
    @pytest.mark.parametrize("array_type2", [np.array])
    def test_parity_residual_plots(self, mock_show, array_type1, array_type2):
        original_data_input = array_type1(self.test_data)
        regression_data_input = array_type2(self.sample_points)
        PolyClass = PolynomialRegression(
            original_data_input, regression_data_input, maximum_polynomial_order=3
        )
        PolyClass.get_feature_vector()
        PolyClass.training()

        PolyClass.parity_residual_plots()

    @pytest.fixture(scope="module")
    @pytest.mark.unit
    def test_confint_regression_01(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
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
        expected_stderror = np.array([4.541936, 0.850783])
        data_feed = PolynomialRegression(
            input_array, input_array, maximum_polynomial_order=1, overwrite=True
        )
        p = data_feed.get_feature_vector()
        res = data_feed.training()
        opt_wts = res.optimal_weights_array.reshape(expected_stderror.shape)
        output = data_feed.confint_regression(confidence=0.99)
        tval = 3.2498355440153697
        np.testing.assert_allclose(
            output["Std. error"].values, expected_stderror, atol=1e-3
        )
        np.testing.assert_allclose(
            output["Conf. int. lower"].values,
            opt_wts - tval * expected_stderror,
            atol=1e-3,
        )
        np.testing.assert_allclose(
            output["Conf. int. upper"].values,
            opt_wts + tval * expected_stderror,
            atol=1e-3,
        )

    @pytest.fixture(scope="module")
    @pytest.mark.unit
    def test_confint_regression_02(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
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
        expected_stderror = np.array([4.541936, 0.850783])
        data_feed = PolynomialRegression(
            input_array, input_array, maximum_polynomial_order=1, overwrite=True
        )
        p = data_feed.get_feature_vector()
        res = data_feed.training()
        opt_wts = res.optimal_weights_array.reshape(expected_stderror.shape)
        output = data_feed.confint_regression(confidence=0.9)
        tval = 1.8331129326536335
        np.testing.assert_allclose(
            output["Std. error"].values, expected_stderror, atol=1e-3
        )
        np.testing.assert_allclose(
            output["Conf. int. lower"].values,
            opt_wts - tval * expected_stderror,
            atol=1e-3,
        )
        np.testing.assert_allclose(
            output["Conf. int. upper"].values,
            opt_wts + tval * expected_stderror,
            atol=1e-3,
        )

    @pytest.fixture(scope="module")
    @pytest.mark.unit
    def test_confint_regression_03(self):
        # Create x vector for ax2 + bx + c: x data supplied in x_vector
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
        expected_stderror = np.array([0, 0, 0])
        data_feed = PolynomialRegression(
            input_array, input_array, maximum_polynomial_order=2, overwrite=True
        )
        p = data_feed.get_feature_vector()
        res = data_feed.training()
        opt_wts = res.optimal_weights_array.reshape(expected_stderror.shape)
        output = data_feed.confint_regression(confidence=0.99)
        tval = 3.2498355440153697
        np.testing.assert_allclose(
            output["Std. error"].values, expected_stderror, atol=1e-3
        )
        np.testing.assert_allclose(
            output["Conf. int. lower"].values,
            opt_wts - tval * expected_stderror,
            atol=1e-3,
        )
        np.testing.assert_allclose(
            output["Conf. int. upper"].values,
            opt_wts + tval * expected_stderror,
            atol=1e-3,
        )


if __name__ == "__main__":
    pytest.main()
