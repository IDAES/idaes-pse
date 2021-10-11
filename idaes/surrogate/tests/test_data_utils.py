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
"""
Tests for data_utils module
"""
import pytest
import pandas as pd
from idaes.surrogate.data_utils import split_training_testing, split_training_testing_validation


class TestDataUtils:
    @pytest.mark.unit
    def test_split_training_testing(self):
        d = {'a':[1,2,3,4,5], 'b':[2,3,4,5,6]}
        df = pd.DataFrame(d)

        # test no seed
        df_training, df_testing = split_training_testing(
            dataframe=df,
            training_fraction=0.7)
        assert len(df_training) == 3
        assert len(df_testing) == 2

        # test seed
        expected_training_df = pd.DataFrame(
            {'a': [2, 5, 3], 'b': [3, 6, 4]})
        expected_testing_df = pd.DataFrame(
            {'a': [1, 4], 'b': [2, 5]})

        df_training, df_testing = split_training_testing(
            dataframe=df,
            training_fraction=0.7,
            seed=42)
        assert len(df_training) == 3
        assert len(df_testing) == 2
        pd.testing.assert_frame_equal(df_training.reset_index(drop=True), expected_training_df)
        pd.testing.assert_frame_equal(df_testing.reset_index(drop=True), expected_testing_df)

    @pytest.mark.unit
    def test_split_training_testing_validation(self):
        d = {'a':[1,2,3,4,5,6,7,8,9], 'b':[2,3,4,5,6,7,8,9,10]}
        df = pd.DataFrame(d)

        # test no seed
        df_training, df_testing, df_validation = \
            split_training_testing_validation(dataframe=df,
                                              training_fraction=0.6,
                                              testing_fraction=0.35)
        assert len(df_training) == 5
        assert len(df_testing) == 3
        assert len(df_validation) == 1

        # test seed
        expected_training_df = pd.DataFrame(
            {'a': [8, 2, 6, 1, 9], 'b': [9, 3, 7, 2, 10]})
        expected_testing_df = pd.DataFrame(
            {'a': [3, 5, 4], 'b': [4, 6, 5]})
        expected_validation_df = pd.DataFrame(
            {'a': [7], 'b': [8]})

        df_training, df_testing, df_validation = \
            split_training_testing_validation(dataframe=df,
                                              training_fraction=0.6,
                                              testing_fraction=0.35,
                                              seed=42)
        
        assert len(df_training) == 5
        assert len(df_testing) == 3
        assert len(df_validation) == 1
        pd.testing.assert_frame_equal(df_training.reset_index(drop=True), expected_training_df)
        pd.testing.assert_frame_equal(df_testing.reset_index(drop=True), expected_testing_df)
        pd.testing.assert_frame_equal(df_validation.reset_index(drop=True), expected_validation_df)

if __name__ == '__main__':
    TestDataUtils().test_split_training_testing()
    TestDataUtils().test_split_training_testing_validation()
    
