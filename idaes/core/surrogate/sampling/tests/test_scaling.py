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
Tests for surrogates/sampling/scaling module
"""
import pytest
import pandas as pd
from idaes.core.surrogate.sampling.scaling import OffsetScaler


class TestScaling:
    @pytest.mark.unit
    def test_OffsetScaler(self):
        df = pd.DataFrame(
            {"A": [1.0, 2.0, 3.0, 4.0, 5.0], "B": [2.0, 5.0, 6.0, 10.0, 12.0]}
        )

        scaler = OffsetScaler.create_from_mean_std(df)
        pd.testing.assert_series_equal(scaler._offset, pd.Series({"A": 3.0, "B": 7.0}))
        pd.testing.assert_series_equal(
            scaler._factor, pd.Series({"A": 1.58113, "B": 4.0})
        )

        scaled_df = scaler.scale(df)
        expected_scaled_df = pd.DataFrame(
            {
                "A": [-1.264911064, -0.632455532, 0, 0.632455532, 1.264911064],
                "B": [-1.25, -0.5, -0.25, 0.75, 1.25],
            }
        )
        pd.testing.assert_frame_equal(scaled_df, expected_scaled_df)

        unscaled_df = scaler.unscale(scaled_df)
        pd.testing.assert_frame_equal(df, unscaled_df)

        # test that error is thrown if columns do not match
        with pytest.raises(ValueError):
            df2 = pd.DataFrame({"A": [1, 2, 3, 4, 5], "B2": [2, 5, 6, 10, 12]})
            scaled_df2 = scaler.scale(df2)

        with pytest.raises(ValueError) as excinfo:
            scaler = OffsetScaler(
                expected_columns=["A", "B"],
                offset_series=pd.Series({"A": 1, "C": 2}),
                factor_series=pd.Series({"A": 3, "B": 4}),
            )
        assert (
            str(excinfo.value)
            == "OffsetScaler was passed an offset series with an index that"
            " does not match expected_columns. Please make sure these labels match."
        )

        with pytest.raises(ValueError) as excinfo:
            scaler = OffsetScaler(
                expected_columns=["A", "B"],
                offset_series=pd.Series({"A": 1, "B": 2}),
                factor_series=pd.Series({"A": 3, "C": 4}),
            )
        assert (
            str(excinfo.value)
            == "OffsetScaler was passed a factor series with an index that"
            " does not match expected_columns. Please make sure these labels match."
        )

    @pytest.mark.unit
    def test_OffsetScaler_normalize(self):
        df = pd.DataFrame(
            {"A": [1.0, 2.0, 3.0, 4.0, 5.0], "B": [2.0, 5.0, 6.0, 10.0, 12.0]}
        )

        scaler = OffsetScaler.create_normalizing_scaler(df)
        pd.testing.assert_series_equal(scaler._offset, pd.Series({"A": 1.0, "B": 2.0}))
        pd.testing.assert_series_equal(scaler._factor, pd.Series({"A": 4.0, "B": 10.0}))

        scaled_df = scaler.scale(df)
        expected_scaled_df = pd.DataFrame(
            {"A": [0.0, 0.25, 0.5, 0.75, 1.0], "B": [0.0, 0.3, 0.4, 0.8, 1.0]}
        )
        pd.testing.assert_frame_equal(scaled_df, expected_scaled_df)

        df = pd.DataFrame({"A": [0.0, 6.0], "B": [3.0, 7.0]})
        scaled_df = scaler.scale(df)
        expected_scaled_df = pd.DataFrame({"A": [-0.25, 1.25], "B": [0.1, 0.5]})
        pd.testing.assert_frame_equal(scaled_df, expected_scaled_df)

        unscaled_df = scaler.unscale(scaled_df)
        pd.testing.assert_frame_equal(df, unscaled_df)

    @pytest.mark.unit
    def test_OffsetScaler_serialization(self):
        df = pd.DataFrame(
            {"A": [4.0, 10.0, 12.0, 20.0, 24.0], "B": [2.0, 5.0, 6.0, 10.0, 12.0]}
        )
        scaler = OffsetScaler.create_from_mean_std(df)

        d = scaler.to_dict()
        assert d["expected_columns"] == ["A", "B"]
        assert d["offset"]["A"] == pytest.approx(14)
        assert d["offset"]["B"] == pytest.approx(7)
        assert d["factor"]["A"] == pytest.approx(8)
        assert d["factor"]["B"] == pytest.approx(4)

        scaler = OffsetScaler.from_dict(d)
        assert scaler._expected_columns == ["A", "B"]
        pd.testing.assert_series_equal(scaler._offset, pd.Series({"A": 14.0, "B": 7.0}))
        pd.testing.assert_series_equal(scaler._factor, pd.Series({"A": 8.0, "B": 4.0}))


if __name__ == "__main__":
    pass
