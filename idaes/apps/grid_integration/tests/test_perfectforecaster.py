#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pytest
from pyomo.common import unittest as pyo_unittest
from idaes.apps.grid_integration.forecaster import PerfectForecaster
import idaes.logger as idaeslog
import pandas as pd
import numpy as np

@pytest.fixture
def wind_df():
    start_year = 2020
    start_mon = 1
    start_day = 1
    start_date = pd.Timestamp(f"{start_year}-{start_mon:02d}-{start_day:02d} 00:00:00")
    ix = pd.date_range(start=start_date, 
                        end=start_date
                        + pd.offsets.DateOffset(days=1)
                        - pd.offsets.DateOffset(hours=1),
                        freq='1H')
    df = pd.DataFrame(index = ix)
    df["303_WIND_1-DACF"] = list(i/100 for i in range(24))
    df["303_WIND_1-RTCF"] = list((i+1)/100 for i in range(24))
    df["Caesar-DALMP"] = list(i for i in range(24))
    df["Caesar-RTLMP"] = list(i+1 for i in range(24))
    return df

@pytest.fixture
def base_perfectforecaster(wind_df):
    return PerfectForecaster(wind_df)

@pytest.mark.unit
def test_create_perfectforecaster(wind_df):
    perfectforecaster = PerfectForecaster(data_path_or_df=wind_df)
    assert perfectforecaster.data is wind_df


@pytest.mark.unit
@pytest.mark.parametrize("value", [np.array([1,2,3,4,5]), [1,2,3,4,5]])
def test_create_perfectforecaster_with_ndarray_and_list(value):
    with pytest.raises(ValueError):
        perfectforecaster = PerfectForecaster(value)

@pytest.mark.unit
def test_get_column_from_data(base_perfectforecaster):
    date = "2020-01-01"
    hour = 0
    horizon = 24
    col = '303_WIND_1-DACF'
    expected_forecast = [i/100 for i in range(24)]
    result_forecast = base_perfectforecaster.get_column_from_data(date, hour, horizon, col)

    pyo_unittest.assertStructuredAlmostEqual(
        first=result_forecast.tolist(), second=expected_forecast
    )
    return

@pytest.mark.unit
def test_forecast_day_ahead_prices(base_perfectforecaster):
    date = "2020-01-01"
    hour = 0
    horizon = 24
    bus = 'Caesar'
    nsp = 0     # this nsp is not used in the self.forecast_day_ahead_prices(), but we have it to make this func consistent with the stochastic bidder and coordinator.
    result_forecast = base_perfectforecaster.forecast_day_ahead_prices(date, hour, bus, horizon, nsp)
    expected_forecast = [i for i in range(24)]
    pyo_unittest.assertStructuredAlmostEqual(
        first=result_forecast.tolist(), second=expected_forecast
    )
    return

@pytest.mark.unit
def test_forecast_real_time_prices(base_perfectforecaster):
    date = "2020-01-01"
    hour = 0
    horizon = 4
    bus = 'Caesar'
    nsp = 0
    result_forecast = base_perfectforecaster.forecast_real_time_prices(date, hour, bus, horizon, nsp)
    expected_forecast = [i+1 for i in range(4)]
    pyo_unittest.assertStructuredAlmostEqual(
        first=result_forecast.tolist(), second=expected_forecast
    )
    return

@pytest.mark.unit
def test_forecast_day_ahead_capacity_factor(base_perfectforecaster):
    date = "2020-01-01"
    hour = 0
    horizon = 24
    gen = '303_WIND_1'
    result_forecast = base_perfectforecaster.forecast_day_ahead_capacity_factor(date, hour, gen, horizon)
    expected_forecast = [i/100 for i in range(24)]
    pyo_unittest.assertStructuredAlmostEqual(
        first=result_forecast.tolist(), second=expected_forecast
    )
    return

@pytest.mark.unit
def test_forecast_real_time_capacity_factor(base_perfectforecaster):
    date = "2020-01-01"
    hour = 0
    horizon = 4
    gen = '303_WIND_1'
    result_forecast = base_perfectforecaster.forecast_real_time_capacity_factor(date, hour, gen, horizon)
    expected_forecast = [(i+1)/100 for i in range(4)]
    pyo_unittest.assertStructuredAlmostEqual(
        first=result_forecast.tolist(), second=expected_forecast
    )
    return