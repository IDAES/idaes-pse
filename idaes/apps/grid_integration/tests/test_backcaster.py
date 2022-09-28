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

import pytest
from pyomo.common import unittest as pyo_unittest
from idaes.apps.grid_integration.forecaster import ForecastError, Backcaster
import warnings
pytestmark = pytest.mark.filterwarnings("error")

@pytest.fixture
def historical_da_prices():
    return {"test_bus": [1] * 24 + [2] * 24 + [3] * 24}


@pytest.fixture
def historical_rt_prices():
    return {"test_bus": [10] * 24 + [20] * 24 + [30] * 24}


@pytest.fixture
def base_backcaster(historical_da_prices, historical_rt_prices):
    return Backcaster(historical_da_prices, historical_rt_prices)


@pytest.mark.unit
def test_create_backcaster(historical_da_prices, historical_rt_prices):
    backcaster = Backcaster(historical_da_prices, historical_rt_prices)
    assert backcaster.historical_da_prices is historical_da_prices
    assert backcaster.historical_rt_prices is historical_rt_prices


@pytest.mark.unit
def test_create_backcaster_with_small_max_historical_days(
    historical_da_prices, historical_rt_prices
):
    max_n_days = 1.0
    with pytest.raises(Warning, match=r".*Dropping the first day's data.*"):
        backcaster = Backcaster(
            historical_da_prices, historical_rt_prices, max_historical_days=max_n_days
        )
        expected_historical_da_prices = {"test_bus": [3] * 24}
        expected_historical_rt_prices = {"test_bus": [30] * 24}

        pyo_unittest.assertStructuredAlmostEqual(
            first=backcaster.historical_da_prices, second=expected_historical_da_prices
        )

        pyo_unittest.assertStructuredAlmostEqual(
            first=backcaster.historical_rt_prices, second=expected_historical_rt_prices
        )


@pytest.mark.unit
@pytest.mark.parametrize("value", ["10", "ten", [10]])
def test_invalid_max_historical_days_type(
    value, historical_da_prices, historical_rt_prices
):
    with pytest.raises(TypeError, match=r".*max_historical_days must be a number.*"):
        Backcaster(
            historical_da_prices, historical_rt_prices, max_historical_days=value
        )


@pytest.mark.unit
@pytest.mark.parametrize("value", [-10, 0, -0.0])
def test_invalid_max_historical_days_value(
    value, historical_da_prices, historical_rt_prices
):
    with pytest.raises(ValueError, match=r".*max_historical_days must be >= 1.*"):
        Backcaster(
            historical_da_prices, historical_rt_prices, max_historical_days=value
        )


@pytest.mark.unit
def test_create_backcaster_with_non_dict(historical_da_prices, historical_rt_prices):
    with pytest.raises(
        TypeError, match=r"Given historical price is not an dictionary object.*"
    ):
        Backcaster(historical_da_prices, [])
        Backcaster(100, historical_rt_prices)


@pytest.mark.unit
def test_create_backcaster_with_empty_dict(historical_da_prices, historical_rt_prices):
    with pytest.raises(ValueError, match=r"Given historical price is empty."):
        Backcaster({}, historical_rt_prices)
        Backcaster(historical_da_prices, {})


@pytest.mark.unit
def test_create_backcaster_with_non_list(historical_da_prices, historical_rt_prices):
    with pytest.raises(TypeError, match=r".*bus test_bus is not a list object.*"):
        Backcaster({"test_bus": {1, 2, 3}}, historical_rt_prices)
        Backcaster(historical_da_prices, {"test_bus": {1, 2, 3}})


@pytest.mark.unit
def test_create_backcaster_with_not_enough_entries(
    historical_da_prices, historical_rt_prices
):
    with pytest.raises(
        ValueError, match=r".*At least a day of the historical prices.*"
    ):
        Backcaster({"test_bus": [1] * 23}, historical_rt_prices)
        Backcaster(historical_da_prices, {"test_bus": [1] * 2})


@pytest.mark.unit
def test_create_backcaster_with_non_24_entries(
    historical_da_prices, historical_rt_prices
):
    with pytest.raises(ValueError, match=r".*should be a multiple of 24.*"):
        Backcaster({"test_bus": [1] * 47}, historical_rt_prices)
        Backcaster(historical_da_prices, {"test_bus": [1] * 46})


# @pytest.mark.unit
# def test_create_backcaster_with_histroical_days_geq_max_historical_days(
#     historical_da_prices, historical_rt_prices
# ):
#     max_h_days = 1.0
#     with pytest.raises(TypeError, match=r".*Drooping the data for the first.*"):
#         Backcaster({"test_bus": [1] * 72}, historical_rt_prices, max_historical_days = max_h_days)
#         Backcaster(historical_da_prices, {"test_bus": [1] * 72}, max_historical_days = max_h_days)


@pytest.mark.unit
def test_forecast_real_time_prices(base_backcaster):

    n_samples = 3
    horizon = 4

    result_forecasts = base_backcaster.forecast_real_time_prices(
        date="2022-05-11", hour=18, bus="test_bus", horizon=horizon, n_samples=n_samples
    )
    expected_forecasts = {
        n_samples - i - 1: [(i + 1) * 10] * horizon for i in range(n_samples)
    }

    pyo_unittest.assertStructuredAlmostEqual(
        first=result_forecasts, second=expected_forecasts
    )


@pytest.mark.unit
def test_forecast_day_ahead_prices(base_backcaster):

    n_samples = 2
    horizon = 48

    result_forecasts = base_backcaster.forecast_day_ahead_prices(
        date="2022-05-11", hour=0, bus="test_bus", horizon=horizon, n_samples=n_samples
    )
    expected_forecasts = {0: [3] * 24 + [1] * 24, 1: [2] * 24 + [3] * 24}

    pyo_unittest.assertStructuredAlmostEqual(
        first=result_forecasts, second=expected_forecasts
    )


@pytest.mark.unit
def test_forecast_day_ahead_and_real_time_prices(base_backcaster):

    da_result_forecasts, rt_forecasts = base_backcaster.forecast_day_ahead_and_real_time_prices(
        date = "2022-05-11", hour = 0, bus = "test_bus", horizon = 48, n_samples = 2)
    expected_da_forecasts = {0: [3] * 24 + [1] * 24, 1: [2] * 24 + [3] * 24}
    
    pyo_unittest.assertStructuredAlmostEqual(
        first=da_result_forecasts, second=expected_da_forecasts
    )

    da_forecasts, rt_result_forecasts = base_backcaster.forecast_day_ahead_and_real_time_prices(
        date = "2022-05-11", hour = 18, bus = "test_bus", horizon = 4, n_samples = 3)

    expected_rt_forecasts = {
        3 - i - 1: [(i + 1) * 10] * 4 for i in range(3)
    }
    pyo_unittest.assertStructuredAlmostEqual(
        first=rt_result_forecasts, second=expected_rt_forecasts
    )


@pytest.mark.unit
def test_forecast_nonexistent_bus_prices(base_backcaster):

    wrong_bus = "test_bussss"

    n_samples = 3
    horizon = 4

    with pytest.raises(
        ForecastError, match=r"No test_bussss real-time price available"
    ):
        base_backcaster.forecast_real_time_prices(
            date="2022-05-11",
            hour=18,
            bus=wrong_bus,
            horizon=horizon,
            n_samples=n_samples,
        )

    n_samples = 2
    horizon = 48

    with pytest.raises(
        ForecastError, match=r"No test_bussss day-ahead price available"
    ):
        result_forecasts = base_backcaster.forecast_day_ahead_prices(
            date="2022-05-11",
            hour=0,
            bus=wrong_bus,
            horizon=horizon,
            n_samples=n_samples,
        )


class MockPrescientHourlyStats:
    def __init__(self, bus_to_prices_dict) -> None:
        self.observed_bus_LMPs = bus_to_prices_dict


@pytest.mark.unit
def test_fetch_hourly_stats_from_prescient(base_backcaster, historical_rt_prices):

    prescient_hourly_stats = MockPrescientHourlyStats({"test_bus":15})
    base_backcaster.fetch_hourly_stats_from_prescient(prescient_hourly_stats)

    expected_current_day_rt_prices = {}
    expected_current_day_rt_prices["test_bus"] = [15]
    
    pyo_unittest.assertStructuredAlmostEqual(
        first=expected_current_day_rt_prices["test_bus"], second=base_backcaster._current_day_rt_prices["test_bus"]
    )

    for i in range(23):
        prescient_hourly_stats.observed_bus_LMPs["test_bus"] = 15
        base_backcaster.fetch_hourly_stats_from_prescient(prescient_hourly_stats)

    expected_current_day_rt_prices_1 = {}
    expected_current_day_rt_prices_1["test_bus"] = []

    expected_historical_rt_prices = [10] * 24 + [20] * 24 + [30] * 24 + [15]*24

    pyo_unittest.assertStructuredAlmostEqual(
        first=expected_current_day_rt_prices_1["test_bus"], second=base_backcaster._current_day_rt_prices["test_bus"]
    )
    pyo_unittest.assertStructuredAlmostEqual(
        first=expected_historical_rt_prices, second=base_backcaster._historical_rt_prices["test_bus"])


class MockRucMarket:
    def __init__(self, da_prices) -> None:
        self.day_ahead_prices = da_prices
    def get(self, bus, t):
        return self.day_ahead_prices[bus][t]

class MockPrescientRucPlan:
    def __init__(self, da_prices) -> None:
        self.ruc_market = MockRucMarket(da_prices)


@pytest.mark.unit
def fetch_day_ahead_stats_from_prescient(base_backcaster, historical_da_prices):

    for i in range(base_backcaster.max_historical_days + 1):
        da_price = {"test_bus": [i] * 24}
        prescient_ruc_plan = MockPrescientRucPlan(da_price)
        base_backcaster.fetch_day_ahead_stats_from_prescient(None, None, prescient_ruc_plan)

    expected_historical_da_prices = {}
    expected_historical_da_prices["test_bus"] = ([j]*24 for j in range(1,base_backcaster.max_historical_days+1))

    pyo_unittest.assertStructuredAlmostEqual(
        first=expected_historical_da_prices, second=base_backcaster._historical_da_prices["test_bus"])


