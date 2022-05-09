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
from abc import ABC, abstractmethod
import numpy as np


class AbstractPriceForecaster(ABC):
    @abstractmethod
    def forecast_day_ahead_and_real_time_prices(
        self, date, hour, bus, horizon, n_samples
    ):
        pass

    @abstractmethod
    def forecast_real_time_prices(self, date, hour, bus, horizon, n_samples):
        pass

    @abstractmethod
    def forecast_day_ahead_prices(self, date, hour, bus, horizon, n_samples):
        pass


class PlaceHolderForecaster(AbstractPriceForecaster):

    """
    This a placeholder for a real price forecaster. This placeholder can takes
    representative daily values and standard deviations for real-time and
    day-ahead prices to forecast prices in the future.
    """

    def __init__(
        self,
        daily_da_price_means: list,
        daily_rt_price_means: list,
        daily_da_price_stds: list,
        daily_rt_price_stds: list,
    ):

        self.daily_da_price_means = daily_da_price_means
        self.daily_rt_price_means = daily_rt_price_means
        self.daily_da_price_stds = daily_da_price_stds
        self.daily_rt_price_stds = daily_rt_price_stds

    def forecast_day_ahead_and_real_time_prices(
        self, date, hour, bus, horizon, n_samples
    ):
        rt_forecast = self.forecast_real_time_prices(
            date, hour, bus, horizon, n_samples
        )
        da_forecast = self.forecast_day_ahead_prices(
            date, hour, bus, horizon, n_samples
        )

        return da_forecast, rt_forecast

    def forecast_day_ahead_prices(self, date, hour, bus, horizon, n_samples):
        return self._forecast(
            means=self.daily_da_price_means,
            stds=self.daily_da_price_stds,
            hour=hour,
            horizon=horizon,
            n_samples=n_samples,
        )

    def forecast_real_time_prices(self, date, hour, bus, horizon, n_samples):
        return self._forecast(
            means=self.daily_rt_price_means,
            stds=self.daily_rt_price_stds,
            hour=hour,
            horizon=horizon,
            n_samples=n_samples,
        )

    def _forecast(self, means, stds, hour, horizon, n_samples):

        corresponding_means = [means[t % 24] for t in range(hour, hour + horizon)]
        corresponding_stds = [stds[t % 24] for t in range(hour, hour + horizon)]

        forecasts_arr = np.random.normal(
            loc=corresponding_means, scale=corresponding_stds, size=(n_samples, horizon)
        )
        forecasts_arr[forecasts_arr < 0] = 0

        return {i: list(forecasts_arr[i]) for i in range(n_samples)}
