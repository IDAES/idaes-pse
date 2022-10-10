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
from numbers import Real
import numpy as np
import idaes.logger as idaeslog


_logger = idaeslog.getLogger(__name__)


class ForecastError(Exception):
    """Error to indicate error with forecasters."""


class AbstractPriceForecaster(ABC):

    """
    The abstract class for price forecaster.
    """

    @abstractmethod
    def forecast_day_ahead_and_real_time_prices(
        self, date, hour, bus, horizon, n_samples
    ):
        """
        Forecast both day-ahead and real-time market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: day-ahead price forecasts

            dict: real-time price forecasts

        """

        pass

    @abstractmethod
    def forecast_real_time_prices(self, date, hour, bus, horizon, n_samples):

        """
        Forecast real-time market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: real-time price forecasts

        """
        pass

    @abstractmethod
    def forecast_day_ahead_prices(self, date, hour, bus, horizon, n_samples):

        """
        Forecast day-ahead market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: day-ahead price forecasts

        """

        pass


class AbstractPrescientPriceForecaster(AbstractPriceForecaster):

    """
    Abstract class for price forecasters that will interface with Prescient.
    """

    @abstractmethod
    def fetch_hourly_stats_from_prescient(self, prescient_hourly_stats):

        """
        This method fetches the hourly stats from Prescient to the price forecaster
        once the hourly stats are published.

        Arguments:
            prescient_hourly_stats: Prescient HourlyStats object.

        Returns:
            None
        """

        pass

    @abstractmethod
    def fetch_day_ahead_stats_from_prescient(self, uc_date, uc_hour, day_ahead_result):

        """
        This method fetches the day-ahead market to the price forecaster after the
        UC is solved from Prescient through the coordinator.

        Arguments:
            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

            day_ahead_result: a Prescient RucPlan object.

        Returns:
            None
        """

        pass


class PlaceHolderForecaster(AbstractPrescientPriceForecaster):

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
        """
        Initialize the PlaceHolderForecaster.

        Arguments:
            daily_da_price_means: list of day-ahead price means

            daily_rt_price_means: list of real-time price means

            daily_da_price_stds: list of price standard deviations

            daily_rt_price_stds: list of real-time price standard deviations

        Returns:
            None
        """

        self.daily_da_price_means = daily_da_price_means
        self.daily_rt_price_means = daily_rt_price_means
        self.daily_da_price_stds = daily_da_price_stds
        self.daily_rt_price_stds = daily_rt_price_stds

    def forecast_day_ahead_and_real_time_prices(
        self, date, hour, bus, horizon, n_samples
    ):
        """
        Forecast both day-ahead and real-time market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: day-ahead price forecasts

            dict: real-time price forecasts

        """

        rt_forecast = self.forecast_real_time_prices(
            date, hour, bus, horizon, n_samples
        )
        da_forecast = self.forecast_day_ahead_prices(
            date, hour, bus, horizon, n_samples
        )

        return da_forecast, rt_forecast

    def forecast_day_ahead_prices(self, date, hour, bus, horizon, n_samples):

        """
        Forecast day-ahead market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: day-ahead price forecasts

        """

        return self._forecast(
            means=self.daily_da_price_means,
            stds=self.daily_da_price_stds,
            hour=hour,
            horizon=horizon,
            n_samples=n_samples,
        )

    def forecast_real_time_prices(self, date, hour, bus, horizon, n_samples):

        """
        Forecast real-time market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: real-time price forecasts

        """

        return self._forecast(
            means=self.daily_rt_price_means,
            stds=self.daily_rt_price_stds,
            hour=hour,
            horizon=horizon,
            n_samples=n_samples,
        )

    def _forecast(self, means, stds, hour, horizon, n_samples):

        """
        Generate price forecasts.

        Arguments:
            means: list of price means

            stds: list of price standard deviations

            hour: intended hour of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: real-time price forecasts
        """

        corresponding_means = [means[t % 24] for t in range(hour, hour + horizon)]
        corresponding_stds = [stds[t % 24] for t in range(hour, hour + horizon)]

        forecasts_arr = np.random.normal(
            loc=corresponding_means, scale=corresponding_stds, size=(n_samples, horizon)
        )
        forecasts_arr[forecasts_arr < 0] = 0

        return {i: list(forecasts_arr[i]) for i in range(n_samples)}

    def fetch_hourly_stats_from_prescient(self, prescient_hourly_stats):

        """
        This method fetches the hourly stats from Prescient to the price forecaster
        once the hourly stats are published.

        Arguments:
            prescient_hourly_stats: Prescient HourlyStats object.

        Returns:
            None
        """

        return

    def fetch_day_ahead_stats_from_prescient(self, uc_date, uc_hour, day_ahead_result):

        """
        This method fetches the day-ahead market to the price forecaster after the
        UC is solved from Prescient through the coordinator.

        Arguments:
            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

            day_ahead_result: a Prescient RucPlan object.

        Returns:
            None
        """

        return


class Backcaster(AbstractPrescientPriceForecaster):

    """
    Generate price forecasts by directly using historical prices.
    """

    def __init__(
        self, historical_da_prices, historical_rt_prices, max_historical_days=10
    ):
        """
        Initialize the Backcaster.

        Arguments:
            historical_da_prices: dictionary of list for historical hourly day-ahead prices

            historical_rt_prices: dictionary of list for historical hourly real-time prices

            max_historical_days: maximum number of days of price data to store on the instance

        Returns:
            None
        """

        self.max_historical_days = max_historical_days
        self.historical_da_prices = historical_da_prices
        self.historical_rt_prices = historical_rt_prices
        self._current_day_rt_prices = {bus: [] for bus in historical_da_prices}

    def _validate_input_historical_price(self, historical_price):

        """
        Validate input historical prices.

        Arguments:
            historical_price: dictionary of list for historical hourly prices

        Returns:
            None
        """

        if not isinstance(historical_price, dict):
            raise TypeError(
                "Given historical price is not an dictionary object. Dictionaries with bus name (str) as keys are expected."
            )

        if len(historical_price) == 0:
            raise ValueError(f"Given historical price is empty.")

        for b, price_list in historical_price.items():
            if not isinstance(price_list, list):
                raise TypeError(
                    f"Given historical price for bus {b} is not a list object. A list of historical prices is expected."
                )

            n_prices = len(price_list)

            if n_prices < 24:
                raise ValueError(
                    f"At least a day of the historical prices (24 entries) is needed. For bus {b}, only {n_prices} are provided."
                )

            if n_prices % 24 != 0:
                raise ValueError(
                    f"The number of prices for each bus should be a multiple of 24. But for bus {b}, {n_prices} are provided."
                )

            num_days = len(historical_price[b]) // 24
            if num_days > self.max_historical_days:
                _logger.warning(
                    f"The number of days in the input historical prices for bus {b} is greater than the max value {self.max_historical_days}."
                    f" Dropping the data for the first {num_days - self.max_historical_days} day(s)."
                )

                historical_price[b] = historical_price[b][
                    ((num_days - self.max_historical_days) * 24) :
                ]

    @property
    def max_historical_days(self):

        """
        Property getter for max_historical_days.

        Returns:
            int: max historical days
        """

        return self._max_historical_days

    @max_historical_days.setter
    def max_historical_days(self, value):

        """
        Property setter for max_historical_days (validate before setting).

        Args:
            value: intended value for max_historical_days

        Returns:
            None
        """

        if not isinstance(value, Real):
            raise TypeError(
                f"max_historical_days must be a number, but {type(value)} is provided."
            )

        value = int(value)

        if value < 1:
            raise ValueError(
                f"max_historical_days must be >= 1, but {value} is provided."
            )

        self._max_historical_days = value

    @property
    def historical_da_prices(self):

        """
        Property getter for historical_da_prices.

        Returns:
            dict: saved historical day-ahead prices
        """

        return self._historical_da_prices

    @historical_da_prices.setter
    def historical_da_prices(self, value):

        """
        Property setter for historical_da_prices (validate before setting).

        Args:
            value: intended value for historical_da_prices

        Returns:
            None
        """

        self._validate_input_historical_price(value)
        self._historical_da_prices = value

    @property
    def historical_rt_prices(self):

        """
        Property getter for historical_rt_prices.

        Returns:
            dict: saved historical real-time prices
        """

        return self._historical_rt_prices

    @historical_rt_prices.setter
    def historical_rt_prices(self, value):

        """
        Property setter for historical_rt_prices (validate before setting).

        Args:
            value: intended value for historical_rt_prices

        Returns:
            None
        """

        self._validate_input_historical_price(value)
        self._historical_rt_prices = value

    def forecast_day_ahead_and_real_time_prices(
        self, date, hour, bus, horizon, n_samples
    ):
        """
        Forecast both day-ahead and real-time market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: day-ahead price forecasts

            dict: real-time price forecasts

        """

        rt_forecast = self.forecast_real_time_prices(
            date, hour, bus, horizon, n_samples
        )
        da_forecast = self.forecast_day_ahead_prices(
            date, hour, bus, horizon, n_samples
        )

        return da_forecast, rt_forecast

    def forecast_real_time_prices(self, date, hour, bus, horizon, n_samples):

        """
        Forecast real-time market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: real-time price forecasts

        """

        return self._forecast(
            historical_price_dict=self.historical_rt_prices,
            market="real-time",
            date=date,
            hour=hour,
            bus=bus,
            horizon=horizon,
            n_samples=n_samples,
        )

    def forecast_day_ahead_prices(self, date, hour, bus, horizon, n_samples):

        """
        Forecast day-ahead market prices.

        Arguments:
            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: day-ahead price forecasts

        """

        return self._forecast(
            historical_price_dict=self.historical_da_prices,
            market="day-ahead",
            date=date,
            hour=0,
            bus=bus,
            horizon=horizon,
            n_samples=n_samples,
        )

    def _forecast(
        self, historical_price_dict, market, date, hour, bus, horizon, n_samples
    ):
        """
        Forecast energy market prices using historical prices.

        Arguments:

            historical_price_dict: the dictionary that holds the intended historical prices

            market: the market that the price forecast is for, e.g., day-ahead

            date: intended date of the forecasts

            hour: intended hour of the forecasts

            bus: intended bus of the forecasts

            horizon: number of the time periods of the forecasts

            n_samples: number of the samples

        Returns:
            dict: price forecasts

        """

        if bus not in historical_price_dict:
            raise ForecastError(f"No {bus} {market} price available.")

        historical_price_len = len(historical_price_dict[bus])
        n_days = historical_price_len // 24

        forecast = {}

        for i in range(n_samples):
            day_idx = n_days - (i % n_days) - 1

            forecast[i] = [
                historical_price_dict[bus][t % historical_price_len]
                for t in range(day_idx * 24 + hour, day_idx * 24 + hour + horizon)
            ]

        return forecast

    def fetch_hourly_stats_from_prescient(self, prescient_hourly_stats):

        """
        This method fetches the hourly real-time prices from Prescient and store
        them on the price forecaster, once they are published. When the stored historical
        data size has exceeded the specified upper bound, drop the oldest data.

        Arguments:
            prescient_hourly_stats: Prescient HourlyStats object.

        Returns:
            None
        """

        # save the newest rt prices
        for b, price_list in self._current_day_rt_prices.items():
            price_list.append(prescient_hourly_stats.observed_bus_LMPs[b])

        # update the historical
        for b in self._current_day_rt_prices:

            # if a full day's data is ready, get them ready for future forecasts
            if len(self._current_day_rt_prices[b]) >= 24:
                self._historical_rt_prices[b] += self._current_day_rt_prices[b]
                self._current_day_rt_prices[b] = []

            # drop oldes historical prices if total stored data exceeded the upper bound
            while len(self._historical_rt_prices[b]) // 24 > self.max_historical_days:
                self._historical_rt_prices[b] = self._historical_rt_prices[b][24:]

        return

    def fetch_day_ahead_stats_from_prescient(self, uc_date, uc_hour, day_ahead_result):

        """
        This method fetches the hourly day-ahead prices from Prescient and store
        them on the price forecaster, once they are published. When the stored historical
        data size has exceeded the specified upper bound, drop the oldest data.

        Arguments:
            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

            day_ahead_result: a Prescient RucPlan object.

        Returns:
            None
        """

        for b in self._historical_da_prices:

            # save the newest da prices
            self._historical_da_prices[b] += [
                day_ahead_result.ruc_market.day_ahead_prices.get((b, t))
                for t in range(24)
            ]

            # drop oldes historical prices if total stored data exceeded the upper bound
            while len(self._historical_da_prices[b]) // 24 > self.max_historical_days:
                self._historical_da_prices[b] = self._historical_da_prices[b][24:]

        return
