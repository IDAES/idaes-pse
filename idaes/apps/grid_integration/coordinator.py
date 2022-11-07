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
from itertools import zip_longest

from pyomo.common.dependencies import attempt_import
from pyomo.common.config import ConfigDict, ConfigValue
import pyomo.environ as pyo
from idaes.apps.grid_integration.utils import convert_marginal_costs_to_actual_costs

prescient, prescient_avail = attempt_import("prescient")


class DoubleLoopCoordinator:

    """
    Coordinate Prescient, tracker and bidder.
    """

    def __init__(self, bidder, tracker, projection_tracker):

        """
        Initializes the DoubleLoopCoordinator object and registers functionalities
        in Prescient's plugin system.

        Arguments:
            bidder: an initialized bidder object

            tracker: an initialized tracker object

            projection_tracker: an initialized tracker object, this object is
                                mimicking the behaviror of the projection SCED in
                                Prescient and to projecting the system states
                                and updating bidder model.

        Returns:
            None
        """

        self.bidder = bidder
        self.tracker = tracker
        self.projection_tracker = projection_tracker

    def register_plugins(self, context, options, plugin_config):

        """
        Register functionalities in Prescient's plugin system.

        Arguments:
            context: Prescient plugin PluginRegistrationContext from prescient.plugins.plugin_registration

            options: Prescient options from prescient.simulator.config

            plugin_config: Prescient plugin config

        Returns:
            None
        """

        self.plugin_config = plugin_config

        context.register_initialization_callback(self.initialize_customized_results)
        context.register_for_hourly_stats(self.push_hourly_stats_to_forecaster)
        context.register_before_ruc_solve_callback(self.pass_static_params_to_DA)
        context.register_before_ruc_solve_callback(self.bid_into_DAM)
        context.register_after_ruc_generation_callback(self.fetch_DA_prices)
        context.register_after_ruc_generation_callback(self.fetch_DA_dispatches)
        context.register_after_ruc_generation_callback(
            self.push_day_ahead_stats_to_forecaster
        )
        context.register_before_operations_solve_callback(self.pass_static_params_to_RT)
        context.register_before_operations_solve_callback(self.bid_into_RTM)
        context.register_after_operations_callback(self.track_sced_signal)
        context.register_update_operations_stats_callback(self.update_observed_dispatch)
        context.register_after_ruc_activation_callback(self.activate_pending_DA_data)
        context.register_finalization_callback(self.write_plugin_results)

        return

    def get_configuration(self, key):

        """
        Register customized commandline options.

        Arguments:
            key: plugin name

        Returns:
            config: Prescient config dict
        """

        config = ConfigDict()

        # Add command line options
        config.declare(
            "bidding_generator",
            ConfigValue(
                domain=str,
                description="Specifies the generator we derive bidding strategis for.",
                default=None,
            ),
        ).declare_as_argument("--bidding-generator")

        return config

    def initialize_customized_results(self, options, simulator):

        """
        This method is outdated.
        """

        simulator.data_manager.extensions["customized_results"] = {}
        customized_results = simulator.data_manager.extensions["customized_results"]

        customized_results["Generator"] = []
        customized_results["Date"] = []
        customized_results["Hour"] = []
        customized_results["State"] = []
        customized_results["RUC Schedule"] = []
        customized_results["SCED Schedule"] = []
        customized_results["Power Output"] = []

        return

    def push_hourly_stats_to_forecaster(self, prescient_hourly_stats):

        """
        This method pushes the hourly stats from Prescient to the price forecaster
        once the hourly stats are published.

        Arguments:
            prescient_hourly_stats: Prescient HourlyStats object.

        Returns:
            None
        """

        self.bidder.forecaster.fetch_hourly_stats_from_prescient(prescient_hourly_stats)

    def push_day_ahead_stats_to_forecaster(
        self, options, simulator, day_ahead_result, uc_date, uc_hour
    ):
        """
        This method pushes the day-ahead market to the price forecaster after the
        UC is solved.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            day_ahead_result: a Prescient RucPlan object.

            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

        Returns:
            None
        """

        self.bidder.forecaster.fetch_day_ahead_stats_from_prescient(
            uc_date, uc_hour, day_ahead_result
        )

    def _update_bids(self, gen_dict, bids, start_hour, horizon):

        """
        This method takes bids and pass the information in the bids to generator
        dictionary from Prescient.

        Arguments:
            gen_dict: a dictionary from Prescient's RUC/SCED instance that stores
            the generator parameters.

            bids: the bids we want to pass into the market.

            start_hour: the effective start hour of the bid

            horizon: the length of the bids

        Returns:
            None
        """

        gen_name = self.bidder.bidding_model_object.model_data.gen_name

        def _update_p_cost(gen_dict, bids, param_name, start_hour, horizon):

            # update the "p_cost" element in the generator's dict
            gen_dict["p_cost"] = {
                "data_type": "time_series",
                "values": [
                    {
                        "data_type": "cost_curve",
                        "cost_curve_type": "piecewise",
                        "values": bids[t][gen_name]["p_cost"],
                    }
                    for t in range(start_hour, horizon + start_hour)
                ],
            }

            # because the p_cost is updated, so delete p_fuel
            if "p_fuel" in gen_dict:
                gen_dict.pop("p_fuel")

            return

        def _update_time_series_params(gen_dict, bids, param_name, start_hour, horizon):

            value_list = [
                bids[t][gen_name].get(param_name, None)
                for t in range(start_hour, start_hour + horizon)
            ]
            if param_name in gen_dict:
                gen_dict[param_name] = {
                    "data_type": "time_series",
                    "values": value_list,
                }

            return

        def _update_non_time_series_params(
            gen_dict, bids, param_name, start_hour, horizon
        ):
            if param_name in gen_dict:
                gen_dict[param_name] = bids[start_hour][gen_name].get(param_name, None)

            return

        param_update_func_map = {
            "p_cost": _update_p_cost,
            "p_max": _update_time_series_params,
            "p_min": _update_time_series_params,
            "p_min_agc": _update_time_series_params,
            "p_max_agc": _update_time_series_params,
            "fixed_commitment": _update_time_series_params,
            "min_up_time": _update_non_time_series_params,
            "min_down_time": _update_non_time_series_params,
            "startup_capacity": _update_time_series_params,
            "shutdown_capacity": _update_time_series_params,
            "startup_fuel": _update_non_time_series_params,
            "startup_cost": _update_non_time_series_params,
        }

        for param_name in bids[start_hour][gen_name].keys():
            update_func = param_update_func_map[param_name]
            update_func(gen_dict, bids, param_name, start_hour, horizon)

        return

    def _pass_DA_bid_to_prescient(self, options, ruc_instance, bids):

        """
        This method passes the bids into the RUC model for day-ahead market clearing.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            ruc_instance: Prescient RUC object

            bids: the bids we want to pass into the day-ahead market.

        Returns:
            None
        """

        gen_name = self.bidder.bidding_model_object.model_data.gen_name

        # fetch the generator's parameter dictionary from Prescient UC instance
        gen_dict = ruc_instance.data["elements"]["generator"][gen_name]

        self._update_bids(gen_dict, bids, start_hour=0, horizon=options.ruc_horizon)

        return

    def assemble_project_tracking_signal(self, options, simulator, hour):

        """
        This function assembles the signals for the tracking model to estimate the
        state of the bidding model at the beginning of next RUC.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            hour: the simulation hour.

        Returns:
            market_signals: the market signals to be tracked.
        """

        tracking_horizon = len(self.projection_tracker.time_set)

        market_signals = self._assemble_sced_tracking_market_signals(
            hour=hour,
            sced_dispatch=None,
            tracking_horizon=tracking_horizon,
        )
        return market_signals

    def project_tracking_trajectory(self, options, simulator, ruc_hour):

        """
        This function projects the full power dispatch trajectory from the
        tracking model so we can use it to update the bidding model, i.e. advance
        the time for the bidding model.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            ruc_hour: the hour RUC is being solved

        Returns:
            full_projected_trajectory: the full projected power dispatch trajectory.
        """

        current_date = simulator.time_manager.current_time.date
        current_hour = simulator.time_manager.current_time.hour

        self._clone_tracking_model()

        for hour in range(ruc_hour, 24):

            # assemble market_signals
            market_signals = self.assemble_project_tracking_signal(
                options=options, simulator=simulator, hour=hour
            )
            # solve tracking
            self.projection_tracker.track_market_dispatch(
                market_dispatch=market_signals, date=current_date, hour=current_hour
            )

        # merge the trajectory
        full_projected_trajectory = {}
        for stat in self.tracker.daily_stats:
            full_projected_trajectory[stat] = self.tracker.daily_stats.get(
                stat
            ) + self.projection_tracker.daily_stats.get(stat)

        # clear the projection stats
        self.projection_tracker.daily_stats = None

        return full_projected_trajectory

    def _clone_tracking_model(self):
        """
        Clone the model in tracker and replace that of projection tracker. In this
        way, tracker and projection tracker have the same states before projection.

        Arguments:
            None

        Returns:
            None
        """

        # iterate all the variables and params and clone the values
        objects_list = [pyo.Var, pyo.Param]
        for obj in objects_list:
            for tracker_obj, proj_tracker_obj in zip_longest(
                self.tracker.model.component_objects(
                    obj, sort=pyo.SortComponents.alphabetizeComponentAndIndex
                ),
                self.projection_tracker.model.component_objects(
                    obj, sort=pyo.SortComponents.alphabetizeComponentAndIndex
                ),
            ):
                if tracker_obj.name != proj_tracker_obj.name:
                    raise ValueError(
                        f"Trying to copy the value of {tracker_obj} to {proj_tracker_obj}, but they do not have the same name and possibly not the corresponding objects. Please make sure tracker and projection tracker do not diverge. "
                    )
                for idx in tracker_obj.index_set():
                    if pyo.value(proj_tracker_obj[idx]) != pyo.value(tracker_obj[idx]):
                        proj_tracker_obj[idx] = round(pyo.value(tracker_obj[idx]), 4)

        return

    def _update_static_params(self, gen_dict, first_ruc):

        """
        Update static parameters in the Prescient generator parameter data dictionary depending on generator type.

        For a thermal generator, the p_cost data will be via ( MWh, $ ) pairs.
        For a renewable generator, the p_cost is a single cost.

        Args:
            gen_dict: Prescient generator parameter data dictionary.
            first_ruc: Is this the very first unit commitment problem

        Returns:
            None
        """

        is_thermal = (
            self.bidder.bidding_model_object.model_data.generator_type == "thermal"
        )
        is_renewable = (
            self.bidder.bidding_model_object.model_data.generator_type == "renewable"
        )

        for param, value in self.bidder.bidding_model_object.model_data:
            if param == "gen_name" or value is None:
                continue
            elif param == "p_cost":
                if is_thermal:
                    curve_value = convert_marginal_costs_to_actual_costs(value)
                    gen_dict[param] = {
                        "data_type": "cost_curve",
                        "cost_curve_type": "piecewise",
                        "values": curve_value,
                    }
                elif is_renewable:
                    gen_dict[param] = value
                else:
                    raise NotImplementedError(
                        "generator_type must be either 'thermal' or 'renewable'"
                    )

                if "p_fuel" in gen_dict:
                    gen_dict.pop("p_fuel")
            elif param == "initial_status" or param == "initial_p_output":
                if first_ruc:
                    gen_dict[param] = value
            else:
                gen_dict[param] = value

            if param == "startup_cost" and "startup_fuel" in gen_dict:
                gen_dict.pop("startup_fuel")

    def pass_static_params_to_DA(
        self, options, simulator, ruc_instance, ruc_date, ruc_hour
    ):
        """
        This method pass static generator parameters to RUC model in Prescient
        before it is solved.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            ruc_instance: Prescient RUC object.

            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

        Returns:
            None
        """

        gen_name = self.bidder.bidding_model_object.model_data.gen_name
        gen_dict = ruc_instance.data["elements"]["generator"][gen_name]
        first_ruc = (
            simulator.time_manager.current_time is None
            or simulator.time_manager.current_time
            == simulator.time_manager.get_first_time_step()
        )
        self._update_static_params(gen_dict, first_ruc=first_ruc)

        return

    def bid_into_DAM(self, options, simulator, ruc_instance, ruc_date, ruc_hour):

        """
        This function uses the bidding objects to bid into the day-ahead market
        (DAM).

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            ruc_instance: Prescient RUC object.

            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

        Returns:
            None
        """

        # check if it is first day
        is_first_day = simulator.time_manager.current_time is None

        if not is_first_day:

            # solve rolling horizon to get the trajectory
            full_projected_trajectory = self.project_tracking_trajectory(
                options, simulator, options.ruc_execution_hour
            )
            # update the bidding model
            self.bidder.update_day_ahead_model(**full_projected_trajectory)

        # generate bids
        bids = self.bidder.compute_day_ahead_bids(date=ruc_date)

        if is_first_day:
            self.current_bids = bids
        self.next_bids = bids

        # pass to prescient
        self._pass_DA_bid_to_prescient(options, ruc_instance, bids)

        return

    def fetch_DA_prices(self, options, simulator, result, uc_date, uc_hour):

        """
        This method fetches the day-ahead market prices from unit commitment results,
        and save it as a coordinator attribute.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            result: a Prescient RucPlan object.

            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

        Returns:
            None
        """

        current_bus = self.bidder.bidding_model_object.model_data.bus
        is_first_day = simulator.time_manager.current_time is None

        DA_prices = [
            result.ruc_market.day_ahead_prices.get((current_bus, t)) for t in range(24)
        ]

        if is_first_day:
            self.current_avail_DA_prices = DA_prices
        else:
            self.current_avail_DA_prices = self.current_DA_prices + DA_prices
        self.next_DA_prices = DA_prices

    def fetch_DA_dispatches(self, options, simulator, result, uc_date, uc_hour):

        """
        This method fetches the day-ahead dispatches from unit commitment results,
        and save it as a coordinator attribute.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            result: a Prescient RucPlan object.

            ruc_date: the date of the day-ahead market we bid into.

            ruc_hour: the hour the RUC is being solved in the day before.

        Returns:
            None
        """

        is_first_day = simulator.time_manager.current_time is None
        gen_name = self.bidder.bidding_model_object.model_data.gen_name
        gen_type = self.bidder.bidding_model_object.model_data.generator_type

        gen_type_dispatch_mapping = {
            "thermal": result.ruc_market.thermal_gen_cleared_DA,
            "renewable": result.ruc_market.renewable_gen_cleared_DA,
            "virtual": result.ruc_market.virtual_gen_cleared_DA,
        }

        dispatch_dict = gen_type_dispatch_mapping[gen_type]
        DA_dispatches = [dispatch_dict.get((gen_name, t)) for t in range(24)]

        if is_first_day:
            # self.current_DA_dispatches = DA_dispatches
            self.current_avail_DA_dispatches = DA_dispatches
        else:
            self.current_avail_DA_dispatches = (
                self.current_DA_dispatches + DA_dispatches
            )
        self.next_DA_dispatches = DA_dispatches

    def _pass_RT_bid_to_prescient(self, options, simulator, sced_instance, bids, hour):

        """
        This method passes the bids into the SCED model for real-time market
        clearing.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            sced_instance: Prescient SCED object

            bids: the bids we want to pass into the real-time market.

            hour: the hour of the real-time market.

        Returns:
            None
        """

        gen_name = self.bidder.bidding_model_object.model_data.gen_name

        # fetch generator's parameter dictionary from SCED instance
        gen_dict = sced_instance.data["elements"]["generator"][gen_name]

        self._update_bids(gen_dict, bids, start_hour=hour, horizon=options.sced_horizon)

        return

    def pass_static_params_to_RT(self, options, simulator, sced_instance):

        """
        This method pass static generator parameters to SCED model in Prescient
        before it is solved.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            sced_instance: Prescient SCED object.

        Returns:
            None
        """

        gen_name = self.bidder.bidding_model_object.model_data.gen_name
        gen_dict = sced_instance.data["elements"]["generator"][gen_name]
        self._update_static_params(gen_dict, first_ruc=False)

        return

    def bid_into_RTM(self, options, simulator, sced_instance):

        """
        This function bids into the real-time market. At this moment I just copy the
        corresponding day-ahead bid here.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            sced_instance: Prescient SCED object.

        Returns:
            None
        """

        # fetch the bids
        date = simulator.time_manager.current_time.date
        hour = simulator.time_manager.current_time.hour

        bids = self.bidder.compute_real_time_bids(
            date=date,
            hour=hour,
            realized_day_ahead_prices=self.current_avail_DA_prices,
            realized_day_ahead_dispatches=self.current_avail_DA_dispatches,
        )

        # pass bids into sced model
        self._pass_RT_bid_to_prescient(options, simulator, sced_instance, bids, hour)

        return

    def assemble_sced_tracking_market_signals(
        self, options, simulator, sced_instance, hour
    ):

        """
        This function assembles the signals for the tracking model.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            sced_instance: Prescient SCED object

            hour: the simulation hour.

        Returns:
            market_signals: the market signals to be tracked.
        """

        gen_name = self.bidder.bidding_model_object.model_data.gen_name

        # fecth the sced signals for the generation from sced instance
        sced_dispatch = sced_instance.data["elements"]["generator"][gen_name]["pg"][
            "values"
        ]
        tracking_horizon = len(self.tracker.time_set)

        market_signals = self._assemble_sced_tracking_market_signals(
            hour=hour,
            sced_dispatch=sced_dispatch,
            tracking_horizon=tracking_horizon,
        )

        return market_signals

    def _assemble_sced_tracking_market_signals(
        self, hour, sced_dispatch, tracking_horizon
    ):

        """
        This function assembles the signals for the tracking model.

        Arguments:

            hour: the simulation hour

            sced_dispatch: current sced dispatch (a list)

            tracking_horizon: length of the tracking horizon


        Returns:
            market_signals: the market signals to be tracked.
        """

        market_signals = []

        # append the sced dispatch
        if sced_dispatch is None:
            dispatch = self.current_DA_dispatches[hour]
        else:
            dispatch = sced_dispatch[0]
        market_signals.append(dispatch)

        # append corresponding RUC dispatch
        for t in range(hour + 1, hour + tracking_horizon):
            if t < len(self.current_avail_DA_dispatches):
                market_signals.append(self.current_avail_DA_dispatches[t])

        return market_signals

    def track_sced_signal(self, options, simulator, sced_instance, lmp_sced):

        """
        This methods uses the tracking object to track the current real-time market
        signals.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            sced_instance: Prescient SCED object

            lmp_sced: Prescient SCED LMP object

        Returns:
            None

        """

        current_date = simulator.time_manager.current_time.date
        current_hour = simulator.time_manager.current_time.hour

        # get market signals
        market_signals = self.assemble_sced_tracking_market_signals(
            options=options,
            simulator=simulator,
            sced_instance=sced_instance,
            hour=current_hour,
        )

        # actual tracking
        implemented_profiles = self.tracker.track_market_dispatch(
            market_dispatch=market_signals, date=current_date, hour=current_hour
        )

        # update models
        self.tracker.update_model(**implemented_profiles)
        self.bidder.update_real_time_model(**implemented_profiles)

        return

    def update_observed_dispatch(self, options, simulator, ops_stats):

        """
        This methods extract the actual power delivered by the tracking model and
        inform Prescient, so Prescient can use this data to calculate the settlement
        and etc.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            ops_stats: Prescient operation statitstic object

        Returns:
            None

        """
        g = self.bidder.bidding_model_object.model_data.gen_name

        # store the dictionaries that the observed/delivered power levels
        observed_dispatch_level_dicts = [
            ops_stats.observed_thermal_dispatch_levels,
            ops_stats.observed_renewables_levels,
            ops_stats.observed_virtual_dispatch_levels,
        ]

        for observed_dispatch_level in observed_dispatch_level_dicts:
            if g in observed_dispatch_level:
                observed_dispatch_level[g] = self.tracker.get_last_delivered_power()

        return

    def activate_pending_DA_data(self, options, simulator):

        """
        This function puts the day-ahead data computed in the day before into effect,
        i.e. the data for the next day become the data for the current day.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

        Returns:
            None
        """

        self.current_bids, self.next_bids = self.next_bids, None
        self.current_DA_prices, self.next_DA_prices = self.next_DA_prices, None
        self.current_DA_dispatches, self.next_DA_dispatches = (
            self.next_DA_dispatches,
            None,
        )

        # update the available DA signals
        self.current_avail_DA_prices = self.current_DA_prices
        self.current_avail_DA_dispatches = self.current_DA_dispatches

    def write_plugin_results(self, options, simulator):

        """
        After the simulation is completed, the plugins can write their own customized
        results. Each plugin will have to have a method named 'write_results'.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

        Returns:
            None

        """

        self.bidder.write_results(path=options.output_directory)
        self.tracker.write_results(path=options.output_directory)

        return
