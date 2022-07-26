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

            tracker: an initialized bidder object

            projection_tracker: an initialized bidder object, this object is
                                mimicking the behaviro of the projection SCED in
                                Prescient and to projecting the system states
                                and updating bidder model.

        Returns:
            None
        """

        self.bidder = bidder
        self.tracker = tracker
        self.projection_tracker = projection_tracker

    def register_plugins(
        self,
        context,
        options,
        plugin_config,
    ):

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
        context.register_before_ruc_solve_callback(self.pass_static_params_to_DA)
        context.register_before_ruc_solve_callback(self.bid_into_DAM)
        context.register_before_operations_solve_callback(self.pass_static_params_to_RT)
        context.register_before_operations_solve_callback(self.bid_into_RTM)
        context.register_after_operations_callback(self.track_sced_signal)
        context.register_update_operations_stats_callback(self.update_observed_dispatch)
        context.register_after_ruc_activation_callback(self.activate_DA_bids)
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
                gen_dict[param_name] = bids[0][gen_name].get(param_name, None)

            return

        param_update_func_map = {
            "p_cost": _update_p_cost,
            "p_max": _update_time_series_params,
            "p_min": _update_time_series_params,
            "fixed_commitment": _update_time_series_params,
            "min_up_time": _update_non_time_series_params,
            "min_down_time": _update_non_time_series_params,
            "startup_capacity": _update_time_series_params,
            "shutdown_capacity": _update_time_series_params,
            "startup_fuel": _update_non_time_series_params,
            "startup_cost": _update_non_time_series_params,
        }

        for param_name in bids[0][gen_name].keys():
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
        state of the bidding model at the begining of next RUC.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

            hour: the simulation hour.

        Returns:
            market_signals: the market signals to be tracked.
        """

        gen_name = self.bidder.bidding_model_object.model_data.gen_name

        # store the dictionaries that have the current ruc signals
        current_ruc_dispatch_dicts = [
            simulator.data_manager.ruc_market_active.thermal_gen_cleared_DA,
            simulator.data_manager.ruc_market_active.renewable_gen_cleared_DA,
            simulator.data_manager.ruc_market_active.virtual_gen_cleared_DA,
        ]

        tracking_horizon = len(self.projection_tracker.time_set)

        market_signals = self._assemble_sced_tracking_market_signals(
            gen_name=gen_name,
            hour=hour,
            sced_dispatch=None,
            tracking_horizon=tracking_horizon,
            current_ruc_dispatch_dicts=current_ruc_dispatch_dicts,
            next_ruc_dispatch_dicts=None,
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

    def _update_static_params(self, gen_dict):

        """
        Update static parameters in the Prescient generator parameter data dictionary.

        Args:
            gen_dict: Prescient generator parameter data dictionary.

        Returns:
            None
        """

        for param, value in self.bidder.bidding_model_object.model_data:
            if param == "gen_name" or value is None:
                continue
            elif param == "p_cost":
                curve_value = convert_marginal_costs_to_actual_costs(value)
                gen_dict[param] = {
                    "data_type": "cost_curve",
                    "cost_curve_type": "piecewise",
                    "values": curve_value,
                }
                if "p_fuel" in gen_dict:
                    gen_dict.pop("p_fuel")
            else:
                # update generator parameters with value
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
        self._update_static_params(gen_dict)

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
            self.bidder.update_model(**full_projected_trajectory)

        # generate bids
        bids = self.bidder.compute_bids(date=ruc_date)

        if is_first_day:
            self.current_bids = bids
        self.next_bids = bids

        # pass to prescient
        self._pass_DA_bid_to_prescient(options, ruc_instance, bids)

        return

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
        self._update_static_params(gen_dict)

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
        hour = simulator.time_manager.current_time.hour
        bids = self.current_bids

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

        # store the dictionaries that have the current ruc signals
        current_ruc_dispatch_dicts = [
            simulator.data_manager.ruc_market_active.thermal_gen_cleared_DA,
            simulator.data_manager.ruc_market_active.renewable_gen_cleared_DA,
            simulator.data_manager.ruc_market_active.virtual_gen_cleared_DA,
        ]

        if simulator.data_manager.ruc_market_pending is not None:
            next_ruc_dispatch_dicts = [
                simulator.data_manager.ruc_market_pending.thermal_gen_cleared_DA,
                simulator.data_manager.ruc_market_pending.renewable_gen_cleared_DA,
                simulator.data_manager.ruc_market_pending.virtual_gen_cleared_DA,
            ]

        else:
            next_ruc_dispatch_dicts = None

        market_signals = self._assemble_sced_tracking_market_signals(
            gen_name=gen_name,
            hour=hour,
            sced_dispatch=sced_dispatch,
            tracking_horizon=tracking_horizon,
            current_ruc_dispatch_dicts=current_ruc_dispatch_dicts,
            next_ruc_dispatch_dicts=next_ruc_dispatch_dicts,
        )

        return market_signals

    @staticmethod
    def _assemble_sced_tracking_market_signals(
        gen_name,
        hour,
        sced_dispatch,
        tracking_horizon,
        current_ruc_dispatch_dicts,
        next_ruc_dispatch_dicts=None,
    ):

        """
        This function assembles the signals for the tracking model.

        Arguments:
            gen_name: the generator's name

            hour: the simulation hour

            sced_dispatch: current sced dispatch (a list)

            tracking_horizon: length of the tracking horizon

            current_ruc_dispatch_dicts: current day's unit commiment dispatch
            dictionaries, including profiles for thermal, renewable and etc.

            next_ruc_dispatch_dicts: next day's unit commiment dispatch
            dictionaries, including profiles for thermal, renewable and etc.


        Returns:
            market_signals: the market signals to be tracked.
        """

        def get_signals(gen_name, t, ruc_dispatch_dicts, market_signals):

            dispatch = None
            for ruc_dispatch in ruc_dispatch_dicts:
                dispatch = ruc_dispatch.get((gen_name, t), None)
                if dispatch is not None:
                    break

            if dispatch is None and len(market_signals) > 0:
                dispatch = market_signals[-1]
            elif dispatch is None:
                raise ValueError(
                    f"No SCED/RUC signal has been found for generator {gen_name} at hour {t}. No previous signal is available for repeating."
                )

            return dispatch

        market_signals = []
        # append the sced dispatch
        if sced_dispatch is None:
            dispatch = get_signals(
                gen_name, hour, current_ruc_dispatch_dicts, market_signals
            )
        else:
            dispatch = sced_dispatch[0]
        market_signals.append(dispatch)

        # append corresponding RUC dispatch
        for t in range(hour + 1, hour + tracking_horizon):

            # next ruc is available: fetch the signal from next ruc
            if t > 23 and next_ruc_dispatch_dicts:
                t = t % 24
                dispatch = get_signals(
                    gen_name, t, next_ruc_dispatch_dicts, market_signals
                )

            # fetch from the current ruc
            else:
                dispatch = get_signals(
                    gen_name, t, current_ruc_dispatch_dicts, market_signals
                )
            market_signals.append(dispatch)

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
        self.tracker.track_market_dispatch(
            market_dispatch=market_signals, date=current_date, hour=current_hour
        )

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

    def activate_DA_bids(self, options, simulator):

        """
        This function puts the day-ahead bids computed in the day before into effect,
        i.e. the bids for the next day become the bids for the current day.

        Arguments:
            options: Prescient options from prescient.simulator.config.

            simulator: Prescient simulator.

        Returns:
            None
        """

        # change bids
        self.current_bids = self.next_bids
        self.next_bids = None

        return

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
