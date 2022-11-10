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
import pyomo.environ as pyo
from collections import deque
import pandas as pd
from idaes.apps.grid_integration import Tracker
from idaes.apps.grid_integration import Bidder
from idaes.apps.grid_integration import PlaceHolderForecaster
from idaes.apps.grid_integration.model_data import ThermalGeneratorModelData

from pyomo.common.dependencies import attempt_import

prescient_simulator, prescient_avail = attempt_import("prescient.simulator")


class ThermalGenerator:

    """
    Simple thermal generator model (MIP). Equations models are from Gao, Knueven,
    Siirola, Miller, Dowling (2022). Multiscale Simulation of Integrated Energy
    System and Electricity Market Interactions. Applied Energy.
    """

    # Using 4 segments to be consistent with models in RTS-GMLC dataset
    segment_number = 4

    def __init__(
        self,
        rts_gmlc_generator_dataframe,
        rts_gmlc_bus_dataframe,
        generator="102_STEAM_3",
    ):

        """
        Initializes the class object by building the thermal generator model.

        Arguments:
            rts_gmlc_generator_dataframe: the RTS-GMLC generator data in Pandas DataFrame
            horizon: the length of the planning horizon of the model.
            generator: a generator in RTS-GMLC

        Returns:
            None
        """

        self.generator = generator
        rts_gmlc_generator_dataframe = rts_gmlc_generator_dataframe.merge(
            rts_gmlc_bus_dataframe[["Bus ID", "Bus Name"]],
            how="left",
            left_on="Bus ID",
            right_on="Bus ID",
        )
        self._model_data_dict = self.assemble_model_data(
            generator_name=generator, gen_params=rts_gmlc_generator_dataframe
        )
        self.result_list = []

    def assemble_model_data(self, generator_name, gen_params):

        """
        This function assembles the parameter data to build the thermal generator
        model, given a list of generator names and the RTS-GMLC data directory.

        Arguments:
            generator_names: a generator name in RTS-GMLC dataset.
            gen_params: the RTS-GMLC generator data in Pandas DataFrame

        Returns:
            model_data: a dictionary which has this structure
            {data type name: value}.
        """

        gen_params = gen_params.set_index("GEN UID", inplace=False)
        properties = [
            "Bus Name",
            "PMin MW",
            "PMax MW",
            "Min Up Time Hr",
            "Min Down Time Hr",
            "Ramp Rate MW/Min",
            "Start Heat Warm MBTU",
            "Fuel Price $/MMBTU",
            "HR_avg_0",
            "HR_incr_1",
            "HR_incr_2",
            "HR_incr_3",
            "Output_pct_1",
            "Output_pct_2",
            "Output_pct_3",
        ]

        # to dict
        model_data = gen_params.loc[generator_name, properties].to_dict()

        model_data["RU"] = model_data["Ramp Rate MW/Min"] * 60
        model_data["RD"] = model_data["RU"]
        model_data["SU"] = min(model_data["PMin MW"], model_data["RU"])
        model_data["SD"] = min(model_data["PMin MW"], model_data["RD"])
        model_data["SU Cost"] = (
            model_data["Start Heat Warm MBTU"] * model_data["Fuel Price $/MMBTU"]
        )

        model_data["Min Load Cost"] = (
            model_data["HR_avg_0"]
            / 1000
            * model_data["Fuel Price $/MMBTU"]
            * model_data["PMin MW"]
        )

        model_data["Power Segments"] = {}
        model_data["Marginal Costs"] = {}

        model_data["Original Marginal Cost Curve"] = {}
        model_data["Original Marginal Cost Curve"][model_data["PMin MW"]] = (
            model_data["Min Load Cost"] / model_data["PMin MW"]
        )

        for l in range(1, ThermalGenerator.segment_number):
            model_data["Power Segments"][l] = (
                model_data["Output_pct_{}".format(l)] * model_data["PMax MW"]
            )
            model_data["Marginal Costs"][l] = (
                model_data["HR_incr_{}".format(l)]
                / 1000
                * model_data["Fuel Price $/MMBTU"]
            )
            model_data["Original Marginal Cost Curve"][
                model_data["Power Segments"][l]
            ] = model_data["Marginal Costs"][l]

        self._model_data = ThermalGeneratorModelData(
            gen_name=generator_name,
            bus=model_data["Bus Name"],
            p_min=model_data["PMin MW"],
            p_max=model_data["PMax MW"],
            min_down_time=model_data["Min Down Time Hr"],
            min_up_time=model_data["Min Up Time Hr"],
            ramp_up_60min=model_data["RU"],
            ramp_down_60min=model_data["RD"],
            shutdown_capacity=model_data["SD"],
            startup_capacity=model_data["SU"],
            initial_status=-1,
            initial_p_output=0,
            fixed_commitment=None,
            production_cost_bid_pairs=[
                (
                    model_data["PMin MW"],
                    model_data["Min Load Cost"] / model_data["PMin MW"],
                )
            ]
            + [
                (model_data["Power Segments"][l], model_data["Marginal Costs"][l])
                for l in range(1, ThermalGenerator.segment_number)
            ],
            startup_cost_pairs=[
                (model_data["Min Down Time Hr"], model_data["SU Cost"])
            ],
        )

        return model_data

    @property
    def model_data(self):
        return self._model_data

    @staticmethod
    def _add_UT_DT_constraints(b):

        """
        This function adds the minimum up/down time constraints using eq. 4 - 5
        in "On mixed-integer programming formulations for the unit commitment
        problem". INFORMS Journal on Computing, 32(4), pp.857-876. Knueven, B.,
        Ostrowski, J. and Watson, J.P., 2020.

        Arguments:
            b: a pyomo block

        Returns:
            None
        """

        def pre_shut_down_trajectory_set_rule(b):
            return (t for t in range(-pyo.value(b.min_dw_time) + 1, 0))

        b.pre_shut_down_trajectory_set = pyo.Set(
            dimen=1, initialize=pre_shut_down_trajectory_set_rule, ordered=True
        )

        def pre_start_up_trajectory_set_rule(b):
            return (t for t in range(-pyo.value(b.min_up_time) + 1, 0))

        b.pre_start_up_trajectory_set = pyo.Set(
            dimen=1, initialize=pre_start_up_trajectory_set_rule, ordered=True
        )

        b.pre_shut_down_trajectory = pyo.Param(
            b.pre_shut_down_trajectory_set, initialize=0, mutable=True
        )
        b.pre_start_up_trajectory = pyo.Param(
            b.pre_start_up_trajectory_set, initialize=0, mutable=True
        )

        def min_down_time_rule(b, h):
            if h < pyo.value(b.min_dw_time):
                return (
                    sum(
                        b.pre_shut_down_trajectory[h0]
                        for h0 in range(h - pyo.value(b.min_dw_time) + 1, 0)
                    )
                    + sum(b.shut_dw[h0] for h0 in range(h + 1))
                    <= 1 - b.on_off[h]
                )
            else:
                return (
                    sum(
                        b.shut_dw[h0]
                        for h0 in range(h - pyo.value(b.min_dw_time) + 1, h + 1)
                    )
                    <= 1 - b.on_off[h]
                )

        b.min_down_time_con = pyo.Constraint(b.HOUR, rule=min_down_time_rule)

        def min_up_time_rule(b, h):
            if h < pyo.value(b.min_up_time):
                return (
                    sum(
                        b.pre_start_up_trajectory[h0]
                        for h0 in range(h - pyo.value(b.min_up_time) + 1, 0)
                    )
                    + sum(b.start_up[h0] for h0 in range(h + 1))
                    <= b.on_off[h]
                )
            else:
                return (
                    sum(
                        b.start_up[h0]
                        for h0 in range(h - pyo.value(b.min_up_time) + 1, h + 1)
                    )
                    <= b.on_off[h]
                )

        b.min_up_time_con = pyo.Constraint(b.HOUR, rule=min_up_time_rule)

        return

    def populate_model(self, b, horizon):

        """
        This function builds the model for a thermal generator.

        Arguments:
            b: a pyomo block

        Returns:
            b: the constructed block.
        """

        model_data = self._model_data_dict

        ## define the sets
        b.HOUR = pyo.Set(initialize=list(range(horizon)))
        b.SEGMENTS = pyo.Set(initialize=list(range(1, ThermalGenerator.segment_number)))

        ## define the parameters

        b.start_up_cost = pyo.Param(initialize=model_data["SU Cost"], mutable=False)

        # capacity of generators: upper bound (MW)
        b.Pmax = pyo.Param(initialize=model_data["PMax MW"], mutable=False)

        # minimum power of generators: lower bound (MW)
        b.Pmin = pyo.Param(initialize=model_data["PMin MW"], mutable=False)

        b.power_segment_bounds = pyo.Param(
            b.SEGMENTS, initialize=model_data["Power Segments"], mutable=False
        )

        # get the cost slopes
        b.F = pyo.Param(
            b.SEGMENTS, initialize=model_data["Marginal Costs"], mutable=False
        )

        b.min_load_cost = pyo.Param(
            initialize=model_data["Min Load Cost"], mutable=False
        )

        # Ramp up limits (MW/h)
        b.ramp_up = pyo.Param(initialize=model_data["RU"], mutable=False)

        # Ramp down limits (MW/h)
        b.ramp_dw = pyo.Param(initialize=model_data["RD"], mutable=False)

        # start up ramp limit
        b.ramp_start_up = pyo.Param(initialize=model_data["SU"], mutable=False)

        # shut down ramp limit
        b.ramp_shut_dw = pyo.Param(initialize=model_data["SD"], mutable=False)

        # minimum down time [hr]
        b.min_dw_time = pyo.Param(
            initialize=int(model_data["Min Down Time Hr"]), mutable=False
        )

        # minimum up time [hr]
        b.min_up_time = pyo.Param(
            initialize=int(model_data["Min Up Time Hr"]), mutable=False
        )

        # on/off status from previous day
        b.pre_on_off = pyo.Param(within=pyo.Binary, default=1, mutable=True)

        # define a function to initialize the previous power params
        def init_pre_pow_fun(b):
            return b.pre_on_off * b.Pmin

        b.pre_P_T = pyo.Param(initialize=init_pre_pow_fun, mutable=True)

        ## define the variables

        # generator power (MW)

        # power generated by thermal generator
        b.P_T = pyo.Var(
            b.HOUR, initialize=model_data["PMin MW"], within=pyo.NonNegativeReals
        )

        # binary variables indicating on/off
        b.on_off = pyo.Var(b.HOUR, initialize=True, within=pyo.Binary)

        # binary variables indicating  start_up
        b.start_up = pyo.Var(b.HOUR, initialize=False, within=pyo.Binary)

        # binary variables indicating shut down
        b.shut_dw = pyo.Var(b.HOUR, initialize=False, within=pyo.Binary)

        # power produced in each segment
        b.power_segment = pyo.Var(
            b.HOUR, b.SEGMENTS, initialize=0, within=pyo.NonNegativeReals
        )

        ## Constraints

        # bounds on gen_pow
        def lhs_bnd_gen_pow_fun(b, h):
            return b.on_off[h] * b.Pmin <= b.P_T[h]

        b.lhs_bnd_gen_pow = pyo.Constraint(b.HOUR, rule=lhs_bnd_gen_pow_fun)

        def rhs_bnd_gen_pow_fun(b, h):
            return b.P_T[h] <= b.on_off[h] * b.Pmax

        b.rhs_bnd_gen_pow = pyo.Constraint(b.HOUR, rule=rhs_bnd_gen_pow_fun)

        # linearized power
        def linear_power_fun(b, h):
            return (
                b.P_T[h]
                == sum(b.power_segment[h, l] for l in b.SEGMENTS) + b.Pmin * b.on_off[h]
            )

        b.linear_power = pyo.Constraint(b.HOUR, rule=linear_power_fun)

        # bounds on segment power
        def seg_pow_bnd_fun(b, h, l):
            if l == 1:
                return (
                    b.power_segment[h, l]
                    <= (b.power_segment_bounds[l] - b.Pmin) * b.on_off[h]
                )
            else:
                return (
                    b.power_segment[h, l]
                    <= (b.power_segment_bounds[l] - b.power_segment_bounds[l - 1])
                    * b.on_off[h]
                )

        b.seg_pow_bnd = pyo.Constraint(b.HOUR, b.SEGMENTS, rule=seg_pow_bnd_fun)

        # start up and shut down logic (Garver 1962)
        def start_up_shut_dw_fun(b, h):
            if h == 0:
                return b.start_up[h] - b.shut_dw[h] == b.on_off[h] - b.pre_on_off
            else:
                return b.start_up[h] - b.shut_dw[h] == b.on_off[h] - b.on_off[h - 1]

        b.start_up_shut_dw = pyo.Constraint(b.HOUR, rule=start_up_shut_dw_fun)

        # either start up or shut down
        def start_up_or_shut_dw_fun(b, h):
            return b.start_up[h] + b.shut_dw[h] <= 1

        b.start_up_or_shut_dw = pyo.Constraint(b.HOUR, rule=start_up_or_shut_dw_fun)

        # ramp up limits
        def ramp_up_fun(b, h):
            """
            h stand for hour
            """
            if h == 0:
                return (
                    b.P_T[h]
                    <= b.pre_P_T
                    + b.ramp_up * b.pre_on_off
                    + b.ramp_start_up * b.start_up[h]
                )
            else:
                return (
                    b.P_T[h]
                    <= b.P_T[h - 1]
                    + b.ramp_up * b.on_off[h - 1]
                    + b.ramp_start_up * b.start_up[h]
                )

        b.ramp_up_con = pyo.Constraint(b.HOUR, rule=ramp_up_fun)

        # ramp shut down limits
        def ramp_shut_dw_fun(b, h):
            """
            h stand for hour.
            """
            if h == 0:
                return b.pre_P_T <= b.Pmax * b.on_off[h] + b.ramp_shut_dw * b.shut_dw[h]
            else:
                return (
                    b.P_T[h - 1] <= b.Pmax * b.on_off[h] + b.ramp_shut_dw * b.shut_dw[h]
                )

        b.ramp_shut_dw_con = pyo.Constraint(b.HOUR, rule=ramp_shut_dw_fun)

        # ramp down limits
        def ramp_dw_fun(b, h):
            """
            h stand for hour.
            """
            if h == 0:
                return (
                    b.pre_P_T - b.P_T[h]
                    <= b.ramp_dw * b.on_off[h] + b.ramp_shut_dw * b.shut_dw[h]
                )
            else:
                return (
                    b.P_T[h - 1] - b.P_T[h]
                    <= b.ramp_dw * b.on_off[h] + b.ramp_shut_dw * b.shut_dw[h]
                )

        b.ramp_dw_con = pyo.Constraint(b.HOUR, rule=ramp_dw_fun)

        ## add min up and down time constraints
        self._add_UT_DT_constraints(b)

        ## Expression
        def prod_cost_fun(b, h):
            return b.min_load_cost * b.on_off[h] + sum(
                b.F[l] * b.power_segment[h, l] for l in b.SEGMENTS
            )

        b.prod_cost_approx = pyo.Expression(b.HOUR, rule=prod_cost_fun)

        # start up costs
        def start_cost_fun(b, h):
            return b.start_up_cost * b.start_up[h]

        b.start_up_cost_expr = pyo.Expression(b.HOUR, rule=start_cost_fun)

        # total cost
        def tot_cost_fun(b, h):
            return b.prod_cost_approx[h] + b.start_up_cost_expr[h]

        b.tot_cost = pyo.Expression(b.HOUR, rule=tot_cost_fun)

        return

    @staticmethod
    def _update_UT_DT(b, implemented_shut_down, implemented_start_up):
        """
        This method updates the parameters in the minimum up/down time
        constraints based on the implemented shut down and start up events.

        Arguments:
            b: the block that needs to be updated
            implemented_shut_down: realized shut down events: [].
            implemented_start_up: realized start up events: []

        Returns:
            None
        """

        pre_shut_down_trajectory_copy = deque([])
        pre_start_up_trajectory_copy = deque([])

        # copy old trajectory
        for t in b.pre_shut_down_trajectory_set:
            pre_shut_down_trajectory_copy.append(
                round(pyo.value(b.pre_shut_down_trajectory[t]))
            )
        for t in b.pre_start_up_trajectory_set:
            pre_start_up_trajectory_copy.append(
                round(pyo.value(b.pre_start_up_trajectory[t]))
            )

        # add implemented trajectory to the queue
        pre_shut_down_trajectory_copy += deque(implemented_shut_down)
        pre_start_up_trajectory_copy += deque(implemented_start_up)

        # pop out outdated trajectory
        while len(pre_shut_down_trajectory_copy) > pyo.value(b.min_dw_time) - 1:
            pre_shut_down_trajectory_copy.popleft()
        while len(pre_start_up_trajectory_copy) > pyo.value(b.min_up_time) - 1:
            pre_start_up_trajectory_copy.popleft()

        # actual update
        for t in b.pre_shut_down_trajectory_set:
            b.pre_shut_down_trajectory[t] = pre_shut_down_trajectory_copy.popleft()

        for t in b.pre_start_up_trajectory_set:
            b.pre_start_up_trajectory[t] = pre_start_up_trajectory_copy.popleft()

        return

    @staticmethod
    def _update_power(b, implemented_power_output):
        """
        This method updates the parameters in the ramping constraints based on
        the implemented power outputs.

        Arguments:
            b: the block that needs to be updated
            implemented_power_output: realized power outputs: []

         Returns:
             None
        """

        b.pre_P_T = round(implemented_power_output[-1], 2)
        b.pre_on_off = round(int(implemented_power_output[-1] > 1e-3))

        return

    def update_model(
        self, b, implemented_shut_down, implemented_start_up, implemented_power_output
    ):

        """
        This method updates the parameters in the model based on
        the implemented power outputs, shut down and start up events.

        Arguments:
            b: the block that needs to be updated
            implemented_shut_down: realized shut down events: [].
            implemented_start_up: realized start up events: []
            implemented_power_output: realized power outputs: []

         Returns:
             None
        """

        self._update_UT_DT(b, implemented_shut_down, implemented_start_up)
        self._update_power(b, implemented_power_output)

        return

    @staticmethod
    def get_implemented_profile(b, last_implemented_time_step):

        """
        This method gets the implemented variable profiles in the last optimization
        solve.

        Arguments:
            b: the block
            model_var: intended variable name in str
            last_implemented_time_step: time index for the last implemented time
                                        step

         Returns:
             profile: the intended profile, {unit: [...]}
        """

        implemented_shut_down = deque(
            [pyo.value(b.shut_dw[t]) for t in range(last_implemented_time_step + 1)]
        )
        implemented_start_up = deque(
            [pyo.value(b.start_up[t]) for t in range(last_implemented_time_step + 1)]
        )
        implemented_power_output = deque(
            [pyo.value(b.P_T[t]) for t in range(last_implemented_time_step + 1)]
        )

        return {
            "implemented_shut_down": implemented_shut_down,
            "implemented_start_up": implemented_start_up,
            "implemented_power_output": implemented_power_output,
        }

    @staticmethod
    def get_last_delivered_power(b, last_implemented_time_step):

        """
        Returns the last delivered power output.

        Arguments:
            None

        Returns:
            None
        """

        return pyo.value(b.P_T[last_implemented_time_step])

    def record_results(self, b, date=None, hour=None, **kwargs):

        """
        Record the operations stats for the model.

        Arguments:

            date: current simulation date
            hour: current simulation hour

        Returns:
            None

        """

        df_list = []

        for t in b.HOUR:

            result_dict = {}
            result_dict["Generator"] = self.generator
            result_dict["Date"] = date
            result_dict["Hour"] = hour

            # simulation inputs
            result_dict["Horizon [hr]"] = int(t)

            # model vars
            result_dict["Thermal Power Generated [MW]"] = float(
                round(pyo.value(b.P_T[t]), 2)
            )

            result_dict["On/off [bin]"] = int(round(pyo.value(b.on_off[t])))
            result_dict["Start Up [bin]"] = int(round(pyo.value(b.start_up[t])))
            result_dict["Shut Down [bin]"] = int(round(pyo.value(b.shut_dw[t])))

            result_dict["Production Cost [$]"] = float(
                round(pyo.value(b.prod_cost_approx[t]), 2)
            )
            result_dict["Start-up Cost [$]"] = float(
                round(pyo.value(b.start_up_cost_expr[t]), 2)
            )
            result_dict["Total Cost [$]"] = float(round(pyo.value(b.tot_cost[t]), 2))

            # calculate mileage
            if t == 0:
                result_dict["Mileage [MW]"] = float(
                    round(abs(pyo.value(b.P_T[t] - b.pre_P_T)), 2)
                )
            else:
                result_dict["Mileage [MW]"] = float(
                    round(abs(pyo.value(b.P_T[t] - b.P_T[t - 1])), 2)
                )

            for key in kwargs:
                result_dict[key] = kwargs[key]

            result_df = pd.DataFrame.from_dict(result_dict, orient="index")
            df_list.append(result_df.T)

        # save the result to object property
        # wait to be written when simulation ends
        self.result_list.append(pd.concat(df_list))

        return

    def write_results(self, path):

        """
        This methods writes the saved operation stats into an csv file.

        Arguments:
            path: the path to write the results.

        Return:
            None
        """

        pd.concat(self.result_list).to_csv(path, index=False)

    @property
    def power_output(self):
        return "P_T"

    @property
    def total_cost(self):
        return ("tot_cost", 1)


if __name__ == "__main__":

    from idaes.apps.grid_integration.examples.utils import (
        rts_gmlc_generator_dataframe,
        rts_gmlc_bus_dataframe,
        prescient_5bus,
        daily_da_price_means,
        daily_rt_price_means,
        daily_da_price_stds,
        daily_rt_price_stds,
    )

    generator = "10_STEAM"
    horizon = 4
    solver = pyo.SolverFactory("cbc")

    run_tracker = True
    run_bidder = True
    run_prescient = True

    if run_tracker:

        # create a tracker model
        tracking_model_object = ThermalGenerator(
            rts_gmlc_generator_dataframe=rts_gmlc_generator_dataframe,
            rts_gmlc_bus_dataframe=rts_gmlc_bus_dataframe,
            generator=generator,
        )
        # make a tracker
        thermal_tracker = Tracker(
            tracking_model_object=tracking_model_object,
            tracking_horizon=horizon,
            n_tracking_hour=1,
            solver=solver,
        )

        market_dispatch = [30, 40, 50, 70]

        thermal_tracker.track_market_dispatch(
            market_dispatch=market_dispatch, date="2021-07-26", hour="17:00"
        )
        thermal_tracker.write_results(path="./")

    if run_bidder:

        # create a tracker model
        bidding_model_object = ThermalGenerator(
            rts_gmlc_generator_dataframe=rts_gmlc_generator_dataframe,
            rts_gmlc_bus_dataframe=rts_gmlc_bus_dataframe,
            generator=generator,
        )

        # create forecaster
        forecaster = PlaceHolderForecaster(
            daily_da_price_means=daily_da_price_means,
            daily_rt_price_means=daily_rt_price_means,
            daily_da_price_stds=daily_da_price_stds,
            daily_rt_price_stds=daily_rt_price_stds,
        )

        thermal_bidder = Bidder(
            bidding_model_object=bidding_model_object,
            day_ahead_horizon=48,
            real_time_horizon=4,
            n_scenario=10,
            solver=solver,
            forecaster=forecaster,
        )

        date = "2020-07-10"
        hour = 13
        bids = thermal_bidder.compute_day_ahead_bids(date, hour)
        thermal_bidder.write_results(path="./")

    if run_prescient and prescient_avail:

        options = {
            "data_path": prescient_5bus,
            "input_format": "rts-gmlc",
            "simulate_out_of_sample": True,
            "run_sced_with_persistent_forecast_errors": True,
            "output_directory": "bidding_plugin_test_thermal_generator",
            "start_date": "07-10-2020",
            "num_days": 2,
            "sced_horizon": 4,
            "compute_market_settlements": True,
            "day_ahead_pricing": "LMP",
            "ruc_mipgap": 0.01,
            "symbolic_solver_labels": True,
            "reserve_factor": 0.0,
            "deterministic_ruc_solver": "cbc",
            "sced_solver": "cbc",
            "plugin": {
                "doubleloop": {
                    "module": "thermal_generator_prescient_plugin.py",
                    "bidding_generator": generator,
                }
            },
        }

        prescient_simulator.Prescient().simulate(**options)
