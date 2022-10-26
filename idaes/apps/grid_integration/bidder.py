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
import pandas as pd
import pyomo.environ as pyo
from pyomo.opt.base.solvers import OptSolver
import os
from abc import ABC, abstractmethod
from idaes.apps.grid_integration.utils import convert_marginal_costs_to_actual_costs
import datetime
from pyomo.common.dependencies import attempt_import

egret, egret_avail = attempt_import("egret")
if egret_avail:
    from egret.model_library.transmission import tx_utils


class AbstractBidder(ABC):

    """
    The abstract class for all the bidder and self-schedulers.
    """

    @abstractmethod
    def update_day_ahead_model(self, **kwargs):

        """
        Update the day-ahead model (advance timesteps) with necessary parameters in kwargs.

        Arguments:
            kwargs: necessary profiles to update the underlying model. {stat_name: [...]}

        Returns:
            None
        """

        pass

    @abstractmethod
    def update_real_time_model(self, **kwargs):

        """
        Update the real-time model (advance timesteps) with necessary parameters in kwargs.

        Arguments:
            kwargs: necessary profiles to update the underlying model. {stat_name: [...]}

        Returns:
            None
        """

        pass

    @abstractmethod
    def compute_day_ahead_bids(self, date, hour, **kwargs):

        """
        Solve the model to bid/self-schedule into the day-ahead market. After solving,
        record the schedule from the solve.

        Arguments:

            date: current simulation date

            hour: current simulation hour

            **kwargs: other information to record

        Returns:
            None
        """

        pass

    @abstractmethod
    def compute_real_time_bids(self, date, hour, **kwargs):

        """
        Solve the model to bid/self-schedule into the real-time market. After solving,
        record the schedule from the solve.

        Arguments:

            date: current simulation date

            hour: current simulation hour

            **kwargs: other information to record

        Returns:
            None
        """

        pass

    @abstractmethod
    def write_results(self, path):

        """
        This methods writes the saved results into an csv file.

        Arguments:
            path: the path to write the results.

        Return:
            None
        """

        pass

    @abstractmethod
    def formulate_DA_bidding_problem(self):

        """
        Formulate the day-ahead bidding optimization problem by adding necessary
        parameters, constraints, and objective function.

        Arguments:
            None

        Returns:
            None
        """

        pass

    @abstractmethod
    def formulate_RT_bidding_problem(self):

        """
        Formulate the real-time bidding optimization problem by adding necessary
        parameters, constraints, and objective function.

        Arguments:
            None

        Returns:
            None
        """

        pass

    @abstractmethod
    def record_bids(self, bids, model, date, hour):

        """
        This function records the bids (schedule) and the details in the
        underlying bidding model.

        Arguments:
            bids: the obtained bids for this date

            model: the model we obtained bids from

            date: the date we bid into

            hour: the hour we bid into

        Returns:
            None

        """

        pass

    @property
    @abstractmethod
    def generator(self):
        return "AbstractGenerator"

    def _check_inputs(self):

        """
        Check if the inputs to construct the tracker is valid. If not raise errors.
        """

        self._check_bidding_model_object()
        self._check_n_scenario()
        self._check_solver()

    def _check_bidding_model_object(self):

        """
        Check if tracking model object has the necessary methods and attributes.
        """

        method_list = ["populate_model", "update_model"]
        attr_list = ["power_output", "total_cost", "model_data"]
        msg = "Bidding model object does not have a "

        for m in method_list:
            obtained_m = getattr(self.bidding_model_object, m, None)
            if obtained_m is None:
                raise AttributeError(
                    msg
                    + f"{m}() method. The bidder object needs the users to implement this method in their model object."
                )

        for attr in attr_list:
            obtained_attr = getattr(self.bidding_model_object, attr, None)
            if obtained_attr is None:
                raise AttributeError(
                    msg
                    + f"{attr} property. The bidder object needs the users to specify this property in their model object."
                )

    def _check_n_scenario(self):

        """
        Check if the number of LMP scenarios is an integer and greater than 0.
        """

        # check if it is an integer
        if not isinstance(self.n_scenario, int):
            raise TypeError(
                f"The number of LMP scenarios should be an integer, but a {type(self.n_scenario).__name__} was given."
            )

        if self.n_scenario <= 0:
            raise ValueError(
                f"The number of LMP scenarios should be greater than zero, but {self.n_scenario} was given."
            )

    def _check_solver(self):

        """
        Check if provides solver is a valid Pyomo solver object.
        """

        if not isinstance(self.solver, OptSolver):
            raise TypeError(
                f"The provided solver {self.solver} is not a valid Pyomo solver."
            )


class StochasticProgramBidder(AbstractBidder):

    """
    Template class for bidders that use scenario-based stochastic programs.
    """

    def __init__(
        self,
        bidding_model_object,
        day_ahead_horizon,
        real_time_horizon,
        n_scenario,
        solver,
        forecaster,
        real_time_underbid_penalty,
    ):

        """
        Initializes the stochastic bidder object.

        Arguments:
            bidding_model_object: the model object for bidding

            day_ahead_horizon: number of time periods in the day-ahead bidding problem

            real_time_horizon: number of time periods in the real-time bidding problem

            n_scenario: number of uncertain LMP scenarios

            solver: a Pyomo mathematical programming solver object

            forecaster: an initialized LMP forecaster object

            real_time_underbid_penalty: penalty for RT power bid that's less than DA power bid, non-negative

        Returns:
            None
        """

        self.bidding_model_object = bidding_model_object
        self.day_ahead_horizon = day_ahead_horizon
        self.real_time_horizon = real_time_horizon
        self.n_scenario = n_scenario
        self.solver = solver
        self.forecaster = forecaster
        self.real_time_underbid_penalty = real_time_underbid_penalty

        self._check_inputs()

        self.generator = self.bidding_model_object.model_data.gen_name

        # day-ahead model
        self.day_ahead_model = self.formulate_DA_bidding_problem()
        self.real_time_model = self.formulate_RT_bidding_problem()

        # declare a list to store results
        self.bids_result_list = []

    def _set_up_bidding_problem(self, horizon):
        """
        Set up the base stochastic programming bidding problems.

        Arguments:
            horizon: number of time periods in the bidding problem

        Returns:
            pyomo.core.base.PyomoModel.ConcreteModel: base bidding model

        """

        model = pyo.ConcreteModel()

        model.SCENARIOS = pyo.Set(initialize=range(self.n_scenario))

        model.fs = pyo.Block(model.SCENARIOS)
        for i in model.SCENARIOS:
            self.bidding_model_object.populate_model(model.fs[i], horizon)

        self._save_power_outputs(model)

        self._add_bidding_params(model)
        self._add_bidding_vars(model)
        self._add_bidding_objective(model)

        return model

    def formulate_DA_bidding_problem(self):

        """
        Set up the day-ahead stochastic programming bidding problems.

        Returns:
            pyomo.core.base.PyomoModel.ConcreteModel: base bidding model

        """

        model = self._set_up_bidding_problem(self.day_ahead_horizon)
        self._add_DA_bidding_constraints(model)

        # do not relax the DA offering UB
        for i in model.SCENARIOS:
            model.fs[i].real_time_underbid_power.fix(0)

        return model

    def formulate_RT_bidding_problem(self):

        """
        Set up the real-time stochastic programming bidding problems.

        Returns:
            pyomo.core.base.PyomoModel.ConcreteModel: base bidding model

        """

        model = self._set_up_bidding_problem(self.real_time_horizon)
        self._add_RT_bidding_constraints(model)

        # relax the DA offering UB
        for i in model.SCENARIOS:
            model.fs[i].real_time_underbid_power.unfix()

        return model

    def _save_power_outputs(self, model):

        """
        Create references of the power output variable in each price scenario
        block.

        Arguments:
            None

        Returns:
            None
        """

        for i in model.SCENARIOS:
            # get the power output
            power_output_name = self.bidding_model_object.power_output
            model.fs[i].power_output_ref = pyo.Reference(
                getattr(model.fs[i], power_output_name)
            )

        return

    def _add_bidding_params(self, model):

        """
        Add necessary bidding parameters to the model, i.e., market energy price.

        Arguments:
            model: bidding model

        Returns:
            None
        """

        for i in model.SCENARIOS:
            time_index = model.fs[i].power_output_ref.index_set()
            model.fs[i].day_ahead_energy_price = pyo.Param(
                time_index, initialize=0, mutable=True
            )
            model.fs[i].real_time_energy_price = pyo.Param(
                time_index, initialize=0, mutable=True
            )
            model.fs[i].real_time_underbid_penalty = pyo.Param(
                initialize=self.real_time_underbid_penalty, mutable=True
            )

        return

    def _add_bidding_vars(self, model):

        """
        Add necessary bidding parameters to the model, i.e., market energy price.

        Arguments:
            model: bidding model

        Returns:
            None
        """

        def relaxed_day_ahead_power_ub_rule(fs, t):
            return (
                fs.power_output_ref[t] + fs.real_time_underbid_power[t]
                >= fs.day_ahead_power[t]
            )

        for i in model.SCENARIOS:
            time_index = model.fs[i].power_output_ref.index_set()
            model.fs[i].day_ahead_power = pyo.Var(
                time_index, initialize=0, within=pyo.NonNegativeReals
            )

            model.fs[i].real_time_underbid_power = pyo.Var(
                time_index, initialize=0, within=pyo.NonNegativeReals
            )

            model.fs[i].day_ahead_power_ub = pyo.Constraint(
                time_index, rule=relaxed_day_ahead_power_ub_rule
            )

        return

    def _add_bidding_objective(self, model):

        """
        Add objective function to the model, i.e., maximizing the expected profit
        of the energy system.

        Arguments:
            model: bidding model

        Returns:
            None
        """

        # declare an empty objective
        model.obj = pyo.Objective(expr=0, sense=pyo.maximize)

        for k in model.SCENARIOS:
            time_index = model.fs[k].power_output_ref.index_set()

            # currently .total_cost is a tuple of 2 items
            # the first item is the name of the cost expression
            # the second item is the weight for the cost
            cost_name = self.bidding_model_object.total_cost[0]
            cost = getattr(model.fs[k], cost_name)
            weight = self.bidding_model_object.total_cost[1]

            for t in time_index:
                model.obj.expr += (
                    model.fs[k].day_ahead_energy_price[t]
                    * model.fs[k].day_ahead_power[t]
                    + model.fs[k].real_time_energy_price[t]
                    * (model.fs[k].power_output_ref[t] - model.fs[k].day_ahead_power[t])
                    - weight * cost[t]
                    - model.fs[k].real_time_underbid_penalty
                    * model.fs[k].real_time_underbid_power[t]
                )

        return

    def _compute_bids(
        self,
        day_ahead_price,
        real_time_energy_price,
        date,
        hour,
        model,
        power_var_name,
        energy_price_param_name,
        market,
    ):

        """
        Solve the model to bid into the markets. After solving, record the bids from the solve.

        Arguments:

            day_ahead_price: day-ahead price forecasts needed to solve the bidding problem

            real_time_energy_price: real-time price forecasts needed to solve the bidding problem

            date: current simulation date

            hour: current simulation hour

            model: bidding model

            power_var_name: the name of the power output (str)

            energy_price_param_name: the name of the energy price forecast params (str)

            market: the market name (str), e.g., Day-ahead, real-time

        Returns:
            dict: the obtained bids
        """

        # update the price forecasts
        self._pass_price_forecasts(model, day_ahead_price, real_time_energy_price)

        self.solver.solve(model, tee=True)

        bids = self._assemble_bids(
            model,
            power_var_name=power_var_name,
            energy_price_param_name=energy_price_param_name,
            hour=hour,
        )

        self.record_bids(bids, model=model, date=date, hour=hour, market=market)

        return bids

    def compute_day_ahead_bids(self, date, hour=0):

        """
        Solve the model to bid into the day-ahead market. After solving, record
        the bids from the solve.

        Arguments:

            date: current simulation date

            hour: current simulation hour

        Returns:
            dict: the obtained bids
        """

        (
            day_ahead_price,
            real_time_energy_price,
        ) = self.forecaster.forecast_day_ahead_and_real_time_prices(
            date=date,
            hour=hour,
            bus=self.bidding_model_object.model_data.bus,
            horizon=self.day_ahead_horizon,
            n_samples=self.n_scenario,
        )

        return self._compute_bids(
            day_ahead_price=day_ahead_price,
            real_time_energy_price=real_time_energy_price,
            date=date,
            hour=hour,
            model=self.day_ahead_model,
            power_var_name="day_ahead_power",
            energy_price_param_name="day_ahead_energy_price",
            market="Day-ahead",
        )

    def compute_real_time_bids(
        self, date, hour, realized_day_ahead_prices, realized_day_ahead_dispatches
    ):

        """
        Solve the model to bid into the real-time market. After solving, record
        the bids from the solve.

        Arguments:

            date: current simulation date

            hour: current simulation hour

        Returns:
            dict: the obtained bids
        """

        real_time_energy_price = self.forecaster.forecast_real_time_prices(
            date=date,
            hour=hour,
            bus=self.bidding_model_object.model_data.bus,
            horizon=self.real_time_horizon,
            n_samples=self.n_scenario,
        )

        self._pass_realized_day_ahead_dispatches(realized_day_ahead_dispatches, hour)
        self._pass_realized_day_ahead_prices(realized_day_ahead_prices, date, hour)

        return self._compute_bids(
            day_ahead_price=None,
            real_time_energy_price=real_time_energy_price,
            date=date,
            hour=hour,
            model=self.real_time_model,
            power_var_name="power_output_ref",
            energy_price_param_name="real_time_energy_price",
            market="Real-time",
        )

    def _pass_realized_day_ahead_prices(self, realized_day_ahead_prices, date, hour):

        """
        Pass the realized day-ahead prices into model parameters.

        Arguments:
            realized_day_ahead_prices: realized day-ahead prices

            date: current simulation date

            hour: current simulation hour

        Returns:
            None
        """

        time_index = self.real_time_model.fs[0].day_ahead_energy_price.index_set()

        # forecast the day-ahead price, if not enough realized data
        if len(realized_day_ahead_prices[hour:]) < len(time_index):
            forecasts = self.forecaster.forecast_day_ahead_prices(
                date=date + datetime.timedelta(days=1),
                hour=0,
                bus=self.bidding_model_object.model_data.bus,
                horizon=self.day_ahead_horizon,
                n_samples=self.n_scenario,
            )

        for s in self.real_time_model.SCENARIOS:
            for t in time_index:
                try:
                    price = realized_day_ahead_prices[t + hour]
                except IndexError as ex:
                    self.real_time_model.fs[s].day_ahead_energy_price[t] = forecasts[s][
                        (t + hour) - 24
                    ]
                else:
                    self.real_time_model.fs[s].day_ahead_energy_price[t] = price

    def _pass_realized_day_ahead_dispatches(self, realized_day_ahead_dispatches, hour):

        """
        Pass the realized day-ahead dispatches into model and fix the corresponding variables.

        Arguments:
            realized_day_ahead_dispatches: realized day-ahead dispatches

            hour: current simulation hour

        Returns:
            None
        """

        time_index = self.real_time_model.fs[0].day_ahead_power.index_set()
        for s in self.real_time_model.SCENARIOS:
            for t in time_index:
                try:
                    dispatch = realized_day_ahead_dispatches[t + hour]
                except IndexError as ex:
                    self.real_time_model.fs[s].day_ahead_power[t].unfix()
                    # unrelax the DA offering UB
                    self.real_time_model.fs[s].real_time_underbid_power[t].fix(0)
                else:
                    self.real_time_model.fs[s].day_ahead_power[t].fix(dispatch)
                    # relax the DA offering UB
                    self.real_time_model.fs[s].real_time_underbid_power[t].unfix()

    def update_day_ahead_model(self, **kwargs):

        """
        This method updates the parameters in the day-ahead model based on the implemented profiles.

        Arguments:
            kwargs: the newly implemented stats. {stat_name: [...]}

        Returns:
            None
        """

        self._update_model(self.day_ahead_model, **kwargs)

    def update_real_time_model(self, **kwargs):

        """
        This method updates the parameters in the real-time model based on the implemented profiles.

        Arguments:
            kwargs: the newly implemented stats. {stat_name: [...]}

        Returns:
            None
        """

        self._update_model(self.real_time_model, **kwargs)

    def _update_model(self, model, **kwargs):

        """
        Update the flowsheets in all the price scenario blocks to advance time
        step.

        Arguments:

            model: bidding model

            kwargs: necessary profiles to update the underlying model. {stat_name: [...]}

        Returns:
            None
        """

        for i in model.SCENARIOS:
            self.bidding_model_object.update_model(b=model.fs[i], **kwargs)

        return

    def record_bids(self, bids, model, date, hour, market):

        """
        This function records the bids and the details in the underlying bidding model.

        Arguments:
            bids: the obtained bids for this date.

            model: bidding model

            date: the date we bid into

            hour: the hour we bid into

        Returns:
            None

        """

        # record bids
        self._record_bids(bids, date, hour, Market=market)

        # record the details of bidding model
        for i in model.SCENARIOS:
            self.bidding_model_object.record_results(
                model.fs[i], date=date, hour=hour, Scenario=i, Market=market
            )

        return

    def _pass_price_forecasts(self, model, day_ahead_price, real_time_energy_price):

        """
        Pass the price forecasts into model parameters.

        Arguments:
            day_ahead_price: day-ahead price forecasts needed to solve the bidding problem

            real_time_energy_price: real-time price forecasts needed to solve the bidding problem

        Returns:
            None
        """

        for i in model.SCENARIOS:

            time_index = model.fs[i].real_time_energy_price.index_set()

            if day_ahead_price is not None:
                for t, p in zip(time_index, day_ahead_price[i]):
                    model.fs[i].day_ahead_energy_price[t] = p

            for t, p in zip(time_index, real_time_energy_price[i]):
                model.fs[i].real_time_energy_price[t] = p

        return

    def write_results(self, path):
        """
        This methods writes the saved operation stats into an csv file.

        Arguments:
            path: the path to write the results.

        Return:
            None
        """

        print("")
        print("Saving bidding results to disk...")
        pd.concat(self.bids_result_list).to_csv(
            os.path.join(path, "bidder_detail.csv"), index=False
        )
        self.bidding_model_object.write_results(
            path=os.path.join(path, "bidding_model_detail.csv")
        )

        return

    @property
    def generator(self):
        return self._generator

    @generator.setter
    def generator(self, name):
        self._generator = name


class SelfScheduler(StochasticProgramBidder):

    """
    Wrap a model object to self schedule into the market using stochastic programming.
    """

    def __init__(
        self,
        bidding_model_object,
        day_ahead_horizon,
        real_time_horizon,
        n_scenario,
        solver,
        forecaster,
        real_time_underbid_penalty=10000,
        fixed_to_schedule=False,
    ):
        """
        Initializes the stochastic self-scheduler object.

        Arguments:
            bidding_model_object: the model object for bidding

            day_ahead_horizon: number of time periods in the day-ahead bidding problem

            real_time_horizon: number of time periods in the real-time bidding problem

            n_scenario: number of uncertain LMP scenarios

            solver: a Pyomo mathematical programming solver object

            forecaster: an initialized LMP forecaster object

            fixed_to_schedule: If True, forece market simulator to give the same schedule.

        Returns:
            None
        """

        super().__init__(
            bidding_model_object,
            day_ahead_horizon,
            real_time_horizon,
            n_scenario,
            solver,
            forecaster,
            real_time_underbid_penalty,
        )
        self.fixed_to_schedule = fixed_to_schedule

    def _add_DA_bidding_constraints(self, model):

        """
        Add bidding constraints to the model, i.e., power outputs in the first
        stage need to be the same across all the scenarios.

        Arguments:
            model: bidding model

        Returns:
            None
        """

        # nonanticipativity constraints
        def day_ahead_bidding_constraints_rule(model, s1, s2, t):
            if s1 == s2:
                return pyo.Constraint.Skip
            return model.fs[s1].day_ahead_power[t] == model.fs[s2].day_ahead_power[t]

        time_index = model.fs[model.SCENARIOS.first()].power_output_ref.index_set()

        model.day_ahead_bidding_constraints = pyo.Constraint(
            model.SCENARIOS,
            model.SCENARIOS,
            time_index,
            rule=day_ahead_bidding_constraints_rule,
        )

        return

    def _add_RT_bidding_constraints(self, model):

        """
        Add bidding constraints to the model, i.e., power outputs in the first
        stage need to be the same across all the scenarios.

        Arguments:
            model: bidding model

        Returns:
            None
        """

        # nonanticipativity constraints
        def real_time_bidding_constraints_rule(model, s1, s2, t):
            if s1 == s2:
                return pyo.Constraint.Skip
            return model.fs[s1].power_output_ref[t] == model.fs[s2].power_output_ref[t]

        time_index = model.fs[model.SCENARIOS.first()].power_output_ref.index_set()

        model.real_time_bidding_constraints = pyo.Constraint(
            model.SCENARIOS,
            model.SCENARIOS,
            time_index,
            rule=real_time_bidding_constraints_rule,
        )

        return

    def _assemble_bids(self, model, power_var_name, energy_price_param_name, hour):

        """
        This methods extract the bids out of the stochastic programming model and
        organize them into self-schedule bids.

        For thermal generators, startup times and costs are set to 0. And the bid price for power outside of p_min and p_max are 0.

        Arguments:
            model: bidding model

            power_var_name: the name of the power output (str)

            energy_price_param_name: the name of the energy price forecast params (str)

            hour: current simulation hour

        Returns:
            dict: the bid we computed.
        """

        bids = {}
        is_thermal = self.bidding_model_object.model_data.generator_type == "thermal"

        power_output_var = getattr(model.fs[0], power_var_name)
        time_index = power_output_var.index_set()

        for t_idx in time_index:

            t = t_idx + hour

            bids[t] = {}

            bids[t][self.generator] = {
                "p_max": round(pyo.value(power_output_var[t_idx]), 4),
                "p_min": self.bidding_model_object.model_data.p_min,
            }

            if self.fixed_to_schedule:
                bids[t][self.generator]["p_min"] = bids[t][self.generator]["p_max"]
                bids[t][self.generator]["fixed_commitment"] = (
                    1 if bids[t][self.generator]["p_min"] > 0 else 0
                )

                if is_thermal:
                    bids[t][self.generator]["min_up_time"] = 0
                    bids[t][self.generator]["min_down_time"] = 0

                    bids[t][self.generator]["startup_fuel"] = [
                        (bids[t][self.generator]["min_down_time"], 0)
                    ]
                    bids[t][self.generator]["startup_cost"] = [
                        (bids[t][self.generator]["min_down_time"], 0)
                    ]

            if is_thermal:

                bids[t][self.generator]["p_cost"] = [
                    (bids[t][self.generator]["p_min"], 0),
                    (bids[t][self.generator]["p_max"], 0),
                ]

                bids[t][self.generator]["startup_capacity"] = bids[t][self.generator][
                    "p_min"
                ]
                bids[t][self.generator]["shutdown_capacity"] = bids[t][self.generator][
                    "p_min"
                ]

        return bids

    def _record_bids(self, bids, date, hour, **kwargs):

        """
        This function records the bids (schedule) we computed for the given date into a
        DataFrame and temporarily stores the DataFrame in an instance attribute
        list called bids_result_list.

        Arguments:
            bids: the obtained bids (schedule) for this date.

            date: the date we bid into

            hour: the hour we bid into

        Returns:
            None

        """

        df_list = []
        for t in bids:
            for g in bids[t]:

                result_dict = {}
                result_dict["Generator"] = g
                result_dict["Date"] = date
                if hour is not None:
                    result_dict["Hour"] = hour

                result_dict["Horizon"] = t
                result_dict["Bid Power [MW]"] = bids[t][g].get("p_max")
                result_dict["Bid Min Power [MW]"] = bids[t][g].get("p_min")

                for k, v in kwargs.items():
                    result_dict[k] = v

                result_df = pd.DataFrame.from_dict(result_dict, orient="index")
                df_list.append(result_df.T)

        # save the result to object property
        # wait to be written when simulation ends
        self.bids_result_list.append(pd.concat(df_list))


class Bidder(StochasticProgramBidder):

    """
    Wrap a model object to bid into the market using stochastic programming.
    """

    def __init__(
        self,
        bidding_model_object,
        day_ahead_horizon,
        real_time_horizon,
        n_scenario,
        solver,
        forecaster,
        real_time_underbid_penalty=10000,
    ):

        """
        Initializes the bidder object.

        Arguments:
            bidding_model_object: the model object for bidding

            day_ahead_horizon: number of time periods in the day-ahead bidding problem

            real_time_horizon: number of time periods in the real-time bidding problem

            n_scenario: number of uncertain LMP scenarios

            solver: a Pyomo mathematical programming solver object

            forecaster: an initialized LMP forecaster object

            real_time_underbid_penalty: penalty for RT power bid that's less than DA power bid, non-negative

        Returns:
            None
        """

        super().__init__(
            bidding_model_object,
            day_ahead_horizon,
            real_time_horizon,
            n_scenario,
            solver,
            forecaster,
            real_time_underbid_penalty,
        )

    def _add_DA_bidding_constraints(self, model):

        """
        Add bidding constraints to the model, i.e., the bid curves need to be
        nondecreasing.

        Arguments:
            model: bidding model

        Returns:
            None
        """

        def day_ahead_bidding_constraints_rule(model, s1, s2, t):
            if s1 == s2:
                return pyo.Constraint.Skip
            return (
                model.fs[s1].day_ahead_power[t] - model.fs[s2].day_ahead_power[t]
            ) * (
                model.fs[s1].day_ahead_energy_price[t]
                - model.fs[s2].day_ahead_energy_price[t]
            ) >= 0

        time_index = model.fs[model.SCENARIOS.first()].power_output_ref.index_set()

        model.day_ahead_bidding_constraints = pyo.Constraint(
            model.SCENARIOS,
            model.SCENARIOS,
            time_index,
            rule=day_ahead_bidding_constraints_rule,
        )

        return

    def _add_RT_bidding_constraints(self, model):

        """
        Add bidding constraints to the model, i.e., the bid curves need to be
        nondecreasing.

        Arguments:
            model: bidding model

        Returns:
            None
        """

        def real_time_bidding_constraints_rule(model, s1, s2, t):
            if s1 == s2:
                return pyo.Constraint.Skip
            return (
                model.fs[s1].power_output_ref[t] - model.fs[s2].power_output_ref[t]
            ) * (
                model.fs[s1].real_time_energy_price[t]
                - model.fs[s2].real_time_energy_price[t]
            ) >= 0

        time_index = model.fs[model.SCENARIOS.first()].power_output_ref.index_set()

        model.real_time_bidding_constraints = pyo.Constraint(
            model.SCENARIOS,
            model.SCENARIOS,
            time_index,
            rule=real_time_bidding_constraints_rule,
        )

        return

    def _assemble_bids(self, model, power_var_name, energy_price_param_name, hour):

        """
        This methods extract the bids out of the stochastic programming model and
        organize them into ( MWh, $ ) pairs.

        Arguments:
            model: bidding model

            power_var_name: the name of the power output (str)

            energy_price_param_name: the name of the energy price forecast params (str)

            hour: current simulation hour

        Returns:
            bids: the bid we computed.
        """

        bids = {}
        gen = self.generator

        for i in model.SCENARIOS:

            power_output_var = getattr(model.fs[i], power_var_name)
            energy_price_param = getattr(model.fs[i], energy_price_param_name)
            time_index = power_output_var.index_set()

            for t in time_index:

                if t not in bids:
                    bids[t] = {}
                if gen not in bids[t]:
                    bids[t][gen] = {}

                power = round(pyo.value(power_output_var[t]), 2)
                marginal_cost = round(pyo.value(energy_price_param[t]), 2)

                # if power lower than pmin, e.g., power = 0, we need to skip this
                # solution, because Prescient is not expecting any power output lower
                # than pmin in the bids
                if power < self.bidding_model_object.model_data.p_min:
                    continue
                elif power in bids[t][gen]:
                    bids[t][gen][power] = min(bids[t][gen][power], marginal_cost)
                else:
                    bids[t][gen][power] = marginal_cost

        for t in time_index:

            # make sure the orignal points in the bids
            for power, marginal_cost in self.bidding_model_object.model_data.p_cost:
                if round(power, 2) not in bids[t][gen]:
                    bids[t][gen][power] = marginal_cost

            pmin = self.bidding_model_object.model_data.p_min

            # sort the curves by power
            bids[t][gen] = dict(sorted(bids[t][gen].items()))

            # make sure the curve is nondecreasing
            pre_power = pmin
            for power, marginal_cost in bids[t][gen].items():

                # ignore pmin, because min load cost is special
                if pre_power == pmin:
                    pre_power = power
                    continue
                bids[t][gen][power] = max(bids[t][gen][power], bids[t][gen][pre_power])
                pre_power = power

            # calculate the actual cost
            bids[t][gen] = convert_marginal_costs_to_actual_costs(
                list(bids[t][gen].items())
            )

        # check if bids are convex
        for t in bids:
            for gen in bids[t]:
                temp_curve = {
                    "data_type": "cost_curve",
                    "cost_curve_type": "piecewise",
                    "values": bids[t][gen],
                }

                try:
                    tx_utils.validate_and_clean_cost_curve(
                        curve=temp_curve,
                        curve_type="cost_curve",
                        p_min=min([p[0] for p in bids[t][gen]]),
                        p_max=max([p[0] for p in bids[t][gen]]),
                        gen_name=gen,
                        t=t,
                    )
                except NameError as ex:
                    raise RuntimeError(
                        "'egret' must be installed to use this functionality"
                    )

        # create full bids: this includes info in addition to costs
        full_bids = {}

        for t_idx in bids:

            t = t_idx + hour

            full_bids[t] = {}
            for gen in bids[t_idx]:
                full_bids[t][gen] = {}
                full_bids[t][gen]["p_cost"] = bids[t_idx][gen]
                full_bids[t][gen]["p_min"] = min([p[0] for p in bids[t_idx][gen]])
                full_bids[t][gen]["p_max"] = max([p[0] for p in bids[t_idx][gen]])
                full_bids[t][gen]["p_min_agc"] = min([p[0] for p in bids[t_idx][gen]])
                full_bids[t][gen]["p_max_agc"] = max([p[0] for p in bids[t_idx][gen]])
                full_bids[t][gen]["startup_capacity"] = full_bids[t][gen]["p_min"]
                full_bids[t][gen]["shutdown_capacity"] = full_bids[t][gen]["p_min"]

                fixed_commitment = getattr(
                    self.bidding_model_object.model_data, "fixed_commitment", None
                )
                if fixed_commitment is not None:
                    full_bids[t][gen]["fixed_commitment"] = fixed_commitment

        return full_bids

    def _record_bids(self, bids, date, hour, **kwargs):

        """
        This method records the bids we computed for the given date into a
        DataFrame. This DataFrame has the following columns: gen, date, hour,
        power 1, ..., power n, price 1, ..., price n. And concatenate the
        DataFrame into a class property 'bids_result_list'. The methods then
        temporarily stores the DataFrame in an instance attribute list called
        bids_result_list.

        Arguments:
            bids: the obtained bids for this date.

            date: the date we bid into

            hour: the hour we bid into

        Returns:
            None

        """

        df_list = []
        for t in bids:
            for gen in bids[t]:

                result_dict = {}
                result_dict["Generator"] = gen
                result_dict["Date"] = date
                result_dict["Hour"] = t

                for k, v in kwargs.items():
                    result_dict[k] = v

                pair_cnt = len(bids[t][gen]["p_cost"])

                for idx, (power, cost) in enumerate(bids[t][gen]["p_cost"]):
                    result_dict[f"Power {idx} [MW]"] = power
                    result_dict[f"Cost {idx} [$]"] = cost

                # place holder, in case different len of bids
                while pair_cnt < self.n_scenario:
                    result_dict[f"Power {pair_cnt} [MW]"] = None
                    result_dict[f"Cost {pair_cnt} [$]"] = None

                    pair_cnt += 1

                result_df = pd.DataFrame.from_dict(result_dict, orient="index")
                df_list.append(result_df.T)

        # save the result to object property
        # wait to be written when simulation ends
        self.bids_result_list.append(pd.concat(df_list))

        return
