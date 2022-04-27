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

from pyomo.common.dependencies import attempt_import

egret, egret_avail = attempt_import("egret")
if egret_avail:
    from egret.model_library.transmission import tx_utils


class AbstractBidder(ABC):

    """
    The abstract class for all the bidder and self-schedulers.
    """

    @abstractmethod
    def update_model(self, **kwargs):

        """
        Update the flowsheets advance timesteps with necessary parameters in kwargs.

        Arguments:
            kwargs: necessary profiles to update the underlying model. {stat_name: [...]}

        Returns:
            None
        """

        pass

    @abstractmethod
    def compute_bids(self, date, hour, **kwargs):

        """
        Solve the model to bid/self-schedule into the markets. After solving, record
        the schedule from the solve.

        Arguments:

            date: current simulation date

            hour: current simulation hour

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
    def formulate_bidding_problem(self):

        """
        Formulate the bidding optimization problem by adding necessary
        parameters, constraints, and objective function.

        Arguments:
            None

        Returns:
            None
        """

        pass

    @abstractmethod
    def record_bids(self, bids, date, hour):

        """
        This function records the bids (schedule) and the details in the
        underlying bidding model.

        Arguments:
            bids: the obtained bids for this date.

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


class SelfScheduler(AbstractBidder):

    """
    Wrap a model object to self schedule into the market using stochastic programming.
    """

    def __init__(
        self,
        bidding_model_object,
        n_scenario,
        horizon,
        solver,
        forecaster,
        fixed_to_schedule=False,
    ):
        self.bidding_model_object = bidding_model_object
        self.n_scenario = n_scenario
        self.horizon = horizon
        self.solver = solver
        self.forecaster = forecaster
        self.generator = self.bidding_model_object.model_data.gen_name
        self.fixed_to_schedule = fixed_to_schedule

        self._check_inputs()

        # add flowsheets to model
        self.model = pyo.ConcreteModel()

        # declare scenario set
        self.model.SCENARIOS = pyo.Set(initialize=range(self.n_scenario))

        # populate scenario blocks
        self.model.fs = pyo.Block(self.model.SCENARIOS)
        for i in self.model.SCENARIOS:
            self.bidding_model_object.populate_model(self.model.fs[i])

        # save power output variable in the model object
        self._save_power_outputs()

        self.formulate_bidding_problem()

        # declare a list to store results
        self.bids_result_list = []

    def _save_power_outputs(self):

        """
        Create references of the power output variable in each price scenario
        block.

        Arguments:
            None

        Returns:
            None
        """

        for i in self.model.SCENARIOS:
            # get the power output
            power_output_name = self.bidding_model_object.power_output
            self.model.fs[i].power_output_ref = pyo.Reference(
                getattr(self.model.fs[i], power_output_name)
            )

        return

    def formulate_bidding_problem(self):

        """
        Formulate the bidding optimization problem by adding necessary
        parameters, constraints, and objective function.

        Arguments:
            None

        Returns:
            None
        """

        # add bidding params
        for i in self.model.SCENARIOS:
            time_index = self.model.fs[i].power_output_ref.index_set()
            self.model.fs[i].energy_price = pyo.Param(
                time_index, initialize=0, mutable=True
            )

        # nonanticipativity constraints
        def nonanticipativity_constraint_rule(model, s1, s2, t):
            if s1 == s2:
                return pyo.Constraint.Skip
            return model.fs[s1].power_output_ref[t] == model.fs[s2].power_output_ref[t]

        time_index = self.model.fs[
            self.model.SCENARIOS.first()
        ].power_output_ref.index_set()

        self.model.NonanticipativityConstraints = pyo.Constraint(
            self.model.SCENARIOS,
            self.model.SCENARIOS,
            time_index,
            rule=nonanticipativity_constraint_rule,
        )

        # declare an empty objective
        self.model.obj = pyo.Objective(expr=0, sense=pyo.maximize)
        for k in self.model.SCENARIOS:
            time_index = self.model.fs[k].power_output_ref.index_set()

            cost_name = self.bidding_model_object.total_cost[0]
            cost = getattr(self.model.fs[k], cost_name)
            weight = self.bidding_model_object.total_cost[1]

            for t in time_index:
                self.model.obj.expr += (
                    self.model.fs[k].power_output_ref[t]
                    * self.model.fs[k].energy_price[t]
                    - weight * cost[t]
                )

        return

    def compute_bids(self, date, hour=None, **kwargs):

        """
        Solve the model to self-schedule into the markets. After solving, record
        the schedule from the solve.

        Arguments:

            date: current simulation date

            hour: current simulation hour

        Returns:
            None
        """

        price_forecasts = self.forecaster.forecast(date=date, hour=hour, **kwargs)

        # update the price forecasts
        self._pass_price_forecasts(price_forecasts)
        self.solver.solve(self.model, tee=True)
        bids = self._assemble_bids()
        self.record_bids(bids, date=date, hour=hour)

        return bids

    def _assemble_bids(self):

        """
        This methods extract the bids out of the stochastic programming model and
        organize them.

        Arguments:

        Returns:
            bids: the bid we computed.
        """

        bids = {}

        for t in range(self.horizon):

            bids[t] = {}

            bids[t][self.generator] = {
                "p_max": round(pyo.value(self.model.fs[0].power_output_ref[t]), 4),
                "p_min": self.bidding_model_object.model_data.p_min,
            }

            if self.fixed_to_schedule:
                bids[t][self.generator]["p_min"] = bids[t][self.generator]["p_max"]
                bids[t][self.generator]["min_up_time"] = 0
                bids[t][self.generator]["min_down_time"] = 0
                bids[t][self.generator]["fixed_commitment"] = (
                    1 if bids[t][self.generator]["p_min"] > 0 else 0
                )
                bids[t][self.generator]["startup_fuel"] = [
                    (bids[t][self.generator]["min_down_time"], 0)
                ]
                bids[t][self.generator]["startup_cost"] = [
                    (bids[t][self.generator]["min_down_time"], 0)
                ]

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

    def _pass_price_forecasts(self, price_forecasts):

        """
        Pass the price forecasts into model parameters.

        Arguments:
            price_forecasts: price forecasts needed to solve the bidding problem. {LMP scenario: [forecast timeseries] }

        Returns:
            None
        """

        for i in self.model.SCENARIOS:
            time_index = self.model.fs[i].energy_price.index_set()
            for t, p in zip(time_index, price_forecasts[i]):
                self.model.fs[i].energy_price[t] = p

        return

    def _record_bids(self, bids, date, hour):

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

                result_df = pd.DataFrame.from_dict(result_dict, orient="index")
                df_list.append(result_df.T)

        # save the result to object property
        # wait to be written when simulation ends
        self.bids_result_list.append(pd.concat(df_list))

    def record_bids(self, bids, date, hour):

        """
        This function records the bids (schedule) and the details in the
        underlying bidding model.

        Arguments:
            bids: the obtained bids for this date.

            date: the date we bid into

            hour: the hour we bid into

        Returns:
            None

        """

        # record bids
        self._record_bids(bids, date, hour)

        # record the details of bidding model
        for i in self.model.SCENARIOS:
            self.bidding_model_object.record_results(
                self.model.fs[i], date=date, hour=hour, Scenario=i
            )

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

    def update_model(self, **kwargs):

        """
        Update the flowsheets in all the price scenario blocks to advance time
        step.

        Arguments:
            kwargs: necessary profiles to update the underlying model. {stat_name: [...]}

        Returns:
            None
        """

        for i in self.model.SCENARIOS:
            self.bidding_model_object.update_model(b=self.model.fs[i], **kwargs)

    @property
    def generator(self):
        return self._generator

    @generator.setter
    def generator(self, name):
        self._generator = name


class Bidder(AbstractBidder):

    """
    Wrap a model object to bid into the market using stochastic programming.
    """

    def __init__(self, bidding_model_object, n_scenario, solver, forecaster):

        """
        Initializes the bidder object.

        Arguments:
            bidding_model_object: the model object for tracking for bidding

            n_scenario: number of LMP scenarios

            solver: a Pyomo mathematical programming solver object

            forecaster: an initialized LMP forecaster object


        Returns:
            None
        """

        # copy the inputs
        self.bidding_model_object = bidding_model_object
        self.n_scenario = n_scenario
        self.solver = solver
        self.forecaster = forecaster

        self._check_inputs()

        # get the generator name
        self.generator = self.bidding_model_object.model_data.gen_name

        # add flowsheets to model
        self.model = pyo.ConcreteModel()

        # declare scenario set
        self.model.SCENARIOS = pyo.Set(initialize=range(self.n_scenario))

        # populate scenario blocks
        self.model.fs = pyo.Block(self.model.SCENARIOS)
        for i in self.model.SCENARIOS:
            self.bidding_model_object.populate_model(self.model.fs[i])

        # save power output variable in the model object
        self._save_power_outputs()

        self.formulate_bidding_problem()

        # declare a list to store results
        self.bids_result_list = []

    def _save_power_outputs(self):

        """
        Create references of the power output variable in each price scenario
        block.

        Arguments:
            None

        Returns:
            None
        """

        for i in self.model.SCENARIOS:
            # get the power output
            power_output_name = self.bidding_model_object.power_output
            self.model.fs[i].power_output_ref = pyo.Reference(
                getattr(self.model.fs[i], power_output_name)
            )

        return

    def formulate_bidding_problem(self):

        """
        Formulate the bidding optimization problem by adding necessary
        parameters, constraints, and objective function.

        Arguments:
            None

        Returns:
            None
        """

        self._add_bidding_params()
        self._add_bidding_constraints()
        self._add_bidding_objective()

        return

    def _add_bidding_params(self):

        """
        Add necessary bidding parameters to the model, i.e., market energy price.

        Arguments:
            None

        Returns:
            None
        """

        for i in self.model.SCENARIOS:
            time_index = self.model.fs[i].power_output_ref.index_set()
            self.model.fs[i].energy_price = pyo.Param(
                time_index, initialize=0, mutable=True
            )
        return

    def _add_bidding_constraints(self):

        """
        Add bidding constraints to the model, i.e., the bid curves need to be
        nondecreasing.

        Arguments:
            None

        Returns:
            None
        """

        def bidding_constraint_rule(model, s1, s2, t):
            if s1 == s2:
                return pyo.Constraint.Skip
            return (
                model.fs[s1].power_output_ref[t] - model.fs[s2].power_output_ref[t]
            ) * (model.fs[s1].energy_price[t] - model.fs[s2].energy_price[t]) >= 0

        time_index = self.model.fs[
            self.model.SCENARIOS.first()
        ].power_output_ref.index_set()

        self.model.bidding_constraints = pyo.Constraint(
            self.model.SCENARIOS,
            self.model.SCENARIOS,
            time_index,
            rule=bidding_constraint_rule,
        )

        return

    def _add_bidding_objective(self):

        """
        Add objective function to the model, i.e., maximizing the expected profit
        of the energy system.

        Arguments:
            None

        Returns:
            None
        """

        # declare an empty objective
        self.model.obj = pyo.Objective(expr=0, sense=pyo.maximize)

        for k in self.model.SCENARIOS:
            time_index = self.model.fs[k].power_output_ref.index_set()

            # currently .total_cost is a tuple of 2 items
            # the first item is the name of the cost expression
            # the second item is the weight for the cost
            cost_name = self.bidding_model_object.total_cost[0]
            cost = getattr(self.model.fs[k], cost_name)
            weight = self.bidding_model_object.total_cost[1]

            for t in time_index:
                self.model.obj.expr += (
                    self.model.fs[k].power_output_ref[t]
                    * self.model.fs[k].energy_price[t]
                    - weight * cost[t]
                )

    def compute_bids(self, date, hour=None, **kwargs):

        """
        Solve the model to bid into the markets. After solving, record the bids
        from the solve.

        Arguments:

            date: current simulation date

            hour: current simulation hour

        Returns:
            None
        """

        price_forecasts = self.forecaster.forecast(date=date, hour=hour, **kwargs)

        # update the price forecasts
        self._pass_price_forecasts(price_forecasts)
        self.solver.solve(self.model, tee=True)
        bids = self._assemble_bids()
        self.record_bids(bids, date=date, hour=hour)

        return bids

    def update_model(self, **kwargs):

        """
        Update the flowsheets in all the price scenario blocks to advance time
        step.

        Arguments:
            kwargs: necessary profiles to update the underlying model. {stat_name: [...]}

        Returns:
            None
        """

        for i in self.model.SCENARIOS:
            self.bidding_model_object.update_model(b=self.model.fs[i], **kwargs)

    def _pass_price_forecasts(self, price_forecasts):

        """
        Pass the price forecasts into model parameters.

        Arguments:
            price_forecasts: price forecasts needed to solve the bidding problem. {LMP scenario: [forecast timeseries] }

        Returns:
            None
        """

        for i in self.model.SCENARIOS:
            time_index = self.model.fs[i].energy_price.index_set()
            for t, p in zip(time_index, price_forecasts[i]):

                # update the price param (mutable) in the pyomo model
                self.model.fs[i].energy_price[t] = p

        return

    def _assemble_bids(self):

        """
        This methods extract the bids out of the stochastic programming model and
        organize them.

        Arguments:

        Returns:
            bids: the bid we computed.
        """

        bids = {}
        gen = self.generator

        for i in self.model.SCENARIOS:
            time_index = self.model.fs[i].energy_price.index_set()
            for t in time_index:

                if t not in bids:
                    bids[t] = {}
                if gen not in bids[t]:
                    bids[t][gen] = {}

                power = round(pyo.value(self.model.fs[i].power_output_ref[t]), 2)
                marginal_cost = round(pyo.value(self.model.fs[i].energy_price[t]), 2)

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
                if power not in bids[t][gen]:
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
                tx_utils.validate_and_clean_cost_curve(
                    curve=temp_curve,
                    curve_type="cost_curve",
                    p_min=min([p[0] for p in bids[t][gen]]),
                    p_max=max([p[0] for p in bids[t][gen]]),
                    gen_name=gen,
                    t=t,
                )

        # create full bids: this includes info in addition to costs
        full_bids = {}

        for t in bids:
            full_bids[t] = {}
            for gen in bids[t]:
                full_bids[t][gen] = {}
                full_bids[t][gen]["p_cost"] = bids[t][gen]
                full_bids[t][gen]["p_min"] = min([p[0] for p in bids[t][gen]])
                full_bids[t][gen]["p_max"] = max([p[0] for p in bids[t][gen]])
                full_bids[t][gen]["startup_capacity"] = full_bids[t][gen]["p_min"]
                full_bids[t][gen]["shutdown_capacity"] = full_bids[t][gen]["p_min"]

                fixed_commitment = getattr(
                    self.bidding_model_object.model_data, "fixed_commitment", None
                )
                if fixed_commitment is not None:
                    full_bids[t][gen]["fixed_commitment"] = fixed_commitment

        return full_bids

    def _record_bids(self, bids, date, hour):

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

    def record_bids(self, bids, date, hour):

        """
        This function records the bids and the details in the underlying bidding model.

        Arguments:
            bids: the obtained bids for this date.

            date: the date we bid into

            hour: the hour we bid into

        Returns:
            None

        """

        # record bids
        self._record_bids(bids, date, hour)

        # record the details of bidding model
        for i in self.model.SCENARIOS:
            self.bidding_model_object.record_results(
                self.model.fs[i], date=date, hour=hour, Scenario=i
            )

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

    @property
    def generator(self):
        return self._generator

    @generator.setter
    def generator(self, name):
        self._generator = name
