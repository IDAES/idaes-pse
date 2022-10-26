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


class Tracker:

    """
    Wrap a model object to track the market dispatch signals. This class interfaces
    with the DoubleLoopCoordinator.
    """

    def __init__(
        self, tracking_model_object, tracking_horizon, n_tracking_hour, solver
    ):

        """
        Initializes the tracker object.

        Arguments:
            tracking_model_object: the model object for tracking

            tracking_horizon: number of time periods in the tracking problem

            n_tracking_hour: number of implemented hours after each solve

            solver: a Pyomo mathematical programming solver object

        Returns:
            None
        """

        # copy and check model object
        self.tracking_model_object = tracking_model_object
        self.tracking_horizon = tracking_horizon
        self.n_tracking_hour = n_tracking_hour
        self.solver = solver
        self._check_inputs()

        # add flowsheet to model
        self.model = pyo.ConcreteModel()
        self.model.fs = pyo.Block()
        self.tracking_model_object.populate_model(self.model.fs, self.tracking_horizon)

        # get the power output
        power_output_name = self.tracking_model_object.power_output
        self.power_output = getattr(self.model.fs, power_output_name)

        # get the time index set
        self.time_set = self.power_output.index_set()

        self.formulate_tracking_problem()

        self.daily_stats = None
        self.projection = None

        self.result_list = []

    def _check_inputs(self):

        """
        Check if the inputs to construct the tracker is valid. If not raise errors.
        """

        self._check_tracking_model_object()
        self._check_n_tracking_hour()
        self._check_solver()

    def _check_tracking_model_object(self):

        """
        Check if tracking model object has the necessary methods and attributes.
        """

        method_list = [
            "populate_model",
            "get_implemented_profile",
            "update_model",
            "get_last_delivered_power",
            "record_results",
            "write_results",
        ]
        attr_list = ["power_output", "total_cost"]
        msg = "Tracking model object does not have a "

        for m in method_list:
            obtained_m = getattr(self.tracking_model_object, m, None)
            if obtained_m is None:
                raise AttributeError(
                    msg
                    + m
                    + "() method. "
                    + "The tracker object needs the users to "
                    + "implement this method in their model object."
                )

        for attr in attr_list:
            obtained_attr = getattr(self.tracking_model_object, attr, None)
            if obtained_attr is None:
                raise AttributeError(
                    msg
                    + attr
                    + " property. "
                    + "The tracker object needs the users to "
                    + "specify this property in their model object."
                )

    def _check_n_tracking_hour(self):

        """
        Check if the number of hour for tracking is an integer and greater than 0.
        """

        # check if it is an integer
        if not isinstance(self.n_tracking_hour, int):
            raise TypeError(
                "The number of hour for tracking should be an integer, "
                + "but a {} was given.".format(type(self.n_tracking_hour).__name__)
            )

        if self.n_tracking_hour <= 0:
            raise ValueError(
                "The number of hour for tracking should be greater than zero, "
                + "but {} was given.".format(self.n_tracking_hour)
            )

    def _check_solver(self):

        """
        Check if provides solver is a valid Pyomo solver object.
        """

        if not isinstance(self.solver, OptSolver):
            raise TypeError(
                "The provided solver {} is not a valid Pyomo solver.".format(
                    self.solver
                )
            )

    def formulate_tracking_problem(self):

        """
        Formulate the tracking optimization problem by adding necessary
        parameters, constraints, and objective function.

        Arguments:
            None

        Returns:
            None
        """

        self._add_tracking_params()
        self._add_tracking_vars()
        self._add_tracking_constraints()
        self._add_tracking_objective()

        return

    def _add_tracking_vars(self):

        """
        Add necessary tracking variables to the model, i.e., power under and over
        delivered.

        Arguments:
            None

        Returns:
            None
        """

        self.model.power_underdelivered = pyo.Var(
            self.time_set, initialize=0, within=pyo.NonNegativeReals
        )
        self.model.power_overdelivered = pyo.Var(
            self.time_set, initialize=0, within=pyo.NonNegativeReals
        )

        return

    def _add_tracking_params(self):

        """
        Add necessary tracking parameters to the model, i.e., market dispatch
        signal.

        Arguments:
            None

        Returns:
            None
        """

        # add params to the model
        self.model.power_dispatch = pyo.Param(
            self.time_set, initialize=0, within=pyo.Reals, mutable=True
        )

        large_penalty = 10000
        penalty_init = {}
        for t in self.time_set:
            if t < self.n_tracking_hour:
                penalty_init[t] = large_penalty
            else:
                penalty_init[t] = large_penalty / (
                    self.tracking_horizon - self.n_tracking_hour
                )
        self.model.deviation_penalty = pyo.Param(
            self.time_set, initialize=penalty_init, mutable=False
        )

        return

    def _add_tracking_constraints(self):

        """
        Add necessary tracking constraints to the model, e.g., power output needs
        to follow market dispatch signals.

        Arguments:
            None

        Returns:
            None
        """

        # declare a constraint list
        def tracking_dispatch_constraint_rule(m, t):
            return (
                self.power_output[t] + self.model.power_underdelivered[t]
                == self.model.power_dispatch[t] + self.model.power_overdelivered[t]
            )

        self.model.tracking_dispatch_constraints = pyo.Constraint(
            self.time_set, rule=tracking_dispatch_constraint_rule
        )

        return

    def _add_tracking_objective(self):

        """
        Add EMPC objective function to the model, i.e., minimizing different costs
        of the energy system.

        Arguments:
            None

        Returns:
            None
        """

        # declare an empty objective
        self.model.obj = pyo.Objective(expr=0, sense=pyo.minimize)

        cost_name = self.tracking_model_object.total_cost[0]
        cost = getattr(self.model.fs, cost_name)
        weight = self.tracking_model_object.total_cost[1]

        for t in self.time_set:
            self.model.obj.expr += weight * cost[t] + self.model.deviation_penalty[
                t
            ] * (self.model.power_underdelivered[t] + self.model.power_overdelivered[t])

        return

    def update_model(self, **profiles):

        """
        This method updates the parameters in the model based on the implemented profiles.

        Arguments:
            profiles: the newly implemented stats. {stat_name: [...]}

        Returns:
            None
        """

        self.tracking_model_object.update_model(self.model.fs, **profiles)

    def track_market_dispatch(self, market_dispatch, date, hour):

        """
        Solve the model to track the market dispatch signals. After solving,
        record the results from the solve and update the model.

        Arguments:
            market_dispatch: a list that contains the market dispatch signals

            date: current simulation date

            hour: current simulation hour

        Returns:
            None
        """

        self._pass_market_dispatch(market_dispatch)

        # solve the model
        self.solver.solve(self.model, tee=False)

        self.record_results(date=date, hour=hour)

        # update the model
        profiles = self.tracking_model_object.get_implemented_profile(
            b=self.model.fs, last_implemented_time_step=self.n_tracking_hour - 1
        )

        self._record_daily_stats(profiles)

        return profiles

    def _record_daily_stats(self, profiles):

        """
        Record the stats that are used to update the model in the past 24 hours.

        Arguments:
            profiles: the newly implemented stats. {stat_name: [...]}

        Returns:
            None
        """

        if self.daily_stats is None:
            self.daily_stats = profiles
        else:
            for k in self.daily_stats:
                self.daily_stats[k] += profiles[k]

        for v in self.daily_stats.values():
            while len(v) >= 24:
                v.popleft()

        return

    def _pass_market_dispatch(self, market_dispatch):

        """
        Pass the received market signals into model parameters.

        Arguments:
            market_dispatch: a list that contains the market dispatch signals

        Returns:
            None
        """

        for t in self.time_set:

            try:
                dispatch = market_dispatch[t]
            except IndexError as ex:
                self.model.tracking_dispatch_constraints[t].deactivate()
            else:
                self.model.power_dispatch[t] = dispatch
                self.model.tracking_dispatch_constraints[t].activate()

        return

    def get_last_delivered_power(self):

        """
        Returns the last delivered power output.

        Arguments:
            None

        Returns:
            None
        """
        return self.tracking_model_object.get_last_delivered_power(
            b=self.model.fs, last_implemented_time_step=self.n_tracking_hour - 1
        )

    def _record_tracker_results(self, **kwargs):

        """
        Record the tracker stats.

        Arguments:
            kwargs: key word arguments that can be passed into tracking model object's record result function.

        Returns:
            None

        """

        df_list = []
        for t in self.time_set:

            result_dict = {}

            result_dict["Date"] = kwargs["date"]
            result_dict["Hour"] = kwargs["hour"]

            result_dict["Horizon [hr]"] = int(t)

            result_dict["Power Dispatch [MW]"] = round(
                pyo.value(self.model.power_dispatch[t]), 2
            )
            result_dict["Power Output [MW]"] = round(pyo.value(self.power_output[t]), 2)
            result_dict["Power Underdelivered [MW]"] = round(
                pyo.value(self.model.power_underdelivered[t]), 2
            )
            result_dict["Power Overdelivered [MW]"] = round(
                pyo.value(self.model.power_overdelivered[t]), 2
            )

            result_df = pd.DataFrame.from_dict(result_dict, orient="index")
            df_list.append(result_df.T)

        # append to result list
        self.result_list.append(pd.concat(df_list))

    def record_results(self, **kwargs):

        """
        Record the operations stats for the model.

        Arguments:
            kwargs: key word arguments that can be passed into tracking model object's record result function.

        Returns:
            None

        """

        # record tracker details
        self._record_tracker_results(**kwargs)

        # tracking model details
        self.tracking_model_object.record_results(self.model.fs, **kwargs)

    def write_results(self, path):
        """
        This methods writes the saved operation stats into an csv file.

        Arguments:
            path: the path to write the results.

        Return:
            None
        """

        print("")
        print("Saving tracking results to disk...")

        pd.concat(self.result_list).to_csv(
            os.path.join(path, "tracker_detail.csv"), index=False
        )
        self.tracking_model_object.write_results(
            path=os.path.join(path, "tracking_model_detail.csv")
        )
