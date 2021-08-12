import pandas as pd
import pyomo.environ as pyo
from collections import deque
import os

class Tracker:

    def __init__(self, tracking_model_object, n_tracking_hour, solver):

        '''
        Initializes the tracker object.

        Arguments:
            tracking_model_object: the initialized model object for tracking
            n_tracking_hour: number of implemented hours after each solve
            solver: a Pyomo mathematical programming solver object

        Returns:
            None
        '''

        self.tracking_model_object = tracking_model_object
        self.n_tracking_hour = n_tracking_hour
        self.solver = solver
        self.model_dict = self.tracking_model_object.model_dict
        self.formulate_tracking_problem()

    def formulate_tracking_problem(self):

        '''
        Formulate the tracking optimization problem by adding necessary
        parameters, constraints, and objective function.

        Arguments:
            None

        Returns:
            None
        '''

        self._get_time_set()

        for g,m in self.model_dict.items():
            self._add_tracking_params(g)
            self._add_tracking_constraints(g)
            self._add_tracking_objective(g)

        return

    def _get_time_set(self):

        '''
        Get the time index sets and store them in a dictionary as a class property.

        Arguments:
            None

        Returns:
            None
        '''

        self.time_set = {}
        for g in self.model_dict:
            self.time_set[g] = self.tracking_model_object.power_output[g].index_set()

        return

    def _add_tracking_params(self,g):

        '''
        Add necessary tracking parameters to the model, i.e., market dispatch
        signal.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        # add params to the model
        self.model_dict[g].power_dispatch = pyo.Param(self.time_set[g], \
                                                      initialize = 0, \
                                                      mutable = True)
        return

    def _add_tracking_constraints(self,g):

        '''
        Add necessary tracking constraints to the model, e.g., power output needs
        to follow market dispatch signals.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        self._add_tracking_dispatch_constraints(g)
        return

    def _add_tracking_dispatch_constraints(self,g):

        '''
        Add tracking constraints to the model, i.e., power output needs
        to follow market dispatch signals.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        # declare a constraint list
        self.model_dict[g].tracking_dispatch_constraints = pyo.ConstraintList()
        for t in self.time_set[g]:
            self.model_dict[g].tracking_dispatch_constraints.add(self.tracking_model_object.power_output[g][t] == self.model_dict[g].power_dispatch[t])

        return

    def _add_tracking_objective(self,g):

        '''
        Add EMPC objective function to the model, i.e., minimizing different costs
        of the energy system.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        # declare an empty objective
        self.model_dict[g].obj = pyo.Objective(expr = 0, sense = pyo.minimize)
        for c, w in self.tracking_model_object.total_cost[g]:
            for t in self.time_set[g]:
                self.model_dict[g].obj.expr += w * c[t]

        return

    def track_market_dispatch(self, market_dispatch, date, hour):

        '''
        Solve the model to track the market dispatch signals. After solving,
        record the results from the solve and update the model.

        Arguments:
            market_dispatch: a dictionary that contains the market dispatch signals
                             that we want to track. {generator name: [float]}
            date: current simulation date
            hour: current simulation hour

        Returns:
            None
        '''

        self._pass_market_dispatch(market_dispatch)

        # solve the model
        for g,m in self.model_dict.items():
            self.solver.solve(m,tee=True)

        self.record_results(date = date, hour = hour)

        # update the model
        self.tracking_model_object.update_model(last_implemented_time_step = self.n_tracking_hour - 1)

    def _pass_market_dispatch(self, market_dispatch):

        '''
        Pass the received market signals into model parameters.

        Arguments:
            market_dispatch: a dictionary that contains the market dispatch signals
                             that we want to track. {generator name: [float]}

        Returns:
            None
        '''

        for g,m in self.model_dict.items():
            for t, dipsatch in zip(self.time_set[g], market_dispatch[g]):
                m.power_dispatch[t] = dipsatch

        return

    def record_results(self, **kwargs):

        '''
        Record the operations stats for the model.

        Arguments:
            kwargs: key word arguments that can be passed into tracking model
                    object's record result function.

        Returns:
            None

        '''

        self.tracking_model_object.record_results(**kwargs)
