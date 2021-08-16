import pandas as pd
import pyomo.environ as pyo
from collections import deque
import os

class Tracker:

    def __init__(self, tracking_model_class, n_tracking_hour, solver, **kwarg):

        '''
        Initializes the tracker object.

        Arguments:
            tracking_model_object: the initialized model object for tracking
            n_tracking_hour: number of implemented hours after each solve
            solver: a Pyomo mathematical programming solver object

        Returns:
            None
        '''

        # create an instance
        self.tracking_model_object = tracking_model_class(**kwarg)

        # add flowsheet to model
        self.model = pyo.ConcreteModel()
        self.model.fs = pyo.Block()
        b = self.tracking_model_object.create_model()
        self.model.fs.transfer_attributes_from(b)

        # get the power output
        power_output_name = self.tracking_model_object.power_output
        self.power_output = getattr(self.model.fs, power_output_name)

        # get the time index set
        self.time_set = self.power_output.index_set()

        self.n_tracking_hour = n_tracking_hour
        self.solver = solver
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

        self._add_tracking_params()
        self._add_tracking_constraints()
        self._add_tracking_objective()

        return

    def _add_tracking_params(self):

        '''
        Add necessary tracking parameters to the model, i.e., market dispatch
        signal.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        # add params to the model
        self.model.power_dispatch = pyo.Param(self.time_set, \
                                              initialize = 0, \
                                              within = pyo.Reals,\
                                              mutable = True)
        return

    def _add_tracking_constraints(self):

        '''
        Add necessary tracking constraints to the model, e.g., power output needs
        to follow market dispatch signals.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        self._add_tracking_dispatch_constraints()
        return

    def _add_tracking_dispatch_constraints(self):

        '''
        Add tracking constraints to the model, i.e., power output needs
        to follow market dispatch signals.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        # declare a constraint list
        self.model.tracking_dispatch_constraints = pyo.ConstraintList()
        for t in self.time_set:
            self.model.tracking_dispatch_constraints.add(self.power_output[t] == self.model.power_dispatch[t])

        return

    def _add_tracking_objective(self):

        '''
        Add EMPC objective function to the model, i.e., minimizing different costs
        of the energy system.

        Arguments:
            g: generator name in str

        Returns:
            None
        '''

        # declare an empty objective
        self.model.obj = pyo.Objective(expr = 0, sense = pyo.minimize)

        cost_name = self.tracking_model_object.total_cost[0]
        cost = getattr(self.model.fs, cost_name)
        weight = self.tracking_model_object.total_cost[1]

        for t in self.time_set:
            self.model.obj.expr += weight * cost[t]

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
        self.solver.solve(self.model,tee=True)

        self.record_results(date = date, hour = hour)

        # update the model
        profiles = self.tracking_model_object.get_implemented_profile(b = self.model.fs, \
                                                                      last_implemented_time_step = self.n_tracking_hour - 1)
        self.tracking_model_object.update_model(self.model.fs, **profiles)

    def _pass_market_dispatch(self, market_dispatch):

        '''
        Pass the received market signals into model parameters.

        Arguments:
            market_dispatch: a dictionary that contains the market dispatch signals
                             that we want to track. {generator name: [float]}

        Returns:
            None
        '''

        for t, dipsatch in zip(self.time_set, market_dispatch):
            self.model.power_dispatch[t] = dipsatch

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

        self.tracking_model_object.record_results(self.model.fs, **kwargs)

    def write_results(self,path):
        '''
        This methods writes the saved operation stats into an csv file.

        Arguments:
            path: the path to write the results.

        Return:
            None
        '''

        print("")
        print('Saving tracking results to disk...')
        self.tracking_model_object.write_results(path = os.path.join(path,'tracking_detail.csv'))
