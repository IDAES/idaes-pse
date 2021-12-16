import pyomo.environ as pyo
from collections import deque
import pandas as pd
# from idaes.apps.grid_integration.tracker import Tracker
import sys
sys.path.append('../')
from tracker import Tracker
from bidder import Bidder
from forecaster import PlaceHolderForecaster

from multiperiod_rankine_cycle import create_multiperiod_rankine_model

class MultiPeriodRankine:
    def __init__(self, horizon=4, generator_data={}):
        self.horizon = horizon
        self.generator_data = generator_data
        self.result_list = []

     #TODO create hooks on blk to integrate interface
    def populate_model(self, blk):
        mp_rankine = create_multiperiod_rankine_model(n_time_points=self.horizon)
        m = mp_rankine.pyomo_model
        blk.m = m
        # blk.Pmax = m.Pmax
        blk.HOUR = pyo.Set(initialize = list(range(self.horizon)))

        #use Expression?
        blk.P_T = pyo.Var(b.HOUR,initialize = model_data['PMin MW'], within = pyo.NonNegativeReals)


    def update_model(self, blk, implemented_power_output):

        '''
        Update `blk` using implemented power output. This effectively advances the time horizon.

        Arguments:
            blk: the block that needs to be updated
            implemented_power_output: realized power outputs: []

         Returns:
             None
        '''
        self._update_power(b, implemented_power_output)

        return

    @staticmethod
    def _update_power(blk, implemented_power_output):
        '''
        Update parameters for ramping constraints based on
        implemented power outputs.

        Arguments:
            blk: the pyomo blk to update
            implemented_power_output: realized power outputs: []

         Returns:
             None
        '''

        b.pre_P_T = round(implemented_power_output[-1],2)
        b.pre_on_off = round(int(implemented_power_output[-1] > 1e-3))

        return

    @staticmethod
    def get_last_delivered_power(blk, last_implemented_time_step):

        '''
        Returns the last delivered power output.

        Arguments:
            blk: the block
            last_implemented_time_step: time index for the last implemented time
                                        step

        Returns:
            Float64: Value of power output in last time step
        '''

        return pyo.value(blk.P_T[last_implemented_time_step])

    @staticmethod
    def get_implemented_profile(blk, last_implemented_time_step):

        '''
        This method gets the implemented variable profiles in the last optimization
        solve.

        Arguments:
            blk: a Pyomo block
            last_implemented_time_step: time index for the last implemented time step

         Returns:
             profile: the intended profile, {unit: [...]}
        '''
        implemented_power_output = deque([pyo.value(blk.P_T[t]) for t in range(last_implemented_time_step + 1)])

        return {'implemented_power_output': implemented_power_output}

    
    def record_results(self, blk, date=None, hour=None, **kwargs):

        '''
        Record the operations stats for the model.

        Arguments:

            date: current simulation date
            hour: current simulation hour

        Returns:
            None

        '''

        df_list = []
        for t in b.HOUR:

            result_dict = {}
            #result_dict['Generator'] = self.generator
            result_dict['Date'] = date
            result_dict['Hour'] = hour

            # simulation inputs
            result_dict['Horizon [hr]'] = int(t)

            # model vars
            result_dict['Thermal Power Generated [MW]'] = float(round(pyo.value(b.P_T[t]),2))

            result_dict['Production Cost [$]'] = float(round(pyo.value(b.prod_cost_approx[t]),2))
            result_dict['Total Cost [$]'] = float(round(pyo.value(b.tot_cost[t]),2))

            # calculate mileage
            if t == 0:
                result_dict['Mileage [MW]'] = float(round(abs(pyo.value(b.P_T[t] - b.pre_P_T)),2))
            else:
                result_dict['Mileage [MW]'] = float(round(abs(pyo.value(b.P_T[t] - b.P_T[t-1])),2))

            for key in kwargs:
                result_dict[key] = kwargs[key]

            result_df = pd.DataFrame.from_dict(result_dict,orient = 'index')
            df_list.append(result_df.T)

        # save the result to object property
        # wait to be written when simulation ends
        self.result_list.append(pd.concat(df_list))

        return

    def write_results(self, path):

        '''
        Write the saved results to a csv file.

        Arguments:
            path: the path to write the results.

        Return:
            None
        '''

        pd.concat(self.result_list).to_csv(path, index = False)

    @property
    def power_output(self):
        return 'P_T'

    @property
    def total_cost(self):
        return ('tot_cost',1)

    @property
    def default_bids(self):
        return self.generator_data['Original Marginal Cost Curve']

    @property
    def pmin(self):
        return self.generator_data['PMin MW']