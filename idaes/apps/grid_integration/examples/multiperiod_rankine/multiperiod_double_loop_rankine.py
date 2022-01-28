import pyomo.environ as pyo
from collections import deque
import pandas as pd
import numpy as np
import sys

from idaes.apps.grid_integration.examples.multiperiod_rankine.multiperiod_rankine_cycle import create_multiperiod_rankine_model

class MultiPeriodRankine:

    def __init__(self, horizon=4, pmin=20, pmax=100, default_bid_curve=None, generator_name="gen"):
        '''
        Arguments:
            horizon::Int64 - number of time points to use for associated multi-period model

        Returns:
            Float64: Value of power output in last time step
        '''
        if default_bid_curve==None:
            self.default_bid_curve = {p: 30 for p in np.linspace(pmin,pmax,5)}
        else:
            self.default_bid_curve = default_bid_curve
        self.horizon = horizon
        self.mp_rankine = None
        self.result_list = []
        self.p_lower = pmin
        self.p_upper = pmax
        self.generator = generator_name

    def populate_model(self, b):
        '''
        Create a rankine-cycle-battery model using the `MultiPeriod` package.

        Arguments:
            blk: this is an empty block passed in from eithe a bidder or tracker

        Returns:
             None
        '''
        blk = b
        if not blk.is_constructed():
            blk.construct()

        mp_rankine = create_multiperiod_rankine_model(n_time_points=self.horizon, pmin=self.p_lower, pmax=self.p_upper)
        blk.rankine = mp_rankine
        blk.rankine_model = mp_rankine.pyomo_model

        active_blks = mp_rankine.get_active_process_blocks()
        active_blks[0].battery.previous_soc.fix(0)
        active_blks[0].rankine.previous_power_output.fix(50.0)
        
        #create expression that references underlying power variables in multi-period rankine
        blk.HOUR = pyo.Set(initialize = range(self.horizon))
        blk.P_T = pyo.Expression(blk.HOUR)
        blk.tot_cost = pyo.Expression(blk.HOUR)
        for (t,b) in enumerate(active_blks):
            blk.P_T[t] = b.rankine.P_to_grid + b.battery.discharge
            blk.tot_cost[t] = b.rankine.fs.operating_cost

        self.mp_rankine = mp_rankine
        return

    def update_model(self, b, implemented_power_output, realized_soc):

        '''
        Update `blk` variables using the actual implemented power output. 

        Arguments:
            blk: the block that needs to be updated
            implemented_battery_charge: list of power ,length n_tracking_horizon
            implemented_battery_discharge: 
            implemented_power_output:

         Returns:
             None
        '''
        blk = b
        mp_rankine = blk.rankine 
        active_blks = mp_rankine.get_active_process_blocks()

        implemented_power = round(implemented_power_output[-1])
        realized_soc = round(realized_soc[-1])

        #update battery and power output based on implemented values
        active_blks[0].rankine.previous_power_output.fix(implemented_power)
        active_blks[0].battery.previous_soc.fix(realized_soc)

        return

    @staticmethod
    def get_last_delivered_power(b, last_implemented_time_step):

        '''
        Returns the last delivered power output.

        Arguments:
            blk: the block
            last_implemented_time_step: time index for the last implemented time
                                        step

        Returns:
            Float64: Value of power output in last time step
        '''
        blk = b
        return pyo.value(blk.P_T[last_implemented_time_step])

    @staticmethod
    def get_implemented_profile(b, last_implemented_time_step):

        '''
        This method gets the implemented variable profiles in the last optimization solve.

        Arguments:
            blk: a Pyomo block
            last_implemented_time_step: time index for the last implemented time step

         Returns:
             profile: the intended profile, {unit: [...]}
        '''
        blk = b
        mp_rankine = blk.rankine
        active_blks = mp_rankine.get_active_process_blocks()
        implemented_power_output = deque([pyo.value(active_blks[t].rankine.power_output) for t in range(last_implemented_time_step + 1)])
        realized_soc = deque([pyo.value(active_blks[t].battery.soc) for t in range(last_implemented_time_step + 1)])

        return {'implemented_power_output':implemented_power_output,'realized_soc':realized_soc}

    def record_results(self, b, date=None, hour=None, **kwargs):

        '''
        Record the operations stats for the model.

        Arguments:
            blk:  pyomo block
            date: current simulation date
            hour: current simulation hour

        Returns:
            None

        '''
        blk = b
        df_list = []
        for t in blk.HOUR:
            result_dict = {}

            #result_dict['Generator'] = self.generator
            result_dict['Date'] = date
            result_dict['Hour'] = hour

            # simulation inputs
            result_dict['Horizon [hr]'] = int(t)

            # model vars
            result_dict['Thermal Power Generated [MW]'] = float(round(pyo.value(blk.P_T[t]),2))
            result_dict['Total Cost [$]'] = float(round(pyo.value(blk.tot_cost[t]),2))

            for key in kwargs:
                result_dict[key] = kwargs[key]

            result_df = pd.DataFrame.from_dict(result_dict,orient = 'index')
            df_list.append(result_df.T)

        # append to result list
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
        return self.default_bid_curve

    @property
    def pmin(self):
        # return self.generator_data['PMin MW']
        return self.p_lower