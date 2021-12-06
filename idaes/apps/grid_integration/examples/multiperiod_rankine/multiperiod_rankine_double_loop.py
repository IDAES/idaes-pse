import pyomo.environ as pyo
from collections import deque
import pandas as pd
import sys

from idaes.apps.grid_integration.examples.multiperiod_rankine.multiperiod_rankine_cycle import create_multiperiod_rankine_model

class MultiPeriodRankine:
    def __init__(self, horizon=4, generator_data={}):
        '''
        Arguments:
            horizon::Int64 - number of time points to use for associated multi-period model
            generator_data::Dict - dictionary of generator parameters (e.g. pmax,pmin,etc...)

        Returns:
            Float64: Value of power output in last time step
        '''
        self.horizon = horizon
        self.generator_data = generator_data
        self.generator = self.generator_data['generator_name']
        self.mp_rankine = None
        self.result_list = []

    def populate_model(self, blk):
        '''
        Create a rankine-cycle-battery model using the `MultiPeriod` package.

        Arguments:
            blk: this is an empty block passed in from eithe a bidder or tracker

         Returns:
             None
        '''
        if not blk.is_constructed():
            blk.construct()

        mp_rankine = create_multiperiod_rankine_model(n_time_points=self.horizon)
        blk.rankine = mp_rankine

        active_blks = mp_rankine.get_active_process_blocks()
        active_blks[0].battery.soc.fix(0) #initial charge is zero
        
        #create expression that references underlying power variables in multi-period rankine
        blk.HOUR = pyo.Set(initialize = range(self.horizon))
        blk.P_T = pyo.Expression(blk.HOUR)
        blk.tot_cost = pyo.Expression(blk.HOUR)
        for (t,b) in enumerate(active_blks):
            blk.P_T[t] = b.rankine.P_to_grid + b.battery.discharge
            blk.tot_cost[t] = b.rankine.fs.operating_cost

        self.mp_rankine = mp_rankine
        return

    #def update_model(self, b, implemented_shut_down,implemented_start_up, implemented_power_output):
    #NOTE: the profiles passed in here come from `get_implemented_profile`
    def update_model(self, blk, implemented_battery_charge, implemented_battery_discharge):

        '''
        Update `blk` variables using the actual implemented power output. 

        Arguments:
            blk: the block that needs to be updated
            implemented_battery_charge: list of power ,length n_tracking_horizon
            implemented_battery_discharge: 

         Returns:
             None
        '''
        mp_rankine = blk.mp_rankine 
        active_blks = mp_rankine.get_active_process_blocks()

        # power = round(implemented_power_output[-1])
        charge = round(implemented_battery_charge[-1]) 
        discharge = round(implemented_battery_discharge[-1]) 
        
        #update battery and power output based on implemented values
        battery_soc = active_blks[0].battery.soc
        new_battery_soc = value(battery_soc) + charge - discharge
        battery_soc.fix(new_battery_soc)

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
        This method gets the implemented variable profiles in the last optimization solve.

        Arguments:
            blk: a Pyomo block
            last_implemented_time_step: time index for the last implemented time step

         Returns:
             profile: the intended profile, {unit: [...]}
        '''
        #QUESTION: what is deque for?
        #implemented_battery_charge = deque([pyo.value(b.P_T[t]) for t in range(last_implemented_time_step + 1)])
        #implemented_battery_discharge = deque([pyo.value(b.P_T[t]) for t in range(last_implemented_time_step + 1)])
        mp_rankine = blk.mp_rankine
        active_blks = mp_rankine.get_active_process_blocks()
        implemented_battery_charge = deque([pyo.value(active_blks[t].rankine.P_to_battery*active_blks[t].battery.efficiency) for t in range(last_implemented_time_step + 1)])
        implemented_battery_discharge = deque([pyo.value(active_blks[t].battery.discharge/active_blks[t].battery.efficiency) for t in range(last_implemented_time_step + 1)])

        # return {'implemented_power_output':implemented_power_output,'implemented_battery_discharge':implemented_battery_discharge}

        return {'implemented_power_output':implemented_battery_charge,'implemented_battery_discharge':implemented_battery_discharge}


    
    def record_results(self, blk, date=None, hour=None, **kwargs):

        '''
        Record the operations stats for the model.

        Arguments:
            blk:  pyomo block
            date: current simulation date
            hour: current simulation hour

        Returns:
            None

        '''

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
            #result_dict['Production Cost [$]'] = float(round(pyo.value(blk.prod_cost_approx[t]),2))
            result_dict['Total Cost [$]'] = float(round(pyo.value(blk.tot_cost[t]),2))

            # calculate mileage
            # if t == 0:
            #     result_dict['Mileage [MW]'] = float(round(abs(pyo.value(blk.P_T[t] - blk.pre_P_T)),2))
            # else:
            #     result_dict['Mileage [MW]'] = float(round(abs(pyo.value(blk.P_T[t] - blk.P_T[t-1])),2))

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
        return self.generator_data['Original Marginal Cost Curve']

    @property
    def pmin(self):
        return self.generator_data['PMin MW']