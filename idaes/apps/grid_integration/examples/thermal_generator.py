import pyomo.environ as pyo
from collections import deque
import pandas as pd
# from idaes.apps.grid_integration.tracker import Tracker
import sys
sys.path.append('../')
from tracker import Tracker
from bidder import Bidder

class ThermalGenerator:

    def __init__(self, rts_gmlc_dataframe, horizon = 48, generators = None):

        '''
        Initializes the class object by building the thermal generator model.

        Arguments:
            rts_gmlc_dataframe: the RTS-GMLC generator data in Pandas DataFrame
            horizon: the length of the planning horizon of the model.
            generators: a list of generators in RTS-GMLC

        Returns:
            None
        '''

        self.generators = generators
        self.model_data = self.assemble_model_data(generator_names = generators, \
                                                   gen_params = rts_gmlc_dataframe)
        self.model_dict = {}
        for g in self.generators:
            self.model_dict[g] = self.build_thermal_generator_model(generator = g,
                                                                    plan_horizon = horizon,
                                                                    segment_number = 4)
        self.result_list = []

    @staticmethod
    def assemble_model_data(generator_names, gen_params):

        '''
        This function assembles the parameter data to build the thermal generator
        model, given a list of generator names and the RTS-GMLC data directory.

        Arguments:
            generator_names: a list of generator names in RTS-GMLC dataset.
            gen_params: the RTS-GMLC generator data in Pandas DataFrame

        Returns:
            model_data: a dictionary which has this structure
            {data type name: {generator name: value}}.
        '''

        gen_params.set_index('GEN UID',inplace = True)
        properties = ['PMin MW', 'PMax MW', 'Min Up Time Hr', 'Min Down Time Hr',\
                      'Ramp Rate MW/Min', 'Start Heat Warm MBTU', 'Fuel Price $/MMBTU',\
                      'HR_avg_0', 'HR_incr_1', 'HR_incr_2', 'HR_incr_3',\
                      'Output_pct_1','Output_pct_2','Output_pct_3']

        # to dict
        model_data = gen_params.loc[generator_names, properties].to_dict('index')

        # customize data
        for g in generator_names:

            model_data[g]['RU'] = model_data[g]['Ramp Rate MW/Min'] * 60
            model_data[g]['RD'] = model_data[g]['RU']
            model_data[g]['SU'] = min(model_data[g]['PMin MW'], model_data[g]['RU'])
            model_data[g]['SD'] = min(model_data[g]['PMin MW'], model_data[g]['RD'])
            model_data[g]['SU Cost'] = model_data[g]['Start Heat Warm MBTU'] * model_data[g]['Fuel Price $/MMBTU']

            model_data[g]['Min Load Cost'] = model_data[g]['HR_avg_0']/1000 * \
                                             model_data[g]['Fuel Price $/MMBTU'] *\
                                             model_data[g]['PMin MW']

            model_data[g]['Power Segments'] = {}
            model_data[g]['Marginal Costs'] = {}

            model_data[g]['Original Marginal Cost Curve'] = {}
            model_data[g]['Original Marginal Cost Curve'][model_data[g]['PMin MW']] = model_data[g]['Min Load Cost']/model_data[g]['PMin MW']

            for l in range(1,4):
                model_data[g]['Power Segments'][l] = model_data[g]['Output_pct_{}'.format(l)] * model_data[g]['PMax MW']
                model_data[g]['Marginal Costs'][l] = model_data[g]['HR_incr_{}'.format(l)]/1000 * model_data[g]['Fuel Price $/MMBTU']
                model_data[g]['Original Marginal Cost Curve'][model_data[g]['Power Segments'][l]] = model_data[g]['Marginal Costs'][l]

        return model_data

    @staticmethod
    def _add_UT_DT_constraints(b):

        '''
        This function adds the minimum up/down time constraints using eq. 4 - 5
        in "On mixed-integer programming formulations for the unit commitment
        problem". INFORMS Journal on Computing, 32(4), pp.857-876. Knueven, B.,
        Ostrowski, J. and Watson, J.P., 2020.

        Arguments:
            b: a pyomo model

        Returns:
            None
        '''

        def pre_shut_down_trajectory_set_rule(b):
            return (t for t in range(-pyo.value(b.min_dw_time) + 1,0))
        b.pre_shut_down_trajectory_set = pyo.Set(dimen = 1,initialize = pre_shut_down_trajectory_set_rule, ordered = True)

        def pre_start_up_trajectory_set_rule(b):
            return (t for t in range(-pyo.value(b.min_up_time) + 1,0))
        b.pre_start_up_trajectory_set = pyo.Set(dimen = 1,initialize = pre_start_up_trajectory_set_rule, ordered = True)

        b.pre_shut_down_trajectory = pyo.Param(b.pre_shut_down_trajectory_set, initialize = 0, mutable = True)
        b.pre_start_up_trajectory = pyo.Param(b.pre_start_up_trajectory_set, initialize = 0, mutable = True)

        def min_down_time_rule(b,h):
            if h < pyo.value(b.min_dw_time):
                return sum(b.pre_shut_down_trajectory[h0] for h0 in range(h - pyo.value(b.min_dw_time) + 1,0)) \
                       + sum(b.shut_dw[h0] for h0 in range(h + 1)) <= 1 - b.on_off[h]
            else:
                return sum(b.shut_dw[h0] for h0 in range(h - pyo.value(b.min_dw_time) + 1, h + 1)) <= 1 - b.on_off[h]
        b.min_down_time_con = pyo.Constraint(b.HOUR,rule = min_down_time_rule)

        def min_up_time_rule(b,h):
            if h < pyo.value(b.min_up_time):
                return sum(b.pre_start_up_trajectory[h0] for h0 in range(h - pyo.value(b.min_up_time) + 1,0)) \
                       + sum(b.start_up[h0] for h0 in range(h + 1)) <= b.on_off[h]
            else:
                return sum(b.start_up[h0] for h0 in range(h - pyo.value(b.min_up_time) + 1, h + 1)) <= b.on_off[h]
        b.min_up_time_con = pyo.Constraint(b.HOUR,rule = min_up_time_rule)

        return

    def build_thermal_generator_model(self,
                                      generator = None,
                                      plan_horizon = 48,
                                      segment_number = 4):

        '''
        This function builds the model for a thermal generator.

        Arguments:
            generator: generator name in str.
            plan_horizon: the length of the planning horizon of the model.
            segment_number: number of segments used in the piecewise linear
            production model.

        Returns:
            b: the constructed model.
        '''

        model_data = self.model_data[generator]
        b = pyo.Block()

        ## define the sets
        b.HOUR = pyo.Set(initialize = range(plan_horizon))
        b.SEGMENTS = pyo.Set(initialize = range(1, segment_number))

        ## define the parameters

        b.start_up_cost = pyo.Param(initialize = model_data['SU Cost'],mutable = False)

        # capacity of generators: upper bound (MW)
        b.Pmax = pyo.Param(initialize = model_data['PMax MW'], mutable = False)

        # minimum power of generators: lower bound (MW)
        b.Pmin = pyo.Param(initialize = model_data['PMin MW'], mutable = False)

        b.power_segment_bounds = pyo.Param(b.SEGMENTS,initialize = model_data['Power Segments'], mutable = False)

        # get the cost slopes
        b.F = pyo.Param(b.SEGMENTS,initialize = model_data['Marginal Costs'], mutable = False)

        b.min_load_cost = pyo.Param(initialize = model_data['Min Load Cost'], mutable = False)

        # Ramp up limits (MW/h)
        b.ramp_up = pyo.Param(initialize = model_data['RU'], mutable = False)

        # Ramp down limits (MW/h)
        b.ramp_dw = pyo.Param(initialize = model_data['RD'], mutable = False)

        # start up ramp limit
        b.ramp_start_up = pyo.Param(initialize = model_data['SU'], mutable = False)

        # shut down ramp limit
        b.ramp_shut_dw = pyo.Param(initialize = model_data['SD'], mutable = False)

        # minimum down time [hr]
        b.min_dw_time = pyo.Param(initialize = int(model_data['Min Down Time Hr']), mutable = False)

        # minimum up time [hr]
        b.min_up_time = pyo.Param(initialize = int(model_data['Min Up Time Hr']), mutable = False)

        # on/off status from previous day
        b.pre_on_off = pyo.Param(within = pyo.Binary,default= 1,mutable = True)

        # define a function to initialize the previous power params
        def init_pre_pow_fun(b):
            return b.pre_on_off*b.Pmin
        b.pre_P_T = pyo.Param(initialize = init_pre_pow_fun, mutable = True)

        ## define the variables

        # generator power (MW)

        # power generated by thermal generator
        b.P_T = pyo.Var(b.HOUR,initialize = model_data['PMin MW'], within = pyo.NonNegativeReals)

        # binary variables indicating on/off
        b.on_off = pyo.Var(b.HOUR,initialize = True, within = pyo.Binary)

        # binary variables indicating  start_up
        b.start_up = pyo.Var(b.HOUR,initialize = False, within = pyo.Binary)

        # binary variables indicating shut down
        b.shut_dw = pyo.Var(b.HOUR,initialize = False, within = pyo.Binary)

        # power produced in each segment
        b.power_segment = pyo.Var(b.HOUR,b.SEGMENTS, within = pyo.NonNegativeReals)

        ## Constraints

        # bounds on gen_pow
        def lhs_bnd_gen_pow_fun(b,h):
            return b.on_off[h] * b.Pmin <= b.P_T[h]
        b.lhs_bnd_gen_pow = pyo.Constraint(b.HOUR,rule = lhs_bnd_gen_pow_fun)

        def rhs_bnd_gen_pow_fun(b,h):
            return b.P_T[h] <= b.on_off[h] * b.Pmax
        b.rhs_bnd_gen_pow = pyo.Constraint(b.HOUR,rule = rhs_bnd_gen_pow_fun)

        # linearized power
        def linear_power_fun(b,h):
            return b.P_T[h] == \
            sum(b.power_segment[h,l] for l in b.SEGMENTS) + b.Pmin*b.on_off[h]
        b.linear_power = pyo.Constraint(b.HOUR,rule = linear_power_fun)

        # bounds on segment power
        def seg_pow_bnd_fun(b,h,l):
            if l == 1:
                return b.power_segment[h,l]<= (b.power_segment_bounds[l] - b.Pmin) * b.on_off[h]
            else:
                return b.power_segment[h,l]<= (b.power_segment_bounds[l] - b.power_segment_bounds[l-1]) * b.on_off[h]
        b.seg_pow_bnd = pyo.Constraint(b.HOUR,b.SEGMENTS,rule = seg_pow_bnd_fun)

        # start up and shut down logic (Arroyo and Conejo 2000)
        def start_up_shut_dw_fun(b,h):
            if h == 0:
                return b.start_up[h] - b.shut_dw[h] == b.on_off[h] - b.pre_on_off
            else:
                return b.start_up[h] - b.shut_dw[h] == b.on_off[h] - b.on_off[h-1]
        b.start_up_shut_dw = pyo.Constraint(b.HOUR,rule = start_up_shut_dw_fun)

        # either start up or shut down
        def start_up_or_shut_dw_fun(b,h):
            return b.start_up[h] + b.shut_dw[h] <= 1
        b.start_up_or_shut_dw = pyo.Constraint(b.HOUR,rule = start_up_or_shut_dw_fun)

        # ramp up limits
        def ramp_up_fun(b,h):
            '''
            h stand for hour
            '''
            if h==0:
                return b.P_T[h] <= b.pre_P_T \
                + b.ramp_up*b.pre_on_off\
                + b.ramp_start_up*b.start_up[h]
            else:
                return b.P_T[h] <= b.P_T[h-1] \
                + b.ramp_up*b.on_off[h-1]\
                + b.ramp_start_up*b.start_up[h]
        b.ramp_up_con = pyo.Constraint(b.HOUR,rule = ramp_up_fun)

        # ramp shut down limits
        def ramp_shut_dw_fun(b,h):
            '''
            h stand for hour.
            '''
            if h==0:
                return b.pre_P_T <= b.Pmax*b.on_off[h] + b.ramp_shut_dw * b.shut_dw[h]
            else:
                return b.P_T[h-1] <= b.Pmax*b.on_off[h] + b.ramp_shut_dw * b.shut_dw[h]
        b.ramp_shut_dw_con = pyo.Constraint(b.HOUR,rule = ramp_shut_dw_fun)

        # ramp down limits
        def ramp_dw_fun(b,h):
            '''
            h stand for hour.
            '''
            if h == 0:
                return b.pre_P_T - b.P_T[h] <= b.ramp_dw * b.on_off[h]\
                + b.ramp_shut_dw * b.shut_dw[h]
            else:
                return b.P_T[h-1] - b.P_T[h] <= b.ramp_dw * b.on_off[h]\
                + b.ramp_shut_dw * b.shut_dw[h]
        b.ramp_dw_con = pyo.Constraint(b.HOUR,rule = ramp_dw_fun)

        ## add min up and down time constraints
        self._add_UT_DT_constraints(b)

        ## Expression
        def prod_cost_fun(b,h):
            return b.min_load_cost * b.on_off[h] \
            + sum(b.F[l]*b.power_segment[h,l] for l in b.SEGMENTS)
        b.prod_cost_approx = pyo.Expression(b.HOUR,rule = prod_cost_fun)

        # start up costs
        def start_cost_fun(b,h):
            return b.start_up_cost * b.start_up[h]
        b.start_up_cost_expr = pyo.Expression(b.HOUR,rule = start_cost_fun)

        # total cost
        def tot_cost_fun(b,h):
            return b.prod_cost_approx[h] + b.start_up_cost_expr[h]
        b.tot_cost = pyo.Expression(b.HOUR,rule = tot_cost_fun)

        return b

    def _update_UT_DT(self, implemented_shut_down, implemented_start_up):
        '''
        This method updates the parameters in the minimum up/down time
        constraints based on the implemented shut down and start up events.

        Arguments:
            implemented_shut_down: realized shut down events: {unit: []}.
            implemented_start_up: realized start up events: {unit: []}

        Returns:
            None
        '''

        # copy to a queue
        pre_shut_down_trajectory_copy = {}
        pre_start_up_trajectory_copy = {}

        for g,b in self.model_dict.items():

            pre_shut_down_trajectory_copy[g] = deque([])
            pre_start_up_trajectory_copy[g] = deque([])

            # copy old trajectory
            for t in b.pre_shut_down_trajectory_set:
                pre_shut_down_trajectory_copy[g].append(round(pyo.value(b.pre_shut_down_trajectory[t])))
            for t in b.pre_start_up_trajectory_set:
                pre_start_up_trajectory_copy[g].append(round(pyo.value(b.pre_start_up_trajectory[t])))

            # add implemented trajectory to the queue
            pre_shut_down_trajectory_copy[g] += deque(implemented_shut_down[g])
            pre_start_up_trajectory_copy[g] += deque(implemented_start_up[g])

            # pop out outdated trajectory
            while len(pre_shut_down_trajectory_copy[g]) > pyo.value(b.min_dw_time) - 1:
                pre_shut_down_trajectory_copy[g].popleft()
            while len(pre_start_up_trajectory_copy[g]) > pyo.value(b.min_up_time) - 1:
                pre_start_up_trajectory_copy[g].popleft()

            # actual update
            for t in b.pre_shut_down_trajectory_set:
                b.pre_shut_down_trajectory[t] = pre_shut_down_trajectory_copy[g].popleft()

            for t in b.pre_start_up_trajectory_set:
                b.pre_start_up_trajectory[t] = pre_start_up_trajectory_copy[g].popleft()

        return

    def _update_power(self, implemented_power_output):
        '''
        This method updates the parameters in the ramping constraints based on
        the implemented power outputs.

        Arguments:
            implemented_power_output: realized power outputs: {unit: []}

         Returns:
             None
        '''

        for g, b in self.model_dict.items():
            b.pre_P_T = round(implemented_power_output[g][-1],2)
            b.pre_on_off = round(int(implemented_power_output[g][-1] > 1e-3))

        return

    def update_model(self,last_implemented_time_step):

        '''
        This method updates the parameters in the model based on
        the implemented power outputs, shut down and start up events.

        Arguments:
            last_implemented_time_step: time index for the last implemented time
                                        step

         Returns:
             None
        '''

        implemented_shut_down = self.get_implemented_profile(model_var = 'shut_dw',\
                                                             last_implemented_time_step = last_implemented_time_step)
        implemented_start_up = self.get_implemented_profile(model_var = 'start_up',\
                                                            last_implemented_time_step = last_implemented_time_step)
        implemented_power_output = self.get_implemented_profile(model_var = 'P_T',\
                                                                last_implemented_time_step = last_implemented_time_step)

        self._update_UT_DT(implemented_shut_down, implemented_start_up)
        self._update_power(implemented_power_output)

        return

    def get_implemented_profile(self, model_var, last_implemented_time_step):

        '''
        This method gets the implemented variable profiles in the last optimization
        solve.

        Arguments:
            model_var: intended variable name in str
            last_implemented_time_step: time index for the last implemented time
                                        step

         Returns:
             profile: the intended profile, {unit: [...]}
        '''

        profile = {}
        for g in self.generators:
            var = getattr(self.model_dict[g], model_var)
            profile[g] = [pyo.value(var[t]) for t in range(last_implemented_time_step + 1)]

        return profile

    def record_results(self, date = None, hour = None, **kwargs):

        '''
        Record the operations stats for the model.

        Arguments:

            date: current simulation date
            hour: current simulation hour

        Returns:
            None

        '''

        df_list = []

        for g, b in self.model_dict.items():
            for t in b.HOUR:

                result_dict = {}
                result_dict['Generator'] = g
                result_dict['Date'] = date
                result_dict['Hour'] = hour

                # simulation inputs
                result_dict['Horizon [hr]'] = int(t)

                # model vars
                result_dict['Thermal Power Generated [MW]'] = float(round(pyo.value(b.P_T[t]),2))

                result_dict['On/off [bin]'] = int(round(pyo.value(b.on_off[t])))
                result_dict['Start Up [bin]'] = int(round(pyo.value(b.start_up[t])))
                result_dict['Shut Down [bin]'] = int(round(pyo.value(b.shut_dw[t])))

                result_dict['Production Cost [$]'] = float(round(pyo.value(b.prod_cost_approx[t]),2))
                result_dict['Start-up Cost [$]'] = float(round(pyo.value(b.start_up_cost_expr[t]),2))
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

    @property
    def power_output(self):
        return {g: b.P_T for g, b in self.model_dict.items()}

    @property
    def total_cost(self):
        return {g: [(b.tot_cost,1)] for g, b in self.model_dict.items()}

    @property
    def default_bids(self):
        return {g: self.model_data[g]['Original Marginal Cost Curve'] for g in self.generators}

    @property
    def pmin(self):
        return {g: self.model_data[g]['PMin MW'] for g in self.generators}

if __name__ == "__main__":

    generator = "102_STEAM_3"
    horizon = 4

    rts_gmlc_dataframe = pd.read_csv('gen.csv')
    # thermal_generator_object = ThermalGenerator(rts_gmlc_dataframe = rts_gmlc_dataframe, \
    #                                             horizon = horizon, \
    #                                             generators = [generator])

    solver = pyo.SolverFactory('cbc')

    run_tracker = True
    run_bidder = False

    if run_tracker:
        # make a tracker
        thermal_tracker = Tracker(tracking_model_class = ThermalGenerator,\
                                  n_tracking_hour = 1, \
                                  solver = solver,\
                                  rts_gmlc_dataframe = rts_gmlc_dataframe,\
                                  horizon = horizon,\
                                  generators = [generator])

        market_dispatch = {generator: [30, 40 , 50, 70]}

        thermal_tracker.track_market_dispatch(market_dispatch = market_dispatch, \
                                              date = "2021-07-26", \
                                              hour = '17:00')

    if run_bidder:

        thermal_bidder = Bidder(bidding_model_object = thermal_generator_object,\
                                 solver = solver)

        price_forecasts = {generator:{0:[25,15,20,10], \
                                      1:[20,12,22,13], \
                                      2:[22,13,28,14]}}
        date = "2021-08-01"
        hour = "13:00"
        bids = thermal_bidder.compute_bids(price_forecasts, date, hour)
