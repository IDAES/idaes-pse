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
    def _add_UT_DT_constraints(m):

        '''
        This function adds the minimum up/down time constraints using eq. 4 - 5
        in "On mixed-integer programming formulations for the unit commitment
        problem". INFORMS Journal on Computing, 32(4), pp.857-876. Knueven, B.,
        Ostrowski, J. and Watson, J.P., 2020.

        Arguments:
            m: a pyomo model

        Returns:
            None
        '''

        def pre_shut_down_trajectory_set_rule(m):
            return (t for t in range(-pyo.value(m.min_dw_time) + 1,0))
        m.pre_shut_down_trajectory_set = pyo.Set(dimen = 1,initialize = pre_shut_down_trajectory_set_rule, ordered = True)

        def pre_start_up_trajectory_set_rule(m):
            return (t for t in range(-pyo.value(m.min_up_time) + 1,0))
        m.pre_start_up_trajectory_set = pyo.Set(dimen = 1,initialize = pre_start_up_trajectory_set_rule, ordered = True)

        m.pre_shut_down_trajectory = pyo.Param(m.pre_shut_down_trajectory_set, initialize = 0, mutable = True)
        m.pre_start_up_trajectory = pyo.Param(m.pre_start_up_trajectory_set, initialize = 0, mutable = True)

        def min_down_time_rule(m,h):
            if h < pyo.value(m.min_dw_time):
                return sum(m.pre_shut_down_trajectory[h0] for h0 in range(h - pyo.value(m.min_dw_time) + 1,0)) \
                       + sum(m.shut_dw[h0] for h0 in range(h + 1)) <= 1 - m.on_off[h]
            else:
                return sum(m.shut_dw[h0] for h0 in range(h - pyo.value(m.min_dw_time) + 1, h + 1)) <= 1 - m.on_off[h]
        m.min_down_time_con = pyo.Constraint(m.HOUR,rule = min_down_time_rule)

        def min_up_time_rule(m,h):
            if h < pyo.value(m.min_up_time):
                return sum(m.pre_start_up_trajectory[h0] for h0 in range(h - pyo.value(m.min_up_time) + 1,0)) \
                       + sum(m.start_up[h0] for h0 in range(h + 1)) <= m.on_off[h]
            else:
                return sum(m.start_up[h0] for h0 in range(h - pyo.value(m.min_up_time) + 1, h + 1)) <= m.on_off[h]
        m.min_up_time_con = pyo.Constraint(m.HOUR,rule = min_up_time_rule)

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
            m: the constructed model.
        '''

        model_data = self.model_data[generator]
        m = pyo.ConcreteModel()

        ## define the sets
        m.HOUR = pyo.Set(initialize = range(plan_horizon))
        m.SEGMENTS = pyo.Set(initialize = range(1, segment_number))

        ## define the parameters

        m.start_up_cost = pyo.Param(initialize = model_data['SU Cost'],mutable = False)

        # capacity of generators: upper bound (MW)
        m.Pmax = pyo.Param(initialize = model_data['PMax MW'], mutable = False)

        # minimum power of generators: lower bound (MW)
        m.Pmin = pyo.Param(initialize = model_data['PMin MW'], mutable = False)

        m.power_segment_bounds = pyo.Param(m.SEGMENTS,initialize = model_data['Power Segments'], mutable = False)

        # get the cost slopes
        m.F = pyo.Param(m.SEGMENTS,initialize = model_data['Marginal Costs'], mutable = False)

        m.min_load_cost = pyo.Param(initialize = model_data['Min Load Cost'], mutable = False)

        # Ramp up limits (MW/h)
        m.ramp_up = pyo.Param(initialize = model_data['RU'], mutable = False)

        # Ramp down limits (MW/h)
        m.ramp_dw = pyo.Param(initialize = model_data['RD'], mutable = False)

        # start up ramp limit
        m.ramp_start_up = pyo.Param(initialize = model_data['SU'], mutable = False)

        # shut down ramp limit
        m.ramp_shut_dw = pyo.Param(initialize = model_data['SD'], mutable = False)

        # minimum down time [hr]
        m.min_dw_time = pyo.Param(initialize = int(model_data['Min Down Time Hr']), mutable = False)

        # minimum up time [hr]
        m.min_up_time = pyo.Param(initialize = int(model_data['Min Up Time Hr']), mutable = False)

        # on/off status from previous day
        m.pre_on_off = pyo.Param(within = pyo.Binary,default= 1,mutable = True)

        # define a function to initialize the previous power params
        def init_pre_pow_fun(m):
            return m.pre_on_off*m.Pmin
        m.pre_P_T = pyo.Param(initialize = init_pre_pow_fun, mutable = True)

        ## define the variables

        # generator power (MW)

        # power generated by thermal generator
        m.P_T = pyo.Var(m.HOUR,initialize = pyo.value(m.Pmin), within = pyo.NonNegativeReals)

        # binary variables indicating on/off
        m.on_off = pyo.Var(m.HOUR,initialize = True, within = pyo.Binary)

        # binary variables indicating  start_up
        m.start_up = pyo.Var(m.HOUR,initialize = False, within = pyo.Binary)

        # binary variables indicating shut down
        m.shut_dw = pyo.Var(m.HOUR,initialize = False, within = pyo.Binary)

        # power produced in each segment
        m.power_segment = pyo.Var(m.HOUR,m.SEGMENTS, within = pyo.NonNegativeReals)

        ## Constraints

        # bounds on gen_pow
        def lhs_bnd_gen_pow_fun(m,h):
            return m.on_off[h] * m.Pmin <= m.P_T[h]
        m.lhs_bnd_gen_pow = pyo.Constraint(m.HOUR,rule = lhs_bnd_gen_pow_fun)

        def rhs_bnd_gen_pow_fun(m,h):
            return m.P_T[h] <= m.on_off[h] * m.Pmax
        m.rhs_bnd_gen_pow = pyo.Constraint(m.HOUR,rule = rhs_bnd_gen_pow_fun)

        # linearized power
        def linear_power_fun(m,h):
            return m.P_T[h] == \
            sum(m.power_segment[h,l] for l in m.SEGMENTS) + m.Pmin*m.on_off[h]
        m.linear_power = pyo.Constraint(m.HOUR,rule = linear_power_fun)

        # bounds on segment power
        def seg_pow_bnd_fun(m,h,l):
            if l == 1:
                return m.power_segment[h,l]<= (m.power_segment_bounds[l] - m.Pmin) * m.on_off[h]
            else:
                return m.power_segment[h,l]<= (m.power_segment_bounds[l] - m.power_segment_bounds[l-1]) * m.on_off[h]
        m.seg_pow_bnd = pyo.Constraint(m.HOUR,m.SEGMENTS,rule = seg_pow_bnd_fun)

        # start up and shut down logic (Arroyo and Conejo 2000)
        def start_up_shut_dw_fun(m,h):
            if h == 0:
                return m.start_up[h] - m.shut_dw[h] == m.on_off[h] - m.pre_on_off
            else:
                return m.start_up[h] - m.shut_dw[h] == m.on_off[h] - m.on_off[h-1]
        m.start_up_shut_dw = pyo.Constraint(m.HOUR,rule = start_up_shut_dw_fun)

        # either start up or shut down
        def start_up_or_shut_dw_fun(m,h):
            return m.start_up[h] + m.shut_dw[h] <= 1
        m.start_up_or_shut_dw = pyo.Constraint(m.HOUR,rule = start_up_or_shut_dw_fun)

        # ramp up limits
        def ramp_up_fun(m,h):
            '''
            h stand for hour
            '''
            if h==0:
                return m.P_T[h] <= m.pre_P_T \
                + m.ramp_up*m.pre_on_off\
                + m.ramp_start_up*m.start_up[h]
            else:
                return m.P_T[h] <= m.P_T[h-1] \
                + m.ramp_up*m.on_off[h-1]\
                + m.ramp_start_up*m.start_up[h]
        m.ramp_up_con = pyo.Constraint(m.HOUR,rule = ramp_up_fun)

        # ramp shut down limits
        def ramp_shut_dw_fun(m,h):
            '''
            h stand for hour.
            '''
            if h==0:
                return m.pre_P_T <= m.Pmax*m.on_off[h] + m.ramp_shut_dw * m.shut_dw[h]
            else:
                return m.P_T[h-1] <= m.Pmax*m.on_off[h] + m.ramp_shut_dw * m.shut_dw[h]
        m.ramp_shut_dw_con = pyo.Constraint(m.HOUR,rule = ramp_shut_dw_fun)

        # ramp down limits
        def ramp_dw_fun(m,h):
            '''
            h stand for hour.
            '''
            if h == 0:
                return m.pre_P_T - m.P_T[h] <= m.ramp_dw * m.on_off[h]\
                + m.ramp_shut_dw * m.shut_dw[h]
            else:
                return m.P_T[h-1] - m.P_T[h] <= m.ramp_dw * m.on_off[h]\
                + m.ramp_shut_dw * m.shut_dw[h]
        m.ramp_dw_con = pyo.Constraint(m.HOUR,rule = ramp_dw_fun)

        ## add min up and down time constraints
        self._add_UT_DT_constraints(m)

        ## Expression
        def prod_cost_fun(m,h):
            return m.min_load_cost * m.on_off[h] \
            + sum(m.F[l]*m.power_segment[h,l] for l in m.SEGMENTS)
        m.prod_cost_approx = pyo.Expression(m.HOUR,rule = prod_cost_fun)

        # start up costs
        def start_cost_fun(m,h):
            return m.start_up_cost * m.start_up[h]
        m.start_up_cost_expr = pyo.Expression(m.HOUR,rule = start_cost_fun)

        # total cost
        def tot_cost_fun(m,h):
            return m.prod_cost_approx[h] + m.start_up_cost_expr[h]
        m.tot_cost = pyo.Expression(m.HOUR,rule = tot_cost_fun)

        return m

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

        for g,m in self.model_dict.items():

            pre_shut_down_trajectory_copy[g] = deque([])
            pre_start_up_trajectory_copy[g] = deque([])

            # copy old trajectory
            for t in m.pre_shut_down_trajectory_set:
                pre_shut_down_trajectory_copy[g].append(round(pyo.value(m.pre_shut_down_trajectory[t])))
            for t in m.pre_start_up_trajectory_set:
                pre_start_up_trajectory_copy[g].append(round(pyo.value(m.pre_start_up_trajectory[t])))

            # add implemented trajectory to the queue
            pre_shut_down_trajectory_copy[g] += deque(implemented_shut_down[g])
            pre_start_up_trajectory_copy[g] += deque(implemented_start_up[g])

            # pop out outdated trajectory
            while len(pre_shut_down_trajectory_copy[g]) > pyo.value(m.min_dw_time) - 1:
                pre_shut_down_trajectory_copy[g].popleft()
            while len(pre_start_up_trajectory_copy[g]) > pyo.value(m.min_up_time) - 1:
                pre_start_up_trajectory_copy[g].popleft()

            # actual update
            for t in m.pre_shut_down_trajectory_set:
                m.pre_shut_down_trajectory[t] = pre_shut_down_trajectory_copy[g].popleft()

            for t in m.pre_start_up_trajectory_set:
                m.pre_start_up_trajectory[t] = pre_start_up_trajectory_copy[g].popleft()

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

        for g, m in self.model_dict.items():
            m.pre_P_T = round(implemented_power_output[g][-1],2)
            m.pre_on_off = round(int(implemented_power_output[g][-1] > 1e-3))

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

        for g, m in self.model_dict.items():
            for t in m.HOUR:

                result_dict = {}
                result_dict['Generator'] = g
                result_dict['Date'] = date
                result_dict['Hour'] = hour

                # simulation inputs
                result_dict['Horizon [hr]'] = int(t)

                # model vars
                result_dict['Thermal Power Generated [MW]'] = float(round(pyo.value(m.P_T[t]),2))

                result_dict['On/off [bin]'] = int(round(pyo.value(m.on_off[t])))
                result_dict['Start Up [bin]'] = int(round(pyo.value(m.start_up[t])))
                result_dict['Shut Down [bin]'] = int(round(pyo.value(m.shut_dw[t])))

                result_dict['Production Cost [$]'] = float(round(pyo.value(m.prod_cost_approx[t]),2))
                result_dict['Start-up Cost [$]'] = float(round(pyo.value(m.start_up_cost_expr[t]),2))
                result_dict['Total Cost [$]'] = float(round(pyo.value(m.tot_cost[t]),2))

                # calculate mileage
                if t == 0:
                    result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[t] - m.pre_P_T)),2))
                else:
                    result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[t] - m.P_T[t-1])),2))

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
        return {g: m.P_T for g, m in self.model_dict.items()}

    @property
    def total_cost(self):
        return {g: [(m.tot_cost,1)] for g, m in self.model_dict.items()}

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
    thermal_generator_object = ThermalGenerator(rts_gmlc_dataframe = rts_gmlc_dataframe, \
                                                horizon = horizon, \
                                                generators = [generator])

    solver = pyo.SolverFactory('cbc')

    run_tracker = True
    run_bidder = False

    if run_tracker:
        # make a tracker
        thermal_tracker = Tracker(tracking_model_object = thermal_generator_object,\
                                  n_tracking_hour = 1, \
                                  solver = solver)

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
