import pyomo.environ as pyo
from collections import deque, OrderedDict
import pandas as pd
# from idaes.apps.grid_integration.tracker import Tracker
import sys
sys.path.append('../')
from tracker import Tracker
from bidder import Bidder

class ThermalGenerator:

    def __init__(self, rts_gmlc_dataframe, horizon = 48,generators = None, n_scenario = 10):

        '''
        Initializes the class object by building the thermal generator model.

        Arguments:
            rts_gmlc_dataframe: the RTS-GMLC generator data in Pandas DataFrame
            horizon: the length of the planning horizon of the model.
            generators: a list of generators in RTS-GMLC
            n_scenario: number of uncertain scenarios.

        Returns:
            None
        '''

        self.n_scenario = n_scenario

        self.model_data = self.assemble_model_data(generator_names = generators, \
                                                   gen_params = rts_gmlc_dataframe)
        self.model = self.build_thermal_generator_model(plan_horizon = horizon,
                                                        segment_number = 4)
        self.result_list = []

    @staticmethod
    def assemble_model_data(generator_names, gen_params, **kwargs):

        '''
        This function assembles the parameter data to build the thermal generator
        model, given a list of generator names and the RTS-GMLC data directory.

        Arguments:
            generator_names: a list of generator names in RTS-GMLC dataset.
            rts_gmlc_data_dir: the RTS-GMLC data directory.

        Returns:
            model_data: a dictionary which has this structure
            {data type name: {generator name: value}}.
        '''

        # read data
        gen_params = gen_params[gen_params['GEN UID'].isin(generator_names)]

        model_data = {}

        # generator names
        model_data['Generator'] = generator_names

        # Pmin [MW]
        model_data['Pmin'] = dict(zip(generator_names,gen_params['PMin MW']))

        # Pmax [MW]
        model_data['Pmax'] = dict(zip(generator_names,gen_params['PMax MW']))

        # minimum up time [MW/hr]
        model_data['UT'] = dict(zip(generator_names,gen_params['Min Up Time Hr'].astype(int)))

        # minimum down time [hr]
        model_data['DT'] = dict(zip(generator_names,gen_params['Min Down Time Hr'].astype(int)))

        ## ramp rates [MW/hr]
        ramp_rates = gen_params['Ramp Rate MW/Min'].values * 60

        # ramp up rate [MW/hr]
        model_data['RU'] = dict(zip(generator_names,ramp_rates))

        # ramp down rate [MW/hr]
        model_data['RD'] = dict(zip(generator_names,ramp_rates))

        # ramp start up [MW/hr]
        model_data['SU'] = {gen: min(model_data['Pmin'][gen],model_data['RU'][gen]) for gen in generator_names}

        # ramp shut down [MW/hr] (will use pmin for now)
        model_data['SD'] = {gen: min(model_data['Pmin'][gen],model_data['RD'][gen]) for gen in generator_names}

        # start up cost [$/SU] (will use the warm start up cost for now)
        start_up_cost = gen_params['Start Heat Warm MBTU'] * gen_params['Fuel Price $/MMBTU']
        model_data['SU Cost'] = dict(zip(generator_names,start_up_cost))

        ## production cost

        # power segments and marginal costs
        model_data['Power Segments'] = {}
        model_data['Marginal Costs'] = {}

        model_data['Min Load Cost'] = dict(zip(generator_names,gen_params['HR_avg_0']/1000 * gen_params['Fuel Price $/MMBTU'] * gen_params['PMin MW']))
        gen_params.set_index('GEN UID',inplace=True)
        for gen in generator_names:
            for l in range(1,4):
                # power segements
                model_data['Power Segments'][(gen,l)] = gen_params.loc[gen,'Output_pct_{}'.format(l)] * gen_params.loc[gen,'PMax MW']

                # segment marginal cost
                model_data['Marginal Costs'][(gen,l)] = gen_params.loc[gen,'HR_incr_{}'.format(l)]/1000 * gen_params.loc[gen,'Fuel Price $/MMBTU']

        gen_params.reset_index(inplace=True)

        # get the original cost curve
        model_data['Original Cost Curve'] = {}
        for gen in generator_names:
            model_data['Original Cost Curve'][gen] = OrderedDict()

            pmin = round(model_data['Pmin'][gen],2)
            model_data['Original Cost Curve'][gen][pmin] = model_data['Min Load Cost'][gen]

            old_p = pmin
            old_cost = model_data['Original Cost Curve'][gen][pmin]
            for l in range(1,4):

                new_p = round(model_data['Power Segments'][(gen,l)],2)
                delta_p = new_p - old_p

                increased_cost = delta_p * model_data['Marginal Costs'][(gen,l)]
                model_data['Original Cost Curve'][gen][new_p] = old_cost + increased_cost

                old_cost += increased_cost
                old_p = new_p

        model_data['Original Marginal Cost Curve'] = {}
        for gen in generator_names:
            model_data['Original Marginal Cost Curve'][gen] = OrderedDict()

            pmin = round(model_data['Pmin'][gen],2)
            model_data['Original Marginal Cost Curve'][gen][pmin] = model_data['Min Load Cost'][gen]/pmin

            for l in range(1,4):
                new_p = round(model_data['Power Segments'][(gen,l)],2)
                model_data['Original Marginal Cost Curve'][gen][new_p] = model_data['Marginal Costs'][(gen,l)]

        for key in kwargs:
            model_data[key] = {gen: kwargs[key] for gen in generator_names}

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
            return ((j,t) for j in m.UNITS for t in range(-pyo.value(m.min_dw_time[j]) + 1,0))
        m.pre_shut_down_trajectory_set = pyo.Set(dimen = 2,initialize = pre_shut_down_trajectory_set_rule, ordered = True)

        def pre_start_up_trajectory_set_rule(m):
            return ((j,t) for j in m.UNITS for t in range(-pyo.value(m.min_up_time[j]) + 1,0))
        m.pre_start_up_trajectory_set = pyo.Set(dimen = 2,initialize = pre_start_up_trajectory_set_rule, ordered = True)

        m.pre_shut_down_trajectory = pyo.Param(m.pre_shut_down_trajectory_set, initialize = 0, mutable = True)
        m.pre_start_up_trajectory = pyo.Param(m.pre_start_up_trajectory_set, initialize = 0, mutable = True)

        def min_down_time_rule(m,j,h,k):
            if h < pyo.value(m.min_dw_time[j]):
                return sum(m.pre_shut_down_trajectory[j,h0] for h0 in range(h - pyo.value(m.min_dw_time[j]) + 1,0)) \
                       + sum(m.shut_dw[j,h0,k] for h0 in range(h + 1)) <= 1 - m.on_off[j,h,k]
            else:
                return sum(m.shut_dw[j,h0,k] for h0 in range(h - pyo.value(m.min_dw_time[j]) + 1, h + 1)) <= 1 - m.on_off[j,h,k]
        m.min_down_time_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = min_down_time_rule)

        def min_up_time_rule(m,j,h,k):
            if h < pyo.value(m.min_up_time[j]):
                return sum(m.pre_start_up_trajectory[j,h0] for h0 in range(h - pyo.value(m.min_up_time[j]) + 1,0)) \
                       + sum(m.start_up[j,h0,k] for h0 in range(h + 1)) <= m.on_off[j,h,k]
            else:
                return sum(m.start_up[j,h0,k] for h0 in range(h - pyo.value(m.min_up_time[j]) + 1, h + 1)) <= m.on_off[j,h,k]
        m.min_up_time_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = min_up_time_rule)

        return

    def build_thermal_generator_model(self,
                                      plan_horizon = 48,
                                      segment_number = 4):

        '''
        This function builds the model for a thermal generator.

        Arguments:
            plan_horizon: the length of the planning horizon of the model.
            segment_number: number of segments used in the piecewise linear
            production model.

        Returns:
            m: the constructed model.
        '''

        model_data = self.model_data
        m = pyo.ConcreteModel()

        ## define the sets
        m.HOUR = pyo.Set(initialize = range(plan_horizon))
        m.SEGMENTS = pyo.Set(initialize = range(1, segment_number))
        m.UNITS = pyo.Set(initialize = model_data['Generator'], ordered = True)
        m.SCENARIOS = pyo.Set(initialize = range(self.n_scenario))

        ## define the parameters

        m.start_up_cost = pyo.Param(m.UNITS,initialize = model_data['SU Cost'],mutable = False)

        # capacity of generators: upper bound (MW)
        m.Pmax = pyo.Param(m.UNITS,initialize = model_data['Pmax'], mutable = False)

        # minimum power of generators: lower bound (MW)
        m.Pmin = pyo.Param(m.UNITS,initialize = model_data['Pmin'], mutable = False)

        m.power_segment_bounds = pyo.Param(m.UNITS,m.SEGMENTS,initialize = model_data['Power Segments'], mutable = False)

        # get the cost slopes
        m.F = pyo.Param(m.UNITS,m.SEGMENTS,initialize = model_data['Marginal Costs'], mutable = False)

        m.min_load_cost = pyo.Param(m.UNITS,initialize = model_data['Min Load Cost'], mutable = False)

        # Ramp up limits (MW/h)
        m.ramp_up = pyo.Param(m.UNITS,initialize = model_data['RU'], mutable = False)

        # Ramp down limits (MW/h)
        m.ramp_dw = pyo.Param(m.UNITS,initialize = model_data['RD'], mutable = False)

        # start up ramp limit
        m.ramp_start_up = pyo.Param(m.UNITS,initialize = model_data['SU'], mutable = False)

        # shut down ramp limit
        m.ramp_shut_dw = pyo.Param(m.UNITS,initialize = model_data['SD'], mutable = False)

        # minimum down time [hr]
        m.min_dw_time = pyo.Param(m.UNITS,initialize = model_data['DT'], mutable = False)

        # minimum up time [hr]
        m.min_up_time = pyo.Param(m.UNITS,initialize = model_data['UT'], mutable = False)

        # power from the previous day (MW)
        # need to assume the power output is at least always at the minimum pow output

        # on/off status from previous day
        m.pre_on_off = pyo.Param(m.UNITS,within = pyo.Binary,default= 1,mutable = True)

        # define a function to initialize the previous power params
        def init_pre_pow_fun(m,j):
            return m.pre_on_off[j]*m.Pmin[j]
        m.pre_P_T = pyo.Param(m.UNITS,initialize = init_pre_pow_fun, mutable = True)

        ## define the variables

        # generator power (MW)

        # power generated by thermal generator
        m.P_T = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.NonNegativeReals)

        # binary variables indicating on/off
        m.on_off = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,initialize = True, within = pyo.Binary)

        # binary variables indicating  start_up
        m.start_up = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,initialize = False, within = pyo.Binary)

        # binary variables indicating shut down
        m.shut_dw = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,initialize = False, within = pyo.Binary)

        # power produced in each segment
        m.power_segment = pyo.Var(m.UNITS,m.HOUR,m.SEGMENTS, m.SCENARIOS, within = pyo.NonNegativeReals)

        ## Constraints

        # bounds on gen_pow
        def lhs_bnd_gen_pow_fun(m,j,h,k):
            return m.on_off[j,h,k] * m.Pmin[j] <= m.P_T[j,h,k]
        m.lhs_bnd_gen_pow = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = lhs_bnd_gen_pow_fun)

        def rhs_bnd_gen_pow_fun(m,j,h,k):
            return m.P_T[j,h,k] <= m.on_off[j,h,k] * m.Pmax[j]
        m.rhs_bnd_gen_pow = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = rhs_bnd_gen_pow_fun)

        # linearized power
        def linear_power_fun(m,j,h,k):
            return m.P_T[j,h,k] == \
            sum(m.power_segment[j,h,l,k] for l in m.SEGMENTS) + m.Pmin[j]*m.on_off[j,h,k]
        m.linear_power = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = linear_power_fun)

        # bounds on segment power
        def seg_pow_bnd_fun(m,j,h,l,k):
            if l == 1:
                return m.power_segment[j,h,l,k]<= (m.power_segment_bounds[j,l] - m.Pmin[j]) * m.on_off[j,h,k]
            else:
                return m.power_segment[j,h,l,k]<= (m.power_segment_bounds[j,l] - m.power_segment_bounds[j,l-1]) * m.on_off[j,h,k]
        m.seg_pow_bnd = pyo.Constraint(m.UNITS,m.HOUR,m.SEGMENTS,m.SCENARIOS,rule = seg_pow_bnd_fun)

        # start up and shut down logic (Arroyo and Conejo 2000)
        def start_up_shut_dw_fun(m,j,h,k):
            if h == 0:
                return m.start_up[j,h,k] - m.shut_dw[j,h,k] == m.on_off[j,h,k] - m.pre_on_off[j]
            else:
                return m.start_up[j,h,k] - m.shut_dw[j,h,k] == m.on_off[j,h,k] - m.on_off[j,h-1,k]
        m.start_up_shut_dw = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = start_up_shut_dw_fun)

        # either start up or shut down
        def start_up_or_shut_dw_fun(m,j,h,k):
            return m.start_up[j,h,k] + m.shut_dw[j,h,k] <= 1
        m.start_up_or_shut_dw = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = start_up_or_shut_dw_fun)

        # ramp up limits
        def ramp_up_fun(m,j,h,k):
            '''
            j,h,k stand for unit, hour,scenario respectively.
            '''
            if h==0:
                return m.P_T[j,h,k] <= m.pre_P_T[j] \
                + m.ramp_up[j]*m.pre_on_off[j]\
                + m.ramp_start_up[j]*m.start_up[j,h,k]
            else:
                return m.P_T[j,h,k] <= m.P_T[j,h-1,k] \
                + m.ramp_up[j]*m.on_off[j,h-1,k]\
                + m.ramp_start_up[j]*m.start_up[j,h,k]
        m.ramp_up_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = ramp_up_fun)

        # ramp shut down limits
        def ramp_shut_dw_fun(m,j,h,k):
            '''
            j,h,k stand for unit, hour,scenario respectively.
            '''
            if h==0:
                return m.pre_P_T[j] <= m.Pmax[j]*m.on_off[j,h,k] + m.ramp_shut_dw[j] * m.shut_dw[j,h,k]
            else:
                return m.P_T[j,h-1,k] <= m.Pmax[j]*m.on_off[j,h,k] + m.ramp_shut_dw[j] * m.shut_dw[j,h,k]
        m.ramp_shut_dw_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = ramp_shut_dw_fun)

        # ramp down limits
        def ramp_dw_fun(m,j,h,k):
            '''
            j,h,k stand for unit, hour,scenario respectively.
            '''
            if h == 0:
                return m.pre_P_T[j] - m.P_T[j,h,k] <= m.ramp_dw[j] * m.on_off[j,h,k]\
                + m.ramp_shut_dw[j] * m.shut_dw[j,h,k]
            else:
                return m.P_T[j,h-1,k] - m.P_T[j,h,k] <= m.ramp_dw[j] * m.on_off[j,h,k]\
                + m.ramp_shut_dw[j] * m.shut_dw[j,h,k]
        m.ramp_dw_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = ramp_dw_fun)

        ## add min up and down time constraints
        self._add_UT_DT_constraints(m)

        ## Expression
        def prod_cost_fun(m,j,h,k):
            return m.min_load_cost[j] * m.on_off[j,h,k] \
            + sum(m.F[j,l]*m.power_segment[j,h,l,k] for l in m.SEGMENTS)
        m.prod_cost_approx = pyo.Expression(m.UNITS,m.HOUR,m.SCENARIOS,rule = prod_cost_fun)

        # start up costs
        def start_cost_fun(m,j,h,k):
            return m.start_up_cost[j]*m.start_up[j,h,k]
        m.start_up_cost_expr = pyo.Expression(m.UNITS,m.HOUR,m.SCENARIOS,rule = start_cost_fun)

        # total cost
        def tot_cost_fun(m,j,h,k):
            return m.prod_cost_approx[j,h,k] + m.start_up_cost_expr[j,h,k]
        m.tot_cost = pyo.Expression(m.UNITS,m.HOUR,m.SCENARIOS,rule = tot_cost_fun)

        return m

    @staticmethod
    def _update_UT_DT(m, implemented_shut_down, implemented_start_up):
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

        for unit in m.UNITS:
            pre_shut_down_trajectory_copy[unit] = deque([])
            pre_start_up_trajectory_copy[unit] = deque([])

        for unit,t in m.pre_shut_down_trajectory_set:
            pre_shut_down_trajectory_copy[unit].append(round(pyo.value(m.pre_shut_down_trajectory[unit,t])))
        for unit,t in m.pre_start_up_trajectory_set:
            pre_start_up_trajectory_copy[unit].append(round(pyo.value(m.pre_start_up_trajectory[unit,t])))

        # add implemented trajectory to the queue
        for unit in m.UNITS:
            pre_shut_down_trajectory_copy[unit] += deque(implemented_shut_down[unit])
            pre_start_up_trajectory_copy[unit] += deque(implemented_start_up[unit])

        # pop out outdated trajectory
        for unit in m.UNITS:

            while len(pre_shut_down_trajectory_copy[unit]) > pyo.value(m.min_dw_time[unit]) - 1:
                pre_shut_down_trajectory_copy[unit].popleft()
            while len(pre_start_up_trajectory_copy[unit]) > pyo.value(m.min_up_time[unit]) - 1:
                pre_start_up_trajectory_copy[unit].popleft()

        # actual update
        for unit,t in m.pre_shut_down_trajectory_set:
            m.pre_shut_down_trajectory[unit,t] = pre_shut_down_trajectory_copy[unit].popleft()

        for unit,t in m.pre_start_up_trajectory_set:
            m.pre_start_up_trajectory[unit,t] = pre_start_up_trajectory_copy[unit].popleft()

        return

    @staticmethod
    def _update_power(m,implemented_power_output):
        '''
        This method updates the parameters in the ramping constraints based on
        the implemented power outputs.

        Arguments:
            implemented_power_output: realized power outputs: {unit: []}

         Returns:
             None
        '''

        for unit in m.UNITS:
            m.pre_P_T[unit] = round(implemented_power_output[unit][-1],2)
            m.pre_on_off[unit] = round(int(implemented_power_output[unit][-1] > 1e-3))

        return

    def update_model(self,last_implemented_time_step):

        '''
        This method updates the parameters in the model based on
        the implemented power outputs, shut down and start up events.

        Arguments:
            implemented_power_output: realized power outputs: {unit: []}
            implemented_shut_down: realized shut down events: {unit: []}.
            implemented_start_up: realized start up events: {unit: []}

         Returns:
             None
        '''

        implemented_shut_down = self.get_implemented_profile(model_var = self.model.shut_dw,\
                                                             last_implemented_time_step = last_implemented_time_step)
        implemented_start_up = self.get_implemented_profile(model_var = self.model.start_up,\
                                                            last_implemented_time_step = last_implemented_time_step)
        implemented_power_output = self.get_implemented_profile(model_var = self.model.P_T,\
                                                                last_implemented_time_step = last_implemented_time_step)

        self._update_UT_DT(self.model,implemented_shut_down, implemented_start_up)
        self._update_power(self.model,implemented_power_output)

        return

    def get_implemented_profile(self, model_var, last_implemented_time_step):

        profile = {}
        for g in self.model.UNITS:
            profile[g] = [pyo.value(model_var[g,t,0]) for t in range(last_implemented_time_step + 1)]

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

        m = self.model
        df_list = []
        for generator in m.UNITS:
            for t in m.HOUR:
                for k in m.SCENARIOS:

                    result_dict = {}
                    result_dict['Generator'] = generator
                    result_dict['Date'] = date
                    result_dict['Hour'] = hour

                    # simulation inputs
                    result_dict['Horizon [hr]'] = int(t)
                    result_dict['Scenario'] = int(k)

                    # model vars
                    result_dict['Thermal Power Generated [MW]'] = float(round(pyo.value(m.P_T[generator,t,k]),2))

                    result_dict['On/off [bin]'] = int(round(pyo.value(m.on_off[generator,t,k])))
                    result_dict['Start Up [bin]'] = int(round(pyo.value(m.start_up[generator,t,k])))
                    result_dict['Shut Down [bin]'] = int(round(pyo.value(m.shut_dw[generator,t,k])))

                    result_dict['Production Cost [$]'] = float(round(pyo.value(m.prod_cost_approx[generator,t,k]),2))
                    result_dict['Start-up Cost [$]'] = float(round(pyo.value(m.start_up_cost_expr[generator,t,k]),2))
                    result_dict['Total Cost [$]'] = float(round(pyo.value(m.tot_cost[generator,t,k]),2))

                    # calculate mileage
                    if t == 0:
                        result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[generator,t,k] - m.pre_P_T[generator])),2))
                    else:
                        result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[generator,t,k] - m.P_T[generator,t-1,k])),2))

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
        return self.model.P_T

    @property
    def total_cost(self):
        return {self.model.tot_cost: 1}

    @property
    def indices(self):
        return {self.model.UNITS: 'Generators', self.model.HOUR: 'Time', self.model.SCENARIOS: 'LMP Scenarios'}

    @property
    def default_bids(self):
        return {gen: self.model_data['Original Marginal Cost Curve'][gen] for gen in self.model.UNITS}

    @property
    def pmin(self):
        return {gen: self.model_data['Pmin'][gen] for gen in self.model.UNITS}

if __name__ == "__main__":

    generator = "102_STEAM_3"
    horizon = 4
    n_scenario = 3

    rts_gmlc_dataframe = pd.read_csv('gen.csv')
    thermal_generator_object = ThermalGenerator(rts_gmlc_dataframe = rts_gmlc_dataframe, \
                                                horizon = horizon, \
                                                generators = [generator], \
                                                n_scenario = n_scenario)

    solver = pyo.SolverFactory('cbc')

    run_bidder = True
    run_tracker = False

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
