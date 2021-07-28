import pyomo.environ as pyo
import pandas as pd
from collections import deque
from thermal_generator import ThermalGenerator

class ThermalGeneratorStorageIES:

    def __init__(self, \
                 rts_gmlc_dataframe, \
                 horizon = 48, \
                 n_scenario = 10, \
                 generators = None, \
                 storage_pmax_ratio = None,\
                 storage_size_hour = None,\
                 round_trip_efficiency = None,\
                 IES_mode = '1a'):

        '''
        Initializes the class object by building the IES model.

        Arguments:
            rts_gmlc_dataframe: the RTS-GMLC generator data in Pandas DataFrame
            horizon: the length of the planning horizon of the model.
            n_scenario: number of uncertain scenarios.
            generators: a list of generators in RTS-GMLC
            storage_pmax_ratio: a list of ratios between storage Pmax and thermal
                                generator Pmax.
            storage_size_hour: a list of storage size in hour.
            round_trip_efficiency: a list of storage system round-trip efficiency
            IES_mode: the operation mode of the IES, in {'1a': hide the storage,
                      i.e., cannot charge/discharge from/to market,
                      '1b': hide the storage, i.e., it can only charge from the
                      market, '2': storage and generators are independent,
                      '3a': storage cannot charge from the market}.

        Returns:
            None
        '''

        # create thermal generator object
        self.thermal_generator_object = ThermalGenerator(rts_gmlc_dataframe = rts_gmlc_dataframe,\
                                                         horizon = horizon,\
                                                         generators = generators,\
                                                         n_scenario = n_scenario)

        self.model = self.thermal_generator_object.model

        # add IES variable & constraints
        self.build_thermal_generator_storage_IES_model(generators = generators,\
                                                       storage_pmax_ratio = storage_pmax_ratio,\
                                                       storage_size_hour = storage_size_hour, \
                                                       round_trip_efficiency = round_trip_efficiency,\
                                                       IES_mode = IES_mode)

        self.result_list = []

    def _add_IES_params(self,model_data):
        '''
        This function adds the IES/storage constraints to the existing thermal
        generator model.

        Arguments:
            model_data: a dictionary which has this structure
            {data type name: {generator name: value}}. It contains the IES-related
            parameter values.

        Returns:
            None
        '''

        m = self.model

        # storage size in MWH
        m.storage_size = pyo.Param(m.UNITS,initialize=model_data['Storage Size'],mutable = False)

        # storage size in hr
        m.storage_size_hour = pyo.Param(m.UNITS,initialize=model_data['Storage Size Hour'],mutable = False)

        # initial soc of the storage
        def pre_soc_init(m,j):
            return m.storage_size[j]/2
        m.pre_SOC = pyo.Param(m.UNITS,initialize = pre_soc_init, mutable = True)

        # pmax of storage
        m.pmax_storage = pyo.Param(m.UNITS,initialize = model_data['Storage Pmax'], mutable = False)

        # storage efficiency
        m.round_trip_eff = pyo.Param(m.UNITS,initialize = model_data['Storage Round-trip Efficiency'], mutable = False)

        return

    def _add_IES_vars(self):

        '''
        This function adds the IES/storage variables to the existing thermal
        generator model.

        Arguments:
            None

        Returns:
            None
        '''

        m = self.model

        # define a rule for storage size
        def storage_power_bnd_rule(m,j,h,k):
            return (0,m.pmax_storage[j])

        def storage_size_bnd_rule(m,j,h,k):
            return (0,m.storage_size[j])

        # power to charge the storage by generator
        m.P_E = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,bounds = storage_power_bnd_rule)

        # power to the market by generator
        m.P_G = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.NonNegativeReals)

        # storage power (MW) to merge with thermal power
        m.P_S = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,bounds = storage_power_bnd_rule) #discharge

        # storage power (MW) interacting with the market
        m.P_D = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,bounds = storage_power_bnd_rule) #discharge
        m.P_C = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,bounds = storage_power_bnd_rule) #charge

        # storage SOC
        m.S_SOC = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,bounds = storage_size_bnd_rule)

        # total power to the grid from the generator side
        m.P_R = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.NonNegativeReals)

        # total power to the grid
        m.P_total = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.NonNegativeReals)

        # charge or discharge binary var
        m.y_S = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.Binary)
        m.y_E = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.Binary)
        m.y_D = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.Binary)
        m.y_C = pyo.Var(m.UNITS,m.HOUR,m.SCENARIOS,within = pyo.Binary)

        return

    def _add_IES_constraints(self):

        '''
        This function adds the IES/storage constraints to the existing thermal
        generator model.

        Arguments:
            None

        Returns:
            None
        '''

        m = self.model

        ## convex hull constraints on storage power
        def convex_hull_on_storage_pow_con_rules(m):
            for j in m.UNITS:
                for h in m.HOUR:
                    for k in m.SCENARIOS:
                        yield m.P_S[j,h,k]<= m.pmax_storage[j] * m.y_S[j,h,k]
                        yield m.P_E[j,h,k]<= m.pmax_storage[j] * m.y_E[j,h,k]
                        yield m.P_D[j,h,k]<= m.pmax_storage[j] * m.y_D[j,h,k]
                        yield m.P_C[j,h,k]<= m.pmax_storage[j] * m.y_C[j,h,k]
        m.UB_on_storage_pow_con = pyo.ConstraintList(rule = convex_hull_on_storage_pow_con_rules)

        def charge_or_discharge_fun3(m,j,h,k):
            return m.y_E[j,h,k] + m.y_S[j,h,k]<=1
        m.discharge_con3 = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = charge_or_discharge_fun3)

        def charge_or_discharge_fun4(m,j,h,k):
            return m.y_D[j,h,k] + m.y_C[j,h,k]<=1
        m.discharge_con4 = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = charge_or_discharge_fun4)

        def discharge_only_when_unit_on_fun(m,j,h,k):
            return 1-m.y_S[j,h,k] + m.on_off[j,h,k] >= 1
        m.discharge_only_when_unit_on = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = discharge_only_when_unit_on_fun)
        m.discharge_only_when_unit_on.deactivate()

        def charge_rate_exp(m,j,h,k):
            return m.P_E[j,h,k] + m.P_C[j,h,k]
        m.charge_rate = pyo.Expression(m.UNITS,m.HOUR,m.SCENARIOS,rule = charge_rate_exp)

        def discharge_rate_exp(m,j,h,k):
            return m.P_S[j,h,k] + m.P_D[j,h,k]
        m.discharge_rate = pyo.Expression(m.UNITS,m.HOUR,m.SCENARIOS,rule = discharge_rate_exp)

        # charging rate Constraints
        def charge_rate_fun(m,j,h,k):
            return m.charge_rate[j,h,k] <= m.pmax_storage[j]
        m.charge_rate_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = charge_rate_fun)

        def discharge_rate_fun(m,j,h,k):
            return m.discharge_rate[j,h,k] <= m.pmax_storage[j]
        m.discharge_rate_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = discharge_rate_fun)

        # thermal generator power balance
        def therm_pow_balance(m,j,h,k):
            return m.P_T[j,h,k] == m.P_E[j,h,k] + m.P_G[j,h,k]
        m.thermal_pow_balance_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = therm_pow_balance)

        # total power to the grid
        def total_gen_pow_G_fun(m,j,h,k):
            return m.P_R[j,h,k] == m.P_G[j,h,k] + m.P_S[j,h,k]
        m.total_pow_G_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = total_gen_pow_G_fun)

        # storage energy balance
        def EnergyBalance(m,j,h,k):
            if h == 0 :
                return m.S_SOC[j,h,k] == m.pre_SOC[j] + m.charge_rate[j,h,k]\
                *m.round_trip_eff[j]**0.5 - m.discharge_rate[j,h,k]/m.round_trip_eff[j]**0.5
            else :
                return m.S_SOC[j,h,k] == m.S_SOC[j,h-1,k] + m.charge_rate[j,h,k]\
                *m.round_trip_eff[j]**0.5 - m.discharge_rate[j,h,k]/m.round_trip_eff[j]**0.5
        m.EnergyBalance_Con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = EnergyBalance)

        # constrain state of charge at end time
        def SOC_endtime_fun(m,j,h,k):
            if h == len(m.HOUR) - 1:
                return  m.S_SOC[j,h,k] == m.pre_SOC[j]
            else:
                return pyo.Constraint.Skip
        m.SOC_endtime_con = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = SOC_endtime_fun)
        if len(m.HOUR)%24 != 0:
            m.SOC_endtime_con.deactivate()

        # energy generated by all the units
        def total_power_fun(m,j,h,k):
            return m.P_total[j,h,k] == m.P_R[j,h,k] + m.P_D[j,h,k] - m. P_C[j,h,k]
        m.tot_power = pyo.Constraint(m.UNITS,m.HOUR,m.SCENARIOS,rule = total_power_fun)

        return

    def _assemble_IES_data(self,
                           generators, \
                           storage_pmax_ratio, \
                           storage_size_hour, \
                           round_trip_efficiency):
        '''
        This function assembles the parameter data to build the IES model given
        a list of generator names.

        Arguments:
            generators: a list of generator names in RTS-GMLC dataset.
            storage_pmax_ratio: a list of ratios between storage Pmax and thermal
                                generator Pmax.
            storage_size_hour: a list of storage size in hour.
            round_trip_efficiency: a list of storage system round-trip efficiency
            IES_mode: the operation mode of the IES, in {'1a': hide the storage,
                      i.e., cannot charge/discharge from/to market,
                      '1b': hide the storage, i.e., it can only charge from the
                      market, '2': storage and generators are independent,
                      '3a': storage cannot charge from the market}.

        Returns:
            model_data: a dictionary which has this structure
            {data type name: {generator name: value}}.
        '''

        model_data = {}

        model_data['Storage Pmax Ratio'] = dict(zip(generators,storage_pmax_ratio))
        model_data['Storage Size Hour'] = dict(zip(generators,storage_size_hour))
        model_data['Storage Round-trip Efficiency'] = dict(zip(generators,round_trip_efficiency))

        model_data['Storage Pmax'] = {gen: pyo.value(self.model.Pmax[gen]) * model_data['Storage Pmax Ratio'][gen] for gen in generators}
        model_data['Storage Size'] = {gen: model_data['Storage Size Hour'][gen] * model_data['Storage Pmax'][gen] for gen in generators}
        model_data['Storage Size Ratio'] = {gen: model_data['Storage Size'][gen] / pyo.value(self.model.Pmax[gen]) for gen in generators}

        return model_data

    def switch_model_mode(self,mode = '1a'):

        '''
        Switch the operation mode of the IES.

        Arguments:
            mode: the operation mode of the IES, in {'1a': hide the storage,
                  i.e., cannot charge/discharge from/to market,
                  '1b': hide the storage, i.e., it can only charge from the
                  market, '2': storage and generators are independent,
                  '3a': storage cannot charge from the market}.

        Returns:
            None
        '''

        # tag the model by the mode
        m = self.model

        # first unfix all the vars we care, in case the model is switched before
        m.P_D.unfix()
        m.P_C.unfix()
        m.y_D.unfix()
        m.y_C.unfix()

        m.P_E.unfix()
        m.P_S.unfix()
        m.y_E.unfix()
        m.y_S.unfix()

        # hide the storage: cannot charge/discharge from/to market
        if mode == '1a' or mode == 'static':
            m.P_D.fix(0)
            m.P_C.fix(0)
            m.y_D.fix(0)
            m.y_C.fix(0)

        # hide the storage: it can only charge from the market
        elif mode == '1b':
            m.P_D.fix(0)
            m.y_D.fix(0)

        # storage and generators are independent
        elif mode == '2':
            m.P_E.fix(0)
            m.P_S.fix(0)
            m.y_E.fix(0)
            m.y_S.fix(0)

        # storage cannot charge from the market
        elif mode == '3a':
            m.P_C.fix(0)
            m.y_C.fix(0)

        return

    def build_thermal_generator_storage_IES_model(self, \
                                                  generators, \
                                                  storage_pmax_ratio, \
                                                  storage_size_hour, \
                                                  round_trip_efficiency,\
                                                  IES_mode):

        '''
        This function augments the given thermal generators to build the IES
        model.

        Arguments:
            generators: a list of generators in RTS-GMLC
            storage_pmax_ratio: a list of ratios between storage Pmax and thermal
                                generator Pmax.
            storage_size_hour: a list of storage size in hour.
            round_trip_efficiency: a list of storage system round-trip efficiency
            IES_mode: the operation mode of the IES, in {'1a': hide the storage,
                      i.e., cannot charge/discharge from/to market,
                      '1b': hide the storage, i.e., it can only charge from the
                      market, '2': storage and generators are independent,
                      '3a': storage cannot charge from the market}.

        Returns:
            None
        '''

        # prepare model data
        IES_data = self._assemble_IES_data(generators, \
                                           storage_pmax_ratio, \
                                           storage_size_hour, \
                                           round_trip_efficiency)

        self._add_IES_params(IES_data)
        self._add_IES_vars()
        self._add_IES_constraints()
        self.switch_model_mode(mode = IES_mode)

        return

    @staticmethod
    def _update_SOC(m,implemented_SOC):

        '''
        This method updates the parameters in the energy balance constraint
        based on the realized state-of-charge profile.

        Arguments:
            implemented_SOC: realized state-of-charge profile: {unit: []}

         Returns:
             None
        '''

        for unit in m.UNITS:
            m.pre_SOC[unit] = round(implemented_SOC[unit][-1],2)

        return

    def update_model(self, \
                     implemented_shut_down, \
                     implemented_start_up, \
                     implemented_power_output, \
                     implemented_SOC):

        '''
        This method updates the parameters in the model based on
        the implemented power outputs, state-of-charge, shut down and start up
        events.

        Arguments:
            implemented_power_output: realized power outputs: {unit: []}
            implemented_shut_down: realized shut down events: {unit: []}.
            implemented_start_up: realized start up events: {unit: []}
            implemented_SOC: realized state-of-charge profile: {unit: []}

         Returns:
             None
        '''

        m = self.model
        self.thermal_generator_object._update_UT_DT(m, implemented_shut_down, implemented_start_up)
        self.thermal_generator_object._update_power(m, implemented_power_output)
        self._update_SOC(m, implemented_SOC)

        return

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
                    result_dict['Scenario'] = int(k)

                    # model vars
                    result_dict['Total Power Output [MW]'] = float(round(pyo.value(m.P_total[generator,t,k]),2))
                    result_dict['Thermal Power Generated [MW]'] = float(round(pyo.value(m.P_T[generator,t,k]),2))
                    result_dict['Thermal Power to Storage [MW]'] = float(round(pyo.value(m.P_E[generator,t,k]),2))
                    result_dict['Thermal Power to Market [MW]'] = float(round(pyo.value(m.P_G[generator,t,k]),2))
                    result_dict['Storage Power to Thermal [MW]'] = float(round(pyo.value(m.P_S[generator,t,k]),2))
                    result_dict['Total Thermal Side Power to Market [MW]'] = float(round(pyo.value(m.P_R[generator,t,k]),2))
                    result_dict['Charge Power [MW]'] = float(round(pyo.value(m.P_C[generator,t,k]),2))
                    result_dict['Disharge Power [MW]'] = float(round(pyo.value(m.P_D[generator,t,k]),2))
                    result_dict['State of Charge [MWh]'] = float(round(pyo.value(m.S_SOC[generator,t,k]),2))

                    result_dict['On/off [bin]'] = int(round(pyo.value(m.on_off[generator,t,k])))
                    result_dict['Start Up [bin]'] = int(round(pyo.value(m.start_up[generator,t,k])))
                    result_dict['Shut Down [bin]'] = int(round(pyo.value(m.shut_dw[generator,t,k])))
                    result_dict['Storage to Thermal [bin]'] = int(round(pyo.value(m.y_S[generator,t,k])))
                    result_dict['Thermal to Storage [bin]'] = int(round(pyo.value(m.y_E[generator,t,k])))
                    result_dict['Charge [bin]'] = int(round(pyo.value(m.y_C[generator,t,k])))
                    result_dict['Dicharge [bin]'] = int(round(pyo.value(m.y_D[generator,t,k])))

                    result_dict['Production Cost [$]'] = float(round(pyo.value(m.prod_cost_approx[generator,t,k]),2))
                    result_dict['Start-up Cost [$]'] = float(round(pyo.value(m.start_up_cost_expr[generator,t,k]),2))
                    result_dict['Total Cost [$]'] = float(round(pyo.value(m.tot_cost[generator,t,k]),2))

                    # result_dict['Periodic Boundary Slack [MWh]'] = float(round(pyo.value(m.pbc_slack[generator,t]),2))
                    # result_dict['Power Output Slack [MW]'] = float(round(pyo.value(m.slack_var_power[generator,t]),2))

                    # calculate mileage
                    if t == 0:
                        result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[generator,t] - m.pre_P_T[generator])),2))
                    else:
                        result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[generator,t] - m.P_T[generator,t-1])),2))

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
        return self.model.P_total

    @property
    def production_cost(self):
        return self.model.tot_cost

    @property
    def indices(self):
        return {self.model.UNITS: 'Generators', self.model.HOUR: 'Time', self.model.SCENARIOS: 'LMP Scenarios'}

if __name__ == "__main__":

    rts_gmlc_dataframe = pd.read_csv('gen.csv')
    IES_object = ThermalGeneratorStorageIES(rts_gmlc_dataframe = rts_gmlc_dataframe, \
                                            horizon = 2, \
                                            n_scenario = 1,\
                                            generators = ["102_STEAM_3"], \
                                            storage_pmax_ratio = [0.1],\
                                            storage_size_hour = [4],\
                                            round_trip_efficiency = [0.88])

# def _check_number_of_None(ls):
#
#     cnt = 0
#
#     for item in ls:
#         if item is None:
#             cnt += 1
#
#     return cnt
#
# def assemble_model_data(generator_names, \
#                         storage_size_ratio, \
#                         round_trip_eff, \
#                         storage_size_hour,\
#                         storage_pmax_ratio,\
#                         rts_gmlc_data_dir,\
#                         **kwargs):
#
#     n_none = _check_number_of_None([storage_size_ratio,storage_size_hour,storage_pmax_ratio])
#
#     assert n_none == 1, \
#     'Please provide EXACTLY 2 of the 3 elements: storage capacity (storage_size_ratio), storage size in hour, and storage Pmax!'
#
#     # read data
#     gen_params = pd.read_csv(rts_gmlc_data_dir + 'gen.csv')
#     gen_params = gen_params[gen_params['GEN UID'].isin(generator_names)]
#
#     model_data = {}
#
#     # generator names
#     model_data['Generator'] = generator_names
#
#     # Pmin [MW]
#     model_data['Pmin'] = df_col_to_dict(df = gen_params,\
#                                         key_col_name = 'GEN UID',\
#                                         value_col_name = 'PMin MW')
#
#     # Pmax [MW]
#     model_data['Pmax'] = df_col_to_dict(df = gen_params,\
#                                         key_col_name = 'GEN UID',\
#                                         value_col_name = 'PMax MW')
#
#     # minimum up time [MW/hr]
#     model_data['UT'] = dict(zip(generator_names,gen_params['Min Up Time Hr'].astype(int)))
#
#     # minimum down time [hr]
#     model_data['DT'] = dict(zip(generator_names,gen_params['Min Down Time Hr'].astype(int)))
#
#     ## ramp rates [MW/hr]
#     ramp_rates = gen_params['Ramp Rate MW/Min'].values * 60
#
#     # ramp up rate [MW/hr]
#     model_data['RU'] = dict(zip(generator_names,ramp_rates))
#
#     # ramp down rate [MW/hr]
#     model_data['RD'] = dict(zip(generator_names,ramp_rates))
#
#     # ramp start up [MW/hr]
#     model_data['SU'] = {gen: min(model_data['Pmin'][gen],model_data['RU'][gen]) for gen in generator_names}
#
#     # ramp shut down [MW/hr] (will use pmin for now)
#     model_data['SD'] = {gen: min(model_data['Pmin'][gen],model_data['RD'][gen]) for gen in generator_names}
#
#     # start up cost [$/SU] (will use the warm start up cost for now)
#     start_up_cost = gen_params['Start Heat Warm MBTU'] * gen_params['Fuel Price $/MMBTU']
#     model_data['SU Cost'] = dict(zip(generator_names,start_up_cost))
#
#     ## production cost
#
#     # power segments and marginal costs
#     model_data['Power Segments'] = {}
#     model_data['Marginal Costs'] = {}
#     for gen in generator_names:
#         df = get_data_given(df = gen_params, generator = gen)
#         for l in range(1,4):
#             # power segements
#             model_data['Power Segments'][(gen,l)] = float(df['Output_pct_{}'.format(l)] * df['PMax MW'])
#
#             # segment marginal cost
#             model_data['Marginal Costs'][(gen,l)] = float(df['HR_incr_{}'.format(l)]/1000 * df['Fuel Price $/MMBTU'])
#
#     model_data['Min Load Cost'] = dict(zip(generator_names,gen_params['HR_avg_0']/1000 * gen_params['Fuel Price $/MMBTU'] * gen_params['PMin MW']))
#
#     if storage_size_ratio is None:
#         model_data['Storage Pmax Ratio'] = dict(zip(generator_names,storage_pmax_ratio))
#         model_data['Storage Size Hour'] = dict(zip(generator_names,storage_size_hour))
#         model_data['Storage Pmax'] = {gen: model_data['Pmax'][gen] * model_data['Storage Pmax Ratio'][gen] for gen in generator_names}
#         model_data['Storage Size'] = {gen: model_data['Storage Size Hour'][gen] * model_data['Storage Pmax'][gen] for gen in generator_names}
#         model_data['Storage Size Ratio'] = {gen: model_data['Storage Size'][gen] / model_data['Pmax'][gen] for gen in generator_names}
#
#     elif storage_size_hour is None:
#
#         # storage size ratio (1-based)
#         model_data['Storage Size Ratio'] = dict(zip(generator_names,storage_size_ratio))
#         model_data['Storage Pmax Ratio'] = dict(zip(generator_names,storage_pmax_ratio))
#         model_data['Storage Pmax'] = {gen: model_data['Pmax'][gen] * model_data['Storage Pmax Ratio'][gen] for gen in generator_names}
#         # storage size [MWh]
#         model_data['Storage Size'] = {gen: model_data['Pmax'][gen] * model_data['Storage Size Ratio'][gen] for gen in generator_names}
#         model_data['Storage Size Hour'] = {gen: model_data['Storage Size'][gen] / model_data['Storage Pmax'][gen] for gen in generator_names}
#
#     elif storage_pmax_ratio is None:
#
#         model_data['Storage Size Ratio'] = dict(zip(generator_names,storage_size_ratio))
#         # storage size [hr]
#         model_data['Storage Size Hour'] = dict(zip(generator_names,storage_size_hour))
#
#         # storage size [MW]
#         model_data['Storage Size'] = {gen: model_data['Pmax'][gen] * model_data['Storage Size Ratio'][gen] for gen in generator_names}
#         model_data['Storage Pmax'] = {gen: model_data['Storage Size'][gen] / model_data['Storage Size Hour'][gen] for gen in generator_names}
#         model_data['Storage Pmax Ratio'] = {gen: model_data['Storage Pmax'][gen] / model_data['Pmax'][gen] for gen in generator_names}
#
#     # storage size [MW]
#     model_data['Storage Round-trip Efficiency'] = dict(zip(generator_names,round_trip_eff))
#
#     for key in kwargs:
#         model_data[key] = {gen: kwargs[key] for gen in generator_names}
#
#     return model_data
#
# def _add_CA_UT_DT_constraints(m,plan_horizon):
#
#     # number of hours the plant has been on/off in the previous day
#     m.pre_up_hour = pyo.Param(m.UNITS,within = pyo.Integers,initialize = 24,mutable = True)
#     m.pre_dw_hour = pyo.Param(m.UNITS,within = pyo.Integers,initialize = 0, mutable = True)
#
#     # expressions for min up/down time
#     def G_expr(m,j):
#         # temp = max(0,value(m.min_up_time[j]-m.pre_up_hour[j]))
#         # return min(plan_horizon,temp*value(m.pre_on_off[j]))
#         return (m.min_up_time[j]-m.pre_up_hour[j])*m.pre_on_off[j]
#     m.G = pyo.Expression(m.UNITS,rule = G_expr)
#
#     def L_expr(m,j):
#         # temp = max(0,value(m.min_dw_time[j]-m.pre_dw_hour[j]))
#         # return min(plan_horizon,temp*(1-value(m.pre_on_off[j])))
#         return (m.min_dw_time[j]-m.pre_dw_hour[j])*(1-m.pre_on_off[j])
#     m.L = pyo.Expression(m.UNITS,rule = L_expr)
#
#     ## min up time constraints
#     def minimum_on_fun(m,j):
#         if pyo.value(m.G[j]) <= 0:
#             return pyo.Constraint.Skip
#         elif 0<pyo.value(m.G[j])<= plan_horizon:
#             return sum(1-m.on_off[j,h] for h in range(int(pyo.value(m.G[j])))) == 0
#         else:
#             return sum(1-m.on_off[j,h] for h in range(plan_horizon)) == 0
#     m.minimum_on_con = pyo.Constraint(m.UNITS,rule = minimum_on_fun)
#
#     # minimum up time in the subsequent hours
#     def sub_minimum_on_fun(m):
#         for j in m.UNITS:
#             if pyo.value(m.G[j])<=0:
#                 low_bnd = 0
#             else:
#                 low_bnd = int(pyo.value(m.G[j]))
#
#             for h0 in range(low_bnd,plan_horizon-m.min_up_time[j]+1):
#                 yield sum(m.on_off[j,h] for h in range(h0,h0 + m.min_up_time[j]))\
#                 >= m.min_up_time[j]*m.start_up[j,h0]
#     m.sub_minimum_on_con = pyo.ConstraintList(rule = sub_minimum_on_fun)
#
#     # minimum up time in the final hours
#     def fin_minimum_on_fun(m):
#         for j in m.UNITS:
#
#             if plan_horizon < pyo.value(m.min_up_time[j]):
#                 # if goes into here, sub_minimum_on_con must have been
#                 # skipped.
#                 if pyo.value(m.G[j]) <= 0:
#                     # if goes into here, minimum_on_con must have been
#                     # skipped
#                     for h0 in range(0,plan_horizon):
#                         yield sum(m.on_off[j,h] - m.start_up[j,h0] \
#                         for h in range(h0,plan_horizon)) >=0
#
#                 else:
#                     for h0 in range(int(pyo.value(m.G[j])),plan_horizon):
#                         yield sum(m.on_off[j,h] - m.start_up[j,h0] \
#                         for h in range(h0,plan_horizon)) >=0
#
#             else:
#                 for h0 in range(plan_horizon-m.min_up_time[j]+1,plan_horizon):
#                     yield sum(m.on_off[j,h] - m.start_up[j,h0] \
#                     for h in range(h0,plan_horizon)) >=0
#     m.fin_minimum_on_con = pyo.ConstraintList(rule = fin_minimum_on_fun)
#
#     ## min down time constraints
#     def minimum_off_fun(m,j):
#         if pyo.value(m.L[j])<= 0:
#             return pyo.Constraint.Skip
#         elif 0<pyo.value(m.L[j])<= plan_horizon:
#             return sum(m.on_off[j,h] for h in range(int(pyo.value(m.L[j])))) == 0
#         else:
#             return sum(m.on_off[j,h] for h in range(plan_horizon)) == 0
#     m.minimum_off_con = pyo.Constraint(m.UNITS,rule = minimum_off_fun)
#
#     # minimum down time in the subsequent hours
#     def sub_minimum_dw_fun(m):
#         for j in m.UNITS:
#             if pyo.value(m.L[j])<=0:
#                 low_bnd = 0
#             else:
#                 low_bnd = int(pyo.value(m.L[j]))
#             for h0 in range(low_bnd,plan_horizon-m.min_dw_time[j]+1):
#                 yield sum(1-m.on_off[j,h] for h in range(h0,h0 + m.min_dw_time[j]))\
#                 >= m.min_dw_time[j]* m.shut_dw[j,h0]
#     m.sub_minimum_dw_con = pyo.ConstraintList(rule = sub_minimum_dw_fun)
#
#     # minimum down time in the final hours
#     def fin_minimum_dw_fun(m):
#         for j in m.UNITS:
#
#             if plan_horizon < pyo.value(m.min_dw_time[j]):
#                 # if goes into here, sub_minimum_dw_con must have been
#                 # skipped.
#                 if pyo.value(m.L[j]) <= 0:
#                     # if goes into here, minimum_off_con must have been
#                     # skipped
#                     for h0 in range(0,plan_horizon):
#                         yield sum(1-m.on_off[j,h] - m.shut_dw[j,h0]\
#                         for h in range(h0,plan_horizon)) >=0
#                 else:
#                     for h0 in range(int(pyo.value(m.L[j])),plan_horizon):
#                         yield sum(1-m.on_off[j,h] - m.shut_dw[j,h0]\
#                         for h in range(h0,plan_horizon)) >=0
#
#             else:
#                 for h0 in range(plan_horizon-m.min_dw_time[j]+1,plan_horizon):
#                     yield sum(1-m.on_off[j,h] - m.shut_dw[j,h0]\
#                     for h in range(h0,plan_horizon)) >=0
#     m.fin_minimum_dw_con = pyo.ConstraintList(rule = fin_minimum_dw_fun)
#
#     return
#
# def _add_UT_DT_constraints(m):
#
#     def pre_shut_down_trajectory_set_rule(m):
#         return ((j,t) for j in m.UNITS for t in range(-pyo.value(m.min_dw_time[j]) + 1,0))
#     m.pre_shut_down_trajectory_set = pyo.Set(dimen = 2,initialize = pre_shut_down_trajectory_set_rule, ordered = True)
#
#     def pre_start_up_trajectory_set_rule(m):
#         return ((j,t) for j in m.UNITS for t in range(-pyo.value(m.min_up_time[j]) + 1,0))
#     m.pre_start_up_trajectory_set = pyo.Set(dimen = 2,initialize = pre_start_up_trajectory_set_rule, ordered = True)
#
#     m.pre_shut_down_trajectory = pyo.Param(m.pre_shut_down_trajectory_set, initialize = 0, mutable = True)
#     m.pre_start_up_trajectory = pyo.Param(m.pre_start_up_trajectory_set, initialize = 0, mutable = True)
#
#     def min_down_time_rule(m,j,h):
#         if h < pyo.value(m.min_dw_time[j]):
#             return sum(m.pre_shut_down_trajectory[j,h0] for h0 in range(h - pyo.value(m.min_dw_time[j]) + 1,0)) \
#                    + sum(m.shut_dw[j,h0] for h0 in range(h + 1)) <= 1 - m.on_off[j,h]
#         else:
#             return sum(m.shut_dw[j,h0] for h0 in range(h - pyo.value(m.min_dw_time[j]) + 1, h + 1)) <= 1 - m.on_off[j,h]
#     m.min_down_time_con = pyo.Constraint(m.UNITS,m.HOUR,rule = min_down_time_rule)
#
#     def min_up_time_rule(m,j,h):
#         if h < pyo.value(m.min_up_time[j]):
#             return sum(m.pre_start_up_trajectory[j,h0] for h0 in range(h - pyo.value(m.min_up_time[j]) + 1,0)) \
#                    + sum(m.start_up[j,h0] for h0 in range(h + 1)) <= m.on_off[j,h]
#         else:
#             return sum(m.start_up[j,h0] for h0 in range(h - pyo.value(m.min_up_time[j]) + 1, h + 1)) <= m.on_off[j,h]
#     m.min_up_time_con = pyo.Constraint(m.UNITS,m.HOUR,rule = min_up_time_rule)
#
#     return
#
# def _add_slack_vars(m):
#
#     def storage_size_bnd_rule(m,j,h):
#         return (0,m.storage_size[j])
#
#     # slack vars for storage periodic boundary constraint
#     m.pbc_slack = pyo.Var(m.UNITS,m.HOUR, bounds = storage_size_bnd_rule)
#     m.pbc_slack_positive = pyo.Var(m.UNITS,m.HOUR, bounds = storage_size_bnd_rule)
#     m.pbc_slack_negative = pyo.Var(m.UNITS,m.HOUR, bounds = storage_size_bnd_rule)
#     m.pbc_slack.fix(0)
#     m.pbc_slack_positive.fix(0)
#     m.pbc_slack_negative.fix(0)
#
#     # slack variables for power dispatch
#     m.slack_var_power = pyo.Var(m.UNITS,m.HOUR,within = pyo.NonNegativeReals)
#     m.slack_var_power.fix(0)
#
#     return
#
# def _add_min_slack_objective(m):
#
#     def min_slack_fun(m):
#         return sum(m.slack_var_power[j,h] + m.pbc_slack_positive[j,h] + \
#                m.pbc_slack_negative[j,h] for j in m.UNITS for h in m.HOUR)
#     m.min_slack = pyo.Objective(rule = min_slack_fun, sense = pyo.minimize)
#     m.min_slack.deactivate()
#
#     return
#
# def _add_slack_constraints(m):
#
#     plan_horizon = len(m.HOUR)
#
#     ## add constraints for power dispatch
#     def abs_con_fun1(m,j,h):
#         return -m.slack_var_power[j,h]<=m.P_total[j,h] - m.power_dispatch[j,h]
#     m.power_slack_con1 = pyo.Constraint(m.UNITS,m.HOUR,rule = abs_con_fun1)
#     m.power_slack_con1.deactivate()
#
#     def abs_con_fun2(m,j,h):
#         return m.slack_var_power[j,h]>=m.P_total[j,h] - m.power_dispatch[j,h]
#     m.power_slack_con2 = pyo.Constraint(m.UNITS,m.HOUR,rule = abs_con_fun2)
#     m.power_slack_con2.deactivate()
#
#     ## constraints for pbc slacks
#     def pbc_slack_fun(m,j,h):
#         return m.pbc_slack[j,h] == m.pbc_slack_positive[j,h] - m.pbc_slack_negative[j,h]
#     m.pbc_slack_con = pyo.Constraint(m.UNITS,m.HOUR,rule = pbc_slack_fun)
#     m.pbc_slack_con.deactivate()
#
#     def SOC_endtime_slack_fun(m,j,h):
#         if h == plan_horizon - 1:
#             return  m.S_SOC[j,h] == m.pre_SOC[j] + m.pbc_slack[j,h]
#         else:
#             return pyo.Constraint.Skip
#     m.SOC_endtime_slack_con = pyo.Constraint(m.UNITS,m.HOUR,rule = SOC_endtime_slack_fun)
#     m.SOC_endtime_slack_con.deactivate()
#
#     return
#
# def _add_slack_bounding_constraints(m):
#
#     # add minimum pbc slack as a param
#     m.pbc_slack_bound = pyo.Param(m.UNITS, initialize = 0, mutable = True)
#
#     def bound_pbc_slack_fun(m):
#         for j in m.UNITS:
#             yield sum(m.pbc_slack[j,h] for h in m.HOUR) <= 1.01 * m.pbc_slack_bound[j]
#             yield sum(m.pbc_slack[j,h] for h in m.HOUR) >= -1.01 * m.pbc_slack_bound[j]
#     m.bound_pbc_slack_con = pyo.ConstraintList(rule = bound_pbc_slack_fun)
#     m.bound_pbc_slack_con.deactivate()
#
#     # add minimum power output slack as a param
#     m.power_slack_bound = pyo.Param(m.UNITS, m.HOUR, initialize = 0, mutable = True)
#
#     def bound_power_slack_fun(m):
#         for j in m.UNITS:
#             for h in m.HOUR:
#                 yield m.slack_var_power[j,h]  <= 1.01 * m.power_slack_bound[j,h]
#                 yield m.slack_var_power[j,h]  >= -1.01 * m.power_slack_bound[j,h]
#     m.bound_power_slack_con = pyo.ConstraintList(rule = bound_power_slack_fun)
#     m.bound_power_slack_con.deactivate()
#
#     return
#
# def _build_model(model_data,plan_horizon,segment_number=4):
#
#     m = pyo.ConcreteModel()
#
#     ## define the sets
#     m.HOUR = pyo.Set(initialize = range(plan_horizon))
#     m.SEGMENTS = pyo.Set(initialize = range(1, segment_number))
#     m.UNITS = pyo.Set(initialize = model_data['Generator'], ordered = True)
#
#     ## define the parameters
#
#     # add power schedule as a Param
#     m.power_dispatch = pyo.Param(m.UNITS,m.HOUR, initialize = 0, mutable = True)
#
#     m.DAM_price = pyo.Param(m.HOUR,initialize = 20, mutable = True)
#
#     m.start_up_cost = pyo.Param(m.UNITS,initialize = model_data['SU Cost'],mutable = False)
#
#     # capacity of generators: upper bound (MW)
#     m.Pmax = pyo.Param(m.UNITS,initialize = model_data['Pmax'], mutable = False)
#
#     # minimum power of generators: lower bound (MW)
#     m.Pmin = pyo.Param(m.UNITS,initialize = model_data['Pmin'], mutable = False)
#
#     m.power_segment_bounds = pyo.Param(m.UNITS,m.SEGMENTS,initialize = model_data['Power Segments'], mutable = False)
#
#     # get the cost slopes
#     m.F = pyo.Param(m.UNITS,m.SEGMENTS,initialize = model_data['Marginal Costs'], mutable = False)
#
#     m.min_load_cost = pyo.Param(m.UNITS,initialize = model_data['Min Load Cost'], mutable = False)
#
#     # Ramp up limits (MW/h)
#     m.ramp_up = pyo.Param(m.UNITS,initialize = model_data['RU'], mutable = False)
#
#     # Ramp down limits (MW/h)
#     m.ramp_dw = pyo.Param(m.UNITS,initialize = model_data['RD'], mutable = False)
#
#     # start up ramp limit
#     m.ramp_start_up = pyo.Param(m.UNITS,initialize = model_data['SU'], mutable = False)
#
#     # shut down ramp limit
#     m.ramp_shut_dw = pyo.Param(m.UNITS,initialize = model_data['SD'], mutable = False)
#
#     # minimum down time [hr]
#     m.min_dw_time = pyo.Param(m.UNITS,initialize = model_data['DT'], mutable = False)
#
#     # minimum up time [hr]
#     m.min_up_time = pyo.Param(m.UNITS,initialize = model_data['UT'], mutable = False)
#
#     # power from the previous day (MW)
#     # need to assume the power output is at least always at the minimum pow output
#
#     # on/off status from previous day
#     m.pre_on_off = pyo.Param(m.UNITS,within = pyo.Binary,default= 1,mutable = True)
#
#     # define a function to initialize the previous power params
#     def init_pre_pow_fun(m,j):
#         return m.pre_on_off[j]*m.Pmin[j]
#     m.pre_P_T = pyo.Param(m.UNITS,initialize = init_pre_pow_fun, mutable = True)
#
#     # storage size in MWH
#     m.storage_size = pyo.Param(m.UNITS,initialize=model_data['Storage Size'],mutable = False)
#
#     # storage size in hr
#     m.storage_size_hour = pyo.Param(m.UNITS,initialize=model_data['Storage Size Hour'],mutable = False)
#
#     # initial soc of the storage
#     def pre_soc_init(m,j):
#         return m.storage_size[j]/2
#     m.pre_SOC = pyo.Param(m.UNITS,initialize = pre_soc_init, mutable = True)
#
#     # pmax of storage
#     m.pmax_storage = pyo.Param(m.UNITS,initialize = model_data['Storage Pmax'], mutable = False)
#
#     # storage efficiency
#     m.round_trip_eff = pyo.Param(m.UNITS,initialize = model_data['Storage Round-trip Efficiency'], mutable = False)
#
#     ## define the variables
#
#     # define a rule for storage size
#     def storage_power_bnd_rule(m,j,h):
#         return (0,m.pmax_storage[j])
#
#     def storage_size_bnd_rule(m,j,h):
#         return (0,m.storage_size[j])
#
#     # generator power (MW)
#
#     # power generated by thermal generator
#     m.P_T = pyo.Var(m.UNITS,m.HOUR,within = pyo.NonNegativeReals)
#
#     # power to charge the storage by generator
#     m.P_E = pyo.Var(m.UNITS,m.HOUR,bounds = storage_power_bnd_rule)
#
#     # power to the market by generator
#     m.P_G = pyo.Var(m.UNITS,m.HOUR,within = pyo.NonNegativeReals)
#
#     # storage power (MW) to merge with thermal power
#     m.P_S = pyo.Var(m.UNITS,m.HOUR,bounds = storage_power_bnd_rule) #discharge
#
#     # storage power (MW) interacting with the market
#     m.P_D = pyo.Var(m.UNITS,m.HOUR,bounds = storage_power_bnd_rule) #discharge
#     m.P_C = pyo.Var(m.UNITS,m.HOUR,bounds = storage_power_bnd_rule) #charge
#
#     # storage SOC
#     m.S_SOC = pyo.Var(m.UNITS,m.HOUR,bounds = storage_size_bnd_rule)
#
#     # total power to the grid from the generator side
#     m.P_R = pyo.Var(m.UNITS,m.HOUR,within = pyo.NonNegativeReals)
#
#     # total power to the grid
#     m.P_total = pyo.Var(m.UNITS,m.HOUR,within = pyo.NonNegativeReals)
#
#     # charge or discharge binary var
#     m.y_S = pyo.Var(m.UNITS,m.HOUR,within = pyo.Binary)
#     m.y_E = pyo.Var(m.UNITS,m.HOUR,within = pyo.Binary)
#     m.y_D = pyo.Var(m.UNITS,m.HOUR,within = pyo.Binary)
#     m.y_C = pyo.Var(m.UNITS,m.HOUR,within = pyo.Binary)
#
#     # binary variables indicating on/off
#     m.on_off = pyo.Var(m.UNITS,m.HOUR,initialize = True, within = pyo.Binary)
#
#     # binary variables indicating  start_up
#     m.start_up = pyo.Var(m.UNITS,m.HOUR,initialize = False, within = pyo.Binary)
#
#     # binary variables indicating shut down
#     m.shut_dw = pyo.Var(m.UNITS,m.HOUR,initialize = False, within = pyo.Binary)
#
#     # power produced in each segment
#     m.power_segment = pyo.Var(m.UNITS,m.HOUR,m.SEGMENTS, within = pyo.NonNegativeReals)
#
#     ## Constraints
#
#     # fix power output at market dispatch
#     def fix_power_output_fun(m,j,h):
#         return m.P_total[j,h] == m.power_dispatch[j,h]
#     m.fix_power_output_con = pyo.Constraint(m.UNITS,m.HOUR,rule = fix_power_output_fun)
#
#     ## convex hull constraints on storage power
#     def convex_hull_on_storage_pow_con_rules(m):
#         for j in m.UNITS:
#             for h in m.HOUR:
#                 yield m.P_S[j,h]<= m.pmax_storage[j] * m.y_S[j,h]
#                 yield m.P_E[j,h]<= m.pmax_storage[j] * m.y_E[j,h]
#                 yield m.P_D[j,h]<= m.pmax_storage[j] * m.y_D[j,h]
#                 yield m.P_C[j,h]<= m.pmax_storage[j] * m.y_C[j,h]
#     m.UB_on_storage_pow_con = pyo.ConstraintList(rule = convex_hull_on_storage_pow_con_rules)
#
#     def charge_or_discharge_fun3(m,j,h):
#         return m.y_E[j,h] + m.y_S[j,h]<=1
#     m.discharge_con3 = pyo.Constraint(m.UNITS,m.HOUR,rule = charge_or_discharge_fun3)
#
#     def charge_or_discharge_fun4(m,j,h):
#         return m.y_D[j,h] + m.y_C[j,h]<=1
#     m.discharge_con4 = pyo.Constraint(m.UNITS,m.HOUR,rule = charge_or_discharge_fun4)
#
#     def discharge_only_when_unit_on_fun(m,j,h):
#         return 1-m.y_S[j,h] + m.on_off[j,h] >= 1
#     m.discharge_only_when_unit_on = pyo.Constraint(m.UNITS,m.HOUR,rule = discharge_only_when_unit_on_fun)
#     m.discharge_only_when_unit_on.deactivate()
#
#     def charge_rate_exp(m,j,h):
#         return m.P_E[j,h] + m.P_C[j,h]
#     m.charge_rate = pyo.Expression(m.UNITS,m.HOUR,rule =charge_rate_exp)
#
#     def discharge_rate_exp(m,j,h):
#         return m.P_S[j,h] + m.P_D[j,h]
#     m.discharge_rate = pyo.Expression(m.UNITS,m.HOUR,rule =discharge_rate_exp)
#
#     # charging rate Constraints
#     def charge_rate_fun(m,j,h):
#         return m.charge_rate[j,h] <= m.pmax_storage[j]
#     m.charge_rate_con = pyo.Constraint(m.UNITS,m.HOUR,rule =charge_rate_fun)
#
#     def discharge_rate_fun(m,j,h):
#         return m.discharge_rate[j,h] <= m.pmax_storage[j]
#     m.discharge_rate_con = pyo.Constraint(m.UNITS,m.HOUR,rule =discharge_rate_fun)
#
#     # bounds on gen_pow
#     def lhs_bnd_gen_pow_fun(m,j,h):
#         return m.on_off[j,h] * m.Pmin[j] <= m.P_T[j,h]
#     m.lhs_bnd_gen_pow = pyo.Constraint(m.UNITS,m.HOUR,rule = lhs_bnd_gen_pow_fun)
#
#     def rhs_bnd_gen_pow_fun(m,j,h):
#         return m.P_T[j,h] <= m.on_off[j,h] * m.Pmax[j]
#     m.rhs_bnd_gen_pow = pyo.Constraint(m.UNITS,m.HOUR,rule = rhs_bnd_gen_pow_fun)
#
#     # thermal generator power balance
#     def therm_pow_balance(m,j,h):
#         return m.P_T[j,h] == m.P_E[j,h] + m.P_G[j,h]
#     m.thermal_pow_balance_con = pyo.Constraint(m.UNITS,m.HOUR,rule =therm_pow_balance)
#
#     # total power to the grid
#     def total_gen_pow_G_fun(m,j,h):
#         return m.P_R[j,h] == m.P_G[j,h] + m.P_S[j,h]
#     m.total_pow_G_con = pyo.Constraint(m.UNITS,m.HOUR,rule =total_gen_pow_G_fun)
#
#     # storage energy balance
#     def EnergyBalance(m,j,h):
#         if h == 0 :
#             return m.S_SOC[j,h] == m.pre_SOC[j] + m.charge_rate[j,h]\
#             *m.round_trip_eff[j]**0.5 - m.discharge_rate[j,h]/m.round_trip_eff[j]**0.5
#         else :
#             return m.S_SOC[j,h] == m.S_SOC[j,h-1] + m.charge_rate[j,h]\
#             *m.round_trip_eff[j]**0.5 - m.discharge_rate[j,h]/m.round_trip_eff[j]**0.5
#     m.EnergyBalance_Con = pyo.Constraint(m.UNITS,m.HOUR,rule = EnergyBalance)
#
#     # constrain state of charge at end time
#     def SOC_endtime_fun(m,j,h):
#         if h == plan_horizon - 1:
#             return  m.S_SOC[j,h] == m.pre_SOC[j]
#         else:
#             return pyo.Constraint.Skip
#     m.SOC_endtime_con = pyo.Constraint(m.UNITS,m.HOUR,rule = SOC_endtime_fun)
#     if plan_horizon%24 != 0:
#         m.SOC_endtime_con.deactivate()
#
#     # linearized power
#     def linear_power_fun(m,j,h):
#         return m.P_T[j,h] == \
#         sum(m.power_segment[j,h,l] for l in m.SEGMENTS) + m.Pmin[j]*m.on_off[j,h]
#     m.linear_power = pyo.Constraint(m.UNITS,m.HOUR,rule = linear_power_fun)
#
#     # bounds on segment power
#     def seg_pow_bnd_fun(m,j,h,l):
#         if l == 1:
#             return m.power_segment[j,h,l]<= (m.power_segment_bounds[j,l] - m.Pmin[j]) * m.on_off[j,h]
#         else:
#             return m.power_segment[j,h,l]<= (m.power_segment_bounds[j,l] - m.power_segment_bounds[j,l-1]) * m.on_off[j,h]
#     m.seg_pow_bnd = pyo.Constraint(m.UNITS,m.HOUR,m.SEGMENTS,rule = seg_pow_bnd_fun)
#
#     # start up and shut down logic (Arroyo and Conejo 2000)
#     def start_up_shut_dw_fun(m,j,h):
#         if h == 0:
#             return m.start_up[j,h] - m.shut_dw[j,h] == m.on_off[j,h] - m.pre_on_off[j]
#         else:
#             return m.start_up[j,h] - m.shut_dw[j,h] == m.on_off[j,h] - m.on_off[j,h-1]
#     m.start_up_shut_dw = pyo.Constraint(m.UNITS,m.HOUR,rule = start_up_shut_dw_fun)
#
#     # either start up or shut down
#     def start_up_or_shut_dw_fun(m,j,h):
#         return m.start_up[j,h] + m.shut_dw[j,h] <= 1
#     m.start_up_or_shut_dw = pyo.Constraint(m.UNITS,m.HOUR,rule = start_up_or_shut_dw_fun)
#
#     # ramp up limits
#     def ramp_up_fun(m,j,h):
#         '''
#         j,h stand for unit, hour,scenario respectively.
#         '''
#         if h==0:
#             return m.P_T[j,h] <= m.pre_P_T[j] \
#             + m.ramp_up[j]*m.pre_on_off[j]\
#             + m.ramp_start_up[j]*m.start_up[j,h]
#         else:
#             return m.P_T[j,h] <= m.P_T[j,h-1] \
#             + m.ramp_up[j]*m.on_off[j,h-1]\
#             + m.ramp_start_up[j]*m.start_up[j,h]
#     m.ramp_up_con = pyo.Constraint(m.UNITS,m.HOUR,rule = ramp_up_fun)
#
#     # ramp shut down limits
#     def ramp_shut_dw_fun(m,j,h):
#         '''
#         j,h stand for unit, hour,scenario respectively.
#         '''
#         if h==0:
#             return m.pre_P_T[j] <= m.Pmax[j]*m.on_off[j,h] + m.ramp_shut_dw[j] * m.shut_dw[j,h]
#         else:
#             return m.P_T[j,h-1] <= m.Pmax[j]*m.on_off[j,h] + m.ramp_shut_dw[j] * m.shut_dw[j,h]
#     m.ramp_shut_dw_con = pyo.Constraint(m.UNITS,m.HOUR,rule = ramp_shut_dw_fun)
#
#     # ramp down limits
#     def ramp_dw_fun(m,j,h):
#         '''
#         j,h stand for unit, hour,scenario respectively.
#         '''
#         if h == 0:
#             return m.pre_P_T[j] - m.P_T[j,h] <= m.ramp_dw[j] * m.on_off[j,h]\
#             + m.ramp_shut_dw[j] * m.shut_dw[j,h]
#         else:
#             return m.P_T[j,h-1] - m.P_T[j,h] <= m.ramp_dw[j] * m.on_off[j,h]\
#             + m.ramp_shut_dw[j] * m.shut_dw[j,h]
#     m.ramp_dw_con = pyo.Constraint(m.UNITS,m.HOUR,rule = ramp_dw_fun)
#
#     ## add min up and down time constraints
#     # add_CA_UT_DT_constraints(m,plan_horizon)
#     _add_UT_DT_constraints(m)
#
#     ## Expression
#     # energy generated by all the units
#     def total_power_fun(m,j,h):
#         return m.P_total[j,h] == m.P_R[j,h] + m.P_D[j,h] - m. P_C[j,h]
#     m.tot_power = pyo.Constraint(m.UNITS,m.HOUR,rule = total_power_fun)
#
#     def prod_cost_fun(m,j,h):
#         return m.min_load_cost[j] * m.on_off[j,h] \
#         + sum(m.F[j,l]*m.power_segment[j,h,l] for l in m.SEGMENTS)
#     m.prod_cost_approx = pyo.Expression(m.UNITS,m.HOUR,rule = prod_cost_fun)
#
#     # start up costs
#     def start_cost_fun(m,j,h):
#         return m.start_up_cost[j]*m.start_up[j,h]
#     m.start_up_cost_expr = pyo.Expression(m.UNITS,m.HOUR,rule = start_cost_fun)
#
#     # total cost
#     def tot_cost_fun(m,j,h):
#         return m.prod_cost_approx[j,h] + m.start_up_cost_expr[j,h]
#     m.tot_cost = pyo.Expression(m.UNITS,m.HOUR,rule = tot_cost_fun)
#
#     ## Objective
#     def exp_revenue_fun(m):
#         return sum(m.P_total[j,h]*m.DAM_price[h]- m.tot_cost[j,h] for h in m.HOUR for j in m.UNITS)
#     m.exp_revenue = pyo.Objective(rule = exp_revenue_fun,sense = pyo.maximize)
#     m.exp_revenue.deactivate()
#
#     return m
#
# def _add_tracking_constraints_and_objectives(m, ramping_weight = 20):
#
#     # add params for objective weights
#     m.cost_weight = pyo.Param(initialize = 1, mutable = False)
#     m.ramping_weight = pyo.Param(initialize = ramping_weight, mutable = False)
#
#     ## USE L1 norm
#     # add a slack variable for ramping
#     m.slack_var_ramp = pyo.Var(m.UNITS,m.HOUR,within = pyo.NonNegativeReals)
#
#     # add constraints for ramping
#     def abs_ramp_con_fun1(m,j,h):
#         if h == 0:
#             return -m.slack_var_ramp[j,h]<=m.P_T[j,h] - m.pre_P_T[j]
#         else:
#             return -m.slack_var_ramp[j,h]<=m.P_T[j,h] - m.P_T[j,h-1]
#     m.abs_ramp_con1 = pyo.Constraint(m.UNITS,m.HOUR,rule = abs_ramp_con_fun1)
#
#     def abs_ramp_con_fun2(m,j,h):
#         if h == 0:
#             return m.slack_var_ramp[j,h]>=m.P_T[j,h] - m.pre_P_T[j]
#         else:
#             return m.slack_var_ramp[j,h]>=m.P_T[j,h] - m.P_T[j,h-1]
#     m.abs_ramp_con2 = pyo.Constraint(m.UNITS,m.HOUR,rule = abs_ramp_con_fun2)
#
#     def power_abs_least(m):
#         return m.ramping_weight * sum(m.slack_var_ramp[j,h] for h in m.HOUR for j in m.UNITS)\
#                + m.cost_weight * sum(m.tot_cost[j,h] for h in m.HOUR for j in m.UNITS)
#     m.abs_dev = pyo.Objective(rule = power_abs_least,sense = pyo.minimize)
#
#     return
#
# def build_tracking_model(model_data,plan_horizon,ramping_weight = 20,segment_number=4,model_mode = '1a'):
#
#     m = _build_model(model_data,plan_horizon,segment_number)
#
#     _add_tracking_constraints_and_objectives(m, ramping_weight = ramping_weight)
#     _add_slack_vars(m)
#     _add_slack_constraints(m)
#     _add_min_slack_objective(m)
#     _add_slack_bounding_constraints(m)
#
#     switch_model_mode(m,mode = model_mode)
#
#     return m
#
# def switch_model_mode(m,mode = '1a'):
#
#     # tag the model by the mode
#     m.mode = mode
#
#     # first unfix all the vars we care, in case the model is switched before
#     m.P_D.unfix()
#     m.P_C.unfix()
#     m.y_D.unfix()
#     m.y_C.unfix()
#
#     m.P_E.unfix()
#     m.P_S.unfix()
#     m.y_E.unfix()
#     m.y_S.unfix()
#
#     # hide the storage: cannot charge/discharge from/to market
#     if mode == '1a' or mode == 'static':
#         m.P_D.fix(0)
#         m.P_C.fix(0)
#         m.y_D.fix(0)
#         m.y_C.fix(0)
#
#     # hide the storage: it can only charge from the market
#     elif mode == '1b':
#         m.P_D.fix(0)
#         m.y_D.fix(0)
#
#     # storage and generators are independent
#     elif mode == '2':
#         m.P_E.fix(0)
#         m.P_S.fix(0)
#         m.y_E.fix(0)
#         m.y_S.fix(0)
#
#     # storage cannot charge from the market
#     elif mode == '3a':
#         m.P_C.fix(0)
#         m.y_C.fix(0)
#
#     return m
#
# def record_results(m, time, market_signals, **kwargs):
#
#     '''
#     market_signals: {generator: {DA_dispatches: [], RT_dispatches; [], DA_LMP: [], RT_LMP; []}}
#     use queues, so we can pop the first one
#     '''
#
#     df_list = []
#     for generator in m.UNITS:
#
#         # upack things
#         DA_dispatches = market_signals[generator]['DA_dispatches']
#         RT_dispatches = market_signals[generator]['RT_dispatches']
#         DA_LMP = market_signals[generator]['DA_LMP']
#         RT_LMP = market_signals[generator]['RT_LMP']
#
#         for t in m.HOUR:
#
#             result_dict = {}
#             result_dict['Generator'] = generator
#             result_dict['Time [hr]'] = time
#
#             # simulation inputs
#             result_dict['Horizon [hr]'] = int(t)
#             result_dict['DA Dispatch [MW]'] = float(round(DA_dispatches[t],2))
#             result_dict['RT Dispatch [MW]'] = float(round(RT_dispatches[t],2))
#             result_dict['DA Commitment [bin]'] = int(result_dict['DA Dispatch [MW]'] > 0)
#             result_dict['RT Commitment [bin]'] = int(result_dict['RT Dispatch [MW]'] > 0)
#             result_dict['DA LMP [$/MWh]'] = float(round(DA_LMP[t],2))
#             result_dict['RT LMP [$/MWh]'] = float(round(RT_LMP[t],2))
#
#             # model vars
#             result_dict['Total Power Output [MW]'] = float(round(pyo.value(m.P_total[generator,t]),2))
#             result_dict['Thermal Power Generated [MW]'] = float(round(pyo.value(m.P_T[generator,t]),2))
#             result_dict['Thermal Power to Storage [MW]'] = float(round(pyo.value(m.P_E[generator,t]),2))
#             result_dict['Thermal Power to Market [MW]'] = float(round(pyo.value(m.P_G[generator,t]),2))
#             result_dict['Storage Power to Thermal [MW]'] = float(round(pyo.value(m.P_S[generator,t]),2))
#             result_dict['Total Thermal Side Power to Market [MW]'] = float(round(pyo.value(m.P_R[generator,t]),2))
#             result_dict['Charge Power [MW]'] = float(round(pyo.value(m.P_C[generator,t]),2))
#             result_dict['Disharge Power [MW]'] = float(round(pyo.value(m.P_D[generator,t]),2))
#             result_dict['State of Charge [MWh]'] = float(round(pyo.value(m.S_SOC[generator,t]),2))
#
#             result_dict['On/off [bin]'] = int(round(pyo.value(m.on_off[generator,t])))
#             result_dict['Start Up [bin]'] = int(round(pyo.value(m.start_up[generator,t])))
#             result_dict['Shut Down [bin]'] = int(round(pyo.value(m.shut_dw[generator,t])))
#             result_dict['Storage to Thermal [bin]'] = int(round(pyo.value(m.y_S[generator,t])))
#             result_dict['Thermal to Storage [bin]'] = int(round(pyo.value(m.y_E[generator,t])))
#             result_dict['Charge [bin]'] = int(round(pyo.value(m.y_C[generator,t])))
#             result_dict['Dicharge [bin]'] = int(round(pyo.value(m.y_D[generator,t])))
#
#             result_dict['Production Cost [$]'] = float(round(pyo.value(m.prod_cost_approx[generator,t]),2))
#             result_dict['Start-up Cost [$]'] = float(round(pyo.value(m.start_up_cost_expr[generator,t]),2))
#             result_dict['Total Cost [$]'] = float(round(pyo.value(m.tot_cost[generator,t]),2))
#             result_dict['DA Revenue [$]'] = float(round(DA_dispatches[t] * DA_LMP[t],2))
#             result_dict['RT Revenue [$]'] = float(round((RT_dispatches[t] - DA_dispatches[t]) * RT_LMP[t],2))
#             result_dict['Total Revenue [$]'] = float(round(result_dict['DA Revenue [$]'] + result_dict['RT Revenue [$]'],2))
#             result_dict['One-settlement Profit [$]'] = float(round(result_dict['RT Revenue [$]'] - result_dict['Total Cost [$]'],2))
#             result_dict['Two-settlement Profit [$]'] = float(round(result_dict['Total Revenue [$]'] - result_dict['Total Cost [$]'],2))
#
#             result_dict['Periodic Boundary Slack [MWh]'] = float(round(pyo.value(m.pbc_slack[generator,t]),2))
#             result_dict['Power Output Slack [MW]'] = float(round(pyo.value(m.slack_var_power[generator,t]),2))
#
#             # calculate mileage
#             if t == 0:
#                 result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[generator,t] - m.pre_P_T[generator])),2))
#             else:
#                 result_dict['Mileage [MW]'] = float(round(abs(pyo.value(m.P_T[generator,t] - m.P_T[generator,t-1])),2))
#
#             for key in kwargs:
#                 result_dict[key] = kwargs[key]
#
#             result_df = pd.DataFrame.from_dict(result_dict,orient = 'index')
#             df_list.append(result_df.T)
#
#     return pd.concat(df_list)
#
# def _update_UT_DT(m,last_implemented_time_step = 0):
#
#     # copy to a queue
#     pre_shut_down_trajectory_copy = {}
#     pre_start_up_trajectory_copy = {}
#
#     for unit in m.UNITS:
#         pre_shut_down_trajectory_copy[unit] = deque([])
#         pre_start_up_trajectory_copy[unit] = deque([])
#
#     for unit,t in m.pre_shut_down_trajectory_set:
#         pre_shut_down_trajectory_copy[unit].append(round(pyo.value(m.pre_shut_down_trajectory[unit,t])))
#     for unit,t in m.pre_start_up_trajectory_set:
#         pre_start_up_trajectory_copy[unit].append(round(pyo.value(m.pre_start_up_trajectory[unit,t])))
#
#     # add implemented trajectory to the queue
#     for unit in m.UNITS:
#         for t in range(last_implemented_time_step + 1):
#             pre_shut_down_trajectory_copy[unit].append(round(pyo.value(m.shut_dw[unit,t])))
#             pre_start_up_trajectory_copy[unit].append(round(pyo.value(m.start_up[unit,t])))
#
#     # pop out outdated trajectory
#     for unit in m.UNITS:
#
#         while len(pre_shut_down_trajectory_copy[unit]) > pyo.value(m.min_dw_time[unit]) - 1:
#             pre_shut_down_trajectory_copy[unit].popleft()
#         while len(pre_start_up_trajectory_copy[unit]) > pyo.value(m.min_up_time[unit]) - 1:
#             pre_start_up_trajectory_copy[unit].popleft()
#
#     # actual update
#     for unit,t in m.pre_shut_down_trajectory_set:
#         m.pre_shut_down_trajectory[unit,t] = pre_shut_down_trajectory_copy[unit].popleft()
#
#     for unit,t in m.pre_start_up_trajectory_set:
#         m.pre_start_up_trajectory[unit,t] = pre_start_up_trajectory_copy[unit].popleft()
#
#     return
#
# def _update_SOC(m,last_implemented_time_step = 0):
#
#     for unit in m.UNITS:
#         m.pre_SOC[unit] = round(pyo.value(m.S_SOC[unit,last_implemented_time_step]),2)
#
#     return
#
# def _relax_pbc_constriant(m):
#
#     n_timesteps = len(m.HOUR)
#
#     for j in m.UNITS:
#         for t in m.HOUR:
#             if t == n_timesteps - 1:
#                 m.pbc_slack[j,t].unfix()
#                 m.pbc_slack_positive[j,t].unfix()
#                 m.pbc_slack_negative[j,t].unfix()
#
#     return
#
# def _enforce_pbc_constriant(m):
#
#     for j in m.UNITS:
#         for t in m.HOUR:
#             m.pbc_slack[j,t].fix(0)
#             m.pbc_slack_positive[j,t].fix(0)
#             m.pbc_slack_negative[j,t].fix(0)
#
#     return
#
# def _unfix_slack_vars(m):
#
#     # unfix pbc slack vars
#     n_timesteps = len(m.HOUR)
#     for j in m.UNITS:
#         for t in m.HOUR:
#             if t == n_timesteps - 1:
#                 m.pbc_slack[j,t].unfix()
#                 m.pbc_slack_positive[j,t].unfix()
#                 m.pbc_slack_negative[j,t].unfix()
#
#     # slack variables for power dispatch
#     m.slack_var_power.unfix()
#
#     return
#
# def _fix_slack_vars(m):
#
#     m.pbc_slack.fix(0)
#     m.pbc_slack_positive.fix(0)
#     m.pbc_slack_negative.fix(0)
#
#     m.slack_var_power.fix(0)
#
#     return
#
# def _get_min_slack(m):
#
#     # get the min slack values for pbc
#     for j in m.UNITS:
#         total_slack = 0
#         for h in m.HOUR:
#             total_slack += pyo.value(m.pbc_slack_positive[j,h] + m.pbc_slack_negative[j,h])
#         m.pbc_slack_bound[j] = total_slack
#
#     for j in m.UNITS:
#         for h in m.HOUR:
#             m.power_slack_bound[j,h] = pyo.value(m.slack_var_power[j,h])
#
#     return
#
# def activate_relaxed_mode(m):
#
#     print('')
#     print('Activating relaxing mode...')
#
#     _unfix_slack_vars(m)
#
#     # deactivate original constraint and activate min slack constraint
#     m.abs_dev.deactivate()
#     m.min_slack.activate()
#
#     # deactivate original constraints
#     m.fix_power_output_con.deactivate()
#     m.SOC_endtime_con.deactivate()
#
#     # activate the relaxed constraints
#     m.power_slack_con1.activate()
#     m.power_slack_con2.activate()
#     m.pbc_slack_con.activate()
#     m.SOC_endtime_slack_con.activate()
#
#     return
#
# def activate_bounded_slack_mode(m):
#
#     _get_min_slack(m)
#
#     # solve the original problem with bounded slack
#     m.bound_pbc_slack_con.activate()
#     m.bound_power_slack_con.activate()
#
#     m.min_slack.deactivate()
#     m.abs_dev.activate()
#
#     return
#
# def deactivate_relaxed_mode(m):
#
#     print('')
#     print('Deactivating relaxing mode..')
#
#     _fix_slack_vars(m)
#
#     m.min_slack.deactivate()
#     m.abs_dev.activate()
#
#     # deactivate slack constraints
#     m.bound_pbc_slack_con.deactivate()
#     m.bound_power_slack_con.deactivate()
#     m.power_slack_con1.deactivate()
#     m.power_slack_con2.deactivate()
#     m.pbc_slack_con.deactivate()
#     m.SOC_endtime_slack_con.deactivate()
#
#     # activate original cosntraints
#     m.fix_power_output_con.activate()
#     m.SOC_endtime_con.activate()
#
#     return
#
# def solve_relaxed_model(m,solver):
#     '''
#     This function takes a hybrid model and relaxes its periodic boundary constriant(PBC)
#     and minimize the slack of that constraint. Then solves the original tracking
#     problem and constrains -1.01 * min slack <= slack <= 1.01 * min slack.
#     '''
#
#     activate_relaxed_mode(m)
#
#     print("")
#     print('Solving the minimizing slack problem...')
#     solver.solve(m,tee = True)
#
#     activate_bounded_slack_mode(m)
#
#     print("")
#     print('Solving the bounded slack problem...')
#     solver_result = solver.solve(m,tee = True)
#     solver_status = str(solver_result.solver.termination_condition)
#
#     deactivate_relaxed_mode(m)
#
#     return solver_status
#
# def _update_power(m,last_implemented_time_step = 0):
#
#     for unit in m.UNITS:
#         m.pre_P_T[unit] = round(pyo.value(m.P_T[unit,last_implemented_time_step]),2)
#         m.pre_on_off[unit] = round(pyo.value(m.on_off[unit,last_implemented_time_step]))
#
#     return
#
# def update_model(m,last_implemented_time_step = 0):
#
#     _update_UT_DT(m,last_implemented_time_step)
#     _update_power(m,last_implemented_time_step)
#     _update_SOC(m,last_implemented_time_step)
#
#     return
#
# def pass_schedule_to_track(m,market_signals,implement_steps = 1):
#
#     for unit in m.UNITS:
#
#         DA_dispatches = market_signals[unit]['DA_dispatches']
#         RT_dispatches = market_signals[unit]['RT_dispatches']
#
#         for t in m.HOUR:
#             if t < implement_steps:
#                 m.power_dispatch[unit,t] = RT_dispatches[t]
#                 # m.P_total[unit,t].fix(RT_dispatches[t])
#             else:
#                 m.power_dispatch[unit,t] = RT_dispatches[t]
#                 # m.P_total[unit,t].fix(RT_dispatches[t])
#
#     return
#
# if __name__ == "__main__":
#
#     model_data = assemble_model_data(generator_names = ['102_STEAM_3'],\
#                                      storage_size_ratio = [0.1],\
#                                      round_trip_eff = [0.88],\
#                                      rts_gmlc_data_dir = '../RTS-GMLC/RTS_Data/SourceData/gen.csv')
#     print(model_data)
#     m = build_tracking_model(model_data = model_data,plan_horizon = 10)
