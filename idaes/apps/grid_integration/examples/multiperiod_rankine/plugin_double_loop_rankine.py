#custom plugin file for running the rankine cycle with a battery as a
#multi-period model in the Prescient double-loop
import pandas as pd
import pyomo.environ as pyo
from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.common.fileutils import this_file_dir

from idaes.apps.grid_integration.tracker import Tracker
from idaes.apps.grid_integration.bidder import Bidder
from idaes.apps.grid_integration.forecaster import PlaceHolderForecaster
from idaes.apps.grid_integration.coordinator import DoubleLoopCoordinator

from multiperiod_rankine_double_loop import MultiPeriodRankine

def assemble_generator_data(generator_name, gen_params):

    '''
    Builds a dictionary of generator parameteres given an rts-gmlc generator name and a dataframe of all generator data

    Arguments:
        generator_name: a generator name in RTS-GMLC dataset.
        gen_params: the RTS-GMLC generator data in Pandas DataFrame

    Returns:
        model_data: a dictionary which has this structure
        {data type name: value}.
    '''

    gen_params = gen_params.set_index('GEN UID',inplace = False)
    properties = ['PMin MW', 'PMax MW', 'Min Up Time Hr', 'Min Down Time Hr',\
                  'Ramp Rate MW/Min', 'Start Heat Warm MBTU', 'Fuel Price $/MMBTU',\
                  'HR_avg_0', 'HR_incr_1', 'HR_incr_2', 'HR_incr_3',\
                  'Output_pct_1','Output_pct_2','Output_pct_3']

    # to dict
    model_data = gen_params.loc[generator_name, properties].to_dict()
    model_data['generator_name'] = generator_name
    model_data['RU'] = model_data['Ramp Rate MW/Min'] * 60
    model_data['RD'] = model_data['RU']
    model_data['SU'] = min(model_data['PMin MW'], model_data['RU'])
    model_data['SD'] = min(model_data['PMin MW'], model_data['RD'])
    model_data['SU Cost'] = model_data['Start Heat Warm MBTU'] * model_data['Fuel Price $/MMBTU']

    model_data['Min Load Cost'] = model_data['HR_avg_0']/1000 * \
                                     model_data['Fuel Price $/MMBTU'] *\
                                     model_data['PMin MW']

    model_data['Power Segments'] = {}
    model_data['Marginal Costs'] = {}

    model_data['Original Marginal Cost Curve'] = {}
    model_data['Original Marginal Cost Curve'][model_data['PMin MW']] = model_data['Min Load Cost']/model_data['PMin MW']

    for l in range(1,4):
        model_data['Power Segments'][l] = model_data['Output_pct_{}'.format(l)] * model_data['PMax MW']
        model_data['Marginal Costs'][l] = model_data['HR_incr_{}'.format(l)]/1000 * model_data['Fuel Price $/MMBTU']
        model_data['Original Marginal Cost Curve'][model_data['Power Segments'][l]] = model_data['Marginal Costs'][l]

    return model_data

generator = "102_STEAM_3" #bidding generator
tracking_horizon = 4  #hours
bidding_horizon = 48   #hours
n_scenario = 10       #for bidding
n_tracking_hour = 1   #advance n_tracking_hour (i.e. assume we solve every hour)

# read generator parameters from rts-gmlc file
rts_gmlc_dataframe = pd.read_csv(this_file_dir()+'../gen.csv')

# create forecaster. QUESTION: What does a forecaster do?
price_forecasts_df = pd.read_csv(this_file_dir()+'../lmp_forecasts_concat.csv')
forecaster = PlaceHolderForecaster(price_forecasts_df = price_forecasts_df)

# create solver
solver = pyo.SolverFactory('ipopt')

#Generator data 
generator_data = assemble_generator_data(generator_name=generator, gen_params=rts_gmlc_dataframe)


#Setup trackers, bidder, and coordinator
#################################################################################
#Tracker
mp_rankine_tracker = MultiPeriodRankine(horizon = tracking_horizon, generator_data=generator_data)
thermal_tracker = Tracker(tracking_model_object = mp_rankine_tracker,\
                          n_tracking_hour = n_tracking_hour, \
                          solver = solver)

#Projection Tracker
mp_rankine_projection_tracker = MultiPeriodRankine(horizon = tracking_horizon, generator_data=generator_data)
thermal_projection_tracker = Tracker(tracking_model_object = mp_rankine_projection_tracker,\
                                     n_tracking_hour = n_tracking_hour, \
                                     solver = solver)

#Bidder
mp_rankine_bidder= MultiPeriodRankine(horizon = bidding_horizon, generator_data=generator_data)
thermal_bidder = Bidder(bidding_model_object = mp_rankine_bidder,\
                        n_scenario = n_scenario,\
                        solver = solver,\
                        forecaster = forecaster)

#Coordinator
coordinator = DoubleLoopCoordinator(bidder = thermal_bidder,\
                                    tracker = thermal_tracker,\
                                    projection_tracker = thermal_projection_tracker)

## Prescient requires the following functions in this module
get_configuration = coordinator.get_configuration
register_plugins = coordinator.register_plugins
