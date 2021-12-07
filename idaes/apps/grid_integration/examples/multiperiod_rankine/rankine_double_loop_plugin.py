#custom plugin file for running the rankine cycle with a battery as a
#multi-period model in the Prescient double-loop
import pandas as pd
import pyomo.environ as pyo
from pyomo.common.config import ConfigDict, ConfigValue
import sys

from idaes.apps.grid_integration.tracker import Tracker
from idaes.apps.grid_integration.bidder import Bidder
from idaes.apps.grid_integration.forecaster import PlaceHolderForecaster
from idaes.apps.grid_integration.coordinator import DoubleLoopCoordinator

from multiperiod_rankine_double_loop import MultiPeriodRankine
from utils import assemble_generator_data

generator = "102_STEAM_3" #bidding generator
tracking_horizon = 4  #hours
bidding_horizon = 4   #hours
n_scenario = 10       #for bidding
n_tracking_hour = 1   #advance n_tracking_hour (i.e. assume we solve every hour)

# read generator parameters from rts-gmlc file
rts_gmlc_dataframe = pd.read_csv('gen.csv')

# create forecaster. QUESTION: What does a forecaster do?
price_forecasts_df = pd.read_csv('lmp_forecasts_concat.csv')
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
