import pandas as pd
import pyomo.environ as pyo
import sys
sys.path.append('../')
from tracker import Tracker
from bidder import Bidder
from forecaster import PlaceHolderForecaster
from prescient_plugin import DoubleLoopCoordinator
from thermal_generator import ThermalGenerator

generator = "102_STEAM_3"
tracking_horizon = 4
bidding_horizon = 48
n_scenario = 10
n_tracking_hour = 1

# read generator param data
rts_gmlc_dataframe = pd.read_csv('gen.csv')

# create forecaster
price_forecasts_df = pd.read_csv('lmp_forecasts_concat.csv')
forecaster = PlaceHolderForecaster(price_forecasts_df = price_forecasts_df)

# create solver
solver = pyo.SolverFactory('cbc')

## create trackers

# make a tracker
tracking_model_object = ThermalGenerator(rts_gmlc_dataframe = rts_gmlc_dataframe,\
                                         horizon = tracking_horizon, \
                                         generator = generator)
thermal_tracker = Tracker(tracking_model_object = tracking_model_object,\
                          n_tracking_hour = n_tracking_hour, \
                          solver = solver)

# make a projection tracker
projection_tracking_model_object = ThermalGenerator(rts_gmlc_dataframe = rts_gmlc_dataframe,\
                                                    horizon = tracking_horizon, \
                                                    generator = generator)
thermal_projection_tracker = Tracker(tracking_model_object = projection_tracking_model_object,\
                                     n_tracking_hour = n_tracking_hour, \
                                     solver = solver)

# create a bidder
bidding_model_object = ThermalGenerator(rts_gmlc_dataframe = rts_gmlc_dataframe,\
                                        horizon = bidding_horizon, \
                                        generator = generator)
thermal_bidder = Bidder(bidding_model_object = bidding_model_object,\
                        n_scenario = n_scenario,\
                        solver = solver,\
                        forecaster = forecaster)

# create coordinator
coordinator = DoubleLoopCoordinator(bidder = thermal_bidder,\
                                    tracker = thermal_tracker,\
                                    projection_tracker = thermal_projection_tracker)
