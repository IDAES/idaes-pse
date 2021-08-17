import pandas as pd
import pyomo.environ as pyo
import sys
sys.path.append('../')
from tracker import Tracker
from bidder import Bidder
from forecaster import WhiteNoiseForecaster
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
forecaster = WhiteNoiseForecaster(price_forecasts_df = price_forecasts_df)

# create solver
solver = pyo.SolverFactory('cbc')

# create trackers
thermal_tracker = Tracker(tracking_model_class = ThermalGenerator,\
                          n_tracking_hour = n_tracking_hour, \
                          solver = solver,\
                          rts_gmlc_dataframe = rts_gmlc_dataframe,\
                          horizon = tracking_horizon,\
                          generator = generator)

thermal_projection_tracker = Tracker(tracking_model_class = ThermalGenerator,\
                                     n_tracking_hour = n_tracking_hour, \
                                     solver = solver,\
                                     rts_gmlc_dataframe = rts_gmlc_dataframe,\
                                     horizon = tracking_horizon,\
                                     generator = generator)

# create bidder
thermal_bidder = Bidder(bidding_model_class = ThermalGenerator,\
                        n_scenario = n_scenario,\
                        solver = solver,\
                        forecaster = forecaster,\
                        rts_gmlc_dataframe = rts_gmlc_dataframe,\
                        horizon = bidding_horizon,\
                        generator = generator)

# create coordinator
coordinator = DoubleLoopCoordinator(bidder = thermal_bidder,\
                                    tracker = thermal_tracker,\
                                    projection_tracker = thermal_projection_tracker)
