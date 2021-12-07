import pandas as pd
import pyomo.environ as pyo
from idaes.apps.grid_integration.tracker import Tracker
from multiperiod_rankine_double_loop import MultiPeriodRankine
from utils import assemble_generator_data
import pytest

generator = "102_STEAM_3" #bidding generator
tracking_horizon = 4  #hours
bidding_horizon = 4   #hours
n_scenario = 10       #for bidding
n_tracking_hour = 1   #advance n_tracking_hour (i.e. assume we solve every hour)

# read generator parameters from rts-gmlc file
rts_gmlc_dataframe = pd.read_csv('gen.csv')

# create solver
solver = pyo.SolverFactory('ipopt')

#Generator data 
generator_data = assemble_generator_data(generator_name=generator, gen_params=rts_gmlc_dataframe)


#Setup trackers, bidder, and coordinator
#################################################################################

n_tracking_hour = 1
solver = pyo.SolverFactory('ipopt')

# create a tracker model
tracking_model_object = MultiPeriodRankine(horizon = tracking_horizon, generator_data=generator_data)
tracker_object = Tracker(tracking_model_object = tracking_model_object,\
                         n_tracking_hour = n_tracking_hour, \
                         solver = solver)

market_dispatch = [30, 40 , 50, 70]
tracker_object.track_market_dispatch(market_dispatch = market_dispatch, \
                                     date = "2021-07-26", \
                                     hour = '17:00')

for t, dispatch in zip(range(tracking_horizon), market_dispatch):
    assert pytest.approx(pyo.value(tracker_object.power_output[t]),abs = 1e-3) == dispatch

last_delivered_power = market_dispatch[0]
assert pytest.approx(tracker_object.get_last_delivered_power(), abs = 1e-3) == last_delivered_power
