# custom plugin file for running the rankine cycle with a battery as a
# multi-period model in the Prescient double-loop
import pickle
import pandas as pd
import pyomo.environ as pyo
from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.common.fileutils import this_file_dir

from idaes.apps.grid_integration.tracker import Tracker
from idaes.apps.grid_integration.bidder import Bidder
from idaes.apps.grid_integration.forecaster import PlaceHolderForecaster
from idaes.apps.grid_integration.coordinator import DoubleLoopCoordinator

from multiperiod_double_loop_rankine import MultiPeriodRankine
import os.path

with open("generator_data.pkl", "rb") as f:
    gen_data = pickle.load(f)

default_bid_curve = gen_data["Original Marginal Cost Curve"]
pmin = gen_data["PMin MW"]
pmax = gen_data["PMax MW"]
gen_name = gen_data["generator_name"]

tracking_horizon = 4  # hours
bidding_horizon = 48  # hours
n_scenario = 10  # for bidding
n_tracking_hour = 1  # advance n_tracking_hour (i.e. assume we solve every hour)

# create forecaster
price_forecasts_df = pd.read_csv(
    os.path.join(this_file_dir(), "../lmp_forecasts_concat.csv")
)
forecaster = PlaceHolderForecaster(price_forecasts_df=price_forecasts_df)

# create solver
solver = pyo.SolverFactory("ipopt")

# Setup trackers, bidder, and coordinator
#################################################################################
# Tracker
mp_rankine_tracker = MultiPeriodRankine(
    horizon=tracking_horizon,
    pmin=pmin,
    pmax=pmax,
    default_bid_curve=default_bid_curve,
    generator_name=gen_name,
)
thermal_tracker = Tracker(
    tracking_model_object=mp_rankine_tracker,
    n_tracking_hour=n_tracking_hour,
    solver=solver,
)

# Projection Tracker
mp_rankine_projection_tracker = MultiPeriodRankine(
    horizon=tracking_horizon,
    pmin=pmin,
    pmax=pmax,
    default_bid_curve=default_bid_curve,
    generator_name=gen_name,
)
thermal_projection_tracker = Tracker(
    tracking_model_object=mp_rankine_projection_tracker,
    n_tracking_hour=n_tracking_hour,
    solver=solver,
)

# Bidder
mp_rankine_bidder = MultiPeriodRankine(
    horizon=bidding_horizon,
    pmin=pmin,
    pmax=pmax,
    default_bid_curve=default_bid_curve,
    generator_name=gen_name,
)
thermal_bidder = Bidder(
    bidding_model_object=mp_rankine_bidder,
    n_scenario=n_scenario,
    solver=solver,
    forecaster=forecaster,
)

# Coordinator
coordinator = DoubleLoopCoordinator(
    bidder=thermal_bidder,
    tracker=thermal_tracker,
    projection_tracker=thermal_projection_tracker,
)

## Prescient requires the following functions in this module
get_configuration = coordinator.get_configuration
register_plugins = coordinator.register_plugins
