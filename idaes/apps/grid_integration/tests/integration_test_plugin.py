#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
import pandas as pd
import pyomo.environ as pyo
import os
from test_tracker import TestingModel
from idaes.apps.grid_integration import Tracker
from idaes.apps.grid_integration import Bidder
from idaes.apps.grid_integration import PlaceHolderForecaster
from idaes.apps.grid_integration import DoubleLoopCoordinator

this_module_dir = os.path.dirname(__file__)

generator = "10_STEAM"
pmin = 30
pmax = 76
tracking_horizon = 4
bidding_horizon = 48
n_scenario = 10
n_tracking_hour = 1

# create forecaster
price_forecasts_df = pd.read_csv(
    os.path.join(this_module_dir, os.pardir, "examples", "lmp_forecasts_concat.csv")
)
forecaster = PlaceHolderForecaster(price_forecasts_df=price_forecasts_df)

# create solver
solver = pyo.SolverFactory("cbc")

## create trackers

# make a tracker
tracking_model_object = TestingModel(horizon=tracking_horizon, name=generator, pmin=pmin, pmax=pmax)
thermal_tracker = Tracker(
    tracking_model_object=tracking_model_object,
    n_tracking_hour=n_tracking_hour,
    solver=solver,
)

# make a projection tracker
projection_tracking_model_object = TestingModel(horizon=tracking_horizon, name=generator, pmin=pmin, pmax=pmax)
thermal_projection_tracker = Tracker(
    tracking_model_object=projection_tracking_model_object,
    n_tracking_hour=n_tracking_hour,
    solver=solver,
)

# create a bidder
bidding_model_object = TestingModel(horizon=bidding_horizon, name=generator, pmin=pmin, pmax=pmax)
thermal_bidder = Bidder(
    bidding_model_object=bidding_model_object,
    n_scenario=n_scenario,
    solver=solver,
    forecaster=forecaster,
)

# create coordinator
coordinator = DoubleLoopCoordinator(
    bidder=thermal_bidder,
    tracker=thermal_tracker,
    projection_tracker=thermal_projection_tracker,
)

## Prescient requires the following functions in this module
get_configuration = coordinator.get_configuration
register_plugins = coordinator.register_plugins
