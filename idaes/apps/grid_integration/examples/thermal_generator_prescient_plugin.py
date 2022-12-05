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
import pyomo.environ as pyo
import os
from idaes.apps.grid_integration import Tracker
from idaes.apps.grid_integration import Bidder
from idaes.apps.grid_integration import PlaceHolderForecaster
from idaes.apps.grid_integration import DoubleLoopCoordinator
from idaes.apps.grid_integration.examples.thermal_generator import ThermalGenerator
from idaes.apps.grid_integration.examples.utils import (
    rts_gmlc_generator_dataframe,
    rts_gmlc_bus_dataframe,
    daily_da_price_means,
    daily_rt_price_means,
    daily_da_price_stds,
    daily_rt_price_stds,
)

this_module_dir = os.path.dirname(__file__)

generator = "10_STEAM"
tracking_horizon = 4
day_ahead_bidding_horizon = 48
real_time_bidding_horizon = tracking_horizon
n_scenario = 10
n_tracking_hour = 1

# create forecaster
forecaster = PlaceHolderForecaster(
    daily_da_price_means=daily_da_price_means,
    daily_rt_price_means=daily_rt_price_means,
    daily_da_price_stds=daily_da_price_stds,
    daily_rt_price_stds=daily_rt_price_stds,
)

# create solver
solver = pyo.SolverFactory("cbc")

## create trackers

# make a tracker
tracking_model_object = ThermalGenerator(
    rts_gmlc_generator_dataframe=rts_gmlc_generator_dataframe,
    rts_gmlc_bus_dataframe=rts_gmlc_bus_dataframe,
    generator=generator,
)
thermal_tracker = Tracker(
    tracking_model_object=tracking_model_object,
    tracking_horizon=tracking_horizon,
    n_tracking_hour=n_tracking_hour,
    solver=solver,
)

# make a projection tracker
projection_tracking_model_object = ThermalGenerator(
    rts_gmlc_generator_dataframe=rts_gmlc_generator_dataframe,
    rts_gmlc_bus_dataframe=rts_gmlc_bus_dataframe,
    generator=generator,
)
thermal_projection_tracker = Tracker(
    tracking_model_object=projection_tracking_model_object,
    tracking_horizon=tracking_horizon,
    n_tracking_hour=n_tracking_hour,
    solver=solver,
)

# create a bidder
bidding_model_object = ThermalGenerator(
    rts_gmlc_generator_dataframe=rts_gmlc_generator_dataframe,
    rts_gmlc_bus_dataframe=rts_gmlc_bus_dataframe,
    generator=generator,
)
thermal_bidder = Bidder(
    bidding_model_object=bidding_model_object,
    day_ahead_horizon=day_ahead_bidding_horizon,
    real_time_horizon=real_time_bidding_horizon,
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
