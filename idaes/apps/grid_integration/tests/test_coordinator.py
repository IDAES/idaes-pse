import pytest
import pyomo.environ as pyo
from idaes.apps.grid_integration.bidder import Bidder
from idaes.apps.grid_integration.tracker import Tracker
from idaes.apps.grid_integration.coordinator import DoubleLoopCoordinator
from test_tracker import TestingModel
from test_bidder import TestingForecaster

tracking_horizon = 4
bidding_horizon = 48
n_scenario = 3
n_tracking_hour = 1

@pytest.fixture
def coordinator_object():

    # create solver
    solver = pyo.SolverFactory('cbc')

    ## create trackers
    # make a tracker
    tracking_model_object = TestingModel(horizon = tracking_horizon)
    thermal_tracker = Tracker(tracking_model_object = tracking_model_object,\
                              n_tracking_hour = n_tracking_hour, \
                              solver = solver)

    # make a projection tracker
    projection_tracking_model_object = TestingModel(horizon = tracking_horizon)
    thermal_projection_tracker = Tracker(tracking_model_object = projection_tracking_model_object,\
                                         n_tracking_hour = n_tracking_hour, \
                                         solver = solver)

    ## create a bidder
    forecaster = TestingForecaster(horizon = bidding_horizon, n_sample = n_scenario)
    bidding_model_object = TestingModel(horizon = bidding_horizon)
    thermal_bidder = Bidder(bidding_model_object = bidding_model_object,\
                            n_scenario = n_scenario,\
                            solver = solver,\
                            forecaster = forecaster)

    ## create coordinator
    coordinator = DoubleLoopCoordinator(bidder = thermal_bidder,\
                                        tracker = thermal_tracker,\
                                        projection_tracker = thermal_projection_tracker)

    return coordinator_object

@pytest.mark.unit
def test_assemble_project_tracking_signal(coordinator_object):

    gen_name = coordinator_object.bidder.generator
    current_ruc_dispatch = {(gen_name, t): (t+1) * 10 for t in range(24)}
    sced_horizon = tracking_horizon

    hour = 10
    expected_signal = [(hour + 1 + i) * 10 for i in range(sced_horizon)]
    signal = coordinator_object._assemble_project_tracking_signal(gen_name = gen_name,\
                                                                  current_ruc_dispatch = current_ruc_dispatch,\
                                                                  hour = hour,\
                                                                  sced_horizon = sced_horizon)
    assert signal == expected_signal

    hour = 23
    expected_signal = [(hour + 1) * 10 for i in range(sced_horizon)]
    signal = coordinator_object._assemble_project_tracking_signal(gen_name = gen_name,\
                                                                  current_ruc_dispatch = current_ruc_dispatch,\
                                                                  hour = hour,\
                                                                  sced_horizon = sced_horizon)
    assert signal == expected_signal

@pytest.mark.unit
def test_assemble_sced_tracking_market_signals(coordinator_object):

    gen_name = coordinator_object.bidder.generator
    sced_dispatch = [20] * tracking_horizon
    sced_horizon = tracking_horizon
    current_ruc_dispatch = {(gen_name, t): (t+1) * 10 for t in range(24)}

    # case 1
    hour = 10
    next_ruc_dispatch = None
    expected_signal = [sced_dispatch[0]] + \
                      [current_ruc_dispatch[(gen_name, t) for t in range(hour + 1, hour + sced_dispatch)]]
    signal = coordinator_object._assemble_sced_tracking_market_signals(gen_name = gen_name, \
                                                                       hour = hour, \
                                                                       sced_dispatch = sced_dispatch, \
                                                                       sced_horizon = sced_horizon, \
                                                                       current_ruc_dispatch = current_ruc_dispatch, \
                                                                       next_ruc_dispatch = next_ruc_dispatch)
    assert signal == expected_signal

    # case 2
    hour = 23
    next_ruc_dispatch = None
    expected_signal = sced_dispatch
    signal = coordinator_object._assemble_sced_tracking_market_signals(gen_name = gen_name, \
                                                                       hour = hour, \
                                                                       sced_dispatch = sced_dispatch, \
                                                                       sced_horizon = sced_horizon, \
                                                                       current_ruc_dispatch = current_ruc_dispatch, \
                                                                       next_ruc_dispatch = next_ruc_dispatch)
    assert signal == expected_signal

    # case 3
    hour = 23
    next_ruc_dispatch = {(gen_name, t): (t+1) * 10 for t in range(24)}
    expected_signal = [sced_dispatch[0]] + \
                      [next_ruc_dispatch[(gen_name, t) for t in range(0, hour + sced_dispatch - 1)]]
    signal = coordinator_object._assemble_sced_tracking_market_signals(gen_name = gen_name, \
                                                                       hour = hour, \
                                                                       sced_dispatch = sced_dispatch, \
                                                                       sced_horizon = sced_horizon, \
                                                                       current_ruc_dispatch = current_ruc_dispatch, \
                                                                       next_ruc_dispatch = next_ruc_dispatch)
    assert signal == expected_signal
