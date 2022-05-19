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
import pytest

import pyomo.environ as pyo
from idaes.apps.grid_integration.bidder import Bidder
from idaes.apps.grid_integration.tracker import Tracker
from idaes.apps.grid_integration.coordinator import DoubleLoopCoordinator
from idaes.apps.grid_integration.tests.util import (
    TestingModel,
    TestingForecaster,
    testing_model_data,
)
from pyomo.common import unittest as pyo_unittest

tracking_horizon = 4
bidding_horizon = 48
n_scenario = 3
n_tracking_hour = 1


@pytest.fixture
def coordinator_object():

    # create solver
    solver = pyo.SolverFactory("cbc")

    ## create trackers
    # make a tracker
    tracking_model_object = TestingModel(
        model_data=testing_model_data, horizon=tracking_horizon
    )
    thermal_tracker = Tracker(
        tracking_model_object=tracking_model_object,
        n_tracking_hour=n_tracking_hour,
        solver=solver,
    )

    # make a projection tracker
    projection_tracking_model_object = TestingModel(
        model_data=testing_model_data, horizon=tracking_horizon
    )
    thermal_projection_tracker = Tracker(
        tracking_model_object=projection_tracking_model_object,
        n_tracking_hour=n_tracking_hour,
        solver=solver,
    )

    ## create a bidder
    forecaster = TestingForecaster(horizon=bidding_horizon, n_sample=n_scenario)
    bidding_model_object = TestingModel(
        model_data=testing_model_data, horizon=bidding_horizon
    )
    thermal_bidder = Bidder(
        bidding_model_object=bidding_model_object,
        n_scenario=n_scenario,
        solver=solver,
        forecaster=forecaster,
    )

    ## create coordinator
    coordinator_object = DoubleLoopCoordinator(
        bidder=thermal_bidder,
        tracker=thermal_tracker,
        projection_tracker=thermal_projection_tracker,
    )

    return coordinator_object


@pytest.mark.unit
def test_assemble_sced_tracking_market_signals(coordinator_object):

    gen_name = coordinator_object.bidder.generator
    # assumes constant sced dispatch signal in the horizon
    constant_dispatch = 20
    sced_dispatch = [constant_dispatch] * tracking_horizon
    current_ruc_dispatch_dicts = [{(gen_name, t): (t + 1) * 10 for t in range(24)}]

    # test case 1: no ruc signals from next day
    hour = 10
    next_ruc_dispatch_dicts = None

    # expected sced signals are: the sced signal from the coming hour and
    # corresponding ruc signals (current day) for the remaining horizon
    expected_signal = [sced_dispatch[0]] + [
        current_ruc_dispatch_dicts[0][(gen_name, t)]
        for t in range(hour + 1, hour + tracking_horizon)
    ]
    signal = coordinator_object._assemble_sced_tracking_market_signals(
        gen_name=gen_name,
        hour=hour,
        sced_dispatch=sced_dispatch,
        tracking_horizon=tracking_horizon,
        current_ruc_dispatch_dicts=current_ruc_dispatch_dicts,
        next_ruc_dispatch_dicts=next_ruc_dispatch_dicts,
    )
    pyo_unittest.assertStructuredAlmostEqual(first=signal, second=expected_signal)

    # test case 2: no ruc signals, but between 2 days
    hour = 23
    next_ruc_dispatch_dicts = None

    # expected sced signals are: because there is no next ruc signals, the expected
    # signal will be the same as the sced dispatch
    expected_signal = [sced_dispatch[0]] * tracking_horizon
    signal = coordinator_object._assemble_sced_tracking_market_signals(
        gen_name=gen_name,
        hour=hour,
        sced_dispatch=sced_dispatch,
        tracking_horizon=tracking_horizon,
        current_ruc_dispatch_dicts=current_ruc_dispatch_dicts,
        next_ruc_dispatch_dicts=next_ruc_dispatch_dicts,
    )
    pyo_unittest.assertStructuredAlmostEqual(first=signal, second=expected_signal)

    # test case 3: with ruc signals, between 2 days
    hour = 23
    next_ruc_dispatch_dicts = [{(gen_name, t): (t + 1) * 10 for t in range(24)}]

    # expected sced signals are: because there is ruc signals, the expected
    # signal will be the sced signal from the coming hour and
    # corresponding ruc signals (next day) for the remaining horizon
    expected_signal = [sced_dispatch[0]] + [
        next_ruc_dispatch_dicts[0][(gen_name, t)]
        for t in range(0, tracking_horizon - 1)
    ]
    signal = coordinator_object._assemble_sced_tracking_market_signals(
        gen_name=gen_name,
        hour=hour,
        sced_dispatch=sced_dispatch,
        tracking_horizon=tracking_horizon,
        current_ruc_dispatch_dicts=current_ruc_dispatch_dicts,
        next_ruc_dispatch_dicts=next_ruc_dispatch_dicts,
    )
    pyo_unittest.assertStructuredAlmostEqual(first=signal, second=expected_signal)
