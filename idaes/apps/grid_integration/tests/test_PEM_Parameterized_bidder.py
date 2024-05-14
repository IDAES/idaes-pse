#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
import pytest
from pyomo.common import unittest as pyo_unittest
from pyomo.opt.base.solvers import OptSolver
from idaes.apps.grid_integration.bidder import PEMParametrizedBidder
from idaes.apps.grid_integration.forecaster import PerfectForecaster
from idaes.apps.grid_integration.tests.util import (
    TestingModel,
    TestingForecaster,
    testing_model_data,
)
from idaes.apps.grid_integration.coordinator import prescient_avail

day_ahead_horizon = 24
real_time_horizon = 4
# instead of using cbc, use a quasi solver to pass the _check_solver.
solver = OptSolver(type="solver")


@pytest.mark.unit
def test_creat_PEMParametrizedBidder_with_wrong_PEM_power():
    bidding_model_object = TestingModel(model_data=testing_model_data)
    forecaster = TestingForecaster(prediction=30)
    renewable_mw = 200
    pem_mw = 300
    pem_marginal_cost = 30
    with pytest.raises(
        ValueError, match=r".*The power of PEM is greater than the renewabele power.*"
    ):
        PEM_bidder = PEMParametrizedBidder(
            bidding_model_object,
            day_ahead_horizon,
            real_time_horizon,
            solver,
            forecaster,
            renewable_mw,
            pem_marginal_cost,
            pem_mw,
        )


@pytest.fixture
def bidder_object():
    forecaster = TestingForecaster(prediction=30)
    bidding_model_object = TestingModel(model_data=testing_model_data)
    bidder_object = PEMParametrizedBidder(
        bidding_model_object=bidding_model_object,
        day_ahead_horizon=day_ahead_horizon,
        real_time_horizon=real_time_horizon,
        solver=solver,
        forecaster=forecaster,
        renewable_mw=400,
        pem_marginal_cost=30,
        pem_mw=200,
    )
    return bidder_object
