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
from idaes.apps.grid_integration.tests.util import (
    TestingModel,
    TestingForecaster,
    testing_model_data,
)
from pyomo.common import unittest as pyo_unittest

from idaes.apps.grid_integration.coordinator import prescient_avail


class TestMissingModel:

    """
    A class for testing missing methods and attributes.
    """

    method_dict = {k: lambda: None for k in ["populate_model", "update_model"]}
    attr_dict = {
        "power_output": "power_output",
        "total_cost": ("tot_cost", 1),
        "generator": "test",
        "pmin": 20.00,
        "default_bids": {p: 30 for p in [20.00, 40.00, 60.00, 80.00, 100.00]},
    }

    def __init__(self, missing_method=None, missing_attr=None):

        """
        Constructs a model class without the specified missing methods and/or
        missing attributes.
        """

        for m in self.method_dict:
            if m != missing_method:
                setattr(self, m, self.method_dict[m])

        for a in self.attr_dict:
            if a != missing_attr:
                setattr(self, a, self.attr_dict[a])


horizon = 4
n_scenario = 3


@pytest.mark.unit
def test_model_object_missing_methods():

    solver = pyo.SolverFactory("cbc")
    forecaster = TestingForecaster(horizon=horizon, n_sample=n_scenario)

    # By definition, the model object should contain these methods
    method_list = ["populate_model", "update_model"]

    # test if the correct error message is raised if a model misses necessary methods
    for m in method_list:
        bidding_model_object = TestMissingModel(missing_method=m)
        with pytest.raises(AttributeError, match=r".*{}().*".format(m)):
            bidder_object = Bidder(
                bidding_model_object=bidding_model_object,
                n_scenario=n_scenario,
                solver=solver,
                forecaster=forecaster,
            )


@pytest.mark.unit
def test_model_object_missing_attr():

    solver = pyo.SolverFactory("cbc")
    forecaster = TestingForecaster(horizon=horizon, n_sample=n_scenario)

    # By definition, the model object should contain these attributes
    attr_list = ["power_output", "total_cost", "model_data"]

    # test if the correct error message is raised if a model misses necessary attributes
    for attr in attr_list:
        bidding_model_object = TestMissingModel(missing_attr=attr)
        with pytest.raises(AttributeError, match=r".*{}().*".format(attr)):
            bidder_object = Bidder(
                bidding_model_object=bidding_model_object,
                n_scenario=n_scenario,
                solver=solver,
                forecaster=forecaster,
            )


@pytest.mark.unit
def test_n_scenario_checker():

    solver = pyo.SolverFactory("cbc")
    forecaster = TestingForecaster(horizon=horizon, n_sample=n_scenario)
    bidding_model_object = TestingModel(model_data=testing_model_data, horizon=horizon)

    # test if bidder raise error when negative number of scenario is given
    with pytest.raises(ValueError, match=r".*greater than zero.*"):
        bidder_object = Bidder(
            bidding_model_object=bidding_model_object,
            n_scenario=-1,
            solver=solver,
            forecaster=forecaster,
        )

    # test if bidder raise error when float number of scenario is given
    with pytest.raises(TypeError, match=r".*should be an integer.*"):
        bidder_object = Bidder(
            bidding_model_object=bidding_model_object,
            n_scenario=3.0,
            solver=solver,
            forecaster=forecaster,
        )


@pytest.mark.unit
def test_solver_checker():

    forecaster = TestingForecaster(horizon=horizon, n_sample=n_scenario)
    bidding_model_object = TestingModel(model_data=testing_model_data, horizon=horizon)

    # test if bidder raise error when invalid solver is provided
    invalid_solvers = [5, "cbc", "ipopt"]
    for s in invalid_solvers:
        with pytest.raises(TypeError, match=r".*not a valid Pyomo solver.*"):
            bidder_object = Bidder(
                bidding_model_object=bidding_model_object,
                n_scenario=n_scenario,
                solver=s,
                forecaster=forecaster,
            )


@pytest.fixture
def bidder_object():

    solver = pyo.SolverFactory("cbc")
    forecaster = TestingForecaster(horizon=horizon, n_sample=n_scenario)

    # create a bidder model
    bidding_model_object = TestingModel(model_data=testing_model_data, horizon=horizon)
    bidder_object = Bidder(
        bidding_model_object=bidding_model_object,
        n_scenario=n_scenario,
        solver=solver,
        forecaster=forecaster,
    )
    return bidder_object


@pytest.mark.component
@pytest.mark.skipif(
    not prescient_avail, reason="Prescient (optional dependency) not available"
)
def test_compute_bids(bidder_object):

    marginal_cost = bidder_object.bidding_model_object.marginal_cost
    gen = bidder_object.generator
    default_bids = bidder_object.bidding_model_object.model_data.p_cost
    pmin = bidder_object.bidding_model_object.pmin
    pmax = bidder_object.bidding_model_object.pmax
    date = "2021-08-20"

    # test bidding when price forecast lower than marginal cost
    # expect default bids
    shift = 1
    fixed_forecast = marginal_cost - shift
    bids = bidder_object.compute_bids(date=date, hour=None, prediction=fixed_forecast)

    expected_bids = {
        t: {
            gen: {
                "p_min": pmin,
                "p_max": pmax,
                "startup_capacity": pmin,
                "shutdown_capacity": pmin,
                "p_cost": [
                    (p, p * marginal_cost - shift * pmin) for p, _ in default_bids
                ],
            }
        }
        for t in range(horizon)
    }

    pyo_unittest.assertStructuredAlmostEqual(first=expected_bids, second=bids)

    # test bidding when price forecast higher than marginal cost
    shift = 1
    fixed_forecast = marginal_cost + shift
    bids = bidder_object.compute_bids(date=date, hour=None, prediction=fixed_forecast)

    expected_bids = {}
    for t in range(horizon):

        expected_bids[t] = {}
        expected_bids[t][gen] = {
            "p_min": pmin,
            "p_max": pmax,
            "startup_capacity": pmin,
            "shutdown_capacity": pmin,
        }
        p_cost = []

        pre_p = 0
        pre_cost = 0

        for p, _ in default_bids:
            # because the price forecast is higher than the marginal costs
            # to have highest profit, power output will be pmax
            # and the bidding costs will be computed with the price forecasts
            if p == pmax:
                p_cost.append((p, (p - pre_p) * (marginal_cost + shift) + pre_cost))
            else:
                p_cost.append((p, (p - pre_p) * marginal_cost + pre_cost))
            pre_p = p
            pre_cost = p_cost[-1][1]

        expected_bids[t][gen]["p_cost"] = p_cost

    pyo_unittest.assertStructuredAlmostEqual(first=expected_bids, second=bids)
