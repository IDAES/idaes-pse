import pytest
import pyomo.environ as pyo
from idaes.apps.grid_integration.bidder import Bidder
from test_tracker import TestingModel


class TestingForecaster:
    """
    A fake forecaster class for testing.
    """
    def __init__(self, horizon, n_sample):
        self.horizon = horizon
        self.n_sample = n_sample

    def forecast(self, date, hour, prediction):
        return {i: [prediction] * self.horizon for i in range(self.n_sample)}


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
    attr_list = ["power_output", "total_cost", "generator", "pmin", "default_bids"]

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
    bidding_model_object = TestingModel(horizon=horizon)

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
    bidding_model_object = TestingModel(horizon=horizon)

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
    bidding_model_object = TestingModel(horizon=horizon)
    bidder_object = Bidder(
        bidding_model_object=bidding_model_object,
        n_scenario=n_scenario,
        solver=solver,
        forecaster=forecaster,
    )
    return bidder_object


@pytest.mark.component
def test_compute_bids(bidder_object):

    marginal_cost = bidder_object.bidding_model_object.marginal_cost
    gen = bidder_object.generator
    default_bids = bidder_object.bidding_model_object.default_bids
    pmin = bidder_object.bidding_model_object.pmin
    pmax = bidder_object.bidding_model_object.pmax
    date = "2021-08-20"

    # test bidding when price forecast lower than marginal cost
    # expect default bids
    shift = 1
    fixed_forecast = marginal_cost - shift
    bids = bidder_object.compute_bids(date=date, hour=None, prediction=fixed_forecast)

    expected_bids = {
        t: {gen: {p: p * marginal_cost - shift * pmin for p in default_bids}}
        for t in range(horizon)
    }

    assert expected_bids == bids

    # test bidding when price forecast higher than marginal cost
    shift = 1
    fixed_forecast = marginal_cost + shift
    bids = bidder_object.compute_bids(date=date, hour=None, prediction=fixed_forecast)

    expected_bids = {}
    for t in range(horizon):

        expected_bids[t] = {}
        expected_bids[t][gen] = {}

        pre_p = 0
        pre_cost = 0

        for p in default_bids:
            # because the price forecast is higher than the marginal costs
            # to have highest profit, power output will be pmax
            # and the bidding costs will be computed with the price forecasts
            if p == pmax:
                expected_bids[t][gen][p] = (p - pre_p) * (
                    marginal_cost + shift
                ) + pre_cost
            else:
                expected_bids[t][gen][p] = (p - pre_p) * marginal_cost + pre_cost
            pre_p = p
            pre_cost = expected_bids[t][gen][p]

    assert expected_bids == bids


@pytest.mark.component
def test_is_convex_bid(bidder_object):

    # case 1: convex bid
    convex_bid = {20.0: 580.0, 40.0: 1180.0, 60.0: 1790.0, 80.0: 2410.0}
    assert bidder_object._is_convex_bid(convex_bid) == True

    # case 2: nonconvex bid
    nonconvex_bid = {20.0: 580.0, 40.0: 1180.0, 60.0: 1770.0, 80.0: 2380.0}
    assert bidder_object._is_convex_bid(nonconvex_bid) == False
