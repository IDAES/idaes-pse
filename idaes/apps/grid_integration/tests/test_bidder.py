import pytest
import pyomo.environ as pyo
from idaes.apps.grid_integration.bidder import Bidder
from test_tracker import TestingModel

class TestingForecaster:

    def __init__(self, horizon, n_sample):
        self.horizon = horizon
        self.n_sample = n_sample

    def forecast(self, date, hour, prediction):
        return {i: [prediction]*self.horizon for i in range(self.n_sample)}

horizon = 4
n_scenario = 3

@pytest.fixture
def bidder_object():

    solver = pyo.SolverFactory('cbc')
    forecaster = TestingForecaster(horizon = horizon, n_sample = n_scenario)

    # create a tracker model
    bidding_model_object = TestingModel(horizon = horizon)
    bidder_object = Bidder(bidding_model_object = bidding_model_object,\
                           n_scenario = n_scenario,\
                           solver = solver,\
                           forecaster = forecaster)
    return bidder_object

@pytest.mark.component
def test_compute_bids(bidder_object):

    marginal_cost = bidder_object.bidding_model_object.marginal_cost
    gen = bidder_object.generator
    default_bids = bidder_object.bidding_model_object.default_bids
    pmin = bidder_object.bidding_model_object.pmin
    pmax = bidder_object.bidding_model_object.pmax
    date = '2021-08-20'

    # test forecaster lower than marginal cost
    shift = 1
    fixed_forecast = marginal_cost - shift
    bids = bidder_object.compute_bids(date = date, \
                                      hour = None, \
                                      prediction = fixed_forecast)

    expected_bids = {t: {gen: {p: p*marginal_cost - shift*pmin for p in default_bids}} for t in range(horizon)}

    assert expected_bids == bids

    # test forecaster lower than marginal cost
    shift = 1
    fixed_forecast = marginal_cost + shift
    bids = bidder_object.compute_bids(date = date, \
                                      hour = None, \
                                      prediction = fixed_forecast)

    expected_bids = {}
    for t in range(horizon):

        expected_bids[t]= {}
        expected_bids[t][gen] = {}

        pre_p = 0
        pre_cost = 0

        for p in default_bids:
            if p == pmax:
                expected_bids[t][gen][p] = (p - pre_p) * (marginal_cost + shift) + pre_cost
            else:
                expected_bids[t][gen][p] = (p - pre_p) * marginal_cost + pre_cost
            pre_p = p
            pre_cost = expected_bids[t][gen][p]

    assert expected_bids == bids
