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

class TestMissingModel:

    method_dict = {k: lambda: None for k in ['populate_model', 'update_model']}
    attr_dict = {'power_output': 'power_output',\
                 'total_cost': ('tot_cost',1),\
                 'generator': 'test',\
                 'pmin': 20.00,\
                 'default_bids': {p: 30 for p in [20.00, 40.00, 60.00, 80.00, 100.00]} }

    def __init__(self, missing_method = None, missing_attr = None):

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

    solver = pyo.SolverFactory('cbc')
    forecaster = TestingForecaster(horizon = horizon, n_sample = n_scenario)

    method_list = ['populate_model', 'update_model']

    for m in method_list:
        bidding_model_object = TestMissingModel(missing_method = m)
        with pytest.raises(AttributeError, match = r".*{}().*".format(m)):
            bidder_object = Bidder(bidding_model_object = bidding_model_object,\
                                   n_scenario = n_scenario,\
                                   solver = solver,\
                                   forecaster = forecaster)

@pytest.mark.unit
def test_model_object_missing_attr():

    solver = pyo.SolverFactory('cbc')
    forecaster = TestingForecaster(horizon = horizon, n_sample = n_scenario)

    attr_list = ['power_output','total_cost','generator','pmin', 'default_bids']

    for attr in attr_list:
        bidding_model_object = TestMissingModel(missing_attr = attr)
        with pytest.raises(AttributeError, match = r".*{}().*".format(attr)):
            bidder_object = Bidder(bidding_model_object = bidding_model_object,\
                                   n_scenario = n_scenario,\
                                   solver = solver,\
                                   forecaster = forecaster)

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
