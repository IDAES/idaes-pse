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

solver = pyo.SolverFactory('ipopt')
forecaster = TestingForecaster(horizon = horizon, n_sample = n_scenario)

# create a tracker model
bidding_model_object = MultiPeriodRankine(horizon = horizon)
bidder_object = Bidder(bidding_model_object = bidding_model_object,\
                       n_scenario = n_scenario,\
                       solver = solver,\
                       forecaster = forecaster)