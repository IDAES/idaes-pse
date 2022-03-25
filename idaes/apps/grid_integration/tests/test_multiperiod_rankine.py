import pytest
import pyomo.environ as pyo
from idaes.apps.grid_integration.tracker import Tracker
from idaes.apps.grid_integration.bidder import Bidder
from idaes.apps.grid_integration.examples.multiperiod_rankine.multiperiod_double_loop_rankine import (
    MultiPeriodRankine,
)
from test_bidder import TestingForecaster
from pyomo.common.unittest import assertStructuredAlmostEqual

# test tracking with MultiPeriodRankine model
@pytest.mark.component
def test_track_market_dispatch():
    tracking_horizon = 4
    n_tracking_hour = 1
    solver = pyo.SolverFactory("ipopt")

    pmin = 20.0
    pmax = 100.0
    tracking_model_object = MultiPeriodRankine(tracking_horizon, pmin=pmin, pmax=pmax)

    tracker_object = Tracker(
        tracking_model_object=tracking_model_object,
        n_tracking_hour=n_tracking_hour,
        solver=solver,
    )
    market_dispatch = [30, 40, 50, 70]
    tracker_object.track_market_dispatch(
        market_dispatch=market_dispatch, date="2021-07-26", hour="17:00"
    )

    blks = tracking_model_object.mp_rankine.get_active_process_blocks()
    assert len(blks) == tracking_horizon

    for blk in blks:
        assert blk.rankine.fs.eq_min_power.lb / 1e6 == pmin
        assert blk.rankine.fs.eq_max_power.ub / 1e6 == pmax

    for t, dispatch in zip(range(tracking_horizon), market_dispatch):
        assert (
            pytest.approx(pyo.value(tracker_object.power_output[t]), abs=1e-3)
            == dispatch
        )

    last_delivered_power = market_dispatch[0]
    assert (
        pytest.approx(tracker_object.get_last_delivered_power(), abs=1e-3)
        == last_delivered_power
    )


# test bidding with MultiPeriodRankine model
@pytest.mark.component
def test_compute_bids():
    bidding_horizon = 4
    n_scenario = 3

    solver = pyo.SolverFactory("ipopt")
    forecaster = TestingForecaster(horizon=bidding_horizon, n_sample=n_scenario)

    pmin = 20.0
    pmax = 100.0
    bidding_model_object = MultiPeriodRankine(
        horizon=bidding_horizon, pmin=pmin, pmax=pmax
    )

    bidder_object = Bidder(
        bidding_model_object=bidding_model_object,
        n_scenario=n_scenario,
        solver=solver,
        forecaster=forecaster,
    )

    date = "2021-08-20"
    shift = 1
    fixed_forecast = 30.0 - shift
    bids = bidder_object.compute_bids(date=date, hour=None, prediction=fixed_forecast)

    blks = bidding_model_object.mp_rankine.get_active_process_blocks()
    assert len(blks) == bidding_horizon

    # test against known solution with ipopt
    known_solution = {
        0: {
            "gen": {
                20.0: 600.0,
                40.0: 1200.0,
                60.0: 1800.0,
                70.42: 2112.6,
                80.0: 2400.0,
                100.0: 3000.0,
            }
        },
        1: {
            "gen": {
                20.0: 600.0,
                40.0: 1200.0,
                60.0: 1800.0,
                80.0: 2400.0,
                91.34: 2740.2,
                100.0: 3000.0,
            }
        },
        2: {
            "gen": {
                20.0: 600.0,
                40.0: 1200.0,
                60.0: 1800.0,
                80.0: 2400.0,
                100.0: 3000.0,
                108.66: 3259.8,
            }
        },
        3: {
            "gen": {
                20.0: 600.0,
                40.0: 1200.0,
                60.0: 1800.0,
                80.0: 2400.0,
                100.0: 3000.0,
                129.58: 3887.4,
            }
        },
    }

    assertStructuredAlmostEqual(bids, known_solution)
