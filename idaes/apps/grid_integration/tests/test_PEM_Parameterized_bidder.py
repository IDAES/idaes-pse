#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
import pytest
import pandas as pd
from pyomo.common import unittest as pyo_unittest
from pyomo.opt.base.solvers import OptSolver
from idaes.apps.grid_integration.bidder import PEMParametrizedBidder
from idaes.apps.grid_integration.forecaster import PerfectForecaster
from idaes.apps.grid_integration.tests.util import (
    ExampleModel,
    ExampleForecaster,
    testing_renewable_data,
)
from idaes.apps.grid_integration.coordinator import prescient_avail
from idaes.apps.grid_integration.utils import convert_marginal_costs_to_actual_costs

day_ahead_horizon = 6
real_time_horizon = 6
# instead of using cbc, use a quasi solver to pass the _check_solver.
solver = OptSolver(type="solver")


@pytest.mark.unit
def test_creat_PEMParametrizedBidder_with_wrong_PEM_power():
    """
    This is to test when we creat the PEM bidder with a PEM power greater than the renewable power,
    there will be an error.
    """
    bidding_model_object = ExampleModel(model_data=testing_renewable_data)
    forecaster = ExampleForecaster(prediction=30)
    renewable_mw = 200
    pem_mw = 300
    pem_marginal_cost = 30
    with pytest.raises(
        ValueError, match=r".*The power of PEM is greater than the renewable power.*"
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


def wind_df():
    """
    This is to define a dataframe fed to the PerfectForecaster with example DA/RT
    capacity factors and LMPs.
    """
    start_year = 2020
    start_mon = 1
    start_day = 1
    start_date = pd.Timestamp(f"{start_year}-{start_mon:02d}-{start_day:02d} 00:00:00")
    ix = pd.date_range(
        start=start_date,
        end=start_date
        + pd.offsets.DateOffset(days=1)
        - pd.offsets.DateOffset(hours=24 - day_ahead_horizon + 1),
        freq="1H",
    )
    df = pd.DataFrame(index=ix)
    gen = testing_renewable_data.gen_name
    bus = testing_renewable_data.bus
    df[f"{gen}-DACF"] = [0.1, 0.1, 0.1, 0.8, 0.8, 0.8]
    df[f"{gen}-RTCF"] = [1.0, 0.2, 0.2, 0.75, 0.75, 0.75]
    df[f"{bus}-DALMP"] = [10] * day_ahead_horizon
    df[f"{bus}-RTLMP"] = [20] * real_time_horizon
    return df


@pytest.fixture
def bidder_object():
    """
    This is to define a bidder object.
    """
    example_wind_df = wind_df()
    forecaster = PerfectForecaster(example_wind_df)
    bidding_model_object = ExampleModel(model_data=testing_renewable_data)
    bidder_object = PEMParametrizedBidder(
        bidding_model_object=bidding_model_object,
        day_ahead_horizon=day_ahead_horizon,
        real_time_horizon=real_time_horizon,
        solver=solver,
        forecaster=forecaster,
        renewable_mw=200,
        pem_marginal_cost=30,
        pem_mw=100,
    )
    return bidder_object


@pytest.mark.component
@pytest.mark.skipif(
    not prescient_avail, reason="Prescient (optional dependency) not available"
)
def test_compute_DA_bids(bidder_object):
    """
    This is to test the if the compute_DA_bids function works correctly.
    """
    gen = bidder_object.generator
    pmin = bidder_object.bidding_model_object.pmin
    pmax = bidder_object.bidding_model_object.pmax
    pem_pmax = bidder_object.pem_mw
    pem_marginal_cost = bidder_object.pem_marginal_cost
    date = "2020-01-01"

    # test DA bidding
    bids = bidder_object.compute_day_ahead_bids(date=date, hour=0)
    expected_DA_cf = [0.1, 0.1, 0.1, 0.8, 0.8, 0.8]
    expected_bids = {}
    for t in range(day_ahead_horizon):
        expect_da_wind = expected_DA_cf[t] * pmax
        if t <= 2:
            expect_bids_curve = [(0, 0), (20, 30)]
        else:
            expect_bids_curve = [
                (0, 0),
                (60, 0),
                (160, 30),
            ]
        expect_cost_curve = convert_marginal_costs_to_actual_costs(expect_bids_curve)
        expected_bids[t] = {
            gen: {
                "p_cost": expect_cost_curve,
                "p_min": pmin,
                "p_max": expect_da_wind,
                "startup_capacity": expect_da_wind,
                "shutdown_capacity": expect_da_wind,
            }
        }
    pyo_unittest.assertStructuredAlmostEqual(first=expected_bids, second=bids)


@pytest.mark.component
@pytest.mark.skipif(
    not prescient_avail, reason="Prescient (optional dependency) not available"
)
def test_compute_RT_bids(bidder_object):
    """
    This is to test the if the compute_DA_bids function works correctly.
    """
    gen = bidder_object.generator
    pmin = bidder_object.bidding_model_object.pmin
    pmax = bidder_object.bidding_model_object.pmax
    pem_pmax = bidder_object.pem_mw
    pem_marginal_cost = bidder_object.pem_marginal_cost
    date = "2020-01-01"
    # totally 3 cases:
    # rt_wind - realized_day_ahead_dispatches > P_pem_max (t = 0)
    # 0 <= rt_wind - realized_day_ahead_dispatches <= P_pem_max (t = 1, 2, 3, 4)
    # rt_wind <= realized_day_ahead_dispatches (t = 5)
    realized_day_ahead_dispatches = [20, 20, 20, 120, 120, 180]
    realized_day_ahead_prices = None

    # testing RT bidding
    bids = bidder_object.compute_real_time_bids(
        date=date,
        hour=0,
        realized_day_ahead_dispatches=realized_day_ahead_dispatches,
        realized_day_ahead_prices=realized_day_ahead_prices,
    )

    expected_RT_cf = [1.0, 0.2, 0.2, 0.75, 0.75, 0.75]
    expected_bids = {}
    for t in range(real_time_horizon):
        expect_rt_wind = expected_RT_cf[t] * pmax
        if t == 0:
            expect_bids_curve = [
                (0, 0),
                (80, 0),
                (180, 30),
            ]
        elif t in [1, 2]:
            expect_bids_curve = [(0, 0), (20, 30)]
        elif t in [3, 4]:
            expect_bids_curve = [(0, 0), (30, 30)]
        else:
            expect_bids_curve = [(0, 0)]
        expect_cost_curve = convert_marginal_costs_to_actual_costs(expect_bids_curve)
        expected_bids[t] = {
            gen: {
                "p_cost": expect_cost_curve,
                "p_min": pmin,
                "p_max": max([p[0] for p in expect_cost_curve]),
                "startup_capacity": expect_rt_wind,
                "shutdown_capacity": expect_rt_wind,
            }
        }
    pyo_unittest.assertStructuredAlmostEqual(first=expected_bids, second=bids)
