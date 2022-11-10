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
from idaes.apps.grid_integration.tracker import Tracker
from idaes.apps.grid_integration.tests.util import TestingModel, testing_model_data


class TestMissingModel:

    """
    A class for testing missing methods and attributes.
    """

    method_list = [
        "populate_model",
        "get_implemented_profile",
        "update_model",
        "get_last_delivered_power",
        "record_results",
        "write_results",
    ]
    method_dict = {k: lambda: None for k in method_list}

    attr_dict = {"power_output": "power_output", "total_cost": ("tot_cost", 1)}

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


@pytest.mark.unit
def test_model_object_missing_methods():

    n_tracking_hour = 1
    solver = pyo.SolverFactory("cbc")

    # By definition, the model object should contain these methods
    method_list = [
        "populate_model",
        "get_implemented_profile",
        "update_model",
        "get_last_delivered_power",
        "record_results",
        "write_results",
    ]

    # test if the correct error message is raised if a model misses necessary methods
    for m in method_list:
        tracking_model_object = TestMissingModel(missing_method=m)
        with pytest.raises(AttributeError, match=r".*{}().*".format(m)):
            tracker_object = Tracker(
                tracking_model_object=tracking_model_object,
                tracking_horizon=horizon,
                n_tracking_hour=n_tracking_hour,
                solver=solver,
            )


@pytest.mark.unit
def test_model_object_missing_attr():

    n_tracking_hour = 1
    solver = pyo.SolverFactory("cbc")
    # By definition, the model object should contain these attributes
    attr_list = ["power_output", "total_cost"]

    # test if the correct error message is raised if a model misses necessary attributes
    for a in attr_list:
        tracking_model_object = TestMissingModel(missing_attr=a)
        with pytest.raises(AttributeError, match=r".*{}.*".format(a)):
            tracker_object = Tracker(
                tracking_model_object=tracking_model_object,
                tracking_horizon=horizon,
                n_tracking_hour=n_tracking_hour,
                solver=solver,
            )


@pytest.mark.unit
def test_n_tracking_hour_checker():

    solver = pyo.SolverFactory("cbc")
    tracking_model_object = TestingModel(model_data=testing_model_data)

    # test if tracker raise error when negative number of hours is given
    n_tracking_hour = -1
    with pytest.raises(ValueError, match=r".*greater than zero.*"):
        tracker_object = Tracker(
            tracking_model_object=tracking_model_object,
            tracking_horizon=horizon,
            n_tracking_hour=n_tracking_hour,
            solver=solver,
        )

    # test if tracker raise error when floating number of hours is given
    n_tracking_hour = 3.0
    with pytest.raises(TypeError, match=r".*should be an integer.*"):
        tracker_object = Tracker(
            tracking_model_object=tracking_model_object,
            tracking_horizon=horizon,
            n_tracking_hour=n_tracking_hour,
            solver=solver,
        )


@pytest.mark.unit
def test_solver_checker():

    n_tracking_hour = 1
    tracking_model_object = TestingModel(model_data=testing_model_data)

    # test if bidder raise error when invalid solver is provided
    invalid_solvers = [5, "cbc", "ipopt"]
    for s in invalid_solvers:
        with pytest.raises(TypeError, match=r".*not a valid Pyomo solver.*"):
            tracker_object = Tracker(
                tracking_model_object=tracking_model_object,
                tracking_horizon=horizon,
                n_tracking_hour=n_tracking_hour,
                solver=s,
            )


@pytest.fixture
def tracker_object():

    n_tracking_hour = 1
    solver = pyo.SolverFactory("cbc")

    # create a tracker model
    tracking_model_object = TestingModel(model_data=testing_model_data)
    tracker_object = Tracker(
        tracking_model_object=tracking_model_object,
        tracking_horizon=horizon,
        n_tracking_hour=n_tracking_hour,
        solver=solver,
    )
    return tracker_object


@pytest.mark.component
def test_track_market_dispatch(tracker_object):

    market_dispatch = [30, 40, 50, 70]
    tracker_object.track_market_dispatch(
        market_dispatch=market_dispatch, date="2021-07-26", hour="17:00"
    )

    for t, dispatch in zip(range(horizon), market_dispatch):
        assert (
            pytest.approx(pyo.value(tracker_object.power_output[t]), abs=1e-3)
            == dispatch
        )

    last_delivered_power = market_dispatch[0]
    assert (
        pytest.approx(tracker_object.get_last_delivered_power(), abs=1e-3)
        == last_delivered_power
    )


@pytest.mark.component
def test_track_deviation_penalty(tracker_object):

    large_penalty = 10000

    assert pytest.approx(large_penalty) == pyo.value(
        tracker_object.model.deviation_penalty[0]
    )

    for t in range(1, horizon):
        assert pytest.approx(
            large_penalty / (horizon - tracker_object.n_tracking_hour)
        ) == pyo.value(tracker_object.model.deviation_penalty[t])
