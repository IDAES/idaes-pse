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
from pyomo.common import unittest as pyo_unittest
from idaes.apps.grid_integration.model_data import GeneratorModelData


@pytest.fixture
def generator_params():
    return {
        "gen_name": "Testing_Generator",
        "generator_type": "thermal",
        "p_min": 30,
        "p_max": 76,
        "min_down_time": 2,
        "min_up_time": 3,
        "ramp_up_60min": 100,
        "ramp_down_60min": 100,
        "shutdown_capacity": 30,
        "startup_capacity": 30,
        "production_cost_bid_pairs": [(30, 25), (45, 23), (60, 27), (76, 35)],
        "startup_cost_pairs": [(2, 1000), (6, 1500), (10, 2000)],
        "fixed_commitment": None,
    }


@pytest.fixture
def generator_data_object(generator_params):
    return GeneratorModelData(**generator_params)


@pytest.mark.unit
def test_create_model_data_object(generator_params, generator_data_object):

    # test scalar values
    for param_name in generator_params:
        if param_name not in ["production_cost_bid_pairs", "startup_cost_pairs"]:
            param_val = getattr(generator_data_object, param_name)
            assert param_val == generator_params[param_name]

    # test bids
    attr_names = ["p_cost", "startup_cost"]
    param_names = ["production_cost_bid_pairs", "startup_cost_pairs"]

    for a_name, p_name in zip(attr_names, param_names):
        bids = getattr(generator_data_object, a_name)
        expected_bid = [(key, val) for key, val in generator_params[p_name]]
        pyo_unittest.assertStructuredAlmostEqual(first=bids, second=expected_bid)


@pytest.mark.unit
@pytest.mark.parametrize(
    "param_name, value",
    [
        ("p_min", "30"),
        ("p_max", "76"),
        ("min_down_time", "10"),
        ("min_up_time", "100"),
        ("ramp_up_60min", "100"),
        ("ramp_down_60min", "-100"),
        ("shutdown_capacity", "100"),
        ("startup_capacity", "100"),
    ],
)
def test_create_model_data_with_non_real_numbers(param_name, value, generator_params):
    generator_params[param_name] = value
    with pytest.raises(
        TypeError, match=f"Value for {param_name} shoulde be real numbers."
    ):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
@pytest.mark.parametrize(
    "param_name, value",
    [
        ("p_min", -30),
        ("min_down_time", -10),
        ("min_up_time", -10),
        ("ramp_up_60min", -0.01),
        ("ramp_down_60min", -100),
    ],
)
def test_create_model_data_with_negative_values(param_name, value, generator_params):
    generator_params[param_name] = value
    with pytest.raises(ValueError, match="Value should be greater than or equal to 0."):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
@pytest.mark.parametrize(
    "param_name, value",
    [("p_max", 25), ("shutdown_capacity", 10), ("startup_capacity", 15)],
)
def test_create_model_data_with_less_than_pmin_data(
    param_name, value, generator_params
):
    generator_params[param_name] = value
    with pytest.raises(
        ValueError, match=f"Value for {param_name} shoulde be greater or equal to Pmin."
    ):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_fixed_commitment_value(generator_params):
    generator_params["fixed_commitment"] = 5
    with pytest.raises(
        ValueError, match=r"^(Value for generator fixed commitment must be one of)"
    ):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_generator_name(generator_params):
    generator_params["gen_name"] = 102_111
    with pytest.raises(TypeError, match=r".*generator names must be str.*"):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_generator_type(generator_params):
    generator_params["generator_type"] = "nuclear"
    with pytest.raises(
        ValueError, match=r"^(Value for generator types must be one of)"
    ):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
@pytest.mark.parametrize(
    "param_name", ["production_cost_bid_pairs", "startup_cost_pairs"]
)
def test_empty_bid_pairs(param_name, generator_params):
    generator_params[param_name] = None
    with pytest.raises(
        ValueError, match=r"Empty production|startup cost pairs are provided."
    ):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
def test_bid_missing_pmin(generator_params):
    generator_params["production_cost_bid_pairs"] = generator_params[
        "production_cost_bid_pairs"
    ][1:]
    with pytest.raises(
        ValueError, match=r"^(The first power output in the bid should be the Pmin)"
    ):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
def test_bid_missing_pmax(generator_params):
    generator_params["production_cost_bid_pairs"].pop()
    with pytest.raises(
        ValueError, match=r"^(The last power output in the bid should be the Pmax)"
    ):
        GeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_start_up_bid(generator_params):
    generator_params["startup_cost_pairs"] = [(0, 0)]
    with pytest.raises(
        ValueError,
        match=r"^(The first startup lag should be the same as minimum down time)",
    ):
        GeneratorModelData(**generator_params)
