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
from idaes.apps.grid_integration.model_data import (
    ThermalGeneratorModelData,
    RenewableGeneratorModelData,
)


@pytest.fixture
def generator_params():
    return {
        "gen_name": "Testing_Generator",
        "bus": "bus5",
        "p_min": 30,
        "p_max": 76,
        "min_down_time": 2,
        "min_up_time": 3,
        "ramp_up_60min": 100,
        "ramp_down_60min": 100,
        "shutdown_capacity": 30,
        "startup_capacity": 30,
        "initial_status": 4,
        "initial_p_output": 30,
        "production_cost_bid_pairs": [(30, 25), (45, 23), (60, 27), (76, 35)],
        "startup_cost_pairs": [(2, 1000), (6, 1500), (10, 2000)],
        "fixed_commitment": None,
    }


@pytest.fixture
def renewable_generator_params():
    return {
        "gen_name": "Testing_Renewable_Generator",
        "bus": "bus5",
        "p_min": 0,
        "p_max": 200,
        "p_cost": 0,
        "fixed_commitment": None,
    }


@pytest.fixture
def generator_data_object(generator_params):
    return ThermalGeneratorModelData(**generator_params)


@pytest.fixture
def renewable_generator_data_object(renewable_generator_params):
    return RenewableGeneratorModelData(**renewable_generator_params)


@pytest.mark.unit
def test_create_thermal_model_data_object(generator_params, generator_data_object):

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
def test_create_renewable_model_data_object(
    renewable_generator_params, renewable_generator_data_object
):

    for name, value in renewable_generator_data_object:
        if name == "generator_type":
            assert value == "renewable"
        else:
            assert renewable_generator_params[name] == value


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
        ThermalGeneratorModelData(**generator_params)


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
        ThermalGeneratorModelData(**generator_params)


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
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_fixed_commitment_value(generator_params):
    generator_params["fixed_commitment"] = 5
    with pytest.raises(
        ValueError, match=r"^(Value for generator fixed commitment must be one of)"
    ):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_generator_name(generator_params):
    generator_params["gen_name"] = 102_111
    with pytest.raises(TypeError, match=r".*gen_name must be str.*"):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_bus_name(generator_params):
    generator_params["bus"] = 102
    with pytest.raises(TypeError, match=r".*bus must be str.*"):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
@pytest.mark.parametrize(
    "param_name", ["production_cost_bid_pairs", "startup_cost_pairs"]
)
def test_empty_bid_pairs(param_name, generator_params):
    generator_params[param_name] = None
    with pytest.raises(
        ValueError, match=r"Empty production|startup cost pairs are provided."
    ):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
def test_bid_missing_pmin(generator_params):
    generator_params["production_cost_bid_pairs"] = generator_params[
        "production_cost_bid_pairs"
    ][1:]
    with pytest.raises(
        ValueError, match=r"^(The first power output in the bid should be the Pmin)"
    ):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
def test_bid_missing_pmax(generator_params):
    generator_params["production_cost_bid_pairs"].pop()
    with pytest.raises(
        ValueError, match=r"^(The last power output in the bid should be the Pmax)"
    ):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
def test_invalid_start_up_bid(generator_params):
    generator_params["startup_cost_pairs"] = [(0, 0)]
    with pytest.raises(
        ValueError,
        match=r"^(The first startup lag should be the same as minimum down time)",
    ):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
def test_model_data_iterator(generator_data_object):

    expected_param_names = [
        "gen_name",
        "bus",
        "p_min",
        "p_max",
        "min_down_time",
        "min_up_time",
        "ramp_up_60min",
        "ramp_down_60min",
        "shutdown_capacity",
        "startup_capacity",
        "p_cost",
        "startup_cost",
        "fixed_commitment",
        "generator_type",
        "initial_status",
        "initial_p_output",
    ]
    iter_result = [name for name, value in generator_data_object]

    expected_param_names.sort()
    iter_result.sort()

    pyo_unittest.assertStructuredAlmostEqual(
        first=expected_param_names, second=iter_result
    )


@pytest.mark.unit
@pytest.mark.parametrize(
    "value, error, msg",
    [
        ("1", TypeError, "Value for initial_status shoulde be real numbers"),
        (0, ValueError, "Value for initial_status cannot be zero"),
    ],
)
def test_invalid_initial_status(value, error, msg, generator_params):

    generator_params["initial_status"] = value
    with pytest.raises(error, match=msg):
        ThermalGeneratorModelData(**generator_params)


@pytest.mark.unit
@pytest.mark.parametrize(
    "initial_status, initial_p_output, error, msg",
    [
        (9, "1", TypeError, "Value for initial_p_output shoulde be real numbers"),
        (
            9,
            29,
            ValueError,
            r".*so the initial power output should at least be p_min.*",
        ),
        (-9, 30, ValueError, r".*so the initial power output should at 0.*"),
    ],
)
def test_invalid_initial_p_output(
    initial_status, initial_p_output, error, msg, generator_params
):

    generator_params["initial_p_output"] = initial_p_output
    generator_params["initial_status"] = initial_status
    with pytest.raises(error, match=msg):
        ThermalGeneratorModelData(**generator_params)
