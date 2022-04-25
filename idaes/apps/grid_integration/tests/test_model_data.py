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
    attr_names = ["default_bids", "default_startup_bids"]
    param_names = ["production_cost_bid_pairs", "startup_cost_pairs"]

    for a_name, p_name in zip(attr_names, param_names):
        bids = getattr(generator_data_object, a_name)
        expected_bid = {key: val for key, val in generator_params[p_name]}
        pyo_unittest.assertStructuredAlmostEqual(first=bids, second=expected_bid)
