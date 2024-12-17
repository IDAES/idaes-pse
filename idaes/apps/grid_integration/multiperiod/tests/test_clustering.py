#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

import pytest
import pandas as pd

from idaes.apps.grid_integration.multiperiod.clustering import (
    generate_daily_data,
    get_optimal_num_clusters,
    cluster_lmp_data,
)


@pytest.mark.unit
def test_generate_daily_data(caplog):
    """Tests the generate_daily_data function"""
    # Test the ValueError exception
    raw_data = [1, 2, 3, 4]
    with pytest.raises(
        ValueError,
        match="Horizon length 10 exceeds the length of the price signal \\(4\\)",
    ):
        generate_daily_data(raw_data, horizon_length=10)

    # Test the logger warning
    raw_data = list(range(8))
    daily_data = generate_daily_data(raw_data, horizon_length=3)
    assert (
        "Length of the price signal is not an integer multiple of horizon_length.\n"
        "\tIgnoring the last 2 elements in the price signal."
    ) in caplog.text
    assert isinstance(daily_data, pd.DataFrame)
    assert daily_data.shape == (2, 3)

    # Test the output
    raw_data = list(range(9))
    daily_data = generate_daily_data(raw_data, horizon_length=3)
    assert daily_data.shape == (3, 3)
    assert daily_data.equals(
        pd.DataFrame({1: [0, 1, 2], 2: [3, 4, 5], 3: [6, 7, 8]}).transpose()
    )
