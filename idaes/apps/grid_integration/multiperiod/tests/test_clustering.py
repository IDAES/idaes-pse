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
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from idaes.apps.grid_integration.multiperiod.clustering import (
    generate_daily_data,
    get_optimal_num_clusters,
    cluster_lmp_data,
    sklearn_avail,
    kneed_avail,
)

pytest.importorskip("sklearn", reason="sklearn not available")
pytest.importorskip("kneed", reason="kneed not available")


@pytest.fixture
def sample_data():
    DATA_DIR = Path(__file__).parent
    file_path = DATA_DIR / "FLECCS_princeton.csv"
    data = pd.read_csv(file_path)
    return data


@pytest.mark.unit
def test_daily_data_size(sample_data):
    # Generate price data for each hour of every day in the data
    daily_data = generate_daily_data(sample_data["BaseCaseTax"], horizon_length=24)

    # Check that there is a row for each horizon length in a representative day (24 x 365)
    assert len(daily_data) == 24
    assert len(daily_data.transpose()) == 365


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
@pytest.mark.skipif(not kneed_avail, reason="kneed (optional dependency) not available")
def test_determine_optimal_num_clusters(sample_data):

    daily_data = generate_daily_data(sample_data["BaseCaseTax"], horizon_length=24)
    n_clusters, inertia_values = get_optimal_num_clusters(daily_data, kmax=30, seed=20)

    assert n_clusters == 11


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
@pytest.mark.skipif(not kneed_avail, reason="kneed (optional dependency) not available")
def test_elbow_plot(sample_data):

    daily_data = generate_daily_data(sample_data["BaseCaseTax"], horizon_length=24)
    n_clusters, inertia_values = get_optimal_num_clusters(
        daily_data, kmax=30, generate_elbow_plot=True
    )

    # Test that a figure was created
    assert plt.gcf() is not None
    # Test that axes were created
    assert plt.gca() is not None
    # Test that the plot has data
    assert plt.gca().has_data()

    plt.close("all")


@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
@pytest.mark.skipif(not kneed_avail, reason="kneed (optional dependency) not available")
@pytest.mark.unit
def test_cluster_lmp_data(sample_data):
    daily_data = generate_daily_data(sample_data["BaseCaseTax"], horizon_length=24)
    n_clusters, inertia_values = get_optimal_num_clusters(daily_data, kmax=30, seed=20)
    lmp_data, weights = cluster_lmp_data(
        sample_data["BaseCaseTax"],
        n_clusters=n_clusters,
        horizon_length=24,
    )

    sum_of_weights = sum(weights.values())
    assert sum_of_weights == 365

    assert len(lmp_data) == n_clusters
