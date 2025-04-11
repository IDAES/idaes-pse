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

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from idaes.apps.grid_integration.pricetaker.clustering import (
    generate_daily_data,
    get_optimal_num_clusters,
    cluster_lmp_data,
    sklearn_avail,
    _normalize_values,
)
import idaes.logger as idaeslog

pytest.importorskip("sklearn", reason="sklearn not available")


@pytest.fixture(name="sample_data")
def sample_data_fixture():
    """Returns a sample price signal for testing"""
    file_path = "lmp_data.csv"
    data = pd.read_csv(file_path)
    return data


@pytest.fixture(name="dummy_data")
def dummy_data_fixture():
    """Returns dummy data for clustering"""
    # Test Data Set: [0,0], [0,1], [1,1], [1, 0], [10, 0], [10, 1],
    # [11, 1], [11, 0], [5,10], [5, 11], [6, 11], [6, 10]
    # Centroids for this set must be [0.5, 0.5], [10.5, 0.5], [5.5, 10.5]
    data = [
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        0,
        10,
        0,
        10,
        1,
        11,
        1,
        11,
        0,
        5,
        10,
        5,
        11,
        6,
        11,
        6,
        10,
    ]
    return data


@pytest.mark.unit
def test_daily_data_logger_message1():
    """Tests ValueError in generate_daily_data function"""
    with pytest.raises(
        ValueError,
        match="Horizon length 24 exceeds the price signal length of 4",
    ):
        generate_daily_data(raw_data=[1, 2, 3, 4], horizon_length=24)


@pytest.mark.unit
def test_daily_data_logger_message2(caplog):
    """Tests non-integer-multiple length warning in generate_daily_data function"""

    with caplog.at_level(idaeslog.WARNING):
        daily_data = generate_daily_data(raw_data=list(range(8)), horizon_length=3)
        assert (
            "Length of the price signal is not an integer multiple of horizon_length.\n"
            "\tIgnoring the last 2 elements in the price signal."
        ) in caplog.text

    assert isinstance(daily_data, pd.DataFrame)
    assert daily_data.shape == (2, 3)


@pytest.mark.unit
def test_generate_daily_data():
    """Tests the generate_daily_data function"""

    raw_data = list(range(9))
    daily_data = generate_daily_data(raw_data, horizon_length=3)
    assert daily_data.shape == (3, 3)
    assert daily_data.equals(
        pd.DataFrame({1: [0, 1, 2], 2: [3, 4, 5], 3: [6, 7, 8]}).transpose()
    )


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_cluster_lmp_data(dummy_data):
    """Tests the cluster_lmp_data function"""

    data_clusters, weights = cluster_lmp_data(
        raw_data=dummy_data, horizon_length=2, n_clusters=3, seed=42
    )
    assert weights == {1: 4, 2: 4, 3: 4}
    assert data_clusters == {
        1: {1: 10.5, 2: 0.5},
        2: {1: 5.5, 2: 10.5},
        3: {1: 0.5, 2: 0.5},
    }


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_optimal_num_clusters_elbow(sample_data):
    """Tests the get_optimal_num_clusters function"""
    samples = generate_daily_data(
        sample_data["BaseCase_2030"].tolist(), horizon_length=24
    )
    n_clusters = get_optimal_num_clusters(
        samples=samples,
        kmin=4,
        kmax=30,
        method="elbow",
        generate_elbow_plot=True,
        seed=42,
    )

    assert n_clusters == 7

    # Test that a figure was created
    assert plt.gcf() is not None
    # Test that axes were created
    assert plt.gca() is not None

    plt.close("all")


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_optimal_num_clusters_silhouette(dummy_data):
    """Tests the get_optimal_num_clusters function"""
    samples = generate_daily_data(dummy_data, horizon_length=2)
    n_clusters = get_optimal_num_clusters(
        samples=samples,
        kmin=2,
        kmax=7,
        method="silhouette",
        generate_elbow_plot=True,
    )

    assert n_clusters == 2

    # Test that a figure was created
    assert plt.gcf() is not None
    # Test that axes were created
    assert plt.gca() is not None

    plt.close("all")


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_optimal_clusters_logger_message1(dummy_data):
    """Tests the error message associated with kmin >= kmax"""

    samples = generate_daily_data(dummy_data, 2)
    with pytest.raises(
        ValueError,
        match="kmin must be less than kmax, but 15 >= 10",
    ):
        get_optimal_num_clusters(
            samples=samples,
            kmin=15,
            kmax=10,
        )


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_optimal_clusters_logger_message2(dummy_data):
    """Tests the error message associated with kmin < 2"""

    samples = generate_daily_data(dummy_data, 2)
    with pytest.raises(
        ValueError,
        match="kmin must be > 1. Received kmin = 1",
    ):
        get_optimal_num_clusters(
            samples=samples,
            kmin=1,
            kmax=10,
        )


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_optimal_clusters_logger_message3(dummy_data):
    """Tests the error message associated with kmax >= len(samples)"""

    samples = generate_daily_data(dummy_data, 2)
    with pytest.raises(
        ValueError,
        match="kmax must be < len\\(samples\\). Received kmax = 12",
    ):
        get_optimal_num_clusters(
            samples=samples,
            kmin=2,
            kmax=12,
        )


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_optimal_clusters_logger_message4(dummy_data):
    """Tests the error message associated with unrecognized method"""

    samples = generate_daily_data(dummy_data, 2)
    with pytest.raises(
        ValueError,
        match=(
            "Unrecognized method foo for optimal number of clusters."
            "\tSupported methods include elbow and silhouette."
        ),
    ):
        get_optimal_num_clusters(
            samples=samples,
            kmin=2,
            kmax=5,
            method="foo",
        )


@pytest.mark.unit
@pytest.mark.skipif(
    not sklearn_avail, reason="sklearn (optional dependency) not available"
)
def test_optimal_clusters_logger_message5(dummy_data, caplog):
    """Tests the logger warning associated with n_clusters close to kmax"""

    samples = generate_daily_data(dummy_data, 2)
    with caplog.at_level(idaeslog.WARNING):
        get_optimal_num_clusters(samples, kmin=4, kmax=5, seed=20)
        assert (
            "Optimal number of clusters is close to kmax: 5. Consider increasing kmax."
            in caplog.text
        )


@pytest.mark.unit
def test_normalize_values():
    """Tests the normalize_values function"""

    data = [1, 2, 3, 4, 5]
    normalized_data = _normalize_values(data)
    assert normalized_data == [0, 0.25, 0.5, 0.75, 1]
