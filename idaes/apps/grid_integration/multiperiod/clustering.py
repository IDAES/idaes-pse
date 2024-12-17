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

"""
Contains functions for clustering the price signal into representative days/periods.
"""

import numpy as np
import pandas as pd
import scipy.cluster.vq as spc
import matplotlib.pyplot as plt
import idaes.logger as idaeslog

_logger = idaeslog.getLogger(__name__)


def generate_daily_data(raw_data: list, horizon_length: int):
    """
    Generates the daily data in a usable format from the raw data provided.

    Args:
        raw_data: list,
            Columnar data for a given LMP signal

        horizon_length: int,
            Length of each representative day/period

    Returns:
        daily_data: pd.DataFrame,
            Input data arranged in a DataFrame, where columns correspond to hours
            and rows correspond to days.
    """
    if horizon_length > len(raw_data):
        raise ValueError(
            f"Horizon length {horizon_length} exceeds the length of the "
            f"price signal ({len(raw_data)})"
        )

    elements_ignored = len(raw_data) % horizon_length
    if elements_ignored:
        _logger.warning(
            f"Length of the price signal is not an integer multiple of horizon_length.\n"
            f"\tIgnoring the last {elements_ignored} elements in the price signal."
        )

    daily_data = {
        j: raw_data[((j - 1) * horizon_length) : (j * horizon_length)]
        for j in range(1, (len(raw_data) // horizon_length) + 1)
    }

    # DataFrame arranges the data for each day as a column. Since most clustering techniques
    # require different samples as rows, we take the transpose before returning the data.
    return pd.DataFrame(daily_data).transpose()


def cluster_lmp_data(
    raw_data: list, n_clusters: int, horizon_length: int, seed: int = 42
):
    """
    Clusters the given price signal into n_clusters using the k-means clustering
    technique.

    Args:
        raw_data: list,
            Columnar data for a given LMP signal.

        n_clusters: int,
            Number of clusters desired for the data (representative days).

        horizon_length: int,
            Length of each cluster (representative day/period)

        seed: int,
            Seed value for initializing random number generator within Kmeans

    Returns:
        lmp_data_clusters: dict
            A dictionary of representative day LMP data, indices are indexed
                by integers starting at 1. Example: ::

                    {1: {1: 4, 2: 3, 3: 5},
                    2: {1: 1, 2: 7, 3: 3}}


        weights: dict
            A dictionary of weights for each representative day, indexed the
                same way as lmp_data. Example: ::

                    {1: 45, 2: 56}
    """
    # reconfiguring raw data
    daily_data = generate_daily_data(raw_data, horizon_length)

    # KMeans clustering with the optimal number of clusters
    centroids, labels = spc.kmeans2(daily_data, n_clusters, seed=seed)

    # Set any centroid values that are < 1e-4 to 0 to avoid noise
    centroids = centroids * (abs(centroids) >= 1e-4)

    # Create dicts for lmp data and the weight of each cluster
    num_days, num_time_periods = centroids.shape
    lmp_data_clusters = {
        d + 1: {t + 1: centroids[d, t] for t in range(num_time_periods)}
        for d in range(num_days)
    }
    weights = {d + 1: sum(labels == d) for d in range(num_days)}

    return lmp_data_clusters, weights


def _compute_sse(data, centroids, idx):
    """
    Function used to compute the inertia (sum of square errors) for k clusters.

    Args:
        data:      Columnar data for a given LMP signal
        centroids: Array of k centroids
        idx:       Index for data

    Returns:
        inertia: Sum of square errors for k clusters
    """
    inertia = 0
    for i, centroid in enumerate(centroids):
        cluster_points = data[idx == i]
        inertia += np.sum((cluster_points - centroid) ** 2)
    return inertia


def get_optimal_num_clusters(
    daily_data,
    kmin: int = 4,
    kmax: int = 30,
    generate_elblow_plot: bool = False,
    seed: int = 42,
):
    """
    Determines the appropriate number of clusters needed for a
    given price signal.

    Args:
        daily_data: LMP signal grouped by days (output of generate_daily_data function)
        kmin:       minimum number of clusters
        kmax:       maximum number of clusters

    Returns:
        n_clusters:     the optimal number of clusters for the given data
        inertia_values: within-cluster sum-of-squares
    """
    if kmin >= kmax:
        raise ValueError(f"kmin must be less than kmax, but {kmin} >= {kmax}")

    k_values = list(range(kmin, kmax + 1))
    inertia_values = []

    for k in k_values:
        centroids, _ = spc.kmeans2(daily_data, k)
        idx, _ = spc.vq(daily_data, centroids)

        # Compute the inertia (SSE) for k clusters
        inertia = _compute_sse(daily_data, centroids, idx)
        inertia_values.append(inertia)

    # Calculate the second derivative
    first_deriv = np.diff(inertia_values)
    second_deriv = np.diff(first_deriv)

    # Determine the optimal number of clusters
    # The +2 accounts for the dimension being reduced twice by derivatives
    n_clusters = np.argmin(second_deriv) + 2

    if n_clusters is None:
        raise ValueError(
            "Could not find elbow point for the given kmin, kmax. "
            "Consider increasing the range of kmin, kmax."
        )

    _logger.info(f"Optimal number of clusters is: {n_clusters}")

    if int(n_clusters) + 2 >= kmax:
        _logger.warning(
            f"Optimal number of clusters is close to kmax: {kmax}. Consider increasing kmax."
        )

    if generate_elblow_plot:
        plt.plot(k_values, inertia_values)
        plt.axvline(x=n_clusters, color="red", linestyle="--", label="Elbow")
        plt.xlabel("Number of clusters")
        plt.ylabel("Inertia")
        plt.title("Elbow Method")
        plt.xlim(kmin, kmax)
        plt.grid()

    return int(n_clusters), inertia_values
