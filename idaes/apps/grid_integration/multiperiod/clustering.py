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
import matplotlib.pyplot as plt
import idaes.logger as idaeslog
from pyomo.common.dependencies import attempt_import

sklearn, sklearn_avail = attempt_import("sklearn")
kneed, kneed_avail = attempt_import("kneed")

if sklearn_avail:
    from sklearn.cluster import KMeans

if kneed_avail:
    from kneed import KneeLocator

_logger = idaeslog.getLogger(__name__)


def generate_daily_data(raw_data: list, horizon_length: int):
    """
    Function used to generate the daily data in a usable format
    from the raw data provided.

    Args:
        raw_data:   Columnar data for a given LMP signal

    Returns:
        daily_data: Correctly formatted daily LMP data for later use
    """
    if horizon_length > len(raw_data):
        raise ValueError(
            f"Horizon length {horizon_length} exceeds the price signal length of {len(raw_data)}"
        )

    elements_ignored = len(raw_data) % horizon_length
    if elements_ignored:
        _logger.warning(
            f"Length of the price signal is not an integer multiple of horizon_length.\n"
            f"\tIgnoring the last {elements_ignored} elements in the price signal."
        )

    day_list = list(range(1, (len(raw_data) // horizon_length) + 1))

    daily_data = pd.DataFrame(columns=day_list)

    # Extracting data to populate empty dataframe
    i = 0
    j = horizon_length
    day = 1

    while j <= len(raw_data):
        daily_data[day] = raw_data[i:j].reset_index(drop=True)
        i = j
        j = j + horizon_length
        day = day + 1

    return daily_data


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
    kmeans = KMeans(n_clusters=n_clusters).fit(daily_data.transpose())
    centroids = kmeans.cluster_centers_
    labels = kmeans.labels_

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
    generate_elbow_plot: bool = False,
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

    np.random.seed(seed)

    for k in k_values:
        kmeans = KMeans(n_clusters=k).fit(daily_data.transpose())
        inertia_values.append(kmeans.inertia_)

    # Identify the "elbow point"
    elbow_point = KneeLocator(
        k_values, inertia_values, curve="convex", direction="decreasing"
    )
    n_clusters = elbow_point.knee

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

    if generate_elbow_plot:
        plt.plot(k_values, inertia_values)
        plt.axvline(x=n_clusters, color="red", linestyle="--", label="Elbow")
        plt.xlabel("Number of clusters")
        plt.ylabel("Inertia")
        plt.title("Elbow Method")
        plt.xlim(kmin, kmax)
        plt.grid()

    return int(n_clusters), inertia_values
