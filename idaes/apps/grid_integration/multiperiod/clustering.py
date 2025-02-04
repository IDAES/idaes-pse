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

from typing import Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyomo.common.dependencies import attempt_import
import idaes.logger as idaeslog

sklearn, sklearn_avail = attempt_import("sklearn")

if sklearn_avail:
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score

_logger = idaeslog.getLogger(__name__)


def generate_daily_data(raw_data: list, horizon_length: int):
    """
    Function used to generate the daily data in a usable format
    from the raw data provided.

    Args:
        raw_data : list,
            Columnar data for a given LMP signal

        horizon_length : int,
            Length of each representative day/period

    Returns:
        daily_data: pd.DataFrame,
            Correctly formatted daily LMP data for later use
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

    daily_data = {
        j: raw_data[((j - 1) * horizon_length) : (j * horizon_length)]
        for j in range(1, (len(raw_data) // horizon_length) + 1)
    }

    # DataFrame arranges the data for each day as a column. Since most clustering techniques
    # require different samples as rows, we take the transpose before returning the data.
    return pd.DataFrame(daily_data).transpose()


def cluster_lmp_data(
    raw_data: list, horizon_length: int, n_clusters: int, seed: int = 42
):
    """
    Clusters the given price signal into n_clusters using the k-means clustering
    technique.

    Args:
        raw_data: list,
            Columnar data for a given LMP signal.

        horizon_length: int,
            Length of each cluster (representative day/period)

        n_clusters: int,
            Number of clusters desired for the data (representative days).

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
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed).fit(daily_data)
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


def get_optimal_num_clusters(
    samples: Union[pd.DataFrame, np.array],
    kmin: int = 1,
    kmax: int = 30,
    method: str = "elbow",
    generate_elbow_plot: bool = False,
    seed: int = 42,
):
    """
    Determines the appropriate number of clusters needed for a
    given price signal.

    Args:
        samples: pd.DataFrame | np.array,
            Set of points with rows containing samples and columns
            containing features

        kmin : int,
            Minimum number of clusters

        kmax : int,
            Maximum number of clusters

        generate_elbow_plot : bool,
            If True, generates an elbow plot for inertia as a function of
            number of clusters

        seed : int,
            Seed value for random number generator

    Returns:
        n_clusters: int,
            The optimal number of clusters for the given data
    """
    if kmin >= kmax:
        raise ValueError(f"kmin must be less than kmax, but {kmin} >= {kmax}")

    # For silhouette method, kmin must be 2. So, we require kmin >= 2 in general
    if kmin <= 1:
        raise ValueError(f"kmin must be > 1. Received kmin = {kmin}")

    if kmax >= len(samples):
        raise ValueError(f"kmax must be < len(samples). Received kmax = {kmax}")

    k_values = list(range(kmin, kmax + 1))
    inertia_values = []
    mean_silhouette = []

    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=seed).fit(samples)
        inertia_values.append(kmeans.inertia_)

        if method == "silhouette":
            # Calculate the average silhouette score, if the chosen method
            # is silhouette
            mean_silhouette.append(silhouette_score(samples, kmeans.labels_))

    # Identify the optimal number of clusters
    if method == "elbow":
        n_clusters = locate_elbow(k_values, inertia_values)

    elif method == "silhouette":
        max_index = mean_silhouette.index(max(mean_silhouette))
        n_clusters = k_values[max_index]

    else:
        raise ValueError(
            f"Unrecognized method {method} for optimal number of clusters."
            f"\tSupported methods include elbow and silhouette."
        )

    _logger.info(f"Optimal number of clusters is: {n_clusters}")

    if n_clusters + 2 >= kmax:
        _logger.warning(
            f"Optimal number of clusters is close to kmax: {kmax}. "
            f"Consider increasing kmax."
        )

    if generate_elbow_plot:
        plt.plot(k_values, inertia_values)
        plt.axvline(x=n_clusters, color="red", linestyle="--", label="Elbow")
        plt.xlabel("Number of clusters")
        plt.ylabel("Inertia")
        plt.title("Elbow Method")
        plt.xlim(kmin, kmax)
        plt.show()

    return n_clusters


def locate_elbow(x: list, y: list):
    """
    Identifies the elbow/knee for the input/output data

    Args:
        x : list
            List of independent variables
        y : list
            List of dependent variables

    Returns:
        opt_x : float
            Optimal x at which curvature changes significantly
    """
    # The implementation is based on
    # Ville Satopaa, Jeannie Albrecht, David Irwin, Barath Raghavan
    # Finding a “Kneedle” in a Haystack:
    # Detecting Knee Points in System Behavior
    # https://raghavan.usc.edu/papers/kneedle-simplex11.pdf

    return len(np.array(x)) + len(np.array(y))
