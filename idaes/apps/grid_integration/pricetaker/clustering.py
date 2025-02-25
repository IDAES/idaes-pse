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
from scipy.interpolate import splrep, splev
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
    # Converting the data to a list (Needed if the input is not a list, e.g., pd.DataFrame)
    raw_data = list(raw_data)

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
    raw_data: list,
    horizon_length: int,
    n_clusters: int,
    seed: int = 42,
    eps: int = 1e-4,
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

        eps: int,
            Centroid values below this threshold are set to 0 to limit noise in the data

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
    centroids = centroids * (abs(centroids) >= eps)

    # Create dicts for lmp data and the weight of each cluster
    # By default, the data is of type numpy.int or numpy.float.
    # Converting the data to python int/float. Otherwise, Pyomo complains!
    num_days, num_time_periods = centroids.shape
    lmp_data_clusters = {
        d + 1: {t + 1: float(centroids[d, t]) for t in range(num_time_periods)}
        for d in range(num_days)
    }
    weights = {d + 1: int(sum(labels == d)) for d in range(num_days)}

    return lmp_data_clusters, weights


def get_optimal_num_clusters(
    samples: Union[pd.DataFrame, np.array],
    kmin: int = 2,
    kmax: int = 30,
    method: str = "silhouette",
    sensitivity: int = 1,
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

        method : str,
            Method for obtaining the optimal number of clusters.
            Supported methods are elbow and silhouette

        sensitivity : int,
            difficulty of detecting knees, where 0 detects knees the easiest

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

    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=seed).fit(samples)
        inertia_values.append(kmeans.inertia_)

    if method == "silhouette":
        mean_silhouette = []
        mean_silhouette.append(silhouette_score(samples, kmeans.labels_))

        # Identify the optimal number of clusters
        max_index = mean_silhouette.index(max(mean_silhouette))
        n_clusters = k_values[max_index]

    elif method == "elbow":
        # Invert inertia values such that plot is concave down and increasing
        inverted_inertia_values = [-y for y in inertia_values]
        n_clusters = locate_elbow(
            x=k_values,
            y=inverted_inertia_values,
            sensitivity=sensitivity,
            kmin=kmin,
            kmax=kmax,
        )

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


def locate_elbow(
    x: list,
    y: list,
    sensitivity: int = 1,
    kmin: int = 2,
    kmax: int = 30,
):
    """
    Identifies the elbow/knee for the input/output data that is concave down and increasing

    Args:
        x : list,
            List of independent variables (i.e. k_values)

        y : list,
            List of dependent variables (i.e. inertia values)

        sensitivity : int,
            difficulty of detecting knees, where 0 detects knees the easiest

        kmin : int,
            Minimum number of clusters

        kmax : int,
            Maximum number of clusters

    Returns:
        opt_x : float
            Optimal x at which curvature changes significantly
    """
    # The implementation is based on
    # Ville Satopaa, Jeannie Albrecht, David Irwin, Barath Raghavan
    # Finding a “Kneedle” in a Haystack:
    # Detecting Knee Points in System Behavior
    # https://raghavan.usc.edu/papers/kneedle-simplex11.pdf

    # TODO: See if make_splrep works and generates similar results
    # Use a smoothing spline that retains the data's original shape
    spline = splrep(x, y, k=3)
    k_smooth = np.linspace(kmin, kmax + 1, 100)
    inertia_smooth = splev(k_smooth, spline)

    k_norm = _normalize_values(k_smooth)
    inertia_norm = _normalize_values(inertia_smooth)

    # Compute the set of differences (x, y) -> (x, y-x)
    inertia_diff = []
    for i in range(len(inertia_norm)):
        set_of_differences = inertia_norm[i] - k_norm[i]
        inertia_diff.append(set_of_differences)

    # Identify local maxima
    local_maxima = [
        (k_norm[i], inertia_diff[i])
        for i in range(1, len(inertia_diff) - 1)
        if inertia_diff[i] > inertia_diff[i - 1]
        and inertia_diff[i] > inertia_diff[i + 1]
    ]

    # Calculate optimal # of clusters based on the # of local maxima
    threshold_values = []
    threshold_triggered = False
    n_clusters = 0

    if len(local_maxima) == 0:
        n_clusters = 0
        _logger.warning(
            "The optimal number of cluster cannot be determined for this dataset."
        )
    elif len(local_maxima) == 1:
        n_clusters_norm = local_maxima[0][0]
        ind = k_norm.index(n_clusters_norm)
        n_clusters = x[ind]
    else:
        n = len(local_maxima)
        summation = 0
        for i in range(0, n - 1):
            summation += local_maxima[i + 1][0] - local_maxima[i][0]
        # For each local maxima, compute the threshold and determine the index of the current and next local maxima
        for i in range(0, n - 1):
            threshold = local_maxima[i][1] - (sensitivity * summation / (n - 1))
            threshold_values.append(threshold)
            l_max_index = inertia_diff.index(local_maxima[i][1])
            next_l_max_index = inertia_diff.index(local_maxima[i + 1][1])
            for j in range(l_max_index + 1, next_l_max_index):
                if inertia_diff[j] < threshold:
                    threshold_triggered = True
                    normalized_optimal_n_clusters = local_maxima[i][0]
                    index = k_norm.index(normalized_optimal_n_clusters)
                    n_clusters = x[index]
            if threshold_triggered:
                break
        # If optimal # of clusters cannot be identified, use the first local maxima
        if not threshold_triggered:
            normalized_optimal_n_clusters = local_maxima[0][0]
            index = k_norm.index(normalized_optimal_n_clusters)
            n_clusters = x[index]
            _logger.warning(
                "The number of optimal clusters could not be accurately identified. "
                "Consider lowering the sensitivity."
            )

    return n_clusters


def _normalize_values(values):
    normalized_values = []

    for i in values:
        max_value = max(values)
        min_value = min(values)
        normalized_value = (i - min_value) / (max_value - min_value)
        normalized_values.append(normalized_value)

    return normalized_values
