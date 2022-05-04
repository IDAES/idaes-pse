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


def compute_fit_metrics(surrogate, dataframe):
    """
    Compute the fit metrics for the surrogate against the
    data in dataframe. The surrogate will be evaluated
    for the input values in the dataframe, and the results
    will be compared against the output values in the dataframe
    to provide the metrics.

    Args:
       surrogate : surrogate object (derived from SurrogateBase)
          This is the surrogate object we want to evaluate for the comparison
       dataframe : pandas DataFrame
          The dataframe that contains the inputs and outputs we want to use
          in the evaluation.

    Returns:
        dict-of-dicts with outer keys representing output labels and inner keys
        representing metrics for that output.
    """
    y = dataframe[surrogate.output_labels()]
    f = surrogate.evaluate_surrogate(dataframe)
    assert f.columns.to_list() == surrogate.output_labels()

    y_mean = y.mean(axis=0)
    SST = ((y - y_mean) ** 2).sum(axis=0)
    SSE = ((y - f) ** 2).sum(axis=0)

    R2 = 1 - SSE / SST
    MAE = (y - f).abs().mean(axis=0)
    maxAE = (y - f).abs().max(axis=0)
    MSE = ((y - f) ** 2).mean(axis=0)
    RMSE = MSE**0.5

    # Reorder indices to have output first
    metrics = {}
    for o in surrogate.output_labels():
        metrics[o] = {
            "RMSE": RMSE[o],
            "MSE": MSE[o],
            "MAE": MAE[o],
            "maxAE": maxAE[o],
            "SSE": SSE[o],
            "R2": R2[o],
        }

    return metrics
