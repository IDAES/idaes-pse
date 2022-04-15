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
import numpy as np
import math


def split_training_validation(dataframe, training_fraction, seed=None):
    """
    Randomly split the dataframe into training and validation data

    Args:
       dataframe : pandas DataFrame
          The dataframe to split
       training_fraction : float between 0 < 1
          The fraction of the overall dataframe (# rows) to include
          as training data. The rest will be returned as validation data
       seed : None or int
          seed used for the random number generator. Use default behavior
          for DataFrame.sample if seed is None

    Returns:
       tuple : (training_dataframe, validation_dataframe)
    """
    return split_dataframe(dataframe, [training_fraction], seed)


def split_training_validation_testing(
    dataframe, training_fraction, validation_fraction, seed=None
):
    """
    Randomly split the dataframe into training, validation, and testing data

    Args:
       dataframe : pandas DataFrame
          The dataframe to split
       training_fraction : float between 0 < 1
          The fraction of the overall dataframe (# rows) to include
          as training data.
          The sum of training_fraction and validation_fraction must be less than 1
       validation_fraction : float between 0 < 1
          The fraction of the overall dataframe (# rows) to include
          as validation data. The rest will be returned as validation data.
          The sum of training_fraction and validation_fraction must be less than 1
       seed : None or int
          seed used for the random number generator. Use default behavior
          for DataFrame.sample if seed is None

    Returns:
       tuple : (training_dataframe, validation_dataframe, testing_dataframe)
    """
    return split_dataframe(dataframe, [training_fraction, validation_fraction], seed)


def split_dataframe(dataframe, fractions, seed=None):
    """
    Randomly splits the dataframe into multiple dataframes.

    Args:
       dataframe: pandas DataFrame
          The dataframe to split
       fractions: list of floats between 0 < 1
          The fraction of the data to include in each dataframe. The list of fractions
          must sum to < 1. If fractions has length N, then N+1 dataframes will be returned
          where the fraction for the last dataframe is 1-sum(fractions).
       seed : None or int
          seed for the random number generator. If None, generator is not seeded.

    Returns:
       tuple : (DataFrame, ...)
    """
    assert sum(fractions) < 1.0

    # note seed=None is the default value for random_state (e.g., not seeded)
    shuffled_df = dataframe.sample(frac=1, random_state=seed).reset_index(drop=True)

    dfs = np.split(
        shuffled_df, [math.floor(f * len(shuffled_df)) for f in np.cumsum(fractions)]
    )
    # reset all the indices
    for df in dfs:
        df.reset_index(drop=True, inplace=True)
    return dfs
