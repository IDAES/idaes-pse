import numpy as np
import math

def split_training_testing(dataframe, training_fraction, seed=None):
    """
    Randomly split the dataframe into training and testing data

    Args:
       dataframe : pandas DataFrame
          The dataframe to split
       training_fraction : float between 0 < 1
          The fraction of the overall dataframe (# rows) to include
          as training data. The rest will be returned as testing data
       seed : None or int
          seed used for the random number generator. Use default behavior
          for DataFrame.sample if seed is None

    Returns:
       tuple : (training_dataframe, testing_dataframe)
    """
    return _split_dataframe(dataframe, [training_fraction], seed)

def split_training_testing_validation(dataframe, training_fraction, testing_fraction, seed=None):
    """
    Randomly split the dataframe into training and testing data

    Args:
       dataframe : pandas DataFrame
          The dataframe to split
       training_fraction : float between 0 < 1
          The fraction of the overall dataframe (# rows) to include
          as training data.
          The sum of training_fraction and testing_fraction must be less than 1
       testing_fraction : float between 0 < 1
          The fraction of the overall dataframe (# rows) to include
          as testing data. The rest will be returned as validation data.
          The sum of training_fraction and testing_fraction must be less than 1
       seed : None or int
          seed used for the random number generator. Use default behavior
          for DataFrame.sample if seed is None

    Returns:
       tuple : (training_dataframe, testing_dataframe)
    """
    return _split_dataframe(dataframe, [training_fraction, testing_fraction], seed)

def _split_dataframe(dataframe, fractions, seed=None):
    """
    Randomly splits the dataframe into two dataframes.

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
    if seed is not None:
        shuffled_df = dataframe.sample(frac=1, random_state=seed).reset_index(drop=True)
    else:
        shuffled_df = dataframe.sample(frac=1).reset_index(drop=True)

    return np.split(shuffled_df, [math.floor(f*len(shuffled_df)) for f in np.cumsum(fractions)])
