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
import pandas as pd


class OffsetScaler(object):
    @staticmethod
    def create_normalizing_scaler(dataframe):
        """
        Creates a scaling object that normalizes the data between 0 and 1

        Args:
           dataframe: pandas DataFrame
              The dataframe containing the data (usually the training data)
              that will be used to compute the scaling factor and offset
        """
        expected_columns = list(dataframe.columns)
        offset = dataframe.min()
        factor = dataframe.max() - dataframe.min()
        return OffsetScaler(expected_columns, offset, factor)

    @staticmethod
    def create_from_mean_std(dataframe):
        """
        Creates a scaling object using the mean as the offset and the
        standard devation as the as the factor

        Args:
           dataframe: pandas DataFrame
              The dataframe containing the data (usually the training data) that will
              be used to compute the mean and standard devation for the scaler
        """
        expected_columns = list(dataframe.columns)
        offset = dataframe.mean()
        factor = dataframe.std()
        return OffsetScaler(expected_columns, offset, factor)

    def __init__(self, expected_columns, offset_series, factor_series):
        """
        This scaling object shifts by the offset and then scales the result
        by the factor. Typically, one would create this with the static method
        create_from_dataframe.

        scaled_data = (data-offset)/factor

        Args:
           expected_columns: list of str
              list of strings indicating the names of the columns in offset_series,
              factor_series, and the dataframe passed to scale and unscale.
           offset_series: pandas Series
              Series with columns (or labels) the same as expected_columns and
              values that represent the offset to be used when shifting the data
           factor_series: pandas Series
              Series with columns (or labels) the same as expected_columns and
              values that represent the factor to be used to scale the shifted data
        """
        self._expected_columns = expected_columns
        self._offset = offset_series
        self._factor = factor_series
        if list(offset_series.index) != expected_columns:
            raise ValueError(
                "OffsetScaler was passed an offset series with an index that"
                " does not match expected_columns. Please make sure these labels match."
            )
        if list(factor_series.index) != expected_columns:
            raise ValueError(
                "OffsetScaler was passed a factor series with an index that"
                " does not match expected_columns. Please make sure these labels match."
            )

    def _verify_columns_match(self, dataframe):
        if self._expected_columns != list(dataframe.columns):
            raise ValueError(
                "OffsetScaler was passed a dataframe that did not contain"
                " the same column labels as those used to create the scaler."
                " Please make sure the column labels match."
            )

    def scale(self, dataframe):
        """
        Return a new dataframe where the values are scaled according to the
        offset and factor

        Args:
           dataframe: pandas Dataframe
              The dataframe to be scaled

        Returns: pandas DataFrame
        """
        self._verify_columns_match(dataframe)
        df = dataframe - self._offset
        df = df.divide(self._factor)
        return df

    def unscale(self, dataframe):
        """
        Return a new dataframe where the values are unscaled according to the
        offset and factor

        Args:
           dataframe: pandas Dataframe
              The dataframe to be unscaled

        Returns: pandas DataFrame
        """
        self._verify_columns_match(dataframe)
        df = dataframe.multiply(self._factor)
        df = df + self._offset
        return df

    def expected_columns(self):
        """
        Return the expected column names for the scaler series objects
        """
        return self._expected_columns

    def offset_series(self):
        """
        Return the offset for the scaler as a pandas Series object
        """
        return self._offset

    def factor_series(self):
        """
        Return the factors for the scaler as a pandas Series object
        """
        return self._factor

    def to_dict(self):
        """
        Returns a dictionary representation of this scaler
        """
        d = dict()
        d["expected_columns"] = list(self._expected_columns)
        d["offset"] = self._offset.to_dict()
        d["factor"] = self._factor.to_dict()
        return d

    @staticmethod
    def from_dict(d):
        """
        Create an instance of this scaler from a dictionary
        (that was created with to_dict)

        Args:
           d : dict
              The dict created with to_dict

        Returns: new OffsetScaler
        """
        expected_columns = d["expected_columns"]
        offset = pd.Series(d["offset"])
        factor = pd.Series(d["factor"])
        return OffsetScaler(expected_columns, offset, factor)
