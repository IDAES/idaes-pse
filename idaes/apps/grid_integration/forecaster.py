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


class PlaceHolderForecaster:
    """
    This a placeholder for a real price forecaster. This placeholder can
    manipulate the forecasts in a dataframe and feed to the bidder.
    """

    def __init__(self, price_forecasts_df):
        self.price_forecasts_df = price_forecasts_df
        self.price_forecasts_df.set_index(["Date", "Hour"], inplace=True)
        self._rename_columns()

    def _rename_columns(self):
        col_mapping = {"Scenario {}".format(i): i - 1 for i in range(1, 11)}
        self.price_forecasts_df.rename(columns=col_mapping, inplace=True)

        return

    def forecast(self, date, **kwargs):
        date = str(date)
        return self.price_forecasts_df.loc[date].to_dict("list")
