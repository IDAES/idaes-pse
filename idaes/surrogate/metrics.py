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

msg = "Surrogate metrics have not yet been calculated. Please call compute"
"metrics first."


class TrainingMetrics(object):
    def __init__(self, surrogate, dataframe):
        self._surrogate = surrogate
        self._measured_data = dataframe

        self._evaluated_data = None

        # root mean-squared error
        self._RMSE = None

        # mean squared error
        self._MSE = None

        # sum of squared errors
        self._SSE = None

        # R-squared
        self._R2 = None

        # TODO: number of datapoints out of bounds

    @property
    def RMSE(self):
        if self._RMSE is None:
            raise ValueError(msg)
        else:
            return self._RMSE

    @property
    def MSE(self):
        if self._MSE is None:
            raise ValueError(msg)
        else:
            return self._MSE

    @property
    def SSE(self):
        if self._SSE is None:
            raise ValueError(msg)
        else:
            return self._SSE

    @property
    def R2(self):
        if self._R2 is None:
            raise ValueError(msg)
        else:
            return self._R2

    @staticmethod
    def build_metrics(surrogate, dataframe):
        metrics = TrainingMetrics(surrogate, dataframe)
        metrics.compute_metrics()

        return metrics

    def evaluate_surrogate(self):
        self._evaluated_data = self._surrogate.evaluate_surrogate(
            self._measured_data)

    def compute_metrics(self):
        if self._evaluated_data is None:
            self.evaluate_surrogate()

        y = self._measured_data[self._surrogate._output_labels]
        f = self._evaluated_data[self._surrogate._output_labels]
        y_mean = y.mean(axis=0)
        self._SST = ((y-y_mean)**2).sum(axis=0)
        self._SSE = ((y-f)**2).sum(axis=0)

        self._R2 = 1-self._SSE/self._SST
        self._MSE = ((y-f)**2).mean(axis=0)
        self._RMSE = self._MSE**0.5
