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
    """
    y = dataframe[surrogate.output_labels()]
    f = surrogate.evaluate_surrogate(dataframe)
    assert f.columns == surrogate.output_labels()

    y_mean = y.mean(axis=0)
    SST = ((y-y_mean)**2).sum(axis=0)
    SSE = ((y-f)**2).sum(axis=0)

    R2 = 1-SSE/SST
    MSE = ((y-f)**2).mean(axis=0)
    RMSE = MSE**0.5

    return TrainingMetrics(RMSE=RMSE, MSE=MSE, SSE=SSE, R2=R2)


#TODO: Maybe this should just be a dictionary (or munch)
#TODO: think of other metrics utility functions we may want
class TrainingMetrics(object):
    def __init__(self, RMSE, MSE, SSE, R2):
        # root mean-squared error
        self._RMSE = RMSE

        # mean squared error
        self._MSE = MSE

        # sum of squared errors
        self._SSE = SSE

        # R-squared
        self._R2 = R2

        # TODO: number of datapoints out of bounds

    @property
    def RMSE(self):
        return self._RMSE

    @property
    def MSE(self):
        return self._MSE

    @property
    def SSE(self):
        return self._SSE

    @property
    def R2(self):
        return self._R2

    def __str__(self):
        ret = '\n'
        for l in self.RMSE.index:
            ret += 'Training metrics for output: {}\n'.format(l)
            ret += '---------------------------------------------------------\n'
            ret += 'RMSE:          {}\n'.format(float(self.RMSE[l]))
            ret += 'MSE:           {}\n'.format(float(self.MSE[l]))
            ret += 'R2:            {}\n'.format(float(self.R2[l]))
            ret += 'SSE:           {}\n'.format(float(self.SSE[l]))

        return ret

