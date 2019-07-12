##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Convenience plotting functions for time-dependent variables.
"""

__author__ = "John Eslick"

import pyomo.environ as pyo
import matplotlib.pyplot as plt

def plot_time_dependent(time, y, ylabel, xlabel="time (s)", title="",
                        legend=None):
    """
    Plot time dependent varaibles with pyplot.

    Args:
        time (ContinuousSet): Time index set
        y (Var, Expression, Reference, or list of): Quantity to plot
            indexed only by time. If you want to plot something that is not
            indexed as required, you can create a Pyomo Reference with the
            correct indexing.  You can plot multiple quntities by providing a
            list or tuple of things to plot.
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        title (str): Plot title
        legend (list of str): Legend string for each y

    Returns:
        None
    """
    if isinstance(y, (tuple, list)):
        n = len(y)
        yn = []
        for z in y:
            yn.append([pyo.value(z[t]) for t in time])
        y = yn
    else:
        n = 1
        y = ([pyo.value(y[t]) for t in time],)

    for i in range(n):
        plt.plot(time, y[i])
    if legend is not None:
        plt.legend(legend)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.show()
