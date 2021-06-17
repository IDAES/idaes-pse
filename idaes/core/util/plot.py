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
"""
Convenience plotting functions for time-dependent variables.
"""

__author__ = "John Eslick"

import pyomo.environ as pyo
import matplotlib.pyplot as plt

def stitch_dynamic(*args):
    """
    Combine time-indexed Pyomo component values from different models into one
    combined time set. This allows you to use multiple models to simulate
    sections of the time domain, and plot them all together.

    Args:
        Positional arguments (): Multiple Pyomo components indexed by time, or
            time sets

    Returns:
        (list) with the time indexed Pyomo compoent values concatonated for
            plotting
    """
    l = []
    for v in args:
        if isinstance(v, pyo.Set):
            l += [t for t in v]
        else:
            l += [pyo.value(v[t]) for t in v]
    return l

def plot_dynamic(time, y, ylabel, xlabel="time (s)", title=None, legend=None):
    """
    Plot time dependent variables with pyplot.

    Args:
        time (ContinuousSet or list-like): Time index set
        y (list-like of list-likes of Var, Expression, Reference, or float):
            List of quantities to plot (multiple quantities can be plotted). Each
            quantity in the list should be indexed only by time. If you want to
            plot something that is not indexed only by time, you can create a
            Pyomo Reference with the correct indexing.
        ylabel (str): Y-axis label, required
        xlabel (str): X-axis label, default = 'time (s)'
        title (str or None): Plot title, default = None
        legend (list-like of str or None): Legend string for each y,
            default = None

    Returns:
        None
    """
    y2 = [None]*len(y)
    for i, z in enumerate(y):
        if isinstance(z, (list, tuple)):
            y2[i] = y[i]
        else:
            y2[i] = [pyo.value(z[t]) for t in time]
    for q in y2:
        plt.plot(time, q)
    if legend is not None:
        plt.legend(legend)
    if title is not None:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.show()
