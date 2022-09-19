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

import matplotlib.pyplot as plt

import pyomo.environ as pyo
from pyomo.core.base.indexed_component_slice import IndexedComponent_slice


def dynamic_value_list(*args, units=None):
    """
    Create a list of values from a Pyomo component. Multiple models representing
    different time periods can also be combined by providing multiple arguments
    to assemble into one time series.

    Args:
        Positional arguments (): Multiple Pyomo components indexed by time, or
            time sets. Slices are ok.

    Returns:
        (list) with the time indexed Pyomo component values concatenated for
            plotting
    """
    l = []
    for v in args:
        if isinstance(v, pyo.Set):
            l += [t for t in v]
        elif isinstance(v, IndexedComponent_slice):
            if units is None:
                l += [pyo.value(sv) for sv in v]
            else:
                l += [pyo.value(pyo.units.convert(sv, units)) for sv in v]
        else:
            if units is None:
                l += [pyo.value(v[t]) for t in v]
            else:
                l += [pyo.value(pyo.units.convert(v[t], units)) for t in v]
    return l


def plot_grid(
    x, y, xlabel, ylabel, cols=1, rows=1, same_x=False, ylabel_title=True, to_file=None
):
    """Make a grid of plots.

    Args:
        x: values for x-axis. If all plots have the same x-axis the same_x arg
            should be True, and only one x-value list is need, otherwise x should
            be a list of x value lists.
        y: list of lists of y values for plots
        xlabel: if same_x=True, a x-axis label, if same_x=False, list of xaxis
            labels
        ylable: list of y-axis labels
        cols: number of columns in plot grid
        rows: number of rows in plot grid
        same_x: If True use the same x-axis for all plots, otherwise expect a
            list of x-axis data for the x-axis of each plot
        ylabel_title: If True, use the ylabel for the plot title, otherwise label
            the y-axis

    Returns:
        matplotlib.pyplot
    """
    assert cols * rows >= len(y)
    for i, yi in enumerate(y):
        yilabel = ylabel[i]
        if same_x:
            xi = x
            xilabel = xlabel
        else:
            xi = x[i]
            xilabel = xlabel[i]

        plt.subplot(rows, cols, i + 1)
        plt.xlabel(xilabel)
        plt.plot(xi, yi)
        if ylabel_title:
            plt.title(yilabel)
        else:
            plt.ylabel(yilabel)

    plt.tight_layout()
    if to_file:
        plt.savefig(to_file)
    else:
        plt.show()
    return plt


def plot_grid_dynamic(
    x, y, xlabel, ylabel, yunits=None, cols=1, rows=1, ylabel_title=True, to_file=None
):
    """Make a grid of dynamic plots, where the x-axices are time.

    Args:
        x: time set or list of time values
        y: list of lists of y values for plots.  The y value lists can be
            time-indexed pyo components or time slices
        xlabel: time axis label
        ylable: list of y-axis labels
        yunits: list of y-axis units, if units are provided, the y-values will
            be converted to the requested units for plotting
        cols: number of columns in plot grid
        rows: number of rows in plot grid
        ylabel_title: If True, use the ylabel for the plot title, otherwise label
            the y-axis

    Returns:
        matplotlib.pyplot
    """
    y_list = []
    if isinstance(x, pyo.Component):
        time = dynamic_value_list(x)
    else:
        time = x
    for i, yi in enumerate(y):
        if isinstance(yi, (pyo.Var, pyo.Expression, IndexedComponent_slice)):
            if yunits is not None:
                y_list.append(dynamic_value_list(yi, units=yunits[i]))
            else:
                y_list.append(dynamic_value_list(yi))
        else:
            y_list.append(yi)
    return plot_grid(
        x=time,
        y=y_list,
        xlabel=xlabel,
        ylabel=ylabel,
        cols=cols,
        rows=rows,
        same_x=True,
        ylabel_title=ylabel_title,
        to_file=to_file,
    )


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
    y2 = [None] * len(y)
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
