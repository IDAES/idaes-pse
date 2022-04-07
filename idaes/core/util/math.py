# -*- coding: UTF-8 -*-
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
This module contains utility functions for mathematical operators of use in
equation oriented models.
"""

from pyomo.environ import Param, log, sqrt

__author__ = "Andrew Lee"


def smooth_abs(a, eps=1e-4):
    """General function for creating an expression for a smooth minimum or
    maximum.

    .. math:: |a| = sqrt(a^2 + eps^2)

    Args:
        a : term to get absolute value from (Pyomo component, float or int)
        eps : smoothing parameter (Param, float or int) (default=1e-4)

    Returns:
        An expression for the smoothed absolute value operation.
    """
    # Check type of eps
    if not (isinstance(eps, (float, int, Param))):
        raise TypeError(
            "smooth_abs eps argument must be a float, int or " "Pyomo Param"
        )

    # Create expression
    try:
        expr = (a**2 + eps**2) ** 0.5
    except TypeError:
        raise TypeError(
            "Unsupported argument type for smooth_abs. Must be "
            "a Pyomo Var, Param or Expression, or a float or int."
        )

    return expr


def smooth_minmax(a, b, eps=1e-4, sense="max"):
    """General function for creating an expression for a smooth minimum or
    maximum. Uses the smooth_abs operator.

    .. math:: minmax(a, b) = 0.5*(a+b +- |a-b|)

    Args:
        a : first term in mix or max function (Pyomo component, float or int)
        b : second term in min or max function (Pyomo component, float or int)
        eps : smoothing parameter (Param, float or int) (default=1e-4)
        sense : 'mim' or 'max' (default = 'max')

    Returns:
        An expression for the smoothed minimum or maximum operation.
    """
    # Check type of eps
    if not (isinstance(eps, (float, int, Param))):
        raise TypeError(
            "Smooth {} eps argument must be a float, int or "
            "Pyomo Param".format(sense)
        )

    # Set sense of expression
    if sense == "max":
        mm = 1
    elif sense == "min":
        mm = -1
    else:
        raise ValueError(
            "Unrecognised sense argument to smooth_minmax. " "Must be 'min' or 'max'."
        )

    # Create expression
    try:
        expr = 0.5 * (a + b + mm * smooth_abs(a - b, eps))
    except TypeError:
        raise TypeError(
            "Unsupported argument type for smooth_{}. Must be "
            "a Pyomo Var, Param or Expression, or a float or int.".format(sense)
        )

    return expr


def smooth_max(a, b, eps=1e-4):
    """Smooth maximum operator, using smooth_abs operator.

    .. math:: max(a, b) = 0.5*(a+b + |a-b|)

    Args:
        a : first term in max function
        b : second term in max function
        eps : smoothing parameter (Param or float, default = 1e-4)

    Returns:
        An expression for the smoothed maximum operation.
    """
    expr = smooth_minmax(a, b, eps, sense="max")
    return expr


def smooth_min(a, b, eps=1e-4):
    """Smooth minimum operator, using smooth_abs operator.

    .. math:: max(a, b) = 0.5*(a+b - |a-b|)

    Args:
        a : first term in min function
        b : second term in min function
        eps : smoothing parameter (Param or float, default = 1e-4)

    Returns:
        An expression for the smoothed minimum operation.
    """
    expr = smooth_minmax(a, b, eps, sense="min")
    return expr


def smooth_bound(val, lb, ub, eps=1e-4, eps_lb=None, eps_ub=None):
    """Returns a smooth expression that returns approximatly the value of val
    between lb and ub, the value of ub if val > ub and the value of lb if
    val < lb. The expression returned is

    ..math:: smooth_min(smooth_max(val, lb, eps_lb), ub, eps_ub)

    **Example Usage** This function is useful when calculating the value of a
    quantitiy with physical limitations that may not be otherwise captured. A
    good example of this is using a controller to manipulate a valve to maintain,
    for example, a set flow rate. Under some process conditions the controller
    may indicate the valve should be more than fully opened or closed. In this
    example, the smooth_bound expression can be applied to the unbounded
    controller output to get a bounded controller output, ensuring the valve
    opening stays in the feasible region, while maintaining smoothness.

    Args:
        val: an expression for a bounded version of this value is returned
        lb: lower bound
        ub: lower bound
        eps: smoothing parameter
        eps_lb: (optional) lower bound smoothing parameter if not provided use eps
        eps_ub: (optional) upper bound smoothing parameter if not provided use eps

    Returns:
        expression: smooth bounded version of val.
    """
    if eps_lb is None:
        eps_lb = eps
    if eps_ub is None:
        eps_ub = eps
    return smooth_min(smooth_max(val, lb, eps_lb), ub, eps_ub)


def safe_sqrt(a, eps=1e-4):
    """Returns the square root of max(a, 0) using the smooth_max expression.
    This can be used to avoid transient evaluation errors when changing a model
    from one state to another.  This can be used when a at the solution is not
    expected to be near 0.

    Args:
        a: Pyomo expression
        eps: epsilon parameter for smooth max

    Returns:
        approximatly sqrt(max(a, 0))
    """
    return sqrt(smooth_max(a, 0, eps))


def safe_log(a, eps=1e-4):
    """Returns the log of max(a, eps) using the smooth_max expression.
    This can be used to avoid transient evaluation errors when changing a model
    from one state to another.  This can be used when at the solution, a >> eps.

    Args:
        a: Pyomo expression
        eps: epsilon parameter for smooth max

    Returns:
        approximatly log(max(a, eps))
    """
    return log(smooth_max(a, eps, eps=eps))
