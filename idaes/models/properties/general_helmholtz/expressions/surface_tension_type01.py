#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""Predefined expression for Helmholtz EoS functions
"""

__author__ = "John Eslick"

from idaes.core.util.math import smooth_max


def surface_tension_type01(model, parameters):
    """Type01 expression for the surface tension

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for the surface tension
    """
    s = parameters["transport"]["surface_tension"]["s"]
    n = parameters["transport"]["surface_tension"]["n"]
    tc = parameters["transport"]["surface_tension"]["Tc"]
    eps = 1e-4
    # Use smooth max to avoid evaluation errors on negative numbers
    # this means the surface tension will go to zero at and above the
    # critical temperature, rather than cause error over critical.
    omTr = smooth_max(1 - model.T_star / model.tau / tc, 0, eps=eps)
    return sum(s[i] * omTr ** n[i] for i in s)
