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
    return sum(s[i] * (1 - model.T_star / model.tau / tc) ** n[i] for i in s)
