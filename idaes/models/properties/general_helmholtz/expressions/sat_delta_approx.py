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
"""Predefined expressions for Helmholtz EoS functions
"""

__author__ = "John Eslick"

import pyomo.environ as pyo


def sat_delta_type01(model, name, parameters):
    """Type01 expression for the approximate saturated reduced density

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for approximate saturated reduced density
    """
    c = parameters["aux"][name]["c"]
    n = parameters["aux"][name]["n"]
    t = parameters["aux"][name]["t"]
    return c + sum(n[i] * (1 - 1 / model.tau) ** t[i] for i in n)


def sat_delta_type02(model, name, parameters):
    """Type01 expression for the approximate saturated reduced density

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for approximate saturated reduced density
    """
    c = parameters["aux"][name]["c"]
    n = parameters["aux"][name]["n"]
    t = parameters["aux"][name]["t"]
    return c * pyo.exp(sum(n[i] * (1 - 1 / model.tau) ** t[i] for i in n))
