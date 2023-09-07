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

import pyomo.environ as pyo


def phi_ideal_expressions_type03(model, parameters):
    """Type03 expression for the ideal part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for ideal part of Helmholtz free energy
    """
    last_term = parameters["eos"]["last_term_ideal"]
    n0 = parameters["eos"]["n0"]
    g0 = parameters["eos"]["g0"]
    rng1 = range(4, last_term + 1)
    return {
        "phii": pyo.log(model.delta)
        + n0[1]
        + n0[2] * model.tau
        + n0[3] * pyo.log(model.tau)
        + sum(n0[i] * model.tau ** g0[i] for i in rng1),
        "phii_d": 1.0 / model.delta,
        "phii_dd": -1.0 / model.delta**2,
        "phii_t": n0[2]
        + n0[3] / model.tau
        + sum(n0[i] * g0[i] * model.tau ** (g0[i] - 1) for i in rng1),
        "phii_tt": -n0[3] / model.tau**2
        + sum(n0[i] * g0[i] * (g0[i] - 1) * model.tau ** (g0[i] - 2) for i in rng1),
        "phii_dt": 0,
    }
