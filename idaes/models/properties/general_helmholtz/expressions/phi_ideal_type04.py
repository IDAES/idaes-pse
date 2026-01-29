#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""Predefined expression for Helmholtz EoS functions
"""
# Extended from phi_ideal_type01.py by John Eslick
__author__ = "Ben Lincoln and Stephen Burroughs"

import pyomo.environ as pyo


def phi_ideal_expressions_type04(model, parameters):
    """Type01 expression for the ideal part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for ideal part of Helmholtz free energy
    """
    last_term = parameters["eos"]["last_term_ideal"]
    n0 = parameters["eos"]["n0"]
    Tc = parameters["basic"]["Tc"]
    g0 = parameters["eos"]["g0"]
    rng = range(4, last_term + 1)
    return {
        "phii": pyo.log(model.delta)
        + n0[1]
        + n0[2] * model.tau
        + (n0[3] - 1) * pyo.log(model.tau)
        + sum(n0[i] * pyo.log(1 - pyo.exp(-g0[i] * model.tau / Tc)) for i in rng),
        "phii_d": 1.0 / model.delta,
        "phii_dd": -1.0 / model.delta**2,
        "phii_t": n0[2]
        + (n0[3] - 1) / model.tau
        + sum(
            n0[i] * g0[i] / (Tc * (pyo.exp((model.tau * g0[i]) / Tc) - 1)) for i in rng
        ),
        "phii_tt": (1 - n0[3]) / model.tau**2
        - sum(
            (n0[i] * g0[i] ** 2 * pyo.exp((model.tau * g0[i]) / Tc))
            / (Tc**2 * (pyo.exp((model.tau * g0[i]) / Tc) - 1) ** 2)
            for i in rng
        ),
        "phii_dt": 0,
    }
