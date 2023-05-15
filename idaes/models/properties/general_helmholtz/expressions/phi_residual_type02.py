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


def phi_residual_expressions_type02(model, parameters):
    """Type02 expression for the residual part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for residual part of Helmholtz free energy
    """
    last_terms = parameters["eos"]["last_term_residual"]
    n = parameters["eos"]["n"]
    t = parameters["eos"]["t"]
    d = parameters["eos"]["d"]
    c = parameters["eos"]["c"]
    a = parameters["eos"]["a"]
    b = parameters["eos"]["b"]
    e = parameters["eos"]["e"]
    g = parameters["eos"]["g"]
    first_term = 1
    rng = []
    for last_term in last_terms:
        rng.append(range(first_term, last_term + 1))
        first_term = last_term + 1
    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng[0])
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            for i in rng[2]
        ),
        "phir_d": sum(
            n[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** t[i] for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 1)
            * model.tau ** t[i]
            * (d[i] - c[i] * model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (d[i] / model.delta - 2 * a[i] * (model.delta - e[i]))
            for i in rng[2]
        ),
        "phir_dd": sum(
            n[i] * d[i] * (d[i] - 1) * model.delta ** (d[i] - 2) * model.tau ** t[i]
            for i in rng[0]
        )
        + sum(
            n[i]
            * pyo.exp(-model.delta ** c[i])
            * model.delta ** (d[i] - 2)
            * model.tau ** t[i]
            * (
                (d[i] - c[i] * model.delta ** c[i])
                * (d[i] - 1 - c[i] * model.delta ** c[i])
                - c[i] ** 2 * model.delta ** c[i]
            )
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (
                -2.0 * a[i] * model.delta ** d[i]
                + 4 * a[i] ** 2 * model.delta ** d[i] * (model.delta - e[i]) ** 2
                - 4 * d[i] * a[i] * model.delta ** (d[i] - 1) * (model.delta - e[i])
                + d[i] * (d[i] - 1) * model.delta ** (d[i] - 2)
            )
            for i in rng[2]
        ),
        "phir_t": sum(
            n[i] * t[i] * model.delta ** d[i] * model.tau ** (t[i] - 1) for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** d[i]
            * model.tau ** (t[i] - 1)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (t[i] / model.tau - 2 * b[i] * (model.tau - g[i]))
            for i in rng[2]
        ),
        "phir_tt": sum(
            n[i] * t[i] * (t[i] - 1) * model.delta ** d[i] * model.tau ** (t[i] - 2)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * (t[i] - 1)
            * model.delta ** d[i]
            * model.tau ** (t[i] - 2)
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (
                (t[i] / model.tau - 2 * b[i] * (model.tau - g[i])) ** 2
                - t[i] / model.tau**2
                - 2 * b[i]
            )
            for i in rng[2]
        ),
        "phir_dt": sum(
            n[i] * t[i] * d[i] * model.delta ** (d[i] - 1) * model.tau ** (t[i] - 1)
            for i in rng[0]
        )
        + sum(
            n[i]
            * t[i]
            * model.delta ** (d[i] - 1)
            * model.tau ** (t[i] - 1)
            * (d[i] - c[i] * model.delta ** c[i])
            * pyo.exp(-model.delta ** c[i])
            for i in rng[1]
        )
        + sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            * (d[i] / model.delta - 2.0 * a[i] * (model.delta - e[i]))
            * (t[i] / model.tau - 2.0 * b[i] * (model.tau - g[i]))
            for i in rng[2]
        ),
    }
