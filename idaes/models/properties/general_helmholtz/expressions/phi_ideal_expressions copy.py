#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""Predefined ideal expressions for Helmholtz EoS functions"""

__author__ = "Stephen Burroughs"

import pyomo.environ as pyo


def phi_ideal_expressions_lead(model, parameters):
    """lead expression for the ideal part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for lead part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    return {
        "phii": pyo.log(model.delta) + a[0] + a[1] * model.tau,
        "phii_d": 1.0 / model.delta,
        "phii_dd": -1.0 / model.delta**2,
        "phii_t": a[1],
        "phii_tt": 0,
        "phii_dt": 0,
    }


def phi_ideal_expressions_logtau(model, parameters):
    """logtau expression for the ideal part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for logtau part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    return {
        "phii": a * pyo.log(model.tau),
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": a / model.tau,
        "phii_tt": -a / model.tau**2,
        "phii_dt": 0,
    }


def phi_ideal_expressions_planck_einstein1(model, parameters):
    """Type01 expression for the first Planck Einstein part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for first Planck Einstein part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    g = parameters["g"]
    rng = range(0, len(a))
    return {
        "phii": sum(a[i] * pyo.log(1 - pyo.exp(-g[i] * model.tau)) for i in rng),
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": sum(a[i] * g[i] / (pyo.exp(g[i] * model.tau) - 1) for i in rng),
        "phii_tt": -sum(
            a[i]
            * g[i] ** 2
            * pyo.exp(-g[i] * model.tau)
            / (1 - pyo.exp(-g[i] * model.tau)) ** 2
            for i in rng
        ),
        "phii_dt": 0,
    }


def phi_ideal_expressions_planck_einstein2(model, parameters):
    """Second Planck Einstein expression for the ideal part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for second Planck Einstein part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    c = parameters["c"]
    g = parameters["g"]
    return {
        "phii": a * pyo.log(c + pyo.exp(g * model.tau)),
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": a * g * pyo.exp(g * model.tau) / (c * pyo.exp(g * model.tau)),
        "phii_tt": a
        * c
        * g**2
        * pyo.exp(g * model.tau)
        / (c + pyo.exp(g * model.tau) ** 2),
        "phii_dt": 0,
    }


def phi_ideal_expressions_cp_constant(model, parameters):
    """Type01 expression for the cp constant part of ideal dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for cp constant part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    t0 = model.T_Star / model.T_ref
    return {
        "phii": a - a * t0 + a * pyo.log(t0),
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": a * (1 / t0 - 1),
        "phii_tt": -a / t0**2,
        "phii_dt": 0,
    }


def phi_ideal_expressions_power(model, parameters):
    """Type01 expression for the Power part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for Power part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    g = parameters["g"]
    rng = range(0, len(a))
    return {
        "phii": sum(a[i] * model.tau ** g[i] for i in rng),
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": sum(a[i] * g[i] * model.tau ** (g[i] - 1) for i in rng),
        "phii_tt": sum(a[i] * (g[i] - 1) * g[i] * model.tau ** (g[i] - 2) for i in rng),
        "phii_dt": 0,
    }


def phi_ideal_expressions_enth_entr_offset(model, parameters):
    """Type01 expression for the  part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for  part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    return {
        "phii": a[0] + a[1] * model.tau,
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": a[1],
        "phii_tt": 0,
        "phii_dt": 0,
    }


def phi_ideal_expressions_GERG_Sinh(model, parameters):
    """Type01 expression for the  part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for  part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    g = parameters["g"]
    Tci_over_Tr = model.Tc / model.T_star
    rng = range(0, len(a))

    return {
        "phii": sum(
            a[i] * pyo.log(abs(pyo.sinh(g[i] * Tci_over_Tr * model.tau))) for i in rng
        ),
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": sum(
            a[i] * g[i] * Tci_over_Tr / pyo.tanh(g[i] * Tci_over_Tr * model.tau)
            for i in rng
        ),
        "phii_tt": sum(
            -a[i]
            * (g[i] * Tci_over_Tr) ** 2
            / (pyo.sinh(g[i] * Tci_over_Tr * model.tau) ** 2)
            for i in rng
        ),
        "phii_dt": 0,
    }


def phi_ideal_expressions_GERG_Cosh(model, parameters):
    """Type01 expression for the  part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for  part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    g = parameters["g"]
    Tci_over_Tr = model.Tc / model.T_star
    rng = range(0, len(a))

    return {
        "phii": sum(
            -a[i] * pyo.log(abs(pyo.cosh(g[i] * Tci_over_Tr * model.tau))) for i in rng
        ),
        "phii_d": 0,
        "phii_dd": 0,
        "phii_t": sum(
            -a[i] * g[i] * Tci_over_Tr * pyo.tanh(g[i] * Tci_over_Tr * model.tau)
            for i in rng
        ),
        "phii_tt": sum(
            -a[i]
            * (g[i] * Tci_over_Tr) ** 2
            / (pyo.cosh(g[i] * Tci_over_Tr * model.tau) ** 2)
            for i in rng
        ),
        "phii_dt": 0,
    }
