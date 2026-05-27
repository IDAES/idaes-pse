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
    }


def phi_ideal_expressions_planck_einstein1(model, parameters):
    """Type01 expression for the first Planck Einstein part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for first Planck Einstein part of ideal Helmholtz free energy
    """
    n = parameters["n"]
    t = parameters["t"]

    rng = range(0, len(n))
    return {
        "phii": sum(n[i] * pyo.log(1 - pyo.exp(-t[i] * model.tau)) for i in rng),
    }


def phi_ideal_expressions_planck_einstein2(model, parameters):
    """Second Planck Einstein expression for the ideal part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for second Planck Einstein part of ideal Helmholtz free energy
    """
    n = parameters["n"]
    l = parameters["l"]
    t = parameters["t"]
    d = parameters["d"]
    rng = range(0, len(n))

    return {
        "phii": sum(n[i] * pyo.log(l[i] + d[i]*pyo.exp(t[i] * model.tau))for i in rng),
    }

def phi_ideal_expressions_planck_einstein3(model, parameters):
    """Type01 expression for the third Planck Einstein part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for first Planck Einstein part of ideal Helmholtz free energy
    """
    n = parameters["n"]
    t = parameters["t"]
    Tc = model.Tc

    rng = range(0, len(n))
    return {
        "phii": sum(n[i] * pyo.log(1 - pyo.exp(-t[i] * model.tau/Tc)) for i in rng),
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
    }


def phi_ideal_expressions_power(model, parameters):
    """Type01 expression for the Power part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for Power part of ideal Helmholtz free energy
    """
    a = parameters["n"]
    t = parameters["t"]
    rng = range(0, len(a))
    return {
        "phii": sum(a[i] * model.tau ** t[i] for i in rng),
    }


def phi_ideal_expressions_enth_entr_offset(model, parameters):
    """Enthalpy/entropy offset part of dimensionless ideal Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for  part of ideal Helmholtz free energy
    """
    a = parameters["a"]
    return {
        "phii": a[0] + a[1] * model.tau,
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
    }
