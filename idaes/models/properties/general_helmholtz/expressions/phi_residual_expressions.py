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
"""Predefined residual expressions for Helmholtz EoS functions"""

__author__ = "Stephen Burroughs"

import pyomo.environ as pyo


def phi_residual_expressions_power(model, parameters):
    """Power expression for the residual part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for power part of residual Helmholtz free energy
    """

    n = parameters["n"]
    t = parameters["t"]
    d = parameters["d"]
    rng = range(0, len(n))

    return {
        "phir": sum(n[i] * model.delta ** d[i] * model.tau ** t[i] for i in rng),
    }


def phi_residual_expressions_gaussian(model, parameters):
    """Gaussian expression for the residual part of dimensionless Helmholtz free energy

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for Gaussian part of residual Helmholtz free energy
    """
    n = parameters["n"]
    t = parameters["t"]
    d = parameters["d"]
    a = parameters["a"]
    b = parameters["b"]
    e = parameters["e"]
    g = parameters["g"]
    rng = range(0, len(n))
    return {
        "phir": sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]) ** 2
            )
            for i in rng
        ),
    }


def phi_residual_expressions_gaussian_GERG2008(model, parameters):
    """Gaussian expression for the residual part of dimensionless Helmholtz free energy - GERG
    Reference: Kunz, O., & Wagner, W. (2012). The GERG-2008 wide-range equation of state for natural gases and other mixtures: An expansion of GERG-2004. Journal of chemical & engineering data, 57(11), 3032-3091.

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for Gaussian part of residual Helmholtz free energy
    """
    n = parameters["n"]
    t = parameters["t"]
    d = parameters["d"]
    a = parameters["a"]
    b = parameters["b"]
    e = parameters["e"]
    g = parameters["g"]
    rng = range(0, len(n))
    return {
        "phir": sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-a[i] * (model.delta - e[i]) ** 2 - b[i] * (model.tau - g[i]))
            for i in rng
        ),
    }


def phi_residual_expressions_gaob(model, parameters):
    """Expression for associating term of the residual part of dimensionless Helmholtz free energy
    Reference: Gao, K., Wu, J., Bell, I. H., Harvey, A. H., & Lemmon, E. W. (2023). A reference equation of state with an associating term for the thermodynamic properties of ammonia. Journal of Physical and Chemical Reference Data, 52(1).
    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expression for associating term of residual Helmholtz free energy
    """
    n = parameters["n"]
    t = parameters["t"]
    d = parameters["d"]
    a = parameters["a"]
    b = parameters["b"]
    bi = parameters["bi"]
    e = parameters["e"]
    g = parameters["g"]
    rng = range(0, len(n))
    return {
        "phir": sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(
                -a[i] * (model.delta - e[i]) ** 2
                + 1 / (b[i] * (model.tau - g[i]) ** 2 + bi[i])
            )
            for i in rng
        ),
    }


def phi_residual_expressions_exponential_delta_tau(model, parameters):
    """Expression for exponentials in the delta and tau family for the residual part of dimensionless Helmholtz free energy
    Reference: Lemmon, E. W., & Jacobsen, R. T. (2005). A new functional form and new fitting techniques for equations of state with application to pentafluoroethane (HFC-125). Journal of physical and chemical reference data, 34(1), 69-108.


    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for exponential part of residual Helmholtz free energy with regard to delta and tau
    """
    n = parameters["n"]
    t = parameters["t"]
    d = parameters["d"]
    l = parameters["l"]
    m = parameters["m"]
    rng = range(0, len(n))

    return {
        "phir": sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-model.delta ** l[i])
            * pyo.exp(-model.tau ** m[i])
            for i in rng
        ),
    }


def phi_residual_expressions_double_exponential(model, parameters):
    """Expression for double exponentials in the delta and tau family for the residual part of dimensionless Helmholtz free energy
    Reference: De Reuck, K. M., & Craven, R. J. B. (1993). Methanol. International Thermodynamic Tables of the Fluid State, vol. 12. IUPAC, Blackwell Scientific Publications, London.


    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for double exponential part of residual Helmholtz free energy
    """
    n = parameters["n"]
    t = parameters["t"]
    d = parameters["d"]
    gd = parameters["gd"]
    gt = parameters["gt"]
    ld = parameters["ld"]
    lt = parameters["lt"]

    rng = range(0, len(n))

    return {
        "phir": sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-gd[i] * model.delta ** ld[i - gt[i] * model.tau ** lt[i]])
            for i in rng
        ),
    }


def phi_residual_expressions_exponential_reduced_density(model, parameters):
    """Expression for reduced density exponentials in the residual part of dimensionless Helmholtz free energy, often just referred to as the exponential term

    Args:
        model (Block): Pyomo model
        parameters (dict): Main parameters dictionary

    Returns:
        dict: Expressions for exponential part of residual Helmholtz free energy with regard to reduced density
    """
    n = parameters["n"]
    t = parameters["t"]
    d = parameters["d"]
    l = parameters["l"]
    g = parameters["g"]
    rng = range(0, len(n))
    return {
        "phir": sum(
            n[i]
            * model.delta ** d[i]
            * model.tau ** t[i]
            * pyo.exp(-g[i] * model.delta ** l[i])
            for i in rng
        ),
    }
