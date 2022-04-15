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
This module contains functions for constructing CoolProp expressions.
"""
from pyomo.environ import exp, units as pyunits, Var
from pyomo.core.expr.calculus.derivatives import Modes, differentiate

from idaes.core.util.exceptions import ConfigurationError


# TODO : Only have temperature derivative expression for exponential_tau form
# TODO : Add other derivative forms as/if required


def parameters_nt_sum(cobj, prop, nlist, tlist):
    """
    Create parameters for expression forms using n-t parameters

    Args:
        cobj: Component object that will contain the parameters
        prop: name of property parameters are associated with
        nlist: list of values for n-parameter
        tlist: list of values for t-parameter

    Returns:
        None
    """
    if len(nlist) != len(tlist):
        raise ConfigurationError(
            f"{cobj.name} mismatched length between n and t parameters "
            f"for CoolProp exponential form for property {prop}. Please "
            f"ensure the number of n and t parameters are equal."
        )

    # Use multiple Vars, instead of single indexed Var, to have same
    # structure as cases where each parameter value has different units

    for i, nval in enumerate(nlist):
        coeff = Var(
            doc="Multiplying parameter for CoolProp exponential form",
            units=pyunits.dimensionless,
        )
        cobj.add_component(prop + "_coeff_n" + str(i + 1), coeff)
        coeff.fix(nval)

    for i, tval in enumerate(tlist):
        coeff = Var(
            doc="Exponent parameter for CoolProp exponential form",
            units=pyunits.dimensionless,
        )
        cobj.add_component(prop + "_coeff_t" + str(i + 1), coeff)
        coeff.fix(tval)


def _nt_sum(cobj, prop, theta):
    """
    Create sum expressions in n-t forms (sum(n[i]*theta**t[i]))

    Args:
        cobj: Component object that will contain the parameters
        prop: name of property parameters are associated with
        theta: expression or variable to use for theta in expression

    Returns:
        Pyomo expression of sum term
    """
    # Build sum term
    i = 1
    s = 0
    while True:
        try:
            ni = getattr(cobj, f"{prop}_coeff_n{i}")
            ti = getattr(cobj, f"{prop}_coeff_t{i}")
            s += ni * theta**ti
            i += 1
        except AttributeError:
            break
    return s


def expression_exponential(cobj, prop, T, yc, tau=False):
    """
    Create expressions for CoolProp exponential sum forms. This function
    supports both exponential forms used by CoolProp:

    Without tau: y = yc * exp(sum(ni*theta^ti))
    With tau: y = yc * exp((Tc/T) * sum(ni*theta^ti))

    Args:
        cobj: Component object that will contain the parameters
        prop: name of property parameters are associated with
        T: temperature to use in expression
        yc: value of property at critical point
        tau: whether tau=Tc/T should be included in expression (default=False)

    Returns:
        Pyomo expression matching CoolProp exponential sum form
    """
    Tc = cobj.temperature_crit
    theta = 1 - T / Tc

    s = _nt_sum(cobj, prop, theta)

    if tau:
        return yc * exp(Tc / T * s)
    else:
        return yc * exp(s)


def dT_expression_exponential(cobj, prop, T, yc, tau=False):
    """
    Create expressions for temperature derivative of CoolProp exponential sum
    forms with tau

    Args:
        cobj: Component object that will contain the parameters
        prop: name of property parameters are associated with
        T: temperature to use in expression
        yc: value of property at critical point
        tau: whether tau=Tc/T should be included in expression (default=False)

    Returns:
        Pyomo expression for temperature derivative of CoolProp exponential sum
        form with tau
    """
    # y = yc * exp(Tc/T * sum(ni*theta^ti))
    # Need d(y)/dT

    y = expression_exponential(cobj, prop, T, yc, tau)

    return differentiate(expr=y, wrt=T, mode=Modes.reverse_symbolic)


def expression_nonexponential(cobj, prop, T, yc):
    """
    Create expressions for CoolProp non-exponential sum forms

    Args:
        cobj: Component object that will contain the parameters
        prop: name of property parameters are associated with
        T: temperature to use in expression
        yc: value of property at critical point

    Returns:
        Pyomo expression mathcing CoolProp non-exponential sum form
    """
    # y = yc * (1 + sum(ni*theta^ti))
    Tc = cobj.temperature_crit
    theta = 1 - T / Tc

    s = _nt_sum(cobj, prop, theta)

    return yc * (1 + s)


def parameters_polynomial(cobj, prop, prop_units, alist, blist):
    """
    Create parameters for expression forms using A-B parameters (rational
    polynomial forms)

    Args:
        cobj: Component object that will contain the parameters
        prop: name of property parameters are associated with
        prop_units: units of measurement for property
        Alist: list of values for A-parameter
        Blist: list of values for B-parameter

    Returns:
        None
    """
    for i, aval in enumerate(alist):
        if i == 0:
            param_units = prop_units
        else:
            param_units = prop_units / pyunits.K**i

        coeff = Var(doc="A parameter for CoolProp polynomial form", units=param_units)
        cobj.add_component(prop + "_coeff_A" + str(i), coeff)
        coeff.fix(aval)

    for i, bval in enumerate(blist):
        if i == 0:
            param_units = pyunits.dimensionless
        else:
            param_units = pyunits.K**-i

        coeff = Var(doc="B parameter for CoolProp exponential form", units=param_units)
        cobj.add_component(prop + "_coeff_B" + str(i), coeff)
        coeff.fix(bval)


def expression_polynomial(cobj, prop, T):
    """
    Create expressions for CoolProp rational polynomial forms

    Args:
        cobj: Component object that will contain the parameters
        prop: name of property parameters are associated with
        T: temperature to use in expression

    Returns:
        Pyomo expression mathcing CoolProp rational polynomial form
    """
    i = 0
    asum = 0
    try:
        while True:
            Ai = getattr(cobj, f"{prop}_coeff_A{i}")
            asum += Ai * T**i
            i += 1
    except AttributeError:
        pass

    i = 0
    bsum = 0
    try:
        while True:
            Bi = getattr(cobj, f"{prop}_coeff_B{i}")
            bsum += Bi * T**i
            i += 1
    except AttributeError:
        pass

    return asum / bsum
