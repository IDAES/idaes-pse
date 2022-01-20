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
This module contains methods for constructing CoolProp expressions.
"""
from pyomo.environ import exp, units as pyunits, Var

from idaes.core.util.exceptions import ConfigurationError


# TODO : Only have temperature derivative expression for exponential_tau form
# TODO : Add other derivative forms as/if required

def parameters_nt_sum(cobj, prop, nlist, tlist):
    """
    Method for creating parameters for expression forms using n-t parameters

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        nlist - list of values for n-parameter
        tlist - list of values for t-parameter

    Returns:
        None
    """
    if len(nlist) != len(tlist):
        raise ConfigurationError(
            f"{cobj.name} mismatched length between n and t parameters "
            f"for CoolProp exponential form for property {prop}. Please "
            f"ensure the number of n and t parameters are equal.")

    for i in range(0, len(nlist)):
        nval = nlist[i]
        cobj.add_component(
            prop+"_coeff_n"+str(i+1),
            Var(doc="Multiplying parameter for CoolProp exponential form",
                units=pyunits.dimensionless))
        getattr(cobj, prop+"_coeff_n"+str(i+1)).fix(nval)

    for i in range(0, len(tlist)):
        tval = tlist[i]
        cobj.add_component(
            prop+"_coeff_t"+str(i+1),
            Var(doc="Exponent parameter for CoolProp exponential form",
                units=pyunits.dimensionless))
        getattr(cobj, prop+"_coeff_t"+str(i+1)).fix(tval)


def _nt_sum(cobj, prop, theta):
    """
    Method for creating sum expressions in n-t forms (sum(n[i]*theta**t[i]))

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        theta - expression or variable to use for theta in expression

    Returns:
        Pyomo expression of sum term
    """
    # Build sum term
    i = 1
    s = 0
    while True:
        try:
            ni = getattr(cobj, prop+"_coeff_n"+str(i))
            ti = getattr(cobj, prop+"_coeff_t"+str(i))
            s += ni*theta**ti
            i += 1
        except AttributeError:
            break
    return s


def _nt_sum_dT(cobj, prop, T, Tc):
    """
    Method for creating expression of temperature derivative of sum term in n-t
    forms (d/dT(sum(n[i]*theta**t[i])))

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        T - temperature used in expression
        Tc - critical temperature of component

    Returns:
        Pyomo expression of derivative of sum term
    """
    Tr = T/Tc

    # Build sum term
    i = 1
    sdT = 0
    while True:
        try:
            ni = getattr(cobj, prop+"_coeff_n"+str(i))
            ti = getattr(cobj, prop+"_coeff_t"+str(i))
            sdT += -ni*(1-Tr)**(ti-1)*((ti-1)*Tr+1)
            i += 1
        except AttributeError:
            break
    return sdT/Tr**2/Tc


def expression_exponential(cobj, prop, T, yc):
    """
    Method for creating expressions for CoolProp exponential sum forms

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        T - temperature to use in expression
        yc - value of property at critical point

    Returns:
        Pyomo expression matching CoolProp exponential sum form
    """
    # y = yc * exp(Tc/T * sum(ni*theta^ti))
    Tc = cobj.temperature_crit
    theta = 1 - T/Tc

    s = _nt_sum(cobj, prop, theta)

    return yc*exp(s)


def expression_exponential_tau(cobj, prop, T, yc):
    """
    Method for creating expressions for CoolProp exponential sum forms with tau

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        T - temperature to use in expression
        yc - value of property at critical point

    Returns:
        Pyomo expression matching CoolProp exponential sum form with tau
    """
    # y = yc * exp(Tc/T * sum(ni*theta^ti))
    Tc = cobj.temperature_crit
    theta = 1 - T/Tc

    s = _nt_sum(cobj, prop, theta)

    return yc*exp(Tc/T*s)


def dT_expression_exponential_tau(cobj, prop, T, yc):
    """
    Method for creating expressions for temperature derivative of CoolProp
    exponential sum forms with tau

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        T - temperature to use in expression
        yc - value of property at critical point

    Returns:
        Pyomo expression for temperature derivative of CoolProp exponential sum
        form with tau
    """
    # y = yc * exp(Tc/T * sum(ni*theta^ti))
    # Need d(y)/dT

    y = expression_exponential_tau(cobj, prop, T, yc)

    Tc = cobj.temperature_crit

    sdT = _nt_sum_dT(cobj, prop, T, Tc)

    return y*sdT


def expression_nonexponential(cobj, prop, T, yc):
    """
    Method for creating expressions for CoolProp non-exponential sum forms

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        T - temperature to use in expression
        yc - value of property at critical point

    Returns:
        Pyomo expression mathcing CoolProp non-exponential sum form
    """
    # y = yc * (1 + sum(ni*theta^ti))
    Tc = cobj.temperature_crit
    theta = 1 - T/Tc

    s = _nt_sum(cobj, prop, theta)

    return yc*(1+s)


def parameters_polynomial(cobj, prop, prop_units, alist, blist):
    """
    Method for creating parameters for expression forms using A-B parameters
    (rational polynomial forms)

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        prop_units - units of measurement for property
        Alist - list of values for A-parameter
        Blist - list of values for B-parameter

    Returns:
        None
    """
    for i in range(0, len(alist)):
        aval = alist[i]
        if i == 0:
            param_units = prop_units
        else:
            param_units = prop_units/pyunits.K**i
        cobj.add_component(
            prop+"_coeff_A"+str(i),
            Var(doc="A parameter for CoolProp polynomial form",
                units=param_units))
        getattr(cobj, prop+"_coeff_A"+str(i)).fix(aval)

    for i in range(0, len(blist)):
        bval = blist[i]
        if i == 0:
            param_units = pyunits.dimensionless
        else:
            param_units = pyunits.K**-i
        cobj.add_component(
            prop+"_coeff_B"+str(i),
            Var(doc="B parameter for CoolProp exponential form",
                units=param_units))
        getattr(cobj, prop+"_coeff_B"+str(i)).fix(bval)


def expression_polynomial(cobj, prop, T):
    """
    Method for creating expressions for CoolProp rational polynomial forms

    Args:
        cobj - component object that will contain the parameters
        prop - name of property parameters are associated with
        T - temperature to use in expression

    Returns:
        Pyomo expression mathcing CoolProp rational polynomial form
    """
    i = 0
    asum = 0
    while True:
        try:
            Ai = getattr(cobj, prop+"_coeff_A"+str(i))
            asum += Ai*T**i
            i += 1
        except AttributeError:
            break

    i = 0
    bsum = 0
    while True:
        try:
            Bi = getattr(cobj, prop+"_coeff_B"+str(i))
            bsum += Bi*T**i
            i += 1
        except AttributeError:
            break

    return asum/bsum
