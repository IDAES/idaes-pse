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


def parameters_exponential(cobj, prop, nlist, tlist):
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


def _exponential_sum(cobj, prop, theta):
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


def expression_exponential(cobj, prop, T, yc):
    # y = yc * exp(Tc/T * sum(ni*theta^ti))
    Tc = cobj.temperature_crit
    theta = 1 - T/Tc

    s = _exponential_sum(cobj, prop, theta)

    return yc*exp(s)


def expression_exponential_tau(cobj, prop, T, yc):
    # y = yc * exp(Tc/T * sum(ni*theta^ti))
    Tc = cobj.temperature_crit
    theta = 1 - T/Tc

    s = _exponential_sum(cobj, prop, theta)

    return yc*exp(Tc/T*s)


def parameters_polynomial(cobj, prop, prop_units, alist, blist):
    # TODO : Get number of terms to build
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
