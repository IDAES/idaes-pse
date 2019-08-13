##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Methods for calculating heat capacities, enthalpies and entropies of pure
components.
"""
from pyomo.environ import log


# -----------------------------------------------------------------------------
# Shomate Equation
# Parameter indices are based on convention used by NIST Webbook
def shomate_cp_ig(b, j, T):
    # Specific heat capacity (const. P)  via the Shomate equation
    return(b._params.cp_ig[j, "A"] +
           b._params.cp_ig[j, "B"]*T +
           b._params.cp_ig[j, "C"]*T**2 +
           b._params.cp_ig[j, "D"]*T**3 +
           b._params.cp_ig[j, "E"]*T**-2)


def shomate_enth_ig(b, j, T):
    # Specific enthalpy via the Shomate equation
    return(b._params.cp_ig[j, "A"]*(T-b._params.temperature_ref) +
           (b._params.cp_ig[j, "B"]/2) *
           (T**2-b._params.temperature_ref**2) +
           (b._params.cp_ig[j, "C"]/3) *
           (T**3-b._params.temperature_ref**3) +
           (b._params.cp_ig[j, "D"]/4) *
           (T**4-b._params.temperature_ref**4) -
           b._params.cp_ig[j, "E"]*(T-b._params.temperature_ref) +
           b._params.cp_ig[j, "F"] - b._params.cp_ig[j, "H"])


def shomate_entr_ig(b, j, T):
    # Specific entropy via the Shomate equation
    return(b._params.cp_ig[j, "A"]*log(T) +
           b._params.cp_ig[j, "B"]*T +
           (b._params.cp_ig[j, "C"]/2)*T**2 +
           (b._params.cp_ig[j, "D"]/3)*T**3 -
           (b._params.cp_ig[j, "E"]/2)*T**-2 +
           b._params.cp_ig[j, "G"])


# -----------------------------------------------------------------------------
# Equation from The Properties of Gases & Liquids, 4th Edition
# Reid, Prausnitz, Poling, McGraw-Hill
# Parameter names use convention from source text
def RPP_cp_ig(b, j, T):
    # Specific enthalpy
    return (b._params.cp_ig[j, "D"]*T**3 +
            b._params.cp_ig[j, "C"]*T**2 +
            b._params.cp_ig[j, "B"]*T +
            b._params.cp_ig[j, "A"])


def RPP_enth_ig(b, j, T):
    # Specific enthalpy
    return ((b._params.cp_ig[j, "D"]/4) *
            (T**4-b._params.temperature_ref**4) +
            (b._params.cp_ig[j, "C"]/3) *
            (T**3-b._params.temperature_ref**3) +
            (b._params.cp_ig[j, "B"]/2) *
            (T**2-b._params.temperature_ref**2) +
            b._params.cp_ig[j, "A"] *
            (T-b._params.temperature_ref))


def RPP_entr_ig(b, j, T):
    # Specific entropy
    return ((b._params.cp_ig[j, 'D']/3)*T**3 +
            (b._params.cp_ig[j, 'C']/2)*T**2 +
            b._params.cp_ig[j, 'B']*T +
            b._params.cp_ig[j, 'A']*log(T))


# -----------------------------------------------------------------------------
# Equation from Chemical Engineers Handbook, 7th Edition
# Perry, McGraw-Hill
# Parameter names use convention from source text
def Perry_cp_liq(b, j, T):
    # Specific enthalpy
    return (b._params.cp_liq[j, "5"]*T**4 +
            b._params.cp_liq[j, "4"]*T**3 +
            b._params.cp_liq[j, "3"]*T**2 +
            b._params.cp_liq[j, "2"]*T +
            b._params.cp_liq[j, "1"])


def Perry_enth_liq(b, j, T):
    # Specific enthalpy
    return ((b._params.cp_liq[j, "5"]/5) *
            (T**5-b._params.temperature_ref**5) +
            (b._params.cp_liq[j, "4"]/4) *
            (T**4-b._params.temperature_ref**4) +
            (b._params.cp_liq[j, "3"]/3) *
            (T**3-b._params.temperature_ref**3) +
            (b._params.cp_liq[j, "2"]/2) *
            (T**2-b._params.temperature_ref**2) +
            b._params.cp_liq[j, "1"] *
            (T-b._params.temperature_ref))


def Perry_entr_liq(b, j, T):
    # Specific entropy
    return ((b._params.cp_liq[j, '5']/4)*T**4 +
            (b._params.cp_liq[j, '4']/3)*T**3 +
            (b._params.cp_liq[j, '3']/2)*T**2 +
            b._params.cp_liq[j, '2']*T +
            b._params.cp_liq[j, '1']*log(T))
