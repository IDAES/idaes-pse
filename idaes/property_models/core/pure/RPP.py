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
Methods for calculating pure component properties from:

The Properties of Gases & Liquids, 4th Edition
Reid, Prausnitz and Polling, 1987, McGraw-Hill

All parameter indicies based on conventions used by the source
"""
from pyomo.environ import exp, log


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
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
# Saturation pressure
# Note that this equation in not valid beyond the critical temperature
def RPP_Psat1(b, T, j):
    x = 1 - T/b._params.temperature_crit[j]

    return (exp((1-x)**-1 * (b._params.pressure_sat_coeff[j, 'A']*x +
                             b._params.pressure_sat_coeff[j, 'B']*x**1.5 +
                             b._params.pressure_sat_coeff[j, 'C']*x**3 +
                             b._params.pressure_sat_coeff[j, 'D']*x**6)) *
            b._params.pressure_crit[j])
