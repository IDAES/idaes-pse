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
Pure component properties as used by the NIST WebBook

https://webbook.nist.gov/chemistry/

Retrieved: September 13th, 2019

All parameter indicies based on conventions used by the source
"""
from pyomo.environ import log


# -----------------------------------------------------------------------------
# Shomate Equation for heat capacities, enthalpy and entropy
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
# Antoine equation for saturation pressure
def antoine(b, T, j):
    return 10**(b._params.antoine_coeff[j, 'A'] -
                b._params.antoine_coeff[j, 'B'] /
                (T + b._params.antoine_coeff[j, 'C']))
