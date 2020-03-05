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
def cp_mol_ig_comp(b, j, T):
    # Specific heat capacity (const. P)  via the Shomate equation
    t = T/1000
    return(b.params.cp_mol_ig_comp_coeff[j, "A"] +
           b.params.cp_mol_ig_comp_coeff[j, "B"]*t +
           b.params.cp_mol_ig_comp_coeff[j, "C"]*t**2 +
           b.params.cp_mol_ig_comp_coeff[j, "D"]*t**3 +
           b.params.cp_mol_ig_comp_coeff[j, "E"]*t**-2)


def enth_mol_ig_comp(b, j, T):
    # Specific enthalpy via the Shomate equation
    t = T/1000
    tr = b.params.temperature_ref/1000
    return 1e3*(b.params.cp_mol_ig_comp_coeff[j, "A"]*(t-tr) +
                (b.params.cp_mol_ig_comp_coeff[j, "B"]/2) *
                (t**2-tr**2) +
                (b.params.cp_mol_ig_comp_coeff[j, "C"]/3) *
                (t**3-tr**3) +
                (b.params.cp_mol_ig_comp_coeff[j, "D"]/4) *
                (t**4-tr**4) -
                b.params.cp_mol_ig_comp_coeff[j, "E"]*(1/t-1/tr) +
                b.params.cp_mol_ig_comp_coeff[j, "F"] -
                b.params.cp_mol_ig_comp_coeff[j, "H"])


def entr_mol_ig_comp(b, j, T):
    # Specific entropy via the Shomate equation
    t = T/1000
    return(b.params.cp_mol_ig_comp_coeff[j, "A"]*log(t) +
           b.params.cp_mol_ig_comp_coeff[j, "B"]*t +
           (b.params.cp_mol_ig_comp_coeff[j, "C"]/2)*t**2 +
           (b.params.cp_mol_ig_comp_coeff[j, "D"]/3)*t**3 -
           (b.params.cp_mol_ig_comp_coeff[j, "E"]/2)*t**-2 +
           b.params.cp_mol_ig_comp_coeff[j, "G"])


# -----------------------------------------------------------------------------
# Antoine equation for saturation pressure
def pressure_sat_comp(b, j, T):
    return 10**(b.params.pressure_sat_comp_coeff[j, 'A'] -
                b.params.pressure_sat_comp_coeff[j, 'B'] /
                (T + b.params.pressure_sat_comp_coeff[j, 'C']))


def pressure_sat_comp_dT(b, j, T):
    return (pressure_sat_comp(b, j, T) *
            b.params.pressure_sat_comp_coeff[j, 'B'] *
            log(10)/(T + b.params.pressure_sat_comp_coeff[j, 'C'])**2)
