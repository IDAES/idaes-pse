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
def cp_mol_ig_comp(b, cobj, T):
    # Specific heat capacity
    return (cobj.cp_mol_ig_comp_coeff["D"]*T**3 +
            cobj.cp_mol_ig_comp_coeff["C"]*T**2 +
            cobj.cp_mol_ig_comp_coeff["B"]*T +
            cobj.cp_mol_ig_comp_coeff["A"])


def enth_mol_ig_comp(b, cobj, T):
    # Specific enthalpy
    return ((cobj.cp_mol_ig_comp_coeff["D"]/4) *
            (T**4-b.params.temperature_ref**4) +
            (cobj.cp_mol_ig_comp_coeff["C"]/3) *
            (T**3-b.params.temperature_ref**3) +
            (cobj.cp_mol_ig_comp_coeff["B"]/2) *
            (T**2-b.params.temperature_ref**2) +
            cobj.cp_mol_ig_comp_coeff["A"] *
            (T-b.params.temperature_ref) +
            cobj.enth_mol_form_phase_comp_ref['Vap'])


def entr_mol_ig_comp(b, cobj, T):
    # Specific entropy
    return ((cobj.cp_mol_ig_comp_coeff['D']/3)*T**3 +
            (cobj.cp_mol_ig_comp_coeff['C']/2)*T**2 +
            cobj.cp_mol_ig_comp_coeff['B']*T +
            cobj.cp_mol_ig_comp_coeff['A']*log(T) +
            cobj.entr_mol_form_phase_comp_ref['Vap'])


# -----------------------------------------------------------------------------
# Saturation pressure
# Note that this equation in not valid beyond the critical temperature
def pressure_sat_comp(b, cobj, T):
    x = 1 - T/cobj.temperature_crit_comp

    return (exp((1-x)**-1 * (cobj.pressure_sat_comp_coeff['A']*x +
                             cobj.pressure_sat_comp_coeff['B']*x**1.5 +
                             cobj.pressure_sat_comp_coeff['C']*x**3 +
                             cobj.pressure_sat_comp_coeff['D']*x**6)) *
            cobj.pressure_crit_comp)


def pressure_sat_comp_dT(b, cobj, T):
    x = 1 - T/cobj.temperature_crit_comp

    return (-pressure_sat_comp(b, cobj, T) *
            ((cobj.pressure_sat_comp_coeff['A'] +
              1.5*cobj.pressure_sat_comp_coeff['B']*x**0.5 +
              3*cobj.pressure_sat_comp_coeff['C']*x**2 +
              6*cobj.pressure_sat_comp_coeff['D']*x**5)/T +
             (cobj.temperature_crit_comp/T**2) *
             (cobj.pressure_sat_comp_coeff['A']*x +
              cobj.pressure_sat_comp_coeff['B']*x**1.5 +
              cobj.pressure_sat_comp_coeff['C']*x**3 +
              cobj.pressure_sat_comp_coeff['D']*x**6)))
