##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
from pyomo.environ import log, Var


# -----------------------------------------------------------------------------
# Shomate Equation for heat capacities, enthalpy and entropy
class cp_mol_ig_comp():
    def build_parameters(cobj):
        cobj.cp_mol_ig_comp_coeff = Var(
                ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'],
                initialize=cobj.config.parameter_data["cp_mol_ig_comp_coeff"],
                doc="Shomate parameters for ideal gas molar heat capacity")

    def return_expression(b, cobj, T):
        # Specific heat capacity (const. P)  via the Shomate equation
        t = T/1000
        return(cobj.cp_mol_ig_comp_coeff["A"] +
               cobj.cp_mol_ig_comp_coeff["B"]*t +
               cobj.cp_mol_ig_comp_coeff["C"]*t**2 +
               cobj.cp_mol_ig_comp_coeff["D"]*t**3 +
               cobj.cp_mol_ig_comp_coeff["E"]*t**-2)


class enth_mol_ig_comp():
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
            cp_mol_ig_comp.build_parameters(cobj)

    def return_expression(b, cobj, T):
        # Specific enthalpy via the Shomate equation
        t = T/1000
        tr = b.params.temperature_ref/1000
        return 1e3*(cobj.cp_mol_ig_comp_coeff["A"]*(t-tr) +
                    (cobj.cp_mol_ig_comp_coeff["B"]/2) *
                    (t**2-tr**2) +
                    (cobj.cp_mol_ig_comp_coeff["C"]/3) *
                    (t**3-tr**3) +
                    (cobj.cp_mol_ig_comp_coeff["D"]/4) *
                    (t**4-tr**4) -
                    cobj.cp_mol_ig_comp_coeff["E"]*(1/t-1/tr) +
                    cobj.cp_mol_ig_comp_coeff["F"] -
                    cobj.cp_mol_ig_comp_coeff["H"])


class entr_mol_ig_comp():
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
            cp_mol_ig_comp.build_parameters(cobj)

    def return_expression(b, cobj, T):
        # Specific entropy via the Shomate equation
        t = T/1000
        return(cobj.cp_mol_ig_comp_coeff["A"]*log(t) +
               cobj.cp_mol_ig_comp_coeff["B"]*t +
               (cobj.cp_mol_ig_comp_coeff["C"]/2)*t**2 +
               (cobj.cp_mol_ig_comp_coeff["D"]/3)*t**3 -
               (cobj.cp_mol_ig_comp_coeff["E"]/2)*t**-2 +
               cobj.cp_mol_ig_comp_coeff["G"])


# -----------------------------------------------------------------------------
# Antoine equation for saturation pressure
class pressure_sat_comp():
    def build_parameters(cobj):
        cobj.pressure_sat_comp_coeff = Var(
                ['A', 'B', 'C'],
                initialize=cobj.config.parameter_data[
                    "pressure_sat_comp_coeff"],
                doc="Antoine coefficients for calculating Psat")

    def return_expression(b, cobj, T, dT=False):
        if dT:
            return pressure_sat_comp.dT_expression(b, cobj, T)

        return 10**(cobj.pressure_sat_comp_coeff['A'] -
                    cobj.pressure_sat_comp_coeff['B'] /
                    (T + cobj.pressure_sat_comp_coeff['C']))

    def dT_expression(b, cobj, T):
        return (pressure_sat_comp.return_expression(b, cobj, T) *
                cobj.pressure_sat_comp_coeff['B'] *
                log(10)/(T + cobj.pressure_sat_comp_coeff['C'])**2)
