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
Methods for calculating pure component properties from:

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

All parameter indicies based on conventions used by the source
"""
from pyomo.environ import log, Var


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
class cp_mol_liq_comp():
    def build_parameters(cobj):
        cobj.cp_mol_liq_comp_coeff = Var(
                ['1', '2', '3', '4', '5'],
                initialize=cobj.config.parameter_data["cp_mol_liq_comp_coeff"],
                doc="Parameters for liquid phase molar heat capacity")

    def return_expression(b, cobj, T):
        # Specific heat capacity
        return 1e-3*(cobj.cp_mol_liq_comp_coeff["5"]*T**4 +
                     cobj.cp_mol_liq_comp_coeff["4"]*T**3 +
                     cobj.cp_mol_liq_comp_coeff["3"]*T**2 +
                     cobj.cp_mol_liq_comp_coeff["2"]*T +
                     cobj.cp_mol_liq_comp_coeff["1"])


class enth_mol_liq_comp():
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
            cp_mol_liq_comp.build_parameters(cobj)

        cobj.enth_mol_form_liq_comp_ref = Var(
                initialize=cobj.config.parameter_data[
                    "enth_mol_form_liq_comp_ref"],
                doc="Liquid phase molar heat of formation @ Tref")

    def return_expression(b, cobj, T):
        # Specific enthalpy
        return (1e-3*((cobj.cp_mol_liq_comp_coeff["5"]/5) *
                      (T**5-b.params.temperature_ref**5) +
                      (cobj.cp_mol_liq_comp_coeff["4"]/4) *
                      (T**4-b.params.temperature_ref**4) +
                      (cobj.cp_mol_liq_comp_coeff["3"]/3) *
                      (T**3-b.params.temperature_ref**3) +
                      (cobj.cp_mol_liq_comp_coeff["2"]/2) *
                      (T**2-b.params.temperature_ref**2) +
                      cobj.cp_mol_liq_comp_coeff["1"] *
                      (T-b.params.temperature_ref)) +
                cobj.enth_mol_form_liq_comp_ref)


class entr_mol_liq_comp():
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
            cp_mol_liq_comp.build_parameters(cobj)

        cobj.entr_mol_form_liq_comp_ref = Var(
                initialize=cobj.config.parameter_data[
                    "entr_mol_form_liq_comp_ref"],
                doc="Liquid phase molar entropy of formation @ Tref")

    def return_expression(b, cobj, T):
        # Specific entropy
        return (1e-3*((cobj.cp_mol_liq_comp_coeff['5']/4)*T**4 +
                      (cobj.cp_mol_liq_comp_coeff['4']/3)*T**3 +
                      (cobj.cp_mol_liq_comp_coeff['3']/2)*T**2 +
                      cobj.cp_mol_liq_comp_coeff['2']*T +
                      cobj.cp_mol_liq_comp_coeff['1']*log(T)) +
                cobj.entr_mol_form_liq_comp_ref)


# -----------------------------------------------------------------------------
# Densities
class dens_mol_liq_comp():
    def build_parameters(cobj):
        cobj.dens_mol_liq_comp_coeff = Var(
                ['1', '2', '3', '4'],
                initialize=cobj.config.parameter_data[
                    "dens_mol_liq_comp_coeff"],
                doc="Parameters for liquid phase molar density")

    def return_expression(b, cobj, T):
        # pg. 2-98
        return (cobj.dens_mol_liq_comp_coeff['1'] /
                cobj.dens_mol_liq_comp_coeff['2']**(
                        1 + (1-T/cobj.dens_mol_liq_comp_coeff['3']) **
                        cobj.dens_mol_liq_comp_coeff['4']))
