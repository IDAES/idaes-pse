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

All parameter indicies and units based on conventions used by the source
"""
from pyomo.environ import log, Var, units as pyunits


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
class cp_mol_liq_comp():
    def build_parameters(cobj):
        cobj.cp_mol_liq_comp_coeff_1 = Var(
            initialize=cobj.config.parameter_data["cp_mol_liq_comp_coeff"]['1'],
            doc="Parameter 1 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-1)
        cobj.cp_mol_liq_comp_coeff_2 = Var(
            initialize=cobj.config.parameter_data["cp_mol_liq_comp_coeff"]['2'],
            doc="Parameter 2 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-2)
        cobj.cp_mol_liq_comp_coeff_3 = Var(
            initialize=cobj.config.parameter_data["cp_mol_liq_comp_coeff"]['3'],
            doc="Parameter 3 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-3)
        cobj.cp_mol_liq_comp_coeff_4 = Var(
            initialize=cobj.config.parameter_data["cp_mol_liq_comp_coeff"]['4'],
            doc="Parameter 4 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-4)
        cobj.cp_mol_liq_comp_coeff_5 = Var(
            initialize=cobj.config.parameter_data["cp_mol_liq_comp_coeff"]['5'],
            doc="Parameter 5 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-5)

    def return_expression(b, cobj, T):
        # Specific heat capacity
        T = pyunits.convert(T, to_units=pyunits.K)
        cp = (cobj.cp_mol_liq_comp_coeff_5*T**4 +
              cobj.cp_mol_liq_comp_coeff_4*T**3 +
              cobj.cp_mol_liq_comp_coeff_3*T**2 +
              cobj.cp_mol_liq_comp_coeff_2*T +
              cobj.cp_mol_liq_comp_coeff_1)

        base_units = b.params.get_metadata().default_units
        cp_units = (base_units["mass"] *
                    base_units["length"]**2 *
                    base_units["time"]**-2 *
                    base_units["amount"]**-1 *
                    base_units["temperature"]**-1)
        return pyunits.convert(cp, cp_units)


class enth_mol_liq_comp():
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_liq_comp_coeff_1"):
            cp_mol_liq_comp.build_parameters(cobj)

        base_units = cobj.parent_block().get_metadata().default_units
        h_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1)

        cobj.enth_mol_form_liq_comp_ref = Var(
                initialize=cobj.config.parameter_data[
                    "enth_mol_form_liq_comp_ref"],
                doc="Liquid phase molar heat of formation @ Tref",
                units=h_units)

    def return_expression(b, cobj, T):
        # Specific enthalpy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        base_units = b.params.get_metadata().default_units
        h_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1)

        h = (pyunits.convert(
                (cobj.cp_mol_liq_comp_coeff_5/5)*(T**5-Tr**5) +
                (cobj.cp_mol_liq_comp_coeff_4/4)*(T**4-Tr**4) +
                (cobj.cp_mol_liq_comp_coeff_3/3)*(T**3-Tr**3) +
                (cobj.cp_mol_liq_comp_coeff_2/2)*(T**2-Tr**2) +
                cobj.cp_mol_liq_comp_coeff_1*(T-Tr), h_units) +
             cobj.enth_mol_form_liq_comp_ref)

        return h


class entr_mol_liq_comp():
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_liq_comp_coeff_1"):
            cp_mol_liq_comp.build_parameters(cobj)

        base_units = cobj.parent_block().get_metadata().default_units
        s_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1 *
                   base_units["temperature"]**-1)

        cobj.entr_mol_form_liq_comp_ref = Var(
                initialize=cobj.config.parameter_data[
                    "entr_mol_form_liq_comp_ref"],
                doc="Liquid phase molar entropy of formation @ Tref",
                units=s_units)

    def return_expression(b, cobj, T):
        # Specific entropy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        base_units = b.params.get_metadata().default_units
        s_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1 *
                   base_units["temperature"]**-1)

        s = (pyunits.convert(
                (cobj.cp_mol_liq_comp_coeff_5/4)*(T**4-Tr**4) +
                (cobj.cp_mol_liq_comp_coeff_4/3)*(T**3-Tr**3) +
                (cobj.cp_mol_liq_comp_coeff_3/2)*(T**2-Tr**2) +
                cobj.cp_mol_liq_comp_coeff_2*(T-Tr) +
                cobj.cp_mol_liq_comp_coeff_1*log(T/Tr), s_units) +
             cobj.entr_mol_form_liq_comp_ref)

        return s


# -----------------------------------------------------------------------------
# Densities
class dens_mol_liq_comp():
    def build_parameters(cobj):
        cobj.dens_mol_liq_comp_coeff_1 = Var(
                initialize=cobj.config.parameter_data[
                    "dens_mol_liq_comp_coeff"]["1"],
                doc="Parameter 1 for liquid phase molar density",
                units=pyunits.kmol*pyunits.m**-3)
        cobj.dens_mol_liq_comp_coeff_2 = Var(
                initialize=cobj.config.parameter_data[
                    "dens_mol_liq_comp_coeff"]["2"],
                doc="Parameter 2 for liquid phase molar density",
                units=None)
        cobj.dens_mol_liq_comp_coeff_3 = Var(
                initialize=cobj.config.parameter_data[
                    "dens_mol_liq_comp_coeff"]["3"],
                doc="Parameter 3 for liquid phase molar density",
                units=pyunits.K)
        cobj.dens_mol_liq_comp_coeff_4 = Var(
                initialize=cobj.config.parameter_data[
                    "dens_mol_liq_comp_coeff"]["4"],
                doc="Parameter 4 for liquid phase molar density",
                units=None)

    def return_expression(b, cobj, T):
        # pg. 2-98
        T = pyunits.convert(T, to_units=pyunits.K)

        rho = (cobj.dens_mol_liq_comp_coeff_1 /
               cobj.dens_mol_liq_comp_coeff_2**(
                   1 + (1-T/cobj.dens_mol_liq_comp_coeff_3) **
                   cobj.dens_mol_liq_comp_coeff_4))

        base_units = b.params.get_metadata().default_units
        rho_units = (base_units["amount"] *
                     base_units["length"]**-3)

        return pyunits.convert(rho, rho_units)
