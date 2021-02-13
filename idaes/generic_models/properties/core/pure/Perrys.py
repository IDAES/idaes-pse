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

from idaes.core.util.misc import set_param_from_config


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
class cp_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        cobj.cp_mol_liq_comp_coeff_1 = Var(
            doc="Parameter 1 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-1)
        set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="1")

        cobj.cp_mol_liq_comp_coeff_2 = Var(
            doc="Parameter 2 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-2)
        set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="2")

        cobj.cp_mol_liq_comp_coeff_3 = Var(
            doc="Parameter 3 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-3)
        set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="3")

        cobj.cp_mol_liq_comp_coeff_4 = Var(
            doc="Parameter 4 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-4)
        set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="4")

        cobj.cp_mol_liq_comp_coeff_5 = Var(
            doc="Parameter 5 for liquid phase molar heat capacity",
            units=pyunits.J*pyunits.kmol**-1*pyunits.K**-5)
        set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="5")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific heat capacity
        T = pyunits.convert(T, to_units=pyunits.K)
        cp = (cobj.cp_mol_liq_comp_coeff_5*T**4 +
              cobj.cp_mol_liq_comp_coeff_4*T**3 +
              cobj.cp_mol_liq_comp_coeff_3*T**2 +
              cobj.cp_mol_liq_comp_coeff_2*T +
              cobj.cp_mol_liq_comp_coeff_1)

        units = b.params.get_metadata().derived_units
        return pyunits.convert(cp, units["heat_capacity_mole"])


class enth_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_liq_comp_coeff_1"):
            cp_mol_liq_comp.build_parameters(cobj)

        if cobj.parent_block().config.include_enthalpy_of_formation:
            units = cobj.parent_block().get_metadata().derived_units

            cobj.enth_mol_form_liq_comp_ref = Var(
                    doc="Liquid phase molar heat of formation @ Tref",
                    units=units["energy_mole"])
            set_param_from_config(cobj, param="enth_mol_form_liq_comp_ref")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific enthalpy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        units = b.params.get_metadata().derived_units

        h_form = (cobj.enth_mol_form_liq_comp_ref if
                  b.params.config.include_enthalpy_of_formation
                  else 0*units["energy_mole"])

        h = (pyunits.convert(
                (cobj.cp_mol_liq_comp_coeff_5/5)*(T**5-Tr**5) +
                (cobj.cp_mol_liq_comp_coeff_4/4)*(T**4-Tr**4) +
                (cobj.cp_mol_liq_comp_coeff_3/3)*(T**3-Tr**3) +
                (cobj.cp_mol_liq_comp_coeff_2/2)*(T**2-Tr**2) +
                cobj.cp_mol_liq_comp_coeff_1*(T-Tr), units["energy_mole"]) +
             h_form)

        return h


class entr_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_liq_comp_coeff_1"):
            cp_mol_liq_comp.build_parameters(cobj)

        units = cobj.parent_block().get_metadata().derived_units

        cobj.entr_mol_form_liq_comp_ref = Var(
                doc="Liquid phase molar entropy of formation @ Tref",
                units=units["entropy_mole"])
        set_param_from_config(cobj, param="entr_mol_form_liq_comp_ref")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific entropy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        units = b.params.get_metadata().derived_units

        s = (pyunits.convert(
                (cobj.cp_mol_liq_comp_coeff_5/4)*(T**4-Tr**4) +
                (cobj.cp_mol_liq_comp_coeff_4/3)*(T**3-Tr**3) +
                (cobj.cp_mol_liq_comp_coeff_3/2)*(T**2-Tr**2) +
                cobj.cp_mol_liq_comp_coeff_2*(T-Tr) +
                cobj.cp_mol_liq_comp_coeff_1*log(T/Tr),
                units["entropy_mole"]) +
             cobj.entr_mol_form_liq_comp_ref)

        return s


# -----------------------------------------------------------------------------
# Densities
class dens_mol_liq_comp():

    @staticmethod
    def build_parameters(cobj):
        cobj.dens_mol_liq_comp_coeff_1 = Var(
                doc="Parameter 1 for liquid phase molar density",
                units=pyunits.kmol*pyunits.m**-3)
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="1")

        cobj.dens_mol_liq_comp_coeff_2 = Var(
                doc="Parameter 2 for liquid phase molar density",
                units=pyunits.dimensionless)
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="2")

        cobj.dens_mol_liq_comp_coeff_3 = Var(
                doc="Parameter 3 for liquid phase molar density",
                units=pyunits.K)
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="3")

        cobj.dens_mol_liq_comp_coeff_4 = Var(
                doc="Parameter 4 for liquid phase molar density",
                units=pyunits.dimensionless)
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="4")

    @staticmethod
    def return_expression(b, cobj, T):
        # pg. 2-98
        T = pyunits.convert(T, to_units=pyunits.K)

        rho = (cobj.dens_mol_liq_comp_coeff_1 /
               cobj.dens_mol_liq_comp_coeff_2**(
                   1 + (1-T/cobj.dens_mol_liq_comp_coeff_3) **
                   cobj.dens_mol_liq_comp_coeff_4))

        units = b.params.get_metadata().derived_units

        return pyunits.convert(rho, units["density_mole"])


# -----------------------------------------------------------------------------
class Perrys(object):
    cp_mol_liq_comp = cp_mol_liq_comp
    enth_mol_liq_comp = enth_mol_liq_comp
    entr_mol_liq_comp = entr_mol_liq_comp
    dens_mol_liq_comp = dens_mol_liq_comp
