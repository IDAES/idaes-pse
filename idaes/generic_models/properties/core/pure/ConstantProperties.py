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
Method to set constant pure component properties:

All parameter indicies and units based on conventions used by the source
"""
from pyomo.environ import log, Var, units as pyunits

from idaes.core.util.misc import set_param_from_config


# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
class Constant(object):
    
    class cp_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            cobj.cp_mol_liq_comp_coeff = Var(
                doc="Parameter for liquid phase molar heat capacity",
                units=pyunits.J*pyunits.kmol**-1*pyunits.K**-1)
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            cp = (cobj.cp_mol_liq_comp_coeff)
            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units["heat_capacity_mole"])


    class enth_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            # if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
            #     cp_mol_liq_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_liq_comp_ref = Var(
                        doc="Liquid phase molar heat of formation @ Tref",
                        units=units["energy_mole"])
                set_param_from_config(cobj, param="enth_mol_form_liq_comp_ref")

            cobj.enth_mol_form_liq_comp_coeff = Var(
                doc="Liquid phase molar heat of formation",
                units=pyunits.kJ/pyunits.mol)
            set_param_from_config(cobj, param="enth_mol_form_liq_comp_coeff")

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
                    (cobj.enth_mol_form_liq_comp_coeff), units["energy_mole"]) +
                h_form)

            return h


    class entr_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            # if not hasattr(cobj, "cp_mol_liq_comp_coeff_1"):
            #     cp_mol_liq_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_liq_comp_ref = Var(
                    doc="Liquid phase molar entropy of formation @ Tref",
                    units=units["entropy_mole"])
            set_param_from_config(cobj, param="entr_mol_form_liq_comp_ref")
           
            cobj.entr_mol_form_liq_comp_coeff = Var(
                    doc="Liquid phase molar entropy of formation",
                    units=units["entropy_mole"])
            set_param_from_config(cobj, param="entr_mol_form_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            s = (pyunits.convert(
                    (cobj.entr_mol_form_liq_comp_coeff),
                    units["entropy_mole"]) +
                cobj.entr_mol_form_liq_comp_ref)

            return s

    # -----------------------------------------------------------------------------
    # Densities
    class dens_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            cobj.dens_mol_liq_comp_coeff = Var(
                    doc="Parameter for liquid phase molar density",
                    units=pyunits.kmol*pyunits.m**-3)
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):

            rho = (cobj.dens_mol_liq_comp_coeff)
            units = b.params.get_metadata().derived_units
            return pyunits.convert(rho, units["density_mole"])