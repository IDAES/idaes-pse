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

"""
from pyomo.environ import log, Var, Param, units as pyunits

from idaes.core.util.misc import set_param_from_config

# -----------------------------------------------------------------------------
# Heat capacities, enthalpies and entropies
class Constant(object):
    
    # Ideal liquid properties methods
    class cp_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.cp_mol_liq_comp_coeff = Var(
                doc="Parameter for liquid phase molar heat capacity",
                units=units["heat_capacity_mole"])
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            cp = cobj.cp_mol_liq_comp_coeff
            return cp

    class enth_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
                Constant.cp_mol_liq_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_liq_comp_ref = Var(
                        doc="Liquid phase molar heat of formation @ Tref",
                        units=units["energy_mole"])
                set_param_from_config(cobj, param="enth_mol_form_liq_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            h_form = (cobj.enth_mol_form_liq_comp_ref if
                    b.params.config.include_enthalpy_of_formation
                    else 0*units["energy_mole"])

            h = cobj.cp_mol_liq_comp_coeff*(T-Tr) + h_form

            return h

    class entr_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
                Constant.cp_mol_liq_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_liq_comp_ref = Var(
                    doc="Liquid phase molar entropy of formation @ Tref",
                    units=units["entropy_mole"])
            set_param_from_config(cobj, param="entr_mol_form_liq_comp_ref")
           
        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            s = cobj.cp_mol_liq_comp_coeff*log(T/Tr) + cobj.entr_mol_form_liq_comp_ref

            return s

    class dens_mol_liq_comp(object):

        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.dens_mol_liq_comp_coeff = Var(
                    doc="Parameter for liquid phase molar density",
                    units=units["density_mole"])
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Molar density
            rho = cobj.dens_mol_liq_comp_coeff
            return rho
    

    # Ideal gas properties methods
    class cp_mol_ig_comp(object):

        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.cp_mol_ig_comp_coeff = Var(
                doc="Parameter for ideal gas molar heat capacity",
                units=units["heat_capacity_mole"])
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            cp = cobj.cp_mol_ig_comp_coeff
            return cp

    class enth_mol_ig_comp(object):

        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
                Constant.cp_mol_ig_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_ig_comp_ref = Var(
                        doc="Ideal gas molar heat of formation @ Tref",
                        units=units["energy_mole"])
                set_param_from_config(cobj, param="enth_mol_form_ig_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            h_form = (cobj.enth_mol_form_ig_comp_ref if
                    b.params.config.include_enthalpy_of_formation
                    else 0*units["energy_mole"])

            h = cobj.cp_mol_ig_comp_coeff*(T-Tr) + h_form

            return h

    class entr_mol_ig_comp(object):

        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
                Constant.cp_mol_ig_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_ig_comp_ref = Var(
                    doc="Ideal gas molar entropy of formation @ Tref",
                    units=units["entropy_mole"])
            set_param_from_config(cobj, param="entr_mol_form_ig_comp_ref")
           
        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            s = cobj.cp_mol_ig_comp_coeff*log(T/Tr) + cobj.entr_mol_form_ig_comp_ref

            return s
