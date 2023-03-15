#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Method to set constant pure component properties:

"""
# TODO: Missing doc strings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from pyomo.environ import log, Var

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
                units=units.HEAT_CAPACITY_MOLE,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            return cobj.cp_mol_liq_comp_coeff

    class enth_mol_liq_comp(object):
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
                Constant.cp_mol_liq_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_liq_comp_ref = Var(
                    doc="Liquid phase molar heat of formation @ Tref",
                    units=units.ENERGY_MOLE,
                )
                set_param_from_config(cobj, param="enth_mol_form_liq_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            h_form = (
                cobj.enth_mol_form_liq_comp_ref
                if b.params.config.include_enthalpy_of_formation
                else 0 * units.ENERGY_MOLE
            )

            return cobj.cp_mol_liq_comp_coeff * (T - Tr) + h_form

    class entr_mol_liq_comp(object):
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff"):
                Constant.cp_mol_liq_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_liq_comp_ref = Var(
                doc="Liquid phase molar entropy of formation @ Tref",
                units=units.ENTROPY_MOLE,
            )
            set_param_from_config(cobj, param="entr_mol_form_liq_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            Tr = b.params.temperature_ref

            return (
                cobj.cp_mol_liq_comp_coeff * log(T / Tr)
                + cobj.entr_mol_form_liq_comp_ref
            )

    class dens_mol_liq_comp(object):
        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.dens_mol_liq_comp_coeff = Var(
                doc="Parameter for liquid phase molar density",
                units=units.DENSITY_MOLE,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Molar density
            return cobj.dens_mol_liq_comp_coeff

    # Ideal gas properties methods
    class cp_mol_ig_comp(object):
        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.cp_mol_ig_comp_coeff = Var(
                doc="Parameter for ideal gas molar heat capacity",
                units=units.HEAT_CAPACITY_MOLE,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            return cobj.cp_mol_ig_comp_coeff

    class enth_mol_ig_comp(object):
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
                Constant.cp_mol_ig_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_ig_comp_ref = Var(
                    doc="Ideal gas molar heat of formation @ Tref",
                    units=units.ENERGY_MOLE,
                )
                set_param_from_config(cobj, param="enth_mol_form_ig_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            h_form = (
                cobj.enth_mol_form_ig_comp_ref
                if b.params.config.include_enthalpy_of_formation
                else 0 * units.ENERGY_MOLE
            )

            return cobj.cp_mol_ig_comp_coeff * (T - Tr) + h_form

    class entr_mol_ig_comp(object):
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff"):
                Constant.cp_mol_ig_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_ig_comp_ref = Var(
                doc="Ideal gas molar entropy of formation @ Tref",
                units=units.ENTROPY_MOLE,
            )
            set_param_from_config(cobj, param="entr_mol_form_ig_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            Tr = b.params.temperature_ref

            return (
                cobj.cp_mol_ig_comp_coeff * log(T / Tr) + cobj.entr_mol_form_ig_comp_ref
            )

    # Ideal solid properties methods
    class cp_mol_sol_comp(object):
        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.cp_mol_sol_comp_coeff = Var(
                doc="Parameter for solid phase molar heat capacity",
                units=units.HEAT_CAPACITY_MOLE,
            )
            set_param_from_config(cobj, param="cp_mol_sol_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            cp = cobj.cp_mol_sol_comp_coeff
            return cp

    class enth_mol_sol_comp(object):
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_sol_comp_coeff"):
                Constant.cp_mol_sol_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_sol_comp_ref = Var(
                    doc="Solid phase molar heat of formation @ Tref",
                    units=units.ENERGY_MOLE,
                )
                set_param_from_config(cobj, param="enth_mol_form_sol_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            units = b.params.get_metadata().derived_units
            Tr = b.params.temperature_ref

            h_form = (
                cobj.enth_mol_form_sol_comp_ref
                if b.params.config.include_enthalpy_of_formation
                else 0 * units.ENERGY_MOLE
            )

            return cobj.cp_mol_sol_comp_coeff * (T - Tr) + h_form

    class entr_mol_sol_comp(object):
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_sol_comp_coeff"):
                Constant.cp_mol_sol_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_sol_comp_ref = Var(
                doc="Solid phase molar entropy of formation @ Tref",
                units=units.ENTROPY_MOLE,
            )
            set_param_from_config(cobj, param="entr_mol_form_sol_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            Tr = b.params.temperature_ref

            return (
                cobj.cp_mol_sol_comp_coeff * log(T / Tr)
                + cobj.entr_mol_form_sol_comp_ref
            )

    class dens_mol_sol_comp(object):
        @staticmethod
        def build_parameters(cobj):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.dens_mol_sol_comp_coeff = Var(
                doc="Parameter for solid phase molar density",
                units=units.DENSITY_MOLE,
            )
            set_param_from_config(cobj, param="dens_mol_sol_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, T):
            # Molar density
            return cobj.dens_mol_sol_comp_coeff

    class visc_d_phase_comp(object):
        @staticmethod
        def build_parameters(cobj, p):
            units = cobj.parent_block().get_metadata().derived_units
            # Calling this a "coefficient" doesn't make much sense, but want to be consistent with other methods
            cobj.add_component(
                f"visc_d_{p}_comp_coeff",
                Var(
                    doc=f"Parameter for {p} phase dynamic viscosity",
                    units=units["dynamic_viscosity"],
                ),
            )
            set_param_from_config(cobj, param=f"visc_d_{p}_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, p, T):
            return getattr(cobj, f"visc_d_{p}_comp_coeff")

    class therm_cond_phase_comp(object):
        @staticmethod
        def build_parameters(cobj, p):
            units = cobj.parent_block().get_metadata().derived_units
            cobj.add_component(
                f"therm_cond_{p}_comp_coeff",
                Var(
                    doc=f"Parameter for {p} phase thermal conductivity",
                    units=units["thermal_conductivity"],
                ),
            )
            set_param_from_config(cobj, param=f"therm_cond_{p}_comp_coeff")

        @staticmethod
        def return_expression(b, cobj, p, T):
            return getattr(cobj, f"therm_cond_{p}_comp_coeff")
