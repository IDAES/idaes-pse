#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Methods for calculating pure component properties from:

The Properties of Gases & Liquids, 4th Edition
Reid, Prausnitz and Polling, 1987, McGraw-Hill

All parameter indicies and units based on conventions used by the source
"""
from pyomo.environ import exp, log, Var, units as pyunits

from idaes.core.util.misc import set_param_from_config


# -----------------------------------------------------------------------------
class RPP4(object):
    # -----------------------------------------------------------------------------
    # Heat capacities, enthalpies and entropies
    class cp_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.cp_mol_ig_comp_coeff_A = Var(
                doc="Parameter A for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.mol / pyunits.K,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="A")

            cobj.cp_mol_ig_comp_coeff_B = Var(
                doc="Parameter B for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.mol / pyunits.K**2,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="B")

            cobj.cp_mol_ig_comp_coeff_C = Var(
                doc="Parameter C for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.mol / pyunits.K**3,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="C")

            cobj.cp_mol_ig_comp_coeff_D = Var(
                doc="Parameter D for ideal gas molar heat capacity",
                units=pyunits.J / pyunits.mol / pyunits.K**4,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="D")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            T = pyunits.convert(T, to_units=pyunits.K)

            cp = (
                cobj.cp_mol_ig_comp_coeff_D * T**3
                + cobj.cp_mol_ig_comp_coeff_C * T**2
                + cobj.cp_mol_ig_comp_coeff_B * T
                + cobj.cp_mol_ig_comp_coeff_A
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)

    class enth_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
                RPP4.cp_mol_ig_comp.build_parameters(cobj)

            if cobj.parent_block().config.include_enthalpy_of_formation:
                units = cobj.parent_block().get_metadata().derived_units

                cobj.enth_mol_form_vap_comp_ref = Var(
                    doc="Vapor phase molar heat of formation @ Tref",
                    units=units.ENERGY_MOLE,
                )
                set_param_from_config(cobj, param="enth_mol_form_vap_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            h_form = (
                cobj.enth_mol_form_vap_comp_ref
                if b.params.config.include_enthalpy_of_formation
                else 0 * units.ENERGY_MOLE
            )

            h = (
                pyunits.convert(
                    (cobj.cp_mol_ig_comp_coeff_D / 4) * (T**4 - Tr**4)
                    + (cobj.cp_mol_ig_comp_coeff_C / 3) * (T**3 - Tr**3)
                    + (cobj.cp_mol_ig_comp_coeff_B / 2) * (T**2 - Tr**2)
                    + cobj.cp_mol_ig_comp_coeff_A * (T - Tr),
                    units.ENERGY_MOLE,
                )
                + h_form
            )

            return h

    class entr_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
                RPP4.cp_mol_ig_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_vap_comp_ref = Var(
                doc="Vapor phase molar entropy of formation @ Tref",
                units=units.ENTROPY_MOLE,
            )
            set_param_from_config(cobj, param="entr_mol_form_vap_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            s = (
                pyunits.convert(
                    (cobj.cp_mol_ig_comp_coeff_D / 3) * (T**3 - Tr**3)
                    + (cobj.cp_mol_ig_comp_coeff_C / 2) * (T**2 - Tr**2)
                    + cobj.cp_mol_ig_comp_coeff_B * (T - Tr)
                    + cobj.cp_mol_ig_comp_coeff_A * log(T / Tr),
                    units.ENTROPY_MOLE,
                )
                + cobj.entr_mol_form_vap_comp_ref
            )

            return s

    # -----------------------------------------------------------------------------
    # Saturation pressure
    # Note that this equation in not valid beyond the critical temperature
    class pressure_sat_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.pressure_sat_comp_coeff_A = Var(
                doc="Coefficient A for calculating Psat", units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="A")

            cobj.pressure_sat_comp_coeff_B = Var(
                doc="Coefficient B for calculating Psat", units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="B")

            cobj.pressure_sat_comp_coeff_C = Var(
                doc="Coefficient C for calculating Psat", units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C")

            cobj.pressure_sat_comp_coeff_D = Var(
                doc="Coefficient D for calculating Psat", units=pyunits.dimensionless
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="D")

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            if dT:
                return RPP4.pressure_sat_comp.dT_expression(b, cobj, T)

            x = 1 - T / cobj.temperature_crit

            return (
                exp(
                    (1 - x) ** -1
                    * (
                        cobj.pressure_sat_comp_coeff_A * x
                        + cobj.pressure_sat_comp_coeff_B * x**1.5
                        + cobj.pressure_sat_comp_coeff_C * x**3
                        + cobj.pressure_sat_comp_coeff_D * x**6
                    )
                )
                * cobj.pressure_crit
            )

        @staticmethod
        def dT_expression(b, cobj, T):
            x = 1 - T / cobj.temperature_crit

            return -RPP4.pressure_sat_comp.return_expression(b, cobj, T) * (
                (
                    cobj.pressure_sat_comp_coeff_A
                    + 1.5 * cobj.pressure_sat_comp_coeff_B * x**0.5
                    + 3 * cobj.pressure_sat_comp_coeff_C * x**2
                    + 6 * cobj.pressure_sat_comp_coeff_D * x**5
                )
                / T
                + (cobj.temperature_crit / T**2)
                * (
                    cobj.pressure_sat_comp_coeff_A * x
                    + cobj.pressure_sat_comp_coeff_B * x**1.5
                    + cobj.pressure_sat_comp_coeff_C * x**3
                    + cobj.pressure_sat_comp_coeff_D * x**6
                )
            )
