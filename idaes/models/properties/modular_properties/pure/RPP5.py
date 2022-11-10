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

The Properties of Gases & Liquids, 5th Edition
Reid, Prausnitz and Polling, 2001, McGraw-Hill

All parameter indicies based on conventions used by the source
"""

from pyomo.environ import log, Var, units as pyunits

from idaes.core.util.misc import set_param_from_config

from idaes.core.util.constants import Constants as const


# -----------------------------------------------------------------------------
class RPP5(object):
    # -----------------------------------------------------------------------------
    # Heat capacities, enthalpies and entropies
    class cp_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.cp_mol_ig_comp_coeff_a0 = Var(
                doc="Parameter a0 for ideal gas molar heat capacity", units=None
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="a0")

            cobj.cp_mol_ig_comp_coeff_a1 = Var(
                doc="Parameter a1 for ideal gas molar heat capacity",
                units=pyunits.K**-1,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="a1")

            cobj.cp_mol_ig_comp_coeff_a2 = Var(
                doc="Parameter a2 for ideal gas molar heat capacity",
                units=pyunits.K**-2,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="a2")

            cobj.cp_mol_ig_comp_coeff_a3 = Var(
                doc="Parameter a3 for ideal gas molar heat capacity",
                units=pyunits.K**-3,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="a3")

            cobj.cp_mol_ig_comp_coeff_a4 = Var(
                doc="Parameter a4 for ideal gas molar heat capacity",
                units=pyunits.K**-4,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="a4")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            T = pyunits.convert(T, to_units=pyunits.K)

            cp = (
                cobj.cp_mol_ig_comp_coeff_a4 * T**4
                + cobj.cp_mol_ig_comp_coeff_a3 * T**3
                + cobj.cp_mol_ig_comp_coeff_a2 * T**2
                + cobj.cp_mol_ig_comp_coeff_a1 * T
                + cobj.cp_mol_ig_comp_coeff_a0
            ) * const.gas_constant

            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)

    class enth_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_a0"):
                RPP5.cp_mol_ig_comp.build_parameters(cobj)

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

            h = (
                pyunits.convert(
                    (
                        (cobj.cp_mol_ig_comp_coeff_a4 / 5) * (T**5 - Tr**5)
                        + (cobj.cp_mol_ig_comp_coeff_a3 / 4) * (T**4 - Tr**4)
                        + (cobj.cp_mol_ig_comp_coeff_a2 / 3) * (T**3 - Tr**3)
                        + (cobj.cp_mol_ig_comp_coeff_a1 / 2) * (T**2 - Tr**2)
                        + cobj.cp_mol_ig_comp_coeff_a0 * (T - Tr)
                    )
                    * const.gas_constant,
                    units.ENERGY_MOLE,
                )
                + cobj.enth_mol_form_vap_comp_ref
            )

            return h

    class entr_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_a0"):
                RPP5.cp_mol_ig_comp.build_parameters(cobj)

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
                    (
                        (cobj.cp_mol_ig_comp_coeff_a4 / 4) * (T**4 - Tr**4)
                        + (cobj.cp_mol_ig_comp_coeff_a3 / 3) * (T**3 - Tr**3)
                        + (cobj.cp_mol_ig_comp_coeff_a2 / 2) * (T**2 - Tr**2)
                        + cobj.cp_mol_ig_comp_coeff_a1 * (T - Tr)
                        + cobj.cp_mol_ig_comp_coeff_a0 * log(T / Tr)
                    )
                    * const.gas_constant,
                    units.ENTROPY_MOLE,
                )
                + cobj.entr_mol_form_vap_comp_ref
            )

            return s

    # -----------------------------------------------------------------------------
    # Antoine equation for saturation pressure
    class pressure_sat_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.pressure_sat_comp_coeff_A = Var(
                doc="Antoine A coefficient for calculating Psat",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="A")

            cobj.pressure_sat_comp_coeff_B = Var(
                doc="Antoine B coefficient for calculating Psat", units=pyunits.K
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="B")

            cobj.pressure_sat_comp_coeff_C = Var(
                doc="Antoine C coefficient for calculating Psat", units=pyunits.K
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="C")

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            if dT:
                return RPP5.pressure_sat_comp.dT_expression(b, cobj, T)

            psat = (
                10
                ** (
                    cobj.pressure_sat_comp_coeff_A
                    - (
                        cobj.pressure_sat_comp_coeff_B
                        / (T + cobj.pressure_sat_comp_coeff_C - 273.15 * pyunits.degK)
                    )
                )
                * pyunits.bar
            )

            base_units = b.params.get_metadata().default_units
            p_units = base_units.PRESSURE
            return pyunits.convert(psat, to_units=p_units)

        @staticmethod
        def dT_expression(b, cobj, T):
            p_sat_dT = (
                RPP5.pressure_sat_comp.return_expression(b, cobj, T)
                * cobj.pressure_sat_comp_coeff_B
                * log(10)
                / (T + cobj.pressure_sat_comp_coeff_C - 273.15 * pyunits.degK) ** 2
            )

            base_units = b.params.get_metadata().default_units
            dp_units = base_units.PRESSURE * base_units.TEMPERATURE**-1
            return pyunits.convert(p_sat_dT, to_units=dp_units)
