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
Pure component properties as used by the NIST WebBook

https://webbook.nist.gov/chemistry/

Retrieved: September 13th, 2019

All parameter indicies and units based on conventions used by the source
"""
from pyomo.environ import Expression, log, Var, units as pyunits

from idaes.core.util.misc import set_param_from_config


# -----------------------------------------------------------------------------
class NIST(object):
    # -----------------------------------------------------------------------------
    # Shomate Equation for heat capacities, enthalpy and entropy
    class cp_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.cp_mol_ig_comp_coeff_A = Var(
                doc="Shomate A parameter for ideal gas molar heat capacity",
                units=pyunits.J * pyunits.mol**-1 * pyunits.K**-1,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="A")

            cobj.cp_mol_ig_comp_coeff_B = Var(
                doc="Shomate B parameter for ideal gas molar heat capacity",
                units=pyunits.J
                * pyunits.mol**-1
                * pyunits.K**-1
                * pyunits.kiloK**-1,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="B")

            cobj.cp_mol_ig_comp_coeff_C = Var(
                doc="Shomate C parameter for ideal gas molar heat capacity",
                units=pyunits.J
                * pyunits.mol**-1
                * pyunits.K**-1
                * pyunits.kiloK**-2,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="C")

            cobj.cp_mol_ig_comp_coeff_D = Var(
                doc="Shomate D parameter for ideal gas molar heat capacity",
                units=pyunits.J
                * pyunits.mol**-1
                * pyunits.K**-1
                * pyunits.kiloK**-3,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="D")

            cobj.cp_mol_ig_comp_coeff_E = Var(
                doc="Shomate E parameter for ideal gas molar heat capacity",
                units=pyunits.J
                * pyunits.mol**-1
                * pyunits.K**-1
                * pyunits.kiloK**2,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="E")

            cobj.cp_mol_ig_comp_coeff_F = Var(
                doc="Shomate F parameter for ideal gas molar heat capacity",
                units=pyunits.kJ * pyunits.mol**-1,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="F")

            cobj.cp_mol_ig_comp_coeff_G = Var(
                doc="Shomate G parameter for ideal gas molar heat capacity",
                units=pyunits.J * pyunits.mol**-1 * pyunits.K**-1,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="G")

            cobj.cp_mol_ig_comp_coeff_H = Var(
                doc="Shomate H parameter for ideal gas molar heat capacity",
                units=pyunits.kJ * pyunits.mol**-1,
            )
            set_param_from_config(cobj, param="cp_mol_ig_comp_coeff", index="H")

            # As the H parameter is the specific heat of formation, create an
            # Expression which converts this to the base units with standard name
            units = cobj.parent_block().get_metadata().derived_units
            cobj.enth_mol_form_vap_comp_ref = Expression(
                expr=pyunits.convert(
                    cobj.cp_mol_ig_comp_coeff_H, to_units=units.ENERGY_MOLE
                ),
                doc="Vapor phase molar heat of formation @ Tref",
            )

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity (const. P)  via the Shomate equation
            t = pyunits.convert(T, to_units=pyunits.kiloK)
            cp = (
                cobj.cp_mol_ig_comp_coeff_A
                + cobj.cp_mol_ig_comp_coeff_B * t
                + cobj.cp_mol_ig_comp_coeff_C * t**2
                + cobj.cp_mol_ig_comp_coeff_D * t**3
                + cobj.cp_mol_ig_comp_coeff_E * t**-2
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)

    class enth_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
                NIST.cp_mol_ig_comp.build_parameters(cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific enthalpy via the Shomate equation
            t = pyunits.convert(T, to_units=pyunits.kiloK)

            h_form = b.params.config.include_enthalpy_of_formation
            H = (
                cobj.cp_mol_ig_comp_coeff_H
                if not h_form
                else 0 * pyunits.kJ * pyunits.mol**-1
            )

            h = (
                cobj.cp_mol_ig_comp_coeff_A * t
                + (cobj.cp_mol_ig_comp_coeff_B / 2) * t**2
                + (cobj.cp_mol_ig_comp_coeff_C / 3) * t**3
                + (cobj.cp_mol_ig_comp_coeff_D / 4) * t**4
                - cobj.cp_mol_ig_comp_coeff_E / t
                + cobj.cp_mol_ig_comp_coeff_F
                - H
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(h, units.ENERGY_MOLE)

    class entr_mol_ig_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
                NIST.cp_mol_ig_comp.build_parameters(cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy via the Shomate equation
            t = pyunits.convert(T, to_units=pyunits.kiloK)

            s = (
                cobj.cp_mol_ig_comp_coeff_A * log(t / pyunits.kiloK)
                + cobj.cp_mol_ig_comp_coeff_B * t  # need to make unitless
                + (cobj.cp_mol_ig_comp_coeff_C / 2) * t**2
                + (cobj.cp_mol_ig_comp_coeff_D / 3) * t**3
                - (cobj.cp_mol_ig_comp_coeff_E / 2) * t**-2
                + cobj.cp_mol_ig_comp_coeff_G
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(s, units.ENTROPY_MOLE)

    # -----------------------------------------------------------------------------
    # Antoine equation for saturation pressure
    class pressure_sat_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.pressure_sat_comp_coeff_A = Var(
                doc="Antoine A coefficient for calculating Psat", units=None
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
                return NIST.pressure_sat_comp.dT_expression(b, cobj, T)

            psat = (
                10
                ** (
                    cobj.pressure_sat_comp_coeff_A
                    - cobj.pressure_sat_comp_coeff_B
                    / (
                        pyunits.convert(T, to_units=pyunits.K)
                        + cobj.pressure_sat_comp_coeff_C
                    )
                )
                * pyunits.bar
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(psat, to_units=units.PRESSURE)

        @staticmethod
        def dT_expression(b, cobj, T):
            p_sat_dT = (
                NIST.pressure_sat_comp.return_expression(b, cobj, T)
                * cobj.pressure_sat_comp_coeff_B
                * log(10)
                / (
                    pyunits.convert(T, to_units=pyunits.K)
                    + cobj.pressure_sat_comp_coeff_C
                )
                ** 2
            )

            units = b.params.get_metadata().derived_units
            dp_units = units.PRESSURE / units.TEMPERATURE
            return pyunits.convert(p_sat_dT, to_units=dp_units)
