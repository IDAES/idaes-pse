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

Perry's Chemical Engineers' Handbook, 7th Edition
Perry, Green, Maloney, 1997, McGraw-Hill

All parameter indicies and units based on conventions used by the source
"""
from pyomo.environ import log, Var, Param, units as pyunits

from idaes.core.util.misc import set_param_from_config

from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


# -----------------------------------------------------------------------------
class Perrys(object):
    # -----------------------------------------------------------------------------
    # Heat capacities, enthalpies and entropies
    class cp_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.cp_mol_liq_comp_coeff_1 = Var(
                doc="Parameter 1 for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-1,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="1")

            cobj.cp_mol_liq_comp_coeff_2 = Var(
                doc="Parameter 2 for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-2,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="2")

            cobj.cp_mol_liq_comp_coeff_3 = Var(
                doc="Parameter 3 for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-3,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="3")

            cobj.cp_mol_liq_comp_coeff_4 = Var(
                doc="Parameter 4 for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-4,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="4")

            cobj.cp_mol_liq_comp_coeff_5 = Var(
                doc="Parameter 5 for liquid phase molar heat capacity",
                units=pyunits.J * pyunits.kmol**-1 * pyunits.K**-5,
            )
            set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="5")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific heat capacity
            T = pyunits.convert(T, to_units=pyunits.K)
            cp = (
                cobj.cp_mol_liq_comp_coeff_5 * T**4
                + cobj.cp_mol_liq_comp_coeff_4 * T**3
                + cobj.cp_mol_liq_comp_coeff_3 * T**2
                + cobj.cp_mol_liq_comp_coeff_2 * T
                + cobj.cp_mol_liq_comp_coeff_1
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)

    class enth_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff_1"):
                Perrys.cp_mol_liq_comp.build_parameters(cobj)

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
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            h_form = (
                cobj.enth_mol_form_liq_comp_ref
                if b.params.config.include_enthalpy_of_formation
                else 0 * units.ENERGY_MOLE
            )

            h = (
                pyunits.convert(
                    (cobj.cp_mol_liq_comp_coeff_5 / 5) * (T**5 - Tr**5)
                    + (cobj.cp_mol_liq_comp_coeff_4 / 4) * (T**4 - Tr**4)
                    + (cobj.cp_mol_liq_comp_coeff_3 / 3) * (T**3 - Tr**3)
                    + (cobj.cp_mol_liq_comp_coeff_2 / 2) * (T**2 - Tr**2)
                    + cobj.cp_mol_liq_comp_coeff_1 * (T - Tr),
                    units.ENERGY_MOLE,
                )
                + h_form
            )

            return h

    class entr_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff_1"):
                Perrys.cp_mol_liq_comp.build_parameters(cobj)

            units = cobj.parent_block().get_metadata().derived_units

            cobj.entr_mol_form_liq_comp_ref = Var(
                doc="Liquid phase molar entropy of formation @ Tref",
                units=units.ENTROPY_MOLE,
            )
            set_param_from_config(cobj, param="entr_mol_form_liq_comp_ref")

        @staticmethod
        def return_expression(b, cobj, T):
            # Specific entropy
            T = pyunits.convert(T, to_units=pyunits.K)
            Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            units = b.params.get_metadata().derived_units

            s = (
                pyunits.convert(
                    (cobj.cp_mol_liq_comp_coeff_5 / 4) * (T**4 - Tr**4)
                    + (cobj.cp_mol_liq_comp_coeff_4 / 3) * (T**3 - Tr**3)
                    + (cobj.cp_mol_liq_comp_coeff_3 / 2) * (T**2 - Tr**2)
                    + cobj.cp_mol_liq_comp_coeff_2 * (T - Tr)
                    + cobj.cp_mol_liq_comp_coeff_1 * log(T / Tr),
                    units.ENTROPY_MOLE,
                )
                + cobj.entr_mol_form_liq_comp_ref
            )

            return s

    # -----------------------------------------------------------------------------
    # Densities

    class dens_mol_liq_comp_eqn_1:
        @staticmethod
        def build_parameters(cobj):
            cobj.dens_mol_liq_comp_coeff_1 = Var(
                doc="Parameter 1 for liquid phase molar density",
                units=pyunits.kmol * pyunits.m**-3,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="1")

            cobj.dens_mol_liq_comp_coeff_2 = Var(
                doc="Parameter 2 for liquid phase molar density",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="2")

            cobj.dens_mol_liq_comp_coeff_3 = Var(
                doc="Parameter 3 for liquid phase molar density", units=pyunits.K
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="3")

            cobj.dens_mol_liq_comp_coeff_4 = Var(
                doc="Parameter 4 for liquid phase molar density",
                units=pyunits.dimensionless,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="4")

        @staticmethod
        def return_expression(b, cobj, T):
            # pg. 2-98
            T = pyunits.convert(T, to_units=pyunits.K)

            rho = cobj.dens_mol_liq_comp_coeff_1 / cobj.dens_mol_liq_comp_coeff_2 ** (
                1
                + (1 - T / cobj.dens_mol_liq_comp_coeff_3)
                ** cobj.dens_mol_liq_comp_coeff_4
            )

            units = b.params.get_metadata().derived_units

            return pyunits.convert(rho, units.DENSITY_MOLE)

    class dens_mol_liq_comp_eqn_2:
        @staticmethod
        def build_parameters(cobj):
            cobj.dens_mol_liq_comp_coeff_1 = Var(
                doc="Parameter 1 for liquid phase molar density",
                units=pyunits.kmol / pyunits.m**3,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="1")

            cobj.dens_mol_liq_comp_coeff_2 = Var(
                doc="Parameter 2 for liquid phase molar density",
                units=pyunits.kmol / pyunits.m**3 / pyunits.K,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="2")

            cobj.dens_mol_liq_comp_coeff_3 = Var(
                doc="Parameter 3 for liquid phase molar density",
                units=pyunits.kmol / pyunits.m**3 / pyunits.K**2,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="3")

            cobj.dens_mol_liq_comp_coeff_4 = Var(
                doc="Parameter 4 for liquid phase molar density",
                units=pyunits.kmol / pyunits.m**3 / pyunits.K**3,
            )
            set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="4")

        @staticmethod
        def return_expression(b, cobj, T):
            # pg. 2-98
            T = pyunits.convert(T, to_units=pyunits.K)

            rho = (
                (cobj.dens_mol_liq_comp_coeff_1)
                + (cobj.dens_mol_liq_comp_coeff_2) * T
                + (cobj.dens_mol_liq_comp_coeff_3) * T**2
                + (cobj.dens_mol_liq_comp_coeff_4) * T**3
            )

            units = b.params.get_metadata().derived_units

            return pyunits.convert(rho, units.DENSITY_MOLE)

    class dens_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            cobj.dens_mol_liq_comp_coeff_eqn_type = Param(
                mutable=True, doc="Liquid molar density equation form"
            )

            set_param_from_config(
                cobj, param="dens_mol_liq_comp_coeff", index="eqn_type"
            )

            if cobj.dens_mol_liq_comp_coeff_eqn_type.value == 1:
                Perrys.dens_mol_liq_comp_eqn_1.build_parameters(cobj)
            elif cobj.dens_mol_liq_comp_coeff_eqn_type.value == 2:
                Perrys.dens_mol_liq_comp_eqn_2.build_parameters(cobj)
            else:
                raise ConfigurationError(
                    f"{cobj.name} unrecognized value for "
                    f"dens_mol_liq_comp equation type: "
                    f"{cobj.dens_mol_liq_comp_coeff_eqn_type}"
                )

        @staticmethod
        def return_expression(b, cobj, T):
            if cobj.dens_mol_liq_comp_coeff_eqn_type.value == 1:
                rho = Perrys.dens_mol_liq_comp_eqn_1.return_expression(b, cobj, T)
            elif cobj.dens_mol_liq_comp_coeff_eqn_type.value == 2:
                rho = Perrys.dens_mol_liq_comp_eqn_2().return_expression(b, cobj, T)
            else:
                raise ConfigurationError(
                    "No expression for eqn_type of "
                    "dens_mol_liq_comp_coeff specified,"
                    "please specify valid flag."
                )
            return rho
