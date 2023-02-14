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
Base class for EoS modules.

Raises NotImplementedErrors for all expected methods in case developer misses
some. EoS developers should overload all these methods.
"""
from pyomo.environ import units as pyunits
from idaes.core.util.constants import Constants as const
from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_component_object as cobj,
)
from idaes.core.util.exceptions import ConfigurationError


class EoSBase:
    @staticmethod
    def gas_constant(b):
        # Utility method to convert gas constant to base units
        base_units = b.params.get_metadata().default_units

        return pyunits.convert(const.gas_constant, to_units=base_units.GAS_CONSTANT)

    @staticmethod
    def common(b, pobj):
        raise NotImplementedError(_msg(b, "common"))

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        raise NotImplementedError(_msg(b, "calculate_scaling_factors"))

    @staticmethod
    def build_parameters(b):
        raise NotImplementedError(_msg(b, "build_parameters"))

    @staticmethod
    def get_vol_mol_pure(b, phase, comp, temperature):
        try:
            vol_mol = get_method(b, "vol_mol_" + phase + "_comp", comp)(
                b, cobj(b, comp), temperature
            )
        except (AttributeError, ConfigurationError):
            # vol_mol not defined, try for dens_mol instead
            try:
                vol_mol = 1 / get_method(b, "dens_mol_" + phase + "_comp", comp)(
                    b, cobj(b, comp), temperature
                )
            except (AttributeError, ConfigurationError):
                # Does not have either vol_mol or dens_mol
                suffix = "_" + phase + "_comp"
                raise ConfigurationError(
                    f"{b.name} does not have a method defined to use "
                    f"when calculating molar volume and density for "
                    f"component {comp} in phase {phase}. Each component "
                    f"must define a method for either vol_mol{suffix} or "
                    f"dens_mol{suffix}."
                )
        return vol_mol

    @staticmethod
    def act_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "act_phase_comp"))

    @staticmethod
    def act_phase_comp_true(b, p, j):
        raise NotImplementedError(_msg(b, "act_phase_comp_true"))

    @staticmethod
    def act_phase_comp_appr(b, p, j):
        raise NotImplementedError(_msg(b, "act_phase_comp_appr"))

    @staticmethod
    def act_coeff_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "act_coeff_phase_comp"))

    @staticmethod
    def act_coeff_phase_comp_true(b, p, j):
        raise NotImplementedError(_msg(b, "act_coeff_phase_comp_true"))

    @staticmethod
    def act_coeff_phase_comp_appr(b, p, j):
        raise NotImplementedError(_msg(b, "act_coeff_phase_comp_appr"))

    @staticmethod
    def cp_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "cp_mol_phase"))

    @staticmethod
    def cp_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "cp_mol_phase_comp"))

    @staticmethod
    def cv_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "cv_mol_phase"))

    @staticmethod
    def cv_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "cv_mol_phase_comp"))

    @staticmethod
    def heat_capacity_ratio_phase(b, p):
        return b.cp_mol_phase[p] / b.cv_mol_phase[p]

    @staticmethod
    def cv_mol_ig_comp_pure(b, j):
        # Method for calculating pure component ideal gas cv from cp
        # For ideal gases, cv = cp - R
        units = b.params.get_metadata().derived_units
        R = pyunits.convert(const.gas_constant, to_units=units.HEAT_CAPACITY_MOLE)
        return get_method(b, "cp_mol_ig_comp", j)(b, cobj(b, j), b.temperature) - R

    @staticmethod
    def cv_mol_ls_comp_pure(b, p, j):
        # Method for calculating pure component liquid and solid cv from cp
        # For ideal (incompressible) liquids and solids, cv = cp
        pobj = b.params.get_phase(p)
        if pobj.is_liquid_phase():
            return get_method(b, "cp_mol_liq_comp", j)(b, cobj(b, j), b.temperature)
        elif pobj.is_solid_phase():
            return get_method(b, "cp_mol_sol_comp", j)(b, cobj(b, j), b.temperature)

    @staticmethod
    def dens_mass_phase(b, p):
        raise NotImplementedError(_msg(b, "dens_mass_phase"))

    @staticmethod
    def dens_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "dens_mol_phase"))

    @staticmethod
    def energy_internal_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "energy_internal_mol_phase"))

    @staticmethod
    def energy_internal_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "energy_internal_mol_phase_comp"))

    @staticmethod
    def energy_internal_mol_ig_comp_pure(b, j):
        # Method for calculating pure component U from H for ideal gases
        units = b.params.get_metadata().derived_units
        R = pyunits.convert(const.gas_constant, to_units=units.HEAT_CAPACITY_MOLE)

        if cobj(b, j).parent_block().config.include_enthalpy_of_formation:
            # First, need to determine correction between U_form and H_form
            # U_form = H_form - delta_n*R*T
            ele_comp = cobj(b, j).config.elemental_composition
            if ele_comp is None:
                raise ConfigurationError(
                    "{} calculation of internal energy requires elemental "
                    "composition of all species. Please set this using the "
                    "elemental_composition argument in the component "
                    "declaration ({}).".format(b.name, j)
                )

            delta_n = 0
            for e, s in ele_comp.items():
                # Check for any element which is vapor at standard state
                if e in ["He", "Ne", "Ar", "Kr", "Xe", "Ra"]:
                    delta_n += -s
                elif e in ["F", "Cl", "H", "N", "O"]:
                    delta_n += -s / 2  # These are diatomic at standard state

            delta_n += 1  # One mole of gaseous compound is formed
            dU_form = delta_n * R * b.params.temperature_ref
        else:
            dU_form = 0  # No heat of formation to correct

        # For ideal gases, U = H - R(T-T_ref) + dU_form
        return (
            get_method(b, "enth_mol_ig_comp", j)(b, cobj(b, j), b.temperature)
            - R * (b.temperature - b.params.temperature_ref)
            + dU_form
        )

    @staticmethod
    def energy_internal_mol_ls_comp_pure(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_liquid_phase():
            mthd = get_method(b, "enth_mol_liq_comp", j)
        elif pobj.is_solid_phase():
            mthd = get_method(b, "enth_mol_sol_comp", j)

        # Method for calculating pure component U from H for liquids & solids
        units = b.params.get_metadata().derived_units
        R = pyunits.convert(const.gas_constant, to_units=units.HEAT_CAPACITY_MOLE)

        if cobj(b, j).parent_block().config.include_enthalpy_of_formation:
            # First, need to determine correction between U_form and H_form
            # U_form = H_form - delta_n*R*T
            ele_comp = cobj(b, j).config.elemental_composition
            if ele_comp is None:
                raise ConfigurationError(
                    "{} calculation of internal energy requires elemental "
                    "composition of all species. Please set this using the "
                    "elemental_composition argument in the component "
                    "declaration ({}).".format(b.name, j)
                )

            delta_n = 0
            for e, s in ele_comp.items():
                # Check for any element which is vapor at standard state
                if e in ["He", "Ne", "Ar", "Kr", "Xe", "Ra"]:
                    delta_n += -s
                elif e in ["F", "Cl", "H", "N", "O"]:
                    delta_n += -s / 2  # These are diatomic at standard state
            dU_form = delta_n * R * b.params.temperature_ref

            # For ideal (incompressible) liquids and solids, U = H + dU_form
            return mthd(b, cobj(b, j), b.temperature) + dU_form
        else:
            # If not including heat of formation, U = H
            return mthd(b, cobj(b, j), b.temperature)

    @staticmethod
    def enth_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "enth_mol_phase"))

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "enth_mol_phase_comp"))

    @staticmethod
    def entr_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "entr_mol_phase"))

    @staticmethod
    def entr_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "entr_mol_phase_comp"))

    @staticmethod
    def fug_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "fug_phase_comp"))

    @staticmethod
    def fug_phase_comp_eq(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_eq"))

    @staticmethod
    def log_fug_phase_comp_eq(b, p, j, pp):
        raise NotImplementedError(_msg(b, "log_fug_phase_comp_eq"))

    @staticmethod
    def fug_coeff_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "fug_coeff_phase_comp"))

    @staticmethod
    def fug_coeff_phase_comp_eq(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_coeff_phase_comp_eq"))

    @staticmethod
    def fug_phase_comp_Tbub(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Tbub"))

    @staticmethod
    def fug_phase_comp_Tdew(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Tdew"))

    @staticmethod
    def fug_phase_comp_Pbub(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Pbub"))

    @staticmethod
    def fug_phase_comp_Pdew(b, p, j, pp):
        raise NotImplementedError(_msg(b, "fug_phase_comp_Pdew"))

    @staticmethod
    def gibbs_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "gibbs_mol_phase"))

    @staticmethod
    def gibbs_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "gibbs_mol_phase_comp"))

    @staticmethod
    def isentropic_speed_sound_phase(b, p):
        raise NotImplementedError(_msg(b, "isentropic_speed_sound_phase"))

    @staticmethod
    def isothermal_speed_sound_phase(b, p):
        raise NotImplementedError(_msg(b, "isothermal_speed_sound_phase"))

    def pressure_osm_phase(b, p):
        raise NotImplementedError(_msg(b, "pressure_osm_phase"))

    @staticmethod
    def vol_mol_phase(b, p):
        raise NotImplementedError(_msg(b, "vol_mol_phase"))

    @staticmethod
    def vol_mol_phase_comp(b, p, j):
        raise NotImplementedError(_msg(b, "vol_mol_phase_comp"))


def _msg(b, attr):
    return (
        "{} Equation of State module has not implemented a method for {}. "
        "Please contact the EoS developer or use a different module.".format(
            b.name, attr
        )
    )
