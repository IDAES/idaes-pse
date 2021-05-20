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
Base class for EoS modules.

Raises NotImplementedErrors for all expected methods in case developer misses
some. EoS developers should overload all these methods.
"""
from pyomo.environ import units as pyunits
from idaes.core.util.constants import Constants as const
from idaes.generic_models.properties.core.generic.utility import (
    get_method, get_component_object as cobj)
from idaes.core.util.exceptions import (
    PropertyNotSupportedError, ConfigurationError)


class EoSBase():

    @staticmethod
    def gas_constant(b):
        # Utility method to convert gas constant to base units
        base_units = b.params.get_metadata().default_units

        r_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["temperature"]**-1 *
                   base_units["amount"]**-1 *
                   base_units["time"]**-2)

        return pyunits.convert(const.gas_constant, to_units=r_units)

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
    def cv_mol_ig_comp_pure(b, j):
        # Method for calculating pure component ideal gas cv from cp
        # For ideal gases, cv = cp - R
        units = b.params.get_metadata().derived_units
        R = pyunits.convert(const.gas_constant,
                            to_units=units["heat_capacity_mole"])
        return (get_method(b, "cp_mol_ig_comp", j)(
            b, cobj(b, j), b.temperature) - R)

    @staticmethod
    def cv_mol_ls_comp_pure(b, j):
        # Method for calculating pure component liquid and solid cv from cp
        # For ideal (incompressible) liquids and solids, cv = cp
        return get_method(b, "cp_mol_liq_comp", j)(
            b, cobj(b, j), b.temperature)

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
        R = pyunits.convert(const.gas_constant,
                            to_units=units["heat_capacity_mole"])

        if cobj(b, j).parent_block().config.include_enthalpy_of_formation:
            # First, need to determine correction between U_form and H_form
            # U_form = H_form - delta_n*R*T
            ele_comp = cobj(b, j).config.elemental_composition
            if ele_comp is None:
                raise ConfigurationError(
                    "{} calculation of internal energy requires elemental "
                    "composition of all species. Please set this using the "
                    "elemental_composition argument in the component "
                    "declaration ({}).".format(b.name, j))

            delta_n = 0
            for e, s in ele_comp.items():
                # Check for any element which is vapor at standard state
                if e in ["He", "Ne", "Ar", "Kr", "Xe", "Ra"]:
                    delta_n += -s
                elif e in ["F", "Cl", "H", "N", "O"]:
                    delta_n += -s/2  # These are diatomic at standard state

            delta_n += 1  # One mole of gaseous compound is formed
            dU_form = delta_n*R*b.params.temperature_ref
        else:
            dU_form = 0  # No heat of formation to correct

        # For ideal gases, U = H - R(T-T_ref) + dU_form
        return (get_method(b, "enth_mol_ig_comp", j)(
            b, cobj(b, j), b.temperature) -
            R*(b.temperature-b.params.temperature_ref) +
            dU_form)

    @staticmethod
    def energy_internal_mol_ls_comp_pure(b, j):
        # Method for calculating pure component U from H for liquids & solids
        units = b.params.get_metadata().derived_units
        R = pyunits.convert(const.gas_constant,
                            to_units=units["heat_capacity_mole"])

        if cobj(b, j).parent_block().config.include_enthalpy_of_formation:
            # First, need to determine correction between U_form and H_form
            # U_form = H_form - delta_n*R*T
            ele_comp = cobj(b, j).config.elemental_composition
            if ele_comp is None:
                raise ConfigurationError(
                    "{} calculation of internal energy requires elemental "
                    "composition of all species. Please set this using the "
                    "elemental_composition argument in the component "
                    "declaration ({}).".format(b.name, j))

            delta_n = 0
            for e, s in ele_comp.items():
                # Check for any element which is vapor at standard state
                if e in ["He", "Ne", "Ar", "Kr", "Xe", "Ra"]:
                    delta_n += -s
                elif e in ["F", "Cl", "H", "N", "O"]:
                    delta_n += -s/2  # These are diatomic at standard state
            dU_form = delta_n*R*b.params.temperature_ref

            # For ideal (incompressible) liquids and solids, U = H + dU_form
            return (get_method(b, "enth_mol_liq_comp", j)(
                b, cobj(b, j), b.temperature) +
                dU_form)
        else:
            # If not including heat of formation, U = H
            return get_method(b, "enth_mol_liq_comp", j)(
                b, cobj(b, j), b.temperature)

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


def _msg(b, attr):
    return ("{} Equation of State module has not implemented a method for {}. "
            "Please contact the EoS developer or use a different module."
            .format(b.name, attr))
