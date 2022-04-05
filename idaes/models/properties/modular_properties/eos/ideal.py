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
Methods for ideal equations of state.

Currently only supports liquid and vapor phases
"""
from pyomo.environ import Expression, log

from idaes.core import Apparent
from idaes.core.util.exceptions import ConfigurationError, PropertyNotSupportedError
from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_component_object as cobj,
)
from .eos_base import EoSBase
from idaes.models.properties.modular_properties.phase_equil.henry import (
    henry_pressure,
    log_henry_pressure,
)


# TODO: Add support for ideal solids
class Ideal(EoSBase):
    # Add attribute indicating support for electrolyte systems
    electrolyte_support = True

    @staticmethod
    def common(b, pobj):
        # No common components required for ideal property calculations
        pass

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def build_parameters(b):
        # No EoS specific parameters required
        pass

    @staticmethod
    def act_phase_comp(b, p, j):
        return b.mole_frac_phase_comp[p, j]

    @staticmethod
    def act_phase_comp_true(b, p, j):
        return b.mole_frac_phase_comp_true[p, j]

    @staticmethod
    def act_phase_comp_appr(b, p, j):
        return b.mole_frac_phase_comp_apparent[p, j]

    @staticmethod
    def act_coeff_phase_comp(b, p, j):
        return 1

    @staticmethod
    def act_coeff_phase_comp_true(b, p, j):
        return 1

    @staticmethod
    def act_coeff_phase_comp_appr(b, p, j):
        return 1

    @staticmethod
    def compress_fact_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return 1
        else:
            return 0

    @staticmethod
    def cp_mol_phase(b, p):
        return sum(
            b.get_mole_frac(p)[p, j] * b.cp_mol_phase_comp[p, j]
            for j in b.components_in_phase(p)
        )

    @staticmethod
    def cp_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return get_method(b, "cp_mol_ig_comp", j)(b, cobj(b, j), b.temperature)
        elif pobj.is_liquid_phase():
            return get_method(b, "cp_mol_liq_comp", j)(b, cobj(b, j), b.temperature)
        elif pobj.is_solid_phase():
            return get_method(b, "cp_mol_sol_comp", j)(b, cobj(b, j), b.temperature)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def cv_mol_phase(b, p):
        return sum(
            b.get_mole_frac(p)[p, j] * b.cv_mol_phase_comp[p, j]
            for j in b.components_in_phase(p)
        )

    @staticmethod
    def cv_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return EoSBase.cv_mol_ig_comp_pure(b, j)
        elif pobj.is_liquid_phase() or pobj.is_solid_phase():
            return EoSBase.cv_mol_ls_comp_pure(b, p, j)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def dens_mass_phase(b, p):
        return b.dens_mol_phase[p] * b.mw_phase[p]

    @staticmethod
    def dens_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return b.pressure / (Ideal.gas_constant(b) * b.temperature)
        else:
            return 1 / b.vol_mol_phase[p]

    @staticmethod
    def energy_internal_mol_phase(b, p):
        return sum(
            b.get_mole_frac(p)[p, j] * b.energy_internal_mol_phase_comp[p, j]
            for j in b.components_in_phase(p)
        )

    @staticmethod
    def energy_internal_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return EoSBase.energy_internal_mol_ig_comp_pure(b, j)
        elif pobj.is_liquid_phase() or pobj.is_solid_phase():
            return EoSBase.energy_internal_mol_ls_comp_pure(b, p, j)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def enth_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return sum(
                b.get_mole_frac(p)[p, j] * b.enth_mol_phase_comp[p, j]
                for j in b.components_in_phase(p)
            )
        elif pobj.is_liquid_phase():
            return (
                sum(
                    b.get_mole_frac(p)[p, j]
                    * get_method(b, "enth_mol_liq_comp", j)(
                        b, cobj(b, j), b.temperature
                    )
                    for j in b.components_in_phase(p)
                )
                + (b.pressure - b.params.pressure_ref) / b.dens_mol_phase[p]
            )
        elif pobj.is_solid_phase():
            return (
                sum(
                    b.get_mole_frac(p)[p, j]
                    * get_method(b, "enth_mol_sol_comp", j)(
                        b, cobj(b, j), b.temperature
                    )
                    for j in b.components_in_phase(p)
                )
                + (b.pressure - b.params.pressure_ref) / b.dens_mol_phase[p]
            )
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return get_method(b, "enth_mol_ig_comp", j)(b, cobj(b, j), b.temperature)
        elif pobj.is_liquid_phase():
            return (
                get_method(b, "enth_mol_liq_comp", j)(b, cobj(b, j), b.temperature)
                + (b.pressure - b.params.pressure_ref) / b.dens_mol_phase[p]
            )
        elif pobj.is_solid_phase():
            return (
                get_method(b, "enth_mol_sol_comp", j)(b, cobj(b, j), b.temperature)
                + (b.pressure - b.params.pressure_ref) / b.dens_mol_phase[p]
            )
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def entr_mol_phase(b, p):
        return sum(
            b.get_mole_frac(p)[p, j] * b.entr_mol_phase_comp[p, j]
            for j in b.components_in_phase(p)
        )

    @staticmethod
    def entr_mol_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return get_method(b, "entr_mol_ig_comp", j)(
                b, cobj(b, j), b.temperature
            ) - Ideal.gas_constant(b) * log(
                b.get_mole_frac(p)[p, j] * b.pressure / b.params.pressure_ref
            )
        elif pobj.is_liquid_phase():
            # Assume no pressure/volume dependecy of entropy for ideal liquids
            return get_method(b, "entr_mol_liq_comp", j)(b, cobj(b, j), b.temperature)
        elif pobj.is_solid_phase():
            # Assume no pressure/volume dependecy of entropy for ideal solids
            return get_method(b, "entr_mol_sol_comp", j)(b, cobj(b, j), b.temperature)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def fug_phase_comp(b, p, j):
        return _fug_phase_comp(b, p, j, b.temperature)

    @staticmethod
    def fug_phase_comp_eq(b, p, j, pp):
        return _fug_phase_comp(b, p, j, b._teq[pp])

    @staticmethod
    def log_fug_phase_comp_eq(b, p, j, pp):
        pobj = b.params.get_phase(p)

        if pobj.is_vapor_phase():
            return log(b.get_mole_frac(p)[p, j]) + log(b.pressure)
        elif pobj.is_liquid_phase():
            if (
                cobj(b, j).config.henry_component is not None
                and p in cobj(b, j).config.henry_component
            ):
                # Use Henry's Law
                return log_henry_pressure(b, p, j, b.temperature)
            elif cobj(b, j).config.has_vapor_pressure:
                # Use Raoult's Law
                return log(b.get_mole_frac(p)[p, j]) + log(
                    get_method(b, "pressure_sat_comp", j)(b, cobj(b, j), b.temperature)
                )
            else:
                return Expression.Skip
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def fug_coeff_phase_comp(b, p, j):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return 1

    @staticmethod
    def fug_coeff_phase_comp_eq(b, p, j, pp):
        pobj = b.params.get_phase(p)
        if not (pobj.is_vapor_phase() or pobj.is_liquid_phase()):
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
        return 1

    # TODO Bubble/dew methods broken for Henry's law---see issue #718
    # Need to update to call either (log)_henry_pressure() or
    # or henry_equilibrium_ratio(). Using the former is a problem, because they
    # don't return Henry pressures at bubble/dew mole fractions
    @staticmethod
    def log_fug_phase_comp_Tbub(b, p, j, pp):
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        if pobj.is_vapor_phase():
            return log(b._mole_frac_tbub[pp[0], pp[1], j]) + log(b.pressure)
        elif pobj.is_liquid_phase():
            if (
                cobj.config.henry_component is not None
                and p in cobj.config.henry_component
            ):
                return log(b.mole_frac_comp[j]) + log(
                    get_method(b, "henry_component", j, p)(
                        b, p, j, b.temperature_bubble[pp]
                    )
                )
            else:
                return log(b.mole_frac_comp[j]) + log(
                    get_method(b, "pressure_sat_comp", j)(
                        b, cobj, b.temperature_bubble[pp]
                    )
                )
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def log_fug_phase_comp_Tdew(b, p, j, pp):
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        if pobj.is_vapor_phase():
            return log(b.mole_frac_comp[j]) + log(b.pressure)
        elif pobj.is_liquid_phase():
            if (
                cobj.config.henry_component is not None
                and p in cobj.config.henry_component
            ):
                return log(b._mole_frac_tdew[pp[0], pp[1], j]) + log(
                    get_method(b, "henry_component", j, p)(
                        b, p, j, b.temperature_dew[pp]
                    )
                )
            else:
                return log(b._mole_frac_tdew[pp[0], pp[1], j]) + log(
                    get_method(b, "pressure_sat_comp", j)(
                        b, cobj, b.temperature_dew[pp]
                    )
                )
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def log_fug_phase_comp_Pbub(b, p, j, pp):
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        if pobj.is_vapor_phase():
            return log(b._mole_frac_pbub[pp[0], pp[1], j]) + log(b.pressure_bubble[pp])
        elif pobj.is_liquid_phase():
            if (
                cobj.config.henry_component is not None
                and p in cobj.config.henry_component
            ):
                return log(b.mole_frac_comp[j]) + log(b.henry[p, j])
            else:
                return log(b.mole_frac_comp[j]) + log(b.pressure_sat_comp[j])
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def log_fug_phase_comp_Pdew(b, p, j, pp):
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        if pobj.is_vapor_phase():
            return log(b.mole_frac_comp[j]) + log(b.pressure_dew[pp])
        elif pobj.is_liquid_phase():
            if (
                cobj.config.henry_component is not None
                and p in cobj.config.henry_component
            ):
                return log(b._mole_frac_pdew[pp[0], pp[1], j]) + log(b.henry[p, j])
            else:
                return log(b._mole_frac_pdew[pp[0], pp[1], j]) + log(
                    b.pressure_sat_comp[j]
                )
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def gibbs_mol_phase(b, p):
        return sum(
            b.get_mole_frac(p)[p, j] * b.gibbs_mol_phase_comp[p, j]
            for j in b.components_in_phase(p)
        )

    @staticmethod
    def gibbs_mol_phase_comp(b, p, j):
        return b.enth_mol_phase_comp[p, j] - b.entr_mol_phase_comp[p, j] * b.temperature

    @staticmethod
    def pressure_osm_phase(b, p):
        try:
            solvent_set = b.params.solvent_set
        except AttributeError:
            raise ConfigurationError(
                f"{b.name} called for pressure_osm, but no solvents were "
                f"defined. Osmotic pressure requires at least one component "
                f"to be declared as a solvent."
            )
        C = 0
        for j in b.component_list:
            if (p, j) in b.phase_component_set and j not in solvent_set:
                c_obj = b.params.get_component(j)
                if isinstance(c_obj, Apparent):
                    i = sum(c_obj.config.dissociation_species.values())
                else:
                    i = 1
                C += i * b.conc_mol_phase_comp[p, j]
        return Ideal.gas_constant(b) * b.temperature * C

    @staticmethod
    def vol_mol_phase(b, p):
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return Ideal.gas_constant(b) * b.temperature / b.pressure
        elif pobj.is_liquid_phase():
            ptype = "liq"
        elif pobj.is_solid_phase():
            ptype = "sol"
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

        v_expr = 0
        for j in b.components_in_phase(p):
            # First try to get a method for vol_mol
            v_comp = Ideal.get_vol_mol_pure(b, ptype, j, b.temperature)
            v_expr += b.get_mole_frac(p)[p, j] * v_comp

        return v_expr


def _invalid_phase_msg(name, phase):
    return (
        "{} received unrecognised phase name {}. Ideal property "
        "libray only supports Vap and Liq phases.".format(name, phase)
    )


def _fug_phase_comp(b, p, j, T):
    pobj = b.params.get_phase(p)

    if pobj.is_vapor_phase():
        return b.get_mole_frac(p)[p, j] * b.pressure
    elif pobj.is_liquid_phase():
        if (
            cobj(b, j).config.henry_component is not None
            and p in cobj(b, j).config.henry_component
        ):
            # Use Henry's Law
            return henry_pressure(b, p, j, T)
        elif cobj(b, j).config.has_vapor_pressure:
            # Use Raoult's Law
            return b.get_mole_frac(p)[p, j] * get_method(b, "pressure_sat_comp", j)(
                b, cobj(b, j), T
            )
        else:
            return Expression.Skip
    else:
        raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))
