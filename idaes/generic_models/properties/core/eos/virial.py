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
Methods for virial equations of state (veos) for gases

Currently supports B-truncated volume explicit virial equation of state with validity at low to moderate pressures where the compresiibility factor is approximately a linear function of  pressure. It is most accurate for non-polar species.

Refernces**
[1] J. M Prausnitz, R, N Lichtenthaler, and E. G de Azvedo,
    Molecular thermodynamics of Fluid-phase Equilibruim 3rd ed. Prentice-Hall, Nj, 1998
**Any standard chemical engineering thermodynamic textbook

Author: Akula Paul

"""

from copy import deepcopy
from pyomo.environ import (Expression,
                           log, exp,
                           value, NonNegativeReals,
                           Var, units as pyunits)
from pyomo.common.config import ConfigBlock, ConfigValue, Bool

from idaes.core.util.exceptions import PropertyNotSupportedError
from idaes.generic_models.properties.core.generic.utility import (
    get_method, get_component_object as gcobj)
from idaes.generic_models.properties.core.eos.eos_base import EoSBase
from idaes.core.util.misc import set_param_from_config
from idaes.generic_models.properties.core.pure.RPP4 import entr_mol_ig_comp


VirialConfig = ConfigBlock()

# virial eos options
VirialConfig.declare("use_pseudocritical_rules", ConfigValue(
    default=True,
    domain=Bool,
    description='''Flag indicating whether pseudocritical rules should be used
    when applicable''',
    doc='''Flag indicating whether to compute mixture property based on
    composition and critical properties of pure components:
**default** - True.
**Valid values:** {
**True** - use pseudocritical_rule: mixing rule based on composition and
critical properties of pure components.
**False** - use mixing rule based on composition and corresponding
pure component properties.}'''))


class enth_mol_ig_comp():
    # This fixes the bug in RPP4.Since its just a specific method,recreating it
    # here
    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
            cp_mol_ig_comp.build_parameters(cobj)

        units = cobj.parent_block().get_metadata().derived_units
        cobj.enth_mol_form_vap_comp_ref = Var(
            doc="Vapor phase molar heat of formation @ Tref",
            units=units["energy_mole"])
        if cobj.parent_block().config.include_enthalpy_of_formation:
            set_param_from_config(cobj, param="enth_mol_form_vap_comp_ref")
        else:
            cobj.enth_mol_form_vap_comp_ref.set_value(0.0)

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific enthalpy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tref = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        units = b.params.get_metadata().derived_units

        h = (pyunits.convert(
            (cobj.cp_mol_ig_comp_coeff_D / 4) * (T**4 - Tref**4) +
            (cobj.cp_mol_ig_comp_coeff_C / 3) * (T**3 - Tref**3) +
            (cobj.cp_mol_ig_comp_coeff_B / 2) * (T**2 - Tref**2) +
            cobj.cp_mol_ig_comp_coeff_A * (T - Tref), units["energy_mole"]) +
            cobj.enth_mol_form_vap_comp_ref)

        return h


class entr_mol_ig_comp_G_H_ref(entr_mol_ig_comp):
    # This updates the build_parameter method to avoid setting S_ref
    # So that S_sef is calculated from H_ref and G_ref.
    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_ig_comp_coeff_A"):
            cp_mol_ig_comp.build_parameters(cobj)

        units = cobj.parent_block().get_metadata().derived_units

        cobj.entr_mol_form_vap_comp_ref = Var(
            doc="Vapor phase molar entropy of formation @ Tref",
            units=units["entropy_mole"])


class Virial(EoSBase):

    @staticmethod
    def common(b, pobj):

        # Add expressions for pseudocritical parameters
        def func_p_omega(m):
            return rule_pseudocritical_omega(m)

        b.add_component('pseudo_omega',
                        Expression(rule=func_p_omega,
                                   doc='Pseudocritical omega for gas mixture'))

        def func_p_Tc(m):
            return rule_pseudocritical_temperature(m)

        b.add_component('pseudo_Tc',
                        Expression(rule=func_p_Tc,
                                   doc='Pseudocritical temperature for gas mixture'))

        def func_p_Pc(m):
            return rule_pseudocritical_pressure(m)

        b.add_component('pseudo_Pc',
                        Expression(rule=func_p_Pc,
                                   doc='Pseudocritical Pressure for gas mixture'))

        # add expression for component critical compressibility factor
        for comp in b.params.component_list:
            cobj = b.params.get_component(comp)
            cobj.compress_fact_crit = Expression(expr=cobj.pressure_crit *
                                                 cobj.volume_crit /
                                                 cobj.temperature_crit /
                                                 Virial.gas_constant(b))

        # Add combining rules proposed by Prausnitz et al. [1]
        def rule_omega_ij(m, i, j):
            return 0.5 * (m.params.get_component(i).omega +
                          m.params.get_component(j).omega)
        b.add_component('omega_ij',
                        Expression(b.component_list,
                                   b.component_list,
                                   rule=rule_omega_ij,
                                   doc='Omega combining rule for i-j molecular pair'))

        def rule_Zc_ij(m, i, j):
            return 0.5 * (m.params.get_component(i).compress_fact_crit +
                          m.params.get_component(j).compress_fact_crit)
        b.add_component('Zc_ij',
                        Expression(b.component_list,
                                   b.component_list,
                                   rule=rule_Zc_ij,
                                   doc='Critical compressibility factor '
                                       'combining rule for i-j molecular pair'))

        def rule_Vc_ij(m, i, j):
            return (0.5 * (m.params.get_component(i).volume_crit**(1 / 3) +
                           m.params.get_component(j).volume_crit**(1 / 3)))**3
        b.add_component('Vc_ij',
                        Expression(b.component_list,
                                   b.component_list,
                                   rule=rule_Vc_ij,
                                   doc='Critical volume combining rule for '
                                       'i-j molecular pair'))

        def rule_Tc_ij(m, i, j):
            return ((1 - m.params.kappa[i, j]) *
                    (m.params.get_component(i).temperature_crit *
                     m.params.get_component(j).temperature_crit)**0.5)
        b.add_component('Tc_ij',
                        Expression(b.component_list,
                                   b.component_list,
                                   rule=rule_Tc_ij,
                                   doc='Critical temperature combining rule '
                                       'for i-j molecular pair'))

        def rule_Pc_ij(m, i, j):
            return (m.Zc_ij[i, j] * m.Tc_ij[i, j] / m.Vc_ij[i, j] *
                    Virial.gas_constant(m))
        b.add_component('Pc_ij',
                        Expression(b.component_list,
                                   b.component_list,
                                   rule=rule_Pc_ij,
                                   doc='Critical pressure combining rule '
                                       'for i-j molecular pair'))

        def rule_B_ij(m, i, j):
            return (m.Tc_ij[i, j] / m.Pc_ij[i, j] * Virial.gas_constant(m) *
                    (f_BO(m.temperature / m.Tc_ij[i, j]) +
                     f_B1(m.temperature / m.Tc_ij[i, j]) * m.omega_ij[i, j]))
        b.add_component('B_ij',
                        Expression(b.component_list,
                                   b.component_list,
                                   rule=rule_B_ij,
                                   doc='Second virial coefficient for '
                                       'i-j molecular pair'))

    @staticmethod
    def calculate_scaling_factors(b, pobj):
        pass

    @staticmethod
    def build_parameters(b):
        param_block = b.parent_block()
        # virial eos supports only vapor phase
        if not b.is_vapor_phase():
            raise PropertyNotSupportedError(f"{b.parent_block().name} "
                                            "received unrecognized phase name"
                                            f"{b}. Virial equation of"
                                            "state supports only Vap phase.")

        # get user specified eos options and save on parameter block
        try:
            use_pseudocritical_rules = b.config.equation_of_state_options[
                "use_pseudocritical_rules"]
            b._use_pseudocritical_rules = use_pseudocritical_rules
            setattr(b.parent_block(), "virial_eos_options",
                    deepcopy(VirialConfig))
            ConfigBlock = getattr(b.parent_block(), "virial_eos_options")
            ConfigBlock.set_value(b.config.equation_of_state_options)
        except TypeError:
            use_pseudocritical_rules = True
            b._use_pseudocritical_rules = use_pseudocritical_rules
            setattr(b.parent_block(), "virial_eos_options",
                    deepcopy(VirialConfig))
            ConfigBlock = getattr(b.parent_block(), "virial_eos_options")
            ConfigBlock.set_value(b.config.equation_of_state_options)

        # Get critical volume of components
        units = b.parent_block().get_metadata().derived_units
        for comp in b.parent_block().component_list:
            cobj = b.parent_block().get_component(comp)
            if not hasattr(cobj, 'volume_crit'):
                setattr(cobj, "volume_crit",
                        Var(doc="Pure component critical Volume",
                            units=units["volume"] / pyunits.mol))
                set_param_from_config(cobj, param="volume_crit")

        # kappa binary interaction parameters for i-j molecular pair
        # for  i=j(like molecules), kappa = 0
        # for i not equal to  j, kappa is a small positive number.
        # set to zero if data is unavailable

        try:
            kappa_data = param_block.config.parameter_data["kappa"]
            param_block.add_component(
                'kappa',
                Var(param_block.component_list,
                    param_block.component_list,
                    within=NonNegativeReals,
                    initialize=kappa_data,
                    doc='Binary interaction parameters',
                    units=pyunits.dimensionless))
        except KeyError:
            param_block.add_component(
                'kappa',
                Var(param_block.component_list,
                    param_block.component_list,
                    within=NonNegativeReals,
                    initialize=0.0,
                    doc='Binary interaction parameters',
                    units=pyunits.dimensionless))

        # Since standard property changes of formation are readily available for
        # Enthalpy(H_ref) and Gibbs free energy (G_ref),
        # This builds G_ref and set its value from config block, and  then
        # computes S_ref = (H_ref-G_ref)/T_ref
        for comp in b.parent_block().component_list:
            cobj = b.parent_block().get_component(comp)
            setattr(cobj, "gibbs_mol_form_vap_comp_ref",
                    Var(doc="Vapor phase molar Gibbs free energy"
                        "formation @ Tref",
                        units=units["energy_mole"]))

            if b.parent_block().config.include_enthalpy_of_formation:
                set_param_from_config(
                    cobj, param="gibbs_mol_form_vap_comp_ref")
                S_ref = pyunits.convert(((cobj.enth_mol_form_vap_comp_ref -
                                          cobj.gibbs_mol_form_vap_comp_ref) /
                                         b.parent_block().temperature_ref),
                                        units['entropy_mole'])
                cobj.entr_mol_form_vap_comp_ref.set_value(value(S_ref))
            else:
                cobj.entr_mol_form_vap_comp_ref.set_value(0.0)
                cobj.gibbs_mol_form_vap_comp_ref.set_value(0.0)

    # -------------------------------------------------------------------------
    # METHODS FOR PURE SPECIES

    @staticmethod
    def compress_fact_vap_comp_pure(b, comp):
        # Method for calculating pure component  vapor compressiblity

        cobj = b.params.get_component(comp)
        Tr = b.temperature / cobj.temperature_crit
        Pr = b.pressure / cobj.pressure_crit
        omega = cobj.omega
        return 1 + Pr / Tr * (f_BO(Tr) + omega * f_B1(Tr))

    @staticmethod
    def vol_mol_vap_comp_pure(b, comp):
        # Method for calculating pure component  vapor molar volume
        return (Virial.compress_fact_vap_comp_pure(b, comp) *
                Virial.gas_constant(b) * b.temperature / b.pressure)

    @staticmethod
    def enth_mol_vap_comp_pure(b, comp):
        # Method for calculating pure component  vapor molar enthalpy
        cobj = b.params.get_component(comp)
        Tc = cobj.temperature_crit
        Tr = b.temperature / Tc
        Pr = b.pressure / cobj.pressure_crit
        omega = cobj.omega

        enth_ideal_gas_j = get_method(b, "enth_mol_ig_comp", comp)(
            b, cobj, b.temperature)

        enth_departure_j = Virial.gas_constant(b) * Tc * Pr * (
            f_BO(Tr) - Tr * f_dBO_dTr(Tr) +
            omega * (f_B1(Tr) - Tr * f_dB1_dTr(Tr)))

        return enth_ideal_gas_j + enth_departure_j

    @staticmethod
    def entr_mol_vap_comp_pure(b, comp):
        # Method for calculating pure component  vapor molar entropy
        cobj = b.params.get_component(comp)
        Tr = b.temperature / cobj.temperature_crit
        Pr = b.pressure / cobj.pressure_crit
        omega = cobj.omega

        # assumes entr_mol_ig_comp method has only the cp integral part
        entr_ideal_gas_comp = get_method(b, "entr_mol_ig_comp", comp)(
            b, cobj, b.temperature) - Virial.gas_constant(b) * log(
            b.pressure / b.params.pressure_ref)

        entr_departure_comp = -Virial.gas_constant(b) * Pr * (
            f_dBO_dTr(Tr) + omega * f_dB1_dTr(Tr))

        return entr_ideal_gas_comp + entr_departure_comp

    @staticmethod
    def gibbs_mol_vap_comp_pure(b, comp):
        # Method for calculating pure component  vapor molar Gibbs free energy
        return (Virial.enth_mol_vap_comp_pure(b, comp) -
                Virial.entr_mol_vap_comp_pure(b, comp) * b.temperature)

    @staticmethod
    def energy_internal_mol_vap_comp_pure(b, comp):
        # Method for calculating pure component  vapor molar internal energy
        return (Virial.enth_mol_vap_comp_pure(b, comp) -
                Virial.vol_mol_vap_comp_pure(b, comp) * b.pressure)

    @staticmethod
    def log_fug_coeff_vap_comp_pure(b, comp):
        # Method for natural log of pure component  vapor fugacity coefficient
        cobj = b.params.get_component(comp)
        Tr = b.temperature / cobj.temperature_crit
        Pr = b.pressure / cobj.pressure_crit
        omega = cobj.omega
        return Pr / Tr * (f_BO(Tr) + omega * f_B1(Tr))

    @staticmethod
    def fug_coeff_vap_comp_pure(b, comp):
        # Method for pure component  vapor fugacity coefficient
        return exp(Virial.log_fug_coeff_vap_comp_pure(b, comp))

    @staticmethod
    def fug_vap_comp_pure(b, comp):
        # Method for pure component  vapor fugacity
        return Virial.fug_coeff_vap_comp_pure(b, comp) * b.pressure

    @staticmethod
    def log_fug_coeff_vap_sat_comp(b, comp):
        # Method for natural log of pure component  vapor fugacity coefficient
        # at saturation pressure
        cobj = b.params.get_component(comp)
        Tr = b.temperature / cobj.temperature_crit
        Pr_sat = (get_method(b, "pressure_sat_comp", comp)(
                  b, cobj, b.temperature)) / cobj.pressure_crit
        omega = cobj.omega
        return Pr_sat / Tr * (f_BO(Tr) + omega * f_B1(Tr))

    @staticmethod
    def fug_coeff_vap_sat_comp(b, comp):
        # Pure component  vapor fugacity coefficient at saturation pressure
        return exp(Virial.log_fug_coeff_vap_sat_comp(b, comp))

    @staticmethod
    def fug_vap_sat_comp(b, comp):
        # Method for pure component  vapor fugacity at saturation pressure
        cobj = b.params.get_component(comp)
        return (get_method(b, "pressure_sat_comp", comp)(b, cobj, b.temperature)
                ) * Virial.fug_coeff_vap_sat_comp(b, comp)

    # -------------------------------------------------------------------------
    # METHODS FOR GAS MIXTURES

    @staticmethod
    def second_virial_coeff_vap(b, p):
        # Method for calculating gas mixture second virial coefficient
        # either via pseoudocritical rules or combining rules
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            if pobj._use_pseudocritical_rules:
                Tr = b.temperature / b.pseudo_Tc
                return (Virial.gas_constant(b) * b.pseudo_Tc /
                        b.pseudo_Pc * (f_BO(Tr) + b.pseudo_omega * f_B1(Tr)))
            else:
                return rule_second_virial_coeff(b, p)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def compress_fact_phase(b, p):
        # Method for calculating gas mixture compressibility factor
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return 1 + (Virial.second_virial_coeff_vap(b, p) *
                        b.pressure / b.temperature / Virial.gas_constant(b))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def vol_mol_phase(b, p):
        # Method for calculating gas mixture molar volume
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return (Virial.gas_constant(b) * b.temperature *
                    b.compress_fact_phase[p] / b.pressure)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def dens_mass_phase(b, p):
        # Method for calculating gas mixture mass density
        return b.dens_mol_phase[p] * b.mw_phase[p]

    @staticmethod
    def dens_mol_phase(b, p):
        # Method for calculating gas mixture molar density
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return b.pressure / (Virial.gas_constant(b) * b.temperature
                                 * b.compress_fact_phase[p])
        else:
            return 1 / b.vol_mol_phase[p]

    @staticmethod
    def enth_mol_phase(b, p):
        # Method for calculating gas mixture molar enthalpy
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            enth_ideal_gas = sum(b.mole_frac_phase_comp[p, j] *
                                 get_method(b, "enth_mol_ig_comp", j)(
                                 b, gcobj(b, j), b.temperature)
                                 for j in b.components_in_phase(p))
            if pobj._use_pseudocritical_rules:
                Tr = b.temperature / b.pseudo_Tc
                Pr = b.pressure / b.pseudo_Pc
                enth_departure = Virial.gas_constant(b) * b.pseudo_Tc * Pr * (
                    f_BO(Tr) - Tr * f_dBO_dTr(Tr) +
                    b.pseudo_omega * (f_B1(Tr) - Tr * f_dB1_dTr(Tr)))
            else:
                enth_departure = (Virial.second_virial_coeff_mix(b, p) -
                                  b.temperature * Virial.dBdT(b, p)) * b.pressure
            return enth_ideal_gas + enth_departure
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def entr_mol_phase(b, p):
        # Method for calculating gas mixture molar entropy
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            # Note: this assumes that "entr_mol_ig_comp" method
            # does not include the ideal-gas entropy change of mixing  and
            # it is therefore included here
            entr_ideal_gas = (sum(b.mole_frac_phase_comp[p, j] *
                                  get_method(b, "entr_mol_ig_comp", j)(
                                  b, gcobj(b, j), b.temperature)
                                  for j in b.components_in_phase(p)) -
                              Virial.gas_constant(b) * (
                              log(b.pressure / b.params.pressure_ref) +
                              sum(b.mole_frac_phase_comp[p, k] *
                                  log(b.mole_frac_phase_comp[p, k])
                                  for k in b.components_in_phase(p))))
            if pobj._use_pseudocritical_rules:
                Tr = b.temperature / b.pseudo_Tc
                Pr = b.pressure / b.pseudo_Pc
                entr_departure = ((f_dBO_dTr(Tr) +
                                   b.pseudo_omega * f_dB1_dTr(Tr)) *
                                  -Virial.gas_constant(b) * Pr)
            else:
                entr_departure = -(b.pressure * Virial.dBdT(b, p))
            return entr_ideal_gas + entr_departure
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def energy_internal_mol_phase(b, p):
        # Method for calculating gas mixture molar internal energy
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return b.enth_mol_phase[p] - b.pressure * b.vol_mol_phase[p]
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def gibbs_mol_phase(b, p):
        # Method for calculating gas mixture molar internal energy
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return (b.enth_mol_phase[p] - b.entr_mol_phase[p] * b.temperature)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def log_fug_coeff_phase_comp(b, p, k):
        # Method for natural logarithm of fugacity coefficient of species k
        # in the gas mixture
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return b.pressure / Virial.gas_constant(b) / b.temperature * (
                2 * sum(b.mole_frac_phase_comp[p, i] * b.B_ij[i, k]
                        for i in b.components_in_phase(p)) -
                Virial.second_virial_coeff_mix(b, p))

        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def fug_coeff_phase_comp(b, p, k):
        # Method for fugacity coeffcient of species k in the gas mixture
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return exp(Virial.log_fug_coeff_phase_comp(b, p, k))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def fug_phase_comp(b, p, k):
        # Method for fugacity  of species k in the gas mixture
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return (Virial.fug_coeff_phase_comp(b, p, k) *
                    b.mole_frac_phase_comp[p, k] * b.pressure)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def cp_mol_phase(b, p):
        # Method for calculating gas mixture molar cp
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            cp_ig = sum(b.mole_frac_phase_comp[p, j] *
                        get_method(b, "cp_mol_ig_comp", j)(b, gcobj(b, j),
                                                           b.temperature)
                        for j in b.components_in_phase(p))
            if pobj._use_pseudocritical_rules:
                Tr = b.temperature / b.pseudo_Tc
                cp_residual = -1 * (Virial.gas_constant(b) * b.pressure *
                                    b.temperature / b.pseudo_Tc / b.pseudo_Pc *
                                    (f_dBO_dTr2(Tr) + b.pseudo_omega *
                                     f_dB1_dTr2(Tr)))
            else:
                cp_residual = -b.pressure * b.temperature * Virial.dBdT2(b, p)
            return cp_ig + cp_residual
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def cv_mol_phase(b, p):
        # Method for calculating gas mixture molar cv
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            if pobj._use_pseudocritical_rules:
                Pr = b.pressure / b.pseudo_Pc
                Tr = b.temperature / b.pseudo_Tc
                cv = (b.cp_mol_phase[p] -
                      Virial.gas_constant(b) *
                      (1 + Pr * (f_dBO_dTr(Tr) + b.pseudo_omega *
                                 f_dB1_dTr(Tr)))**2)
            else:
                cv = (b.cp_mol_phase[p] - 1 / Virial.gas_constant(b) *
                      (Virial.gas_constant(b) + b.pressure *
                       Virial.dBdT(b, p))**2)
            return cv
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    # -------------------------------------------------------------------------
    # Partial Molar Properties(PMP) : for consisitency with the
    # derivations of PMP, pseudocritial rules are not used
    @staticmethod
    def gibbs_mol_phase_comp(b, p, j):
        # Method for partial molar Gibbs energy for species j in gas mixture
        pobj = b.params.get_phase(p)
        cobj = b.params.get_component(j)
        enth_ig_comp = get_method(b, "enth_mol_ig_comp", j)(
            b, cobj, b.temperature)
        # assumes entr_mol_ig_comp method has only the cp integral part
        entr_ig_comp = get_method(b, "entr_mol_ig_comp", j)(
            b, cobj, b.temperature) - Virial.gas_constant(b) * log(
            b.pressure / b.params.pressure_ref)

        gibbs_ig_comp = (enth_ig_comp - b.temperature * entr_ig_comp)

        if pobj.is_vapor_phase():
            return (gibbs_ig_comp +
                    Virial.gas_constant(b) *
                    b.temperature * (log(b.mole_frac_phase_comp[p, j]) +
                                     Virial.log_fug_coeff_phase_comp(b, p, j)))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def entr_mol_phase_comp(b, p, k):
        # Method for partial molar entropy for species k in gas mixture
        pobj = b.params.get_phase(p)
        #cobj = b.params.get_component(k)
        if pobj.is_vapor_phase():
            return (get_method(b, "entr_mol_ig_comp", k)(b, gcobj(b, k),
                                                         b.temperature) -
                    Virial.gas_constant(b) *
                    (log(b.mole_frac_phase_comp[p, k]) +
                     Virial.log_fug_coeff_phase_comp(b, p, k) +
                     log(b.pressure / b.params.pressure_ref)) +
                    b.pressure / b.temperature *
                    (2 * sum(b.mole_frac_phase_comp[p, i] * b.B_ij[i, k]
                             for i in b.components_in_phase(p)) -
                     Virial.second_virial_coeff_mix(b, p)) -
                    (2 * Virial.gas_constant(b) *
                     sum(b.mole_frac_phase_comp[p, j] / b.Pc_ij[j, k] *
                         (f_dBO_dTr(b.temperature / b.Tc_ij[j, k]) +
                          b.omega_ij[j, k] *
                          f_dB1_dTr(b.temperature / b.Tc_ij[j, k]))
                         for j in b.components_in_phase(p)) -
                        Virial.dBdT(b, p)) * b.pressure)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def enth_mol_phase_comp(b, p, j):
        # Method for partial molar enthalpy for species j in gas mixture
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return (b.gibbs_mol_phase_comp[p, j] + b.temperature *
                    b.entr_mol_phase_comp[p, j])
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def vol_mol_phase_comp(b, p, j):
        # Method for partial molar volume for species j in gas mixture
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return (Virial.gas_constant(b) * b.temperature / b.pressure -
                    Virial.second_virial_coeff_mix(b, p) +
                    2 * sum(b.mole_frac_phase_comp[p, i] * b.B_ij[i, j]
                            for i in b.components_in_phase(p)))
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def energy_internal_mol_phase_comp(b, p, j):
        # Method for partial molar internal energy for species j in gas mixture
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return (b.enth_mol_phase_comp[p, j] - b.pressure *
                    b.vol_mol_phase_comp[p, j])
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    # -------------------------------------------------------------------------
    # Second virial coefficient and its temperature derivatives  are based on
    # combining rules (even if the option to use pseudocritical rules is True)
    # This situation arises for partial molar properties and fugacity of
    # species in  solution.

    @staticmethod
    def second_virial_coeff_mix(b, p):
        # Method for calculating gas mixture second virial coefficient
        # via combining rule only
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return rule_second_virial_coeff(b, p)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def dBdT(b, p):
        # Method for calculating first temperature derivative of gas mixture
        # second virial coefficient
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return sum(sum(b.mole_frac_phase_comp[p, i] *
                           b.mole_frac_phase_comp[p, j] / b.Pc_ij[i, j] * (
                           f_dBO_dTr(b.temperature / b.Tc_ij[i, j]) +
                           b.omega_ij[i, j] * f_dB1_dTr(b.temperature / b.Tc_ij[i, j]))
                           for i in b.components_in_phase(p))
                       for j in b.components_in_phase(p)) * Virial.gas_constant(b)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

    @staticmethod
    def dBdT2(b, p):
        # Method for calculating second temperature derivative of gas
        # mixture second virial coefficient
        pobj = b.params.get_phase(p)
        if pobj.is_vapor_phase():
            return sum(sum(b.mole_frac_phase_comp[p, i] *
                           b.mole_frac_phase_comp[p, j] / b.Pc_ij[i, j] /
                           b.Tc_ij[i, j] * (
                           f_dBO_dTr2(b.temperature / b.Tc_ij[i, j]) +
                           b.omega_ij[i, j] *
                           f_dB1_dTr2(b.temperature / b.Tc_ij[i, j]))
                           for i in b.components_in_phase(p))
                       for j in b.components_in_phase(p)) * Virial.gas_constant(b)
        else:
            raise PropertyNotSupportedError(_invalid_phase_msg(b.name, p))

# -----------------------------------------------------------------------------
# Combining rule for second virial coeffcient


def rule_second_virial_coeff(m, phase):
    return sum(sum(m.mole_frac_phase_comp[phase, i] *
                   m.mole_frac_phase_comp[phase, j] *
                   m.B_ij[i, j]
                   for i in m.components_in_phase(phase))
               for j in m.components_in_phase(phase))

# Pseudocritical rules


def rule_pseudocritical_omega(m):
    return sum(m.mole_frac_phase_comp['Vap', i] *
               m.params.get_component(i).omega
               for i in m.components_in_phase('Vap'))


def rule_pseudocritical_temperature(m):
    return sum(m.mole_frac_phase_comp['Vap', i] *
               m.params.get_component(i).temperature_crit
               for i in m.components_in_phase('Vap'))


def rule_pseudocritical_pressure(m):
    return sum(m.mole_frac_phase_comp['Vap', i] *
               m.params.get_component(i).pressure_crit
               for i in m.components_in_phase('Vap'))

# -----------------------------------------------------------------------------
# Rules for Abbortt equations and their derivatives


def f_BO(Tr):
    return 0.083 - 0.422 / Tr**1.6


def f_B1(Tr):
    return 0.139 - 0.172 / Tr**4.2


def f_dBO_dTr(Tr):
    return 0.6752 / Tr**2.6


def f_dB1_dTr(Tr):
    return 0.7224 / Tr**5.2


def f_dBO_dTr2(Tr):
    return -1.75552 / Tr**3.6


def f_dB1_dTr2(Tr):
    return -3.75648 / Tr**6.2

#------------------------------------------------------------------------------
# virial eos supports only gas properties: error message if used for phases
# other that 'Vap'.


def _invalid_phase_msg(name, phase):
    return (f"{name} received a non-vapor phase {phase}. Virial equation of "
            "state only supports vapor phase properties")
