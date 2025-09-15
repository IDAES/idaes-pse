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
eNRTL Liquid Phase properties for aqueous MEA solvent with CO2.

The following apparent species are used to represent the mixture, along with
the method used to calculate their vapor pressure:

    Carbon Dioxide (CO2) - Henry's Law
    Monoethanolamine (MEA) - non-volatile,
    Water (H2O) - Raoult's Law

Additionally, the following true ionic species are required for calculating
transport properties:

    MEA+, MEACOO-, HCO3-


"""
# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

from copy import deepcopy
from functools import partial

# Import Pyomo units
from pyomo.environ import Constraint, exp, Expression, log, Reals, units as pyunits, value, Var
import pyomo.environ as pyo
from pyomo.core.expr.calculus.derivatives import Modes, differentiate
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent, check_units_equivalent

# Import IDAES cores
from idaes.core import AqueousPhase, Solvent, Solute, Anion, Cation, Zwitterion
from idaes.core.util.exceptions import BurntToast, ConfigurationError
from idaes.core.util.constants import Constants

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.base.generic_property import StateIndex
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.eos.enrtl import ENRTL, EnthMolPhaseBasis
from idaes.models.properties.modular_properties.eos.enrtl_reference_states import InfiniteDilutionSingleSolvent

from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    log_power_law_equil,
)
from idaes.models.properties.modular_properties.base.utility import ConcentrationForm
from idaes.models.properties.modular_properties.phase_equil.henry import HenryType
from idaes.models.properties.modular_properties.pure import NIST
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import gibbs_energy
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn

from idaes.core.util.misc import set_param_from_config
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)

def initialize_inherent_reactions(indexed_blk):
    m = pyo.ConcreteModel()
    idx0 = indexed_blk.index_set().first()
    params = indexed_blk[idx0].params
    m.idx_set = pyo.Set(initialize=[j for j in indexed_blk.index_set()])
    m.inherent_reaction_idx = pyo.Set(initialize=[j for j in params.inherent_reaction_idx])
    m.true_species_set = pyo.Set(initialize=[j for j in params.true_species_set])
    m.reduced_inherent_reaction_extent = Var(
        m.idx_set,
        m.inherent_reaction_idx,
        initialize=0,
        units=pyunits.dimensionless,
        doc="Apparent extent of inherent reactions",
    )
    m.pseudo_mole_frac_comp_true = Var(
        m.idx_set,
        m.true_species_set,
        initialize=1/len(m.idx_set),
        units=pyunits.dimensionless,
        bounds=[0, 1.001],
        doc="Moles of true species j divided by total number of moles of apparent species",
    )
    @m.Expression(m.idx_set)
    def reduced_total_mol(b, idx):
        return sum(b.pseudo_mole_frac_comp_true[idx, j] for j in m.true_species_set)
    
    @m.Param(m.idx_set, m.true_species_set)
    def mole_frac_comp_initial(b, idx, j):
        # Don't like accessing objects in initialization functions
        # but I don't see a better way without adding a ton of verbosity
        # This needs to be changed if dissociation species are present
        if j in params.apparent_species_set:
            return pyo.value(indexed_blk[idx].mole_frac_comp[j])
        else:
            return 0
    
    from pyomo.util.calc_var_value import calculate_variable_from_constraint

    for idx in m.idx_set:
        for r in m.inherent_reaction_idx:
            calculate_variable_from_constraint(
                indexed_blk[idx].log_k_eq[r],
                indexed_blk[idx].log_k_eq_constraint[r]
            )

    @m.Param(m.idx_set, m.inherent_reaction_idx)
    def log_equil_const(b, idx, r):
        return pyo.value(indexed_blk[idx].log_k_eq[r])

    @m.Constraint(m.idx_set, m.true_species_set)
    def material_balance_eqn(b, idx, j):
        return b.pseudo_mole_frac_comp_true[idx, j] == (
            b.mole_frac_comp_initial[idx, j]
            + sum(
                params.inherent_reaction_stoichiometry[r, "Liq", j]
                * b.reduced_inherent_reaction_extent[idx, r]
                for r in b.inherent_reaction_idx
            )
        )
    @m.Objective()
    def obj(b):
        return sum(
            -sum(
                b.log_equil_const[idx, r]
                * b.reduced_inherent_reaction_extent[idx, r]
                for r in b.inherent_reaction_idx
            )
            - b.reduced_total_mol[idx] * pyo.log(b.reduced_total_mol[idx])
            + sum(
                b.pseudo_mole_frac_comp_true[idx, j]
                * pyo.log(b.pseudo_mole_frac_comp_true[idx, j])
                for j in b.true_species_set
            )
            for idx in b.idx_set
        )
    from idaes.core.util.model_statistics import degrees_of_freedom as dof
    assert dof(m) == len(m.idx_set) * len(m.inherent_reaction_idx)

    from idaes.core.solvers import get_solver
    solver_obj = get_solver(
        "ipopt",
        options={
            "halt_on_ampl_error": "yes",
            "bound_relax_factor": 0,
            # "jac_c_constant": "yes"
        }
    )
    res = solver_obj.solve(m, tee=True)
    pyo.assert_optimal_termination(res)

    for idx in m.idx_set:
        for r in m.inherent_reaction_idx:
            indexed_blk[idx].apparent_inherent_reaction_extent[r].value = pyo.value(
                m.reduced_inherent_reaction_extent[idx, r]
                * indexed_blk[idx].flow_mol
            )

        for j in m.true_species_set:
            indexed_blk[idx].flow_mol_phase_comp_true["Liq", j].value = pyo.value(
                m.pseudo_mole_frac_comp_true[idx, j]
                * indexed_blk[idx].flow_mol
            )

            indexed_blk[idx].mole_frac_phase_comp_true["Liq", j].value = pyo.value(
                m.pseudo_mole_frac_comp_true[idx, j]
                / m.reduced_total_mol[idx]
            )


def constant_density(b, *args, **kwargs):
    """Assume constant density of pure water"""
    return (1000 * pyunits.kg / pyunits.m**3) / b.mw 

def create_heat_capacity_no_inherent_rxns(indexed_blk):
    p = "Liq"
    i0 = indexed_blk.index_set().first()
    params = indexed_blk[i0].params
    pobj = params.get_phase(p)
    pname = pobj.local_name
    
    def rule_cp_phase_excess(b, p):
        return differentiate(
                expr=b.enth_mol_phase_excess[p],
                wrt=b.temperature,
                mode=Modes.reverse_symbolic
            )
    
    def rule_cp_phase(b, p):
        pobj = b.params.get_phase(p)
        pname = pobj.local_name
        cp_excess = getattr(b, f"{pname}_cp_phase_excess")
        return (
            sum(
                b.mole_frac_phase_comp_apparent[p, i] * Ideal.cp_mol_phase_comp(b, p, i)
                for i in b.components_in_phase(pname, true_basis=False)
            )
            + cp_excess
        )
    
    for blk in indexed_blk.values():
        blk.add_component(
            f"{pname}_cp_phase_excess",
            Expression(
                rule=partial(rule_cp_phase_excess, p=p),
                doc=f"Excess heat capacity of phase {pname}"
            )
        )
        blk.add_component(
            f"{pname}_cp",
            Expression(
                rule=partial(rule_cp_phase, p=p),
                doc=f"Heat capacity of phase {pname}"
            )
        )

def create_heat_capacity_eqns(indexed_blk):
    p = "Liq"
    i0 = indexed_blk.index_set().first()
    params = indexed_blk[i0].params
    pobj = params.get_phase(p)
    pname = pobj.local_name
    rxn_set = params.inherent_reaction_idx
    units = params.get_metadata().get_derived_units
    
    def rule_d_flow_mol_phase_comp_true_dT(b, j, p):
        params = b.params
        rxn_set = params.inherent_reaction_idx
        pobj = params.get_phase(p)
        pname = pobj.local_name
        dxi_dT = getattr(b, f"{pname}_d_apparent_inherent_reaction_extent_dT")
        return sum(
            params.inherent_reaction_stoichiometry[r, p, j] * dxi_dT[r] 
            for r in rxn_set
            # Want this conditional to simplify Pyomo expressions because
            # they no longer immediately cancel zeros
            if params.inherent_reaction_stoichiometry[r, p, j] != 0
        )
    
    def rule_d_flow_mol_phase_true_dT(b, p):
        pobj = b.params.get_phase(p)
        pname = pobj.local_name
        dn_dT = getattr(b, f"{pname}_d_flow_mol_phase_comp_true_dT")
        return sum(
            dn_dT[j] for j in b.components_in_phase(pname, true_basis=True)
        )
    
    def rule_d_mole_frac_phase_comp_true_dT(b, j, p):
        pobj = b.params.get_phase(p)
        pname = pobj.local_name
        dn_dT = getattr(b, f"{pname}_d_flow_mol_phase_comp_true_dT")
        dn_tot_dT = getattr(b, f"{pname}_d_flow_mol_phase_true_dT")
        return (
            dn_dT[j] / b.flow_mol_phase_comp_true[p, j]
            - b.mole_frac_phase_comp_true[p, j] * dn_tot_dT
        )

    def rule_d_equilibrium_rxn_dT(b, r, p):
        params = b.params
        pobj = params.get_phase(p)
        pname = pobj.local_name
        dx_true_dT = getattr(b, f"{pname}_d_mole_frac_phase_comp_true_dT")
        d_log_gamma_dT = getattr(b, f"{pname}_d_log_gamma_dT" )
        return (
            sum(
                params.inherent_reaction_stoichiometry[r, p, i] * (
                    dx_true_dT[i] 
                    * exp (-b.log_mole_frac_phase_comp_true[p, i])
                    + d_log_gamma_dT[i]
                )
                for i in b.components_in_phase(pname, true_basis=True)
                # Want this conditional to simplify Pyomo expressions because
                # they no longer immediately cancel zeros
                if params.inherent_reaction_stoichiometry[r, p, i] != 0
            )
            == b.dh_rxn[r] / (ENRTL.gas_constant(b) * b.temperature**2)
        )
    
    def rule_cp_phase_excess(b, p):
        return differentiate(
                expr=b.enth_mol_phase_excess[p],
                wrt=b.temperature,
                mode=Modes.reverse_symbolic
            )
    
    def rule_cp_rxn(b, r):
        return differentiate(
            expr=b.dh_rxn[r],
            wrt=b.temperature,
            mode=Modes.reverse_symbolic
        )

    def rule_cp_phase(b, p):        
        params = b.params
        rxn_set = params.inherent_reaction_idx   
        pobj = b.params.get_phase(p)
        pname = pobj.local_name
        dxi_dT = getattr(b, f"{pname}_d_apparent_inherent_reaction_extent_dT")
        cp_excess = getattr(b, f"{pname}_cp_phase_excess")
        cp_rxn = getattr(b, f"{pname}_cp_rxn")
        return (
            sum(
                b.mole_frac_phase_comp_apparent[p, i] * Ideal.cp_mol_phase_comp(b, p, i)
                for i in b.components_in_phase(pname, true_basis=False)
            )
            + sum(
                b.dh_rxn[r] * dxi_dT[r]
                + b.apparent_inherent_reaction_extent[r]
                * cp_rxn[r]
                for r in rxn_set
            )
            + cp_excess
        )
        

    for blk in indexed_blk.values():
        if not hasattr(blk, f"{pname}_d_log_gamma_dT"):
            ENRTL._create_d_log_gamma_dT(blk, pname)
        dxi_dT = blk.add_component(
            f"{pname}_d_apparent_inherent_reaction_extent_dT",
            Var(
                rxn_set,
                doc="Temperature derivative of inherent reaction extent",
                units=units("amount")/units("temperature"),
                initialize=1
            )
        )
        blk.add_component(
            f"{pname}_d_flow_mol_phase_comp_true_dT",
            Expression(
                blk.components_in_phase(pname, true_basis=True),
                rule=partial(rule_d_flow_mol_phase_comp_true_dT, p=p),
                doc="Temperature derivative of true component phase flow",
            )
        )
        blk.add_component(
            f"{pname}_d_flow_mol_phase_true_dT",
            Expression(
                rule=partial(rule_d_flow_mol_phase_true_dT, p=p),
                doc="Temperature derivative of true total phase flow",
            )
        )
        blk.add_component(
            f"{pname}_d_mole_frac_phase_comp_true_dT",
            Expression(
                blk.components_in_phase(pname, true_basis=True),
                rule=partial(rule_d_mole_frac_phase_comp_true_dT, p=p),
                doc="Temperature derivative of true species mole fraction",
            )
        )
        drxn_dT = blk.add_component(
            f"{pname}_d_equilibrium_rxn_dT",
            Constraint(
                rxn_set,
                rule=partial(rule_d_equilibrium_rxn_dT, p=p),
                doc="Temperature derivative of reaction equilibrium equation"
            )
        )
        blk.add_component(
            f"{pname}_cp_phase_excess",
            Expression(
                rule=partial(rule_cp_phase_excess, p=p),
                doc=f"Excess heat capacity of phase {pname}"
            )
        )
        blk.add_component(
            f"{pname}_cp_rxn",
            Expression(
                rxn_set,
                rule=rule_cp_rxn,
                doc=f"Derivative of heat of reaction"
            )
        )
        blk.add_component(
            f"{pname}_cp",
            Expression(
                rule=partial(rule_cp_phase, p=p),
                doc=f"Heat capacity of phase {pname}"
            )
        )
        blk.cp_vars = [dxi_dT,]
        blk.cp_cons = [drxn_dT,]

def deactivate_heat_capacity_eqns(indexed_blk):
    for blk in indexed_blk.values():
        for var in blk.cp_vars:
            var.fix()
        for con in blk.cp_cons:
            con.deactivate()

def activate_heat_capacity_eqns(indexed_blk):
    for blk in indexed_blk.values():
        for var in blk.cp_vars:
            var.unfix()
        for con in blk.cp_cons:
            con.activate()



def scale_heat_capacity_eqns(indexed_blk):
    p = "Liq"
    i0 = indexed_blk.index_set().first()
    params = indexed_blk[i0].params
    pobj = params.get_phase(p)
    pname = pobj.local_name
    rxn_set = params.inherent_reaction_idx
    for blk in indexed_blk.values():
        pass
        


class LinearCpThermoPure:
    class cp_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff_A"):
                cobj.cp_mol_liq_comp_coeff_A = Var(
                    doc="Parameter A for liquid heat capacity",
                    units=pyunits.J * pyunits.mol**-1 * pyunits.K**-1,
                )
                set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="A")

                cobj.cp_mol_liq_comp_coeff_B = Var(
                    doc="Parameter B for liquid heat capacity",
                    units=pyunits.J
                    * pyunits.mol**-1
                    * pyunits.K**-2,
                )
                set_param_from_config(cobj, param="cp_mol_liq_comp_coeff", index="B")

        @staticmethod
        def return_expression(b, cobj, T):
            t = pyunits.convert(T, pyunits.K)
            cp = (
                cobj.cp_mol_liq_comp_coeff_A
                + cobj.cp_mol_liq_comp_coeff_B * t
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)
        
    class enth_mol_liq_comp:
        @staticmethod
        def build_parameters(cobj):
            if not hasattr(cobj, "cp_mol_liq_comp_coeff_A"):
                LinearCpThermoPure.cp_mol_ig_comp.build_parameters(cobj)

        @staticmethod
        def return_expression(b, cobj, T):
            t = pyunits.convert(T, to_units=pyunits.K)
            t0 = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

            h = (
                cobj.cp_mol_liq_comp_coeff_A * t
                + 0.5 * cobj.cp_mol_liq_comp_coeff_B * t ** 2
                # Set enthalpy to zero at T = T_ref
                - cobj.cp_mol_liq_comp_coeff_A * t0
                - 0.5 * cobj.cp_mol_liq_comp_coeff_B * t0 ** 2
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(h, units.ENERGY_MOLE)

class CpMolSolvent:
    # Method to calculate cp for solvents [1]
    @staticmethod
    def build_parameters(cobj):
        cobj.cp_mass_liq_comp_coeff_1 = Var(
            doc="Parameter 1 for liquid phase mass specific heat capacity",
            units=pyunits.kJ * pyunits.kg**-1 * pyunits.K**-1,
        )
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="1")

        cobj.cp_mass_liq_comp_coeff_2 = Var(
            doc="Parameter 2 for liquid phase mass specific heat capacity",
            units=pyunits.kJ * pyunits.kg**-1 * pyunits.K**-2,
        )
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="2")

        cobj.cp_mass_liq_comp_coeff_3 = Var(
            doc="Parameter 3 for liquid phase mass specific heat capacity",
            units=pyunits.kJ * pyunits.kg**-1 * pyunits.K**-3,
        )
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="3")

        cobj.cp_mass_liq_comp_coeff_4 = Var(
            doc="Parameter 4 for liquid phase mass specific heat capacity",
            units=pyunits.kJ * pyunits.kg**-1 * pyunits.K**-4,
        )
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="4")

        cobj.cp_mass_liq_comp_coeff_5 = Var(
            doc="Parameter 5 for liquid phase mass specific heat capacity",
            units=pyunits.kJ * pyunits.kg**-1 * pyunits.K**-5,
        )
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="5")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific heat capacity
        # Need to convert K to C - fake units for now.
        T = pyunits.convert(T, to_units=pyunits.K) - 273.15 * pyunits.K
        cp = cobj.mw * (
            cobj.cp_mass_liq_comp_coeff_5 * T**4
            + cobj.cp_mass_liq_comp_coeff_4 * T**3
            + cobj.cp_mass_liq_comp_coeff_3 * T**2
            + cobj.cp_mass_liq_comp_coeff_2 * T
            + cobj.cp_mass_liq_comp_coeff_1
        )

        units = b.params.get_metadata().derived_units
        return pyunits.convert(cp, units.HEAT_CAPACITY_MOLE)


class EnthMolSolvent:
    # Method to calculate specific enthalpy of solvents
    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mass_liq_comp_coeff_1"):
            CpMolSolvent.build_parameters(cobj)

        cobj.dh_vap = Var(
            doc="Heat of vaporization of component @ Tref",
            units=pyunits.J / pyunits.mol,
        )
        set_param_from_config(cobj, param="dh_vap")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific enthalpy
        T = pyunits.convert(T, to_units=pyunits.K) - 273.15 * pyunits.K
        Tr = (
            pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)
            - 273.15 * pyunits.K
        )

        units = b.params.get_metadata().derived_units

        h = (
            pyunits.convert(
                cobj.mw
                * (
                    (cobj.cp_mass_liq_comp_coeff_5 / 5) * (T**5 - Tr**5)
                    + (cobj.cp_mass_liq_comp_coeff_4 / 4) * (T**4 - Tr**4)
                    + (cobj.cp_mass_liq_comp_coeff_3 / 3) * (T**3 - Tr**3)
                    + (cobj.cp_mass_liq_comp_coeff_2 / 2) * (T**2 - Tr**2)
                    + cobj.cp_mass_liq_comp_coeff_1 * (T - Tr)
                ),
                units.ENERGY_MOLE,
            )
            - cobj.dh_vap
        )

        return h

class HenryCO2ChenEtAl:
    # Correlation for Henry's constant of CO2 in H2O from:
    # Chen, Chau-Chyun, et al. “Extension and Application of the Pitzer Equation 
    # for Vapor-Liquid Equilibrium of Aqueous Electrolyte Systems with Molecular 
    # Solutes.” AIChE Journal, vol. 25, no. 5, Sept. 1979, pp. 820–31. DOI.org 
    # (Crossref), https://doi.org/10.1002/aic.690250510.

    @staticmethod
    def build_parameters(cobj, phase, h_type):
        if hasattr(cobj, "henry_coeff_1"):
            # Created by enthalpy of vaporization method
            return 
        cobj.henry_coeff_1 = Var(
            doc="Henry's constant coefficient 1",
            units=pyunits.K,
        )
        set_param_from_config(cobj, param="henry_coeff", index="1")

        cobj.henry_coeff_2 = Var(
            doc="Henry's constant coefficient 2",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="henry_coeff", index="2")

        cobj.henry_coeff_3 = Var(
            doc="Henry's constant coefficient 3",
            units=1/pyunits.K,
        )
        set_param_from_config(cobj, param="henry_coeff", index="3")

        cobj.henry_coeff_4 = Var(
            doc="Henry's constant coefficient 4",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="henry_coeff", index="4")

    @staticmethod
    def return_log_expression(b, p, j, T):
        derived_units = b.params.get_metadata().get_derived_units
        cobj = b.params.get_component(j)
        M_H2O = b.params.get_component("H2O").mw # Henry's law in water doesn't make sense without water
        # TODO add Poynting term once we have somewhat reliable volumetric data
        return  (
            cobj.henry_coeff_1 / T
            + cobj.henry_coeff_2 * log(T / pyunits.K)
            + cobj.henry_coeff_3 * T
            + cobj.henry_coeff_4
            # + log(
            #     pyunits.convert(
            #         1 / M_H2O * pyunits.atm * pyunits.kg/pyunits.mol, # Convert from molality basis to mole fraction basis
            #         to_units=derived_units("pressure")
            #     )
            # )
        )

    @staticmethod
    def return_expression(b, p, j, T):
        derived_units = b.params.get_metadata().get_derived_units
        return exp(HenryCO2ChenEtAl.return_log_expression(b, p, j, T)) * derived_units("pressure")

    @staticmethod
    def dT_expression(b, p, j, T=None):
        raise NotImplementedError("No dT method for Henry method")

class EnthMolCO2:
    # Only contribution to enthalpy is heat of absorption
    @staticmethod
    def build_parameters(cobj):
        HenryCO2ChenEtAl.build_parameters(cobj, None, None) # Last two arguments are unused
        NIST.enth_mol_ig_comp.build_parameters(cobj)

    @staticmethod
    def return_expression(b, cobj, T):
        units = b.params.get_metadata().derived_units
        R = pyunits.convert(Constants.gas_constant, to_units=units["gas_constant"])
        # TODO is this sign convention right?
        dH_vap = R * (
            cobj.henry_coeff_1
            - cobj.henry_coeff_2 * T
            - cobj.henry_coeff_3 * T**2
        )
        return dH_vap + NIST.enth_mol_ig_comp.return_expression(b, cobj, T)

class PressureSatSolvent:
    # Method for calculating saturation pressure ofsolvents
    @staticmethod
    def build_parameters(cobj):
        cobj.pressure_sat_comp_coeff_1 = Var(
            doc="Coefficient 1 for calculating Psat", units=pyunits.dimensionless
        )
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="1")

        cobj.pressure_sat_comp_coeff_2 = Var(
            doc="Coefficient 2 for calculating Psat", units=pyunits.K
        )
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="2")

        cobj.pressure_sat_comp_coeff_3 = Var(
            doc="Coefficient 3 for calculating Psat", units=pyunits.dimensionless
        )
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="3")

        cobj.pressure_sat_comp_coeff_4 = Var(
            doc="Coefficient 4 for calculating Psat", units=pyunits.K**-2
        )
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="4")

    @staticmethod
    def return_expression(b, cobj, T, dT=False):
        if dT:
            raise Exception("No dT method for pressure sat")

        return exp(
            PressureSatSolvent.return_log_expression(b, cobj, T, dT=dT)
        ) * pyunits.Pa
    
    @staticmethod
    def return_log_expression(b, cobj, T, dT=False):
        if dT:
            raise Exception("No dT method for pressure sat")

        return (
            cobj.pressure_sat_comp_coeff_1
            + cobj.pressure_sat_comp_coeff_2 / T
            + cobj.pressure_sat_comp_coeff_3 * log(T / pyunits.K)
            + cobj.pressure_sat_comp_coeff_4 * T**2
        )

class PressureSatAntoine:
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

            cobj.pressure_sat_comp_coeff_D = Var(
                doc="Typically 1, used to specify the units of pressure", units=pyunits.Pa
            )
            set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="D")

        @staticmethod
        def return_log_expression(b, cobj, T, dT=False):
            units = b.params.get_metadata().derived_units
            if dT:
                d_logPsat_dT = (
                    cobj.pressure_sat_comp_coeff_B
                    / (
                        pyunits.convert(T, to_units=pyunits.K)
                        + cobj.pressure_sat_comp_coeff_C
                    ) ** 2
                )
                return pyunits.convert(d_logPsat_dT, to_units=1/units.TEMPERATURE)
            
            return (
                cobj.pressure_sat_comp_coeff_A
                - cobj.pressure_sat_comp_coeff_B
                / (
                    pyunits.convert(T, to_units=pyunits.K)
                    + cobj.pressure_sat_comp_coeff_C
                )
                + log(
                    pyunits.convert(
                        cobj.pressure_sat_comp_coeff_D,
                        to_units=units.PRESSURE
                    )
                )
            )

        @staticmethod
        def return_expression(b, cobj, T, dT=False):
            if dT:
                return PressureSatAntoine.pressure_sat_comp.dT_expression(b, cobj, T)

            psat = (
                exp(
                    cobj.pressure_sat_comp_coeff_A
                    - cobj.pressure_sat_comp_coeff_B
                    / (
                        pyunits.convert(T, to_units=pyunits.K)
                        + cobj.pressure_sat_comp_coeff_C
                    )
                )
                * cobj.pressure_sat_comp_coeff_D
            )

            units = b.params.get_metadata().derived_units
            return pyunits.convert(psat, to_units=units.PRESSURE)

        @staticmethod
        def dT_expression(b, cobj, T):
            p_sat_dT = (
                PressureSatAntoine.pressure_sat_comp.return_expression(b, cobj, T)
                * cobj.pressure_sat_comp_coeff_B
                / (
                    pyunits.convert(T, to_units=pyunits.K)
                    + cobj.pressure_sat_comp_coeff_C
                )
                ** 2
            )

            units = b.params.get_metadata().derived_units
            dp_units = units.PRESSURE / units.TEMPERATURE
            return pyunits.convert(p_sat_dT, to_units=dp_units)

class logKwMarshallFranck:
    """Source: Marshall, William L., and E. U. Franck. 
    "Ion Product of Water Substance, 0-1000 °C, 1-10,000 Bars New International Formulation and Its Background."
    Journal of Physical and Chemical Reference Data, vol. 10, no. 2, Apr. 1981, pp. 295–304.
    DOI.org (Crossref), https://doi.org/10.1063/1.555643.
    """
    @staticmethod
    def build_parameters(rblock, config):
        if hasattr(rblock, "k_eq_coeff_A"):
            # Already created by enthalpy of reaction
            return
        for i, idx in enumerate(["A", "B", "C", "D"]):
            coeff = Var(
                doc=f"Equilibrium constant coefficient {idx}", units=pyunits.K ** i
            )
            rblock.add_component(f"k_eq_coeff_{idx}", coeff)
            set_param_from_config(rblock, param="k_eq_coeff", index=idx, config=config)
        for i, idx in enumerate(["E", "F", "G"]):
            coeff = Var(
                doc=f"Equilibrium constant coefficient {idx}", units=pyunits.K ** i
            )
            rblock.add_component(coeff, f"k_eq_coeff_{idx}")
            set_param_from_config(rblock, param="k_eq_coeff", index=idx, config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        return exp(b.log_k_eq[r_idx]) * pyunits.dimensionless

    @staticmethod
    def return_log_expression(b, rblock, r_idx, T):
        return b.log_k_eq[r_idx] == (
            rblock.k_eq_coeff_A
            + rblock.k_eq_coeff_B / T
            + rblock.k_eq_coeff_C / T**2
            + rblock.k_eq_coeff_D / T**3
            + log(
                pyunits.convert(
                    b.dens_mass_phase_comp["Liq","H2O"],
                    to_units=pyunits.g / pyunits.mL
                )
            ) * (
                rblock.k_eq_coeff_E
                + rblock.k_eq_coeff_F / T
                + rblock.k_eq_coeff_G / T**2  
            )
            + 2 * log(0.01802) # Convert from molality to mole fraction
        ) * log(10) # Change base from log 10 to natural log

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1

class enthAutoionizationMarshallFranck:
    """Source: Marshall, William L., and E. U. Franck. 
    "Ion Product of Water Substance, 0-1000 °C, 1-10,000 Bars New International Formulation and Its Background."
    Journal of Physical and Chemical Reference Data, vol. 10, no. 2, Apr. 1981, pp. 295–304.
    DOI.org (Crossref), https://doi.org/10.1063/1.555643.
    """
    @staticmethod
    def build_parameters(rblock, config):
        logKwMarshallFranck.build_parameters(rblock, config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        return -Constants.gas_constant * (
            rblock.k_eq_coeff_B
            + 2 * rblock.k_eq_coeff_C / T
            + 3 * rblock.k_eq_coeff_D / T**2
            # + log(
            #     pyunits.convert(
            #         b.dens_mass_phase_comp["Liq","H2O"],
            #         to_units=pyunits.g / pyunits.mL
            #     )
            # ) * (
            #     rblock.k_eq_coeff_F 
            #     + 2 * rblock.k_eq_coeff_G / T
            # + (
            #     rblock.k_eq_coeff_E
            #     + rblock.k_eq_coeff_F / T
            #     + rblock.k_eq_coeff_G / T**2
            # # This property doesn't exist, but is a placeholder
            # ) * b.coeff_thermal_expansion_phase_comp["Liq", "H2O"] 
        ) * log(10) # Change base from log 10 to natural log


    @staticmethod
    def calculate_scaling_factors(b, rblock):
        v = abs(value(rblock.dh_rxn_ref))

        # Need to make sure dh_rxn is not 0 to avoid division by 0
        if v != 0:
            return 1 / abs(value(rblock.dh_rxn_ref))
        else:
            return 1

# -----------------------------------------------------------------------------
# Equilibrium constant model
class KeqCullinaneRochelle:
    @staticmethod
    def build_parameters(rblock, config):
        if hasattr(rblock, "k_eq_coeff_1"):
            # Created by enthalpy of rxn
            return
        rblock.k_eq_coeff_1 = Var(
            doc="Equilibrium constant coefficient 1", units=pyunits.dimensionless
        )
        set_param_from_config(rblock, param="k_eq_coeff", index="1", config=config)

        rblock.k_eq_coeff_2 = Var(
            doc="Equilibrium constant coefficient 2", units=pyunits.K
        )
        set_param_from_config(rblock, param="k_eq_coeff", index="2", config=config)

        rblock.k_eq_coeff_3 = Var(
            doc="Equilibrium constant coefficient 3", units=pyunits.dimensionless
        )
        set_param_from_config(rblock, param="k_eq_coeff", index="3", config=config)

        rblock.k_eq_coeff_4 = Var(
            doc="Equilibrium constant coefficient 3", units=1/pyunits.K
        )
        set_param_from_config(rblock, param="k_eq_coeff", index="4", config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        return exp(b.log_k_eq[r_idx]) * ((pyunits.m) ** 3 / pyunits.mol)

    @staticmethod
    def return_log_expression(b, rblock, r_idx, T):
        return b.log_k_eq[r_idx] == (
            rblock.k_eq_coeff_1
            + rblock.k_eq_coeff_2 / T
            + rblock.k_eq_coeff_3 * log(T / pyunits.K)
            + rblock.k_eq_coeff_4 * T
        )

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1

class enthRxnCullinaneRochelle:
    @staticmethod
    def build_parameters(rblock, config):
        KeqCullinaneRochelle.build_parameters(rblock, config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        return Constants.gas_constant * (
            -rblock.k_eq_coeff_2
            + rblock.k_eq_coeff_3 * T
            + rblock.k_eq_coeff_4 * T ** 2
        )

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        pass

class relative_permittivity_Akula(object):
    """Methods for relative permittivity."""
    
    @staticmethod
    def build_parameters(cobj):
        param_block = cobj.parent_block()
        units = param_block.get_metadata().derived_units

        cobj.relative_permittivity_liq_comp_coeff_A = Var(
            doc="Relative permittivity parameter A", units=pyunits.dimensionless
        )
        set_param_from_config(cobj, param="relative_permittivity_liq_comp_coeff", index="A")
        cobj.relative_permittivity_liq_comp_coeff_B = Var(
            doc="Relative permittivity parameter B", units=units["temperature"]
        )
        set_param_from_config(cobj, param="relative_permittivity_liq_comp_coeff", index="B")
        cobj.relative_permittivity_liq_comp_coeff_C = Var(
            doc="Relative permittivity parameter C", units=units["temperature"]
        )
        set_param_from_config(cobj, param="relative_permittivity_liq_comp_coeff", index="C")

    @staticmethod
    def return_expression(b, cobj, T):
        return (
            cobj.relative_permittivity_liq_comp_coeff_A
            + cobj.relative_permittivity_liq_comp_coeff_B
            * (
                1 / b.temperature - 1 / cobj.relative_permittivity_liq_comp_coeff_C
            )
        )

class AkulaTau(object):
    """Class for methods assuming tau calculated as in Akula et al. (2023)"""

    @staticmethod
    def build_parameters(b):
        param_block = b.parent_block()
        units = param_block.get_metadata().derived_units

        # Get user provided values for tau (if present)
        try:
            tau_A_data = param_block.config.parameter_data[b.local_name + "_tau_A"]
        except KeyError:
            tau_A_data = {}

        try:
            tau_B_data = param_block.config.parameter_data[b.local_name + "_tau_B"]
        except KeyError:
            tau_B_data = {}

        # Check for unused parameters in tau_data
        for data in tau_A_data, tau_B_data:
            for i, j in data.keys():
                if (i, j) not in b.component_pair_set:
                    raise ConfigurationError(
                        "{} eNRTL tau parameter provided for invalid "
                        "component pair {}. Please check typing and only provide "
                        "parameters for valid species pairs.".format(b.name, (i, j))
                    )

        def tau_A_init(b, i, j):
            try:
                return tau_A_data[(i, j)]
            except KeyError:
                # Default interaction value is 0
                return 0
        
        def tau_B_init(b, i, j):
            try:
                return tau_B_data[(i, j)]
            except KeyError:
                # Default interaction value is 0
                return 0

        b.add_component(
            "tau_A",
            Var(
                b.component_pair_set,
                within=Reals,
                initialize=tau_A_init,
                doc="Binary interaction energy parameters",
                units=pyunits.dimensionless,
            ),
        )

        b.add_component(
            "tau_B",
            Var(
                b.component_pair_set,
                within=Reals,
                initialize=tau_B_init,
                doc="Binary interaction energy parameters",
                units=units["temperature"],
            ),
        )

    @staticmethod
    def return_expression(b, pobj, i, j, T):
        if (i, j) in pobj.tau_A:
            return pobj.tau_A[i, j] + pobj.tau_B[i, j] / T
        elif i == j:
            return 0
        else:
            raise BurntToast(
                "{} tau rule encountered unexpected index {}. Please contact"
                "the IDAES Developers with this bug.".format(b.name, (i, j))
            )

class VolMolSolvent:
    # Weiland Method for calculating molar volume of pure solvents [2]

    @staticmethod
    def build_parameters(cobj):
        cobj.dens_mol_liq_comp_coeff_1 = Var(
            doc="Parameter 1 for liquid phase molar density",
            units=pyunits.g / pyunits.mL / pyunits.K**2,
        )
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="1")

        cobj.dens_mol_liq_comp_coeff_2 = Var(
            doc="Parameter 2 for liquid phase molar density",
            units=pyunits.g / pyunits.mL / pyunits.K,
        )
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="2")

        cobj.dens_mol_liq_comp_coeff_3 = Var(
            doc="Parameter 3 for liquid phase molar density",
            units=pyunits.g / pyunits.mL,
        )
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="3")

    @staticmethod
    def return_expression(b, cobj, T):
        T = pyunits.convert(T, to_units=pyunits.K)

        rho = (
            cobj.dens_mol_liq_comp_coeff_1 * T**2
            + cobj.dens_mol_liq_comp_coeff_2 * T
            + cobj.dens_mol_liq_comp_coeff_3
        )
        vol_mol = cobj.mw / rho

        units = b.params.get_metadata().derived_units

        return pyunits.convert(vol_mol, units.VOLUME / units.AMOUNT)

class VolMolCO2:
    # Weiland Method for calculating molar volume of dissolved CO2 [2]

    @staticmethod
    def build_parameters(cobj):
        cobj.vol_mol_liq_comp_coeff_a = Var(
            doc="Parameter a for liquid phase molar volume",
            units=pyunits.mL / pyunits.mol,
        )
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="a")

        cobj.vol_mol_liq_comp_coeff_b = Var(
            doc="Parameter b for liquid phase molar volume",
            units=pyunits.mL / pyunits.mol,
        )
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="b")

        cobj.vol_mol_liq_comp_coeff_c = Var(
            doc="Parameter c for liquid phase molar volume",
            units=pyunits.mL / pyunits.mol,
        )
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="c")

        cobj.vol_mol_liq_comp_coeff_d = Var(
            doc="Parameter d for liquid phase molar volume",
            units=pyunits.mL / pyunits.mol,
        )
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="d")

        cobj.vol_mol_liq_comp_coeff_e = Var(
            doc="Parameter e for liquid phase molar volume",
            units=pyunits.mL / pyunits.mol,
        )
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="e")

    @staticmethod
    def return_expression(b, cobj, T):
        x = b.mole_frac_comp

        vol_mol = (
            cobj.vol_mol_liq_comp_coeff_a
            + (cobj.vol_mol_liq_comp_coeff_b + cobj.vol_mol_liq_comp_coeff_c * x["MEA"])
            * x["MEA"]
            * x["H2O"]
            / x["CO2"]
            + (cobj.vol_mol_liq_comp_coeff_d + cobj.vol_mol_liq_comp_coeff_e * x["MEA"])
            * x["MEA"]
        )

        units = b.params.get_metadata().derived_units

        return pyunits.convert(vol_mol, units.VOLUME / units.AMOUNT)

_component_params = {
    "H2O": {
        "type": Solvent,
        "cp_mol_liq_comp": NIST,
        #"diffus_phase_comp": {"Liq": DiffusNone},
        "enth_mol_liq_comp": NIST,
        "pressure_sat_comp": PressureSatSolvent,
        "vol_mol_liq_comp": VolMolSolvent,
        # Need pure water density for correlation for log K_w
        "relative_permittivity_liq_comp": relative_permittivity_Akula,
        "parameter_data": {
            "mw": (0.01802, pyunits.kg / pyunits.mol),
            "cp_mol_liq_comp_coeff": {
                "A": -203.6060,
                "B": 1523.290,
                "C": -3196.413,
                "D": 2474.455,
                "E": 3.855326,
                "F": -256.5478,
                "G": -488.7163,
                "H": -285.8304,
            },
            "dens_mol_liq_comp_coeff": {
                "1": (-3.2484e-6, pyunits.g / pyunits.mL / pyunits.K**2),  # [2]
                "2": (0.00165, pyunits.g / pyunits.mL / pyunits.K),
                "3": (0.793, pyunits.g / pyunits.mL),
            },
            "dh_vap": 43.99e3,
            "pressure_sat_comp_coeff": {
                "1": 72.55,
                "2": -7206.70,
                "3": -7.1385,
                "4": 4.05e-6,
            },
            "relative_permittivity_liq_comp_coeff": {
                "A": 78.51,
                "B": 31989,
                "C": 298.15,
            },
            "temperature_crit": (647.13, pyunits.K),
        },
    },
    "MEA": {
            "type": Solvent,
            "cp_mol_liq_comp": NIST,
            # "diffus_phase_comp": {"Liq": DiffusMEA},
            "enth_mol_liq_comp": NIST,
            "pressure_sat_comp": PressureSatSolvent,
            "vol_mol_liq_comp": VolMolSolvent,
            "relative_permittivity_liq_comp": relative_permittivity_Akula,
            "parameter_data": {
                "mw": (0.06108, pyunits.kg / pyunits.mol),
                "cp_mol_liq_comp_coeff": {
                    "A": 78.25,
                    "B": 293.2,
                    "C": 0,
                    "D": 0,
                    "E": 0,
                    "F": -36.362, # Normalize enthalpy to 0 at 298
                    "G": 0,
                    "H": 0,
                },
                "dens_mol_liq_comp_coeff": {
                    "1": (-5.35162e-7, pyunits.g / pyunits.mL / pyunits.K**2),  # [2]
                    "2": (-4.51417e-4, pyunits.g / pyunits.mL / pyunits.K),
                    "3": (1.19451, pyunits.g / pyunits.mL),
                },
                #"dh_vap": 58000,  # [3]
                "diffus_phase_comp_coeff": {
                    "1": -13.275,
                    "2": -2198.3,
                    "3": -7.8142e-5,
                },
                "pressure_sat_comp_coeff": {
                    "1": 172.78,
                    "2": -13492,
                    "3": -21.914,
                    "4": 1.38e-5,
                },
                "relative_permittivity_liq_comp_coeff": {
                    "A": 35.76,
                    "B": 14836,
                    "C": 298.15,
                },
                "temperature_crit": (614.45, pyunits.K),
            },
    },
    "CO2": {
        "type": Solute,
        # "cp_mol_liq_comp": CpMolCO2,
        #"diffus_phase_comp": {"Liq": DiffusCO2},
        "enth_mol_liq_comp": EnthMolCO2,
        "henry_component": {
            "Liq": {
                "method": HenryCO2ChenEtAl,
                "type": HenryType.Kpx,
                "basis": StateIndex.true,
            }
        },
        # "vol_mol_liq_comp": VolMolCO2,
        "dens_mol_liq_comp": constant_density,
        "parameter_data": {
            "mw": (0.04401, pyunits.kg / pyunits.mol),
            "cp_mol_ig_comp_coeff": {
                "A": 24.99735,
                "B": 55.18696,
                "C": -33.69137,
                "D": 7.948387,
                "E": -0.136638,
                "F": -403.6075,
                "G": 228.2431,
                "H": -393.5224,
            },
            "henry_coeff": {
                "1": (-8477.711, pyunits.K),
                "2": (-21.95743, pyunits.dimensionless),
                "3": (0.005780758, 1/pyunits.K),
                # Addition term is to convert from molality and ATM to mole fraction and Pa
                "4": (155.1699 + 15.5426, pyunits.dimensionless),
            },
            # "diffus_phase_comp_coeff": {
            #     "1": 2.35e-6,
            #     "2": 2.9837e-8,
            #     "3": -9.7078e-9,
            #     "4": -2119,
            #     "5": -20.132,
            # },
            # "vol_mol_liq_comp_coeff": {
            #     "a": (10.2074, pyunits.mL / pyunits.mol),  # [2]
            #     "b": (-2.2642, pyunits.mL / pyunits.mol),
            #     "c": (3.0059, pyunits.mL / pyunits.mol),
            #     "d": (207, pyunits.mL / pyunits.mol),
            #     "e": (-563.3701, pyunits.mL / pyunits.mol),
            # },
        },
    },
    "H3O^+": {
        "type": Cation,
        "charge": +1,
        # "diffus_phase_comp": {"Liq": DiffusIons},
        "parameter_data": {
            "born_radius": (3.0, pyunits.angstrom)
        },
    },
    "OH^-": {
        "type": Anion,
        "charge": -1,
        # "diffus_phase_comp": {"Liq": DiffusIons},
        "parameter_data": {
            "born_radius": (3.0, pyunits.angstrom)
        },
    },
    "HCO3^-": {
        "type": Anion,
        "charge": -1,
        # "diffus_phase_comp": {"Liq": DiffusNone},
        "parameter_data":{
            "born_radius": (3.0, pyunits.angstrom)
        }
    },
    "CO3^2-": {
        "type": Anion,
        "charge": -2,
        # "diffus_phase_comp": {"Liq": DiffusNone},
        "parameter_data":{
            "born_radius": (3.0, pyunits.angstrom)
        }
    },
    "MEAH^+": {
        "type": Cation,
        "charge": +1,
        # "diffus_phase_comp": {"Liq": DiffusIons},
        "parameter_data": {
            "diffus_phase_comp_coeff": {"1": -22.64, "2": -1000.0, "3": -0.7},
            "born_radius": (300.0, pyunits.angstrom)
        },
    },
    "MEACOO^-": {
        "type": Anion,
        "charge": -1,
        # "diffus_phase_comp": {"Liq": DiffusIons},
        "parameter_data": {
            "diffus_phase_comp_coeff": {"1": -22.64, "2": -1000.0, "3": -0.7},
            "born_radius": (300.0, pyunits.angstrom)
        },
    },
}

_inherent_rxn_dict = {
    ("H2O",): {
        "H2O_autoionization":{
            "stoichiometry": {
                ("Liq", "H2O"): -2,
                ("Liq", "H3O^+"): 1,
                ("Liq", "OH^-"): 1,
            },
            "equilibrium_constant": KeqCullinaneRochelle,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.activity,
            "heat_of_reaction": enthRxnCullinaneRochelle,
            "parameter_data": {
                "k_eq_coeff": {
                    # "A": (-4.098, pyunits.dimensionless),
                    # "B": (-3245.2, pyunits.K),
                    # "C": (2.2362e5, pyunits.K**2),
                    # "D": (-3.984e7, pyunits.K**3),
                    # "E": (13.957, pyunits.dimensionless),
                    # "F": (-1262.3, pyunits.K),
                    # "G": (8.5641e5, pyunits.K**2),
                    # Cullinane and Rochelle (2005)
                    "1": (132.9, pyunits.dimensionless),
                    "2": (-13446, pyunits.K),
                    "3": (-22.48, pyunits.dimensionless),
                    "4": (0.0, 1/pyunits.K),
                }
            },
        },
    },

    ("H2O", "CO2"): {
        "bicarbonate_formation": {
            "stoichiometry": {
                ("Liq", "CO2"): -1,
                ("Liq", "H2O"): -2,
                ("Liq", "HCO3^-"): 1,
                ("Liq", "H3O^+"): 1,
            },
            "equilibrium_constant": KeqCullinaneRochelle,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.activity,
            "heat_of_reaction": enthRxnCullinaneRochelle,
            "parameter_data": {
                # Cullinane and Rochelle (2005)
                "k_eq_coeff": {
                    "1": (231.4, pyunits.dimensionless),
                    "2": (-12092, pyunits.K),
                    "3": (-36.78, pyunits.dimensionless),
                    "4": (0.0, 1/pyunits.K)
                }
            },
        },
        "carbonate_formation": {
            "stoichiometry": {
                ("Liq", "HCO3^-"): -1,
                ("Liq", "H2O"): -1,
                ("Liq", "CO3^2-"): 1,
                ("Liq", "H3O^+"): 1,
            },
            "equilibrium_constant": KeqCullinaneRochelle,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.activity,
            "heat_of_reaction": enthRxnCullinaneRochelle,
            "parameter_data": {
                # Cullinane and Rochelle (2005)
                "k_eq_coeff": {
                    "1": (216.0, pyunits.dimensionless),
                    "2": (-12432, pyunits.K),
                    "3": (-35.48, pyunits.dimensionless),
                    "4": (0.0, 1/pyunits.K)
                }
            },
        },
    },

    ("H2O", "MEA"): {
        "MEA_protonation": {
            "stoichiometry": {
                ("Liq", "MEA"): -1,
                ("Liq", "H3O^+"): -1,
                ("Liq", "MEAH^+"): 1,
                ("Liq", "H2O"): 1,
            },
            "equilibrium_constant": KeqCullinaneRochelle,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.activity,
            "heat_of_reaction": enthRxnCullinaneRochelle,
            "parameter_data": {
                "k_eq_coeff": {
                    # Fit to aggregated data from Kim et al (2011)
                    # Converted from molality to mole fraction
                    # with symmetric activity coeff for MEA
                    "1": (317.994, pyunits.dimensionless),
                    "2": (0, pyunits.K),
                    "3": (-56.339, pyunits.dimensionless),
                    "4": (0.10199, 1/pyunits.K),
                }
            },
        },
    },

    ("H2O", "CO2", "MEA"): {
        "MEA_carbamate_formation": {
            "stoichiometry": {
                ("Liq", "MEA"): -1,
                ("Liq", "CO2"): -1,
                ("Liq", "H2O"): -1,
                ("Liq", "MEACOO^-"): 1,
                ("Liq", "H3O^+"): 1,
            },
            "equilibrium_constant": KeqCullinaneRochelle,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.activity,
            "heat_of_reaction": enthRxnCullinaneRochelle,
            "parameter_data": {
                "k_eq_coeff": {
                    "1": (296.4, pyunits.dimensionless),
                    "2": (3359.0, pyunits.K),
                    "3": (0.0, pyunits.dimensionless),
                    "4": (0.0, 1/pyunits.K),
                }
            }
        },
    },
}

_parameter_data = {
    "Liq_alpha": {
        # Molecule-molecule interactions have default value
        # of 0.3, all others have default value of 0.2
        ("H2O", "MEA"): 0.2,
        # ("H2O", "CO2"): 0.2,
        # ("MEA", "CO2"): 0.2,
        # ("H2O", "MEA_+, MEACOO_-"): 0.2,
        # ("H2O", "MEA_+, HCO3_-"): 0.2,
        # ("MEA", "MEA_+, MEACOO_-"): 0.2,
        # ("MEA", "MEA_+, HCO3_-"): 0.2,
        # ("CO2", "MEA_+, MEACOO_-"): 0.2,
        # ("CO2", "MEA_+, HCO3_-"): 0.2,
    },
    "Liq_tau_A": {
        # Default value of 0
        ("MEA", "H2O"): 1.5201,
        ("H2O", "MEA"): 0.1559,
        # ("H2O", "MEA_+, MEACOO_-"): 18.588108594566588,
        # ("H2O", "MEA_+, HCO3_-"): 8.5721,
        # ("MEA_+, MEACOO_-", "H2O"): -7.055694714651596,
        # ("MEA_+, HCO3_-", "H2O"): -4.0092,
        # ("MEA", "MEA_+, MEACOO_-"): 8,
        # ("MEA", "MEA_+, HCO3_-"): 8,
        # ("CO2", "MEA_+, MEACOO_-"): 8,
        # ("CO2", "MEA_+, HCO3_-"): 8,
        # ("MEA_+, MEACOO_-", "MEA"): -4,
        # ("MEA_+, HCO3_-", "MEA"): -4,
        # ("MEA_+, MEACOO_-", "CO2"): -4,
        # ("MEA_+, HCO3_-", "CO2"): -4,
    },
    "Liq_tau_B": {
        # Default value of 0
        ("MEA", "H2O"): -910.30,
        ("H2O", "MEA"): 110.80,
        # ("H2O", "MEA_+, MEACOO_-"): -1533.6489230168652,
        # ("H2O", "MEA_+, HCO3_-"): 0,
        # ("MEA_+, MEACOO_-", "H2O"): 1224.0939968361579,
        # ("MEA_+, HCO3_-", "H2O"): 0,
        # ("MEA", "MEA_+, MEACOO_-"): 0,
        # ("MEA", "MEA_+, HCO3_-"): 0,
        # ("CO2", "MEA_+, MEACOO_-"): 0,
        # ("CO2", "MEA_+, HCO3_-"): 0,
        # ("MEA_+, MEACOO_-", "MEA"): 0,
        # ("MEA_+, HCO3_-", "MEA"): 0,
        # ("MEA_+, MEACOO_-", "CO2"): 0,
        # ("MEA_+, HCO3_-", "CO2"): 0,
    },
}

_combined_rxn_template = {
    "stoichiometry": {},
    "equilibrium_constant": KeqCullinaneRochelle,
    "heat_of_reaction": enthRxnCullinaneRochelle,
    "equilibrium_form": log_power_law_equil,
    "concentration_form": ConcentrationForm.activity,
    "parameter_data": {
        "k_eq_coeff": {
            "1": (0.0, pyunits.dimensionless),
            "2": (0.0, pyunits.K),
            "3": (0.0, pyunits.dimensionless),
            "4": (0.0, 1/pyunits.K)
        }
    },
}

# returns a configuration dictionary for the list of specified components
# Note: H2O must be a component
def get_prop_dict(components=None, excluded_rxns=None, rxn_combinations=None):
    """
    excluded_rxns is a Python *set* containing reactions to exclude
    rxn_combinations is a dictionary of the form
    {
    "combination1_name": {
        "reaction1_name": reaction_stoich_coeff,
        "reaction2_name": reaction_stoich_coeff,
        ...
        },
        ...
    }
    All reactions featured in any combination are removed from the list of reactions included
    """
    if components is None:
        components = list(_component_params.keys())
    # Create set for use of subset methods.
    # Cannot use as primary iterable because sets are unordered and so
    # iterating over components might produce different orders when
    # run at different times (been a problem in the past)
    comp_set = set(components) 

    assert "H2O" in comp_set
    assert comp_set.issubset({"H2O", "MEA", "CO2"})
    configuration = {
        "components": {},  # fill in later based on selected components
        "parameter_data": {},
        "phases": {
            "Liq": {
                "type": AqueousPhase,
                "equation_of_state": ENRTL,
                "equation_of_state_options": {
                    "property_basis": "true",
                    "tau_rule": AkulaTau,
                    "reference_state": InfiniteDilutionSingleSolvent,
                    "reference_component": "H2O",
                    "enth_mol_phase_basis": EnthMolPhaseBasis.apparent,
                },
                # "surf_tens_phase": SurfTens,
                # "therm_cond_phase": ThermalCond,
                # "visc_d_phase": Viscosity,
                "parameter_data": {},
            }
        },
        # Set base units of measurement
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        # Specifying state definition
        "state_definition": FTPx,
        "state_bounds": {
            "flow_mol": (0, 1, 1000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 298.15, 500, pyunits.K),
            "pressure": (5e2, 101325, 1e7, pyunits.Pa),
        },
        "state_components": StateIndex.apparent,
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "inherent_reactions": {}
    }
    # Add inherent reactions appropriate to the chosen components
    raw_inherent_reactions = {}
    for key, rxns in _inherent_rxn_dict.items():
        # Keys have to be tuples, so convert them to sets for comparison
        if set(key).issubset(comp_set):
            for rxn_name, rxn_dict in rxns.items():
                raw_inherent_reactions[rxn_name] = deepcopy(rxn_dict)

    if rxn_combinations is not None:
        combined_rxns_dict = {}
        if excluded_rxns is None:
            excluded_rxns = set()
        else:
            excluded_rxns=set(excluded_rxns)
        for combo_name, combo_dict in rxn_combinations.items():
            combined_rxn_dict = deepcopy(_combined_rxn_template)
            for component_rxn, rxn_stoich_coeff in combo_dict.items():
                assert component_rxn in raw_inherent_reactions
                assert isinstance(rxn_stoich_coeff, int)
                assert raw_inherent_reactions[component_rxn]["equilibrium_constant"] is KeqCullinaneRochelle
                assert raw_inherent_reactions[component_rxn]["heat_of_reaction"] is enthRxnCullinaneRochelle
                assert raw_inherent_reactions[component_rxn]["concentration_form"] is ConcentrationForm.activity
                for pc_pair, stoich_coeff in raw_inherent_reactions[component_rxn]["stoichiometry"].items():
                    assert isinstance(stoich_coeff, int)
                    combined_rxn_dict["stoichiometry"][pc_pair] = (
                        combined_rxn_dict["stoichiometry"].get(pc_pair, 0)
                        + rxn_stoich_coeff * stoich_coeff
                    )
                for idx, coeff_tuple in raw_inherent_reactions[component_rxn]["parameter_data"]["k_eq_coeff"].items():
                    assert_units_equivalent(
                        coeff_tuple[1],
                        combined_rxn_dict["parameter_data"]["k_eq_coeff"][idx][1]
                    )
                    coeff_tuple_old = combined_rxn_dict["parameter_data"]["k_eq_coeff"][idx]
                    combined_rxn_dict["parameter_data"]["k_eq_coeff"][idx] = (
                        coeff_tuple_old[0] + rxn_stoich_coeff * coeff_tuple[0],
                        coeff_tuple_old[1]
                    ) 
                excluded_rxns.add(component_rxn)
            # Clean up any stoichiometric coefficients equal to zero
            redundant_comps = set()
            for pc_pair, stoich_coeff in combined_rxn_dict["stoichiometry"].items():
                if stoich_coeff == 0:
                    redundant_comps.add(pc_pair)
            
            for pc_pair in redundant_comps:
                del combined_rxn_dict["stoichiometry"][pc_pair]
            combined_rxns_dict[combo_name] = combined_rxn_dict
            
    
    if excluded_rxns is not None:
        for rxn in excluded_rxns:
            del raw_inherent_reactions[rxn]
    
    if rxn_combinations is None:
        configuration["inherent_reactions"] = raw_inherent_reactions
    else:
        combined_rxns_dict.update(raw_inherent_reactions)
        configuration["inherent_reactions"] = combined_rxns_dict
    

    # Iterate over reaction stoichiometry to find out
    # which components to include
    for comp, comp_dict in _component_params.items():
        include_comp = False
        for rxn_dict in configuration["inherent_reactions"].values():
            for phase_comp_tuple in rxn_dict["stoichiometry"]:
                if phase_comp_tuple[1] == comp:
                    include_comp = True
        if include_comp:
            configuration["components"][comp] = deepcopy(comp_dict)
    included_comps = {key for key in configuration["components"].keys()}

    # In case all reactions involving a component are removed
    for comp in components:
        if comp not in included_comps:
            configuration["components"][comp] = deepcopy(_component_params[comp])
    
    for param, param_dict in _parameter_data.items():
        # If no parameters are filled in, creating an empty dict won't cause harm
        configuration["parameter_data"][param] = {}
        for idx, val in param_dict.items():
            # Deaggregate cation anion pairs by splitting along ,
            relevant_components = idx[0].split(",") + idx[1].split(",")
            # Strip white space from cation anion pairs
            relevant_components = set(comp.strip() for comp in relevant_components)
            if relevant_components.issubset(included_comps):
                configuration["parameter_data"][param][idx] = deepcopy(val)

    return configuration