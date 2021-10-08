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
Ideal Liquid Phase properties for aqueous MEA solvent with CO2.

The following apparent species are used to represent the mixture, along with
the method used to calculate theri vapor pressure:

    Carbon Dioxide (CO2) - Henry's Law
    Monoethanolamine (MEA) - non-volatile,
    Water (H2O) - Raoult's Law

Additionally, the following true ionic species are requried for calculating
transport properties:

    MEA_+, MEACOO_-, HCO3_-

Assumptions:
    * No heats of formation
    * CO2 contribution to enthalpy is due to heat of absorption

References:
    [1] Hilliard thesis (1998)
    [2] Morgan et.al (2015)
"""
# Import Pyomo units
from pyomo.environ import exp, log, units as pyunits, Var

# Import IDAES cores
from idaes.core import AqueousPhase, Solvent, Solute, Anion, Cation

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.core.util.misc import set_param_from_config
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


# -----------------------------------------------------------------------------
# Pure Component Property methods for aqueous MEA
class CpMolCO2():
    # No contribution form CO2, return 0*J/mol.K
    @staticmethod
    def return_expression(*args, **kwargs):
        # Need a very small number to avoid unit consistency probelms
        return 1e-20*pyunits.J/pyunits.mol/pyunits.K


class CpMolSolvent():
    # Method to calculate cp for solvents [1]
    @staticmethod
    def build_parameters(cobj):
        cobj.cp_mass_liq_comp_coeff_1 = Var(
            doc="Parameter 1 for liquid phase mass specific heat capacity",
            units=pyunits.kJ*pyunits.kg**-1*pyunits.K**-1)
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="1")

        cobj.cp_mass_liq_comp_coeff_2 = Var(
            doc="Parameter 2 for liquid phase mass specific heat capacity",
            units=pyunits.kJ*pyunits.kg**-1*pyunits.K**-2)
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="2")

        cobj.cp_mass_liq_comp_coeff_3 = Var(
            doc="Parameter 3 for liquid phase mass specific heat capacity",
            units=pyunits.kJ*pyunits.kg**-1*pyunits.K**-3)
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="3")

        cobj.cp_mass_liq_comp_coeff_4 = Var(
            doc="Parameter 4 for liquid phase mass specific heat capacity",
            units=pyunits.kJ*pyunits.kg**-1*pyunits.K**-4)
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="4")

        cobj.cp_mass_liq_comp_coeff_5 = Var(
            doc="Parameter 5 for liquid phase mass specific heat capacity",
            units=pyunits.kJ*pyunits.kg**-1*pyunits.K**-5)
        set_param_from_config(cobj, param="cp_mass_liq_comp_coeff", index="5")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific heat capacity
        # Need to convert K to C - fake units for now.
        T = pyunits.convert(T, to_units=pyunits.K) - 273.15*pyunits.K
        cp = cobj.mw*(cobj.cp_mass_liq_comp_coeff_5*T**4 +
                      cobj.cp_mass_liq_comp_coeff_4*T**3 +
                      cobj.cp_mass_liq_comp_coeff_3*T**2 +
                      cobj.cp_mass_liq_comp_coeff_2*T +
                      cobj.cp_mass_liq_comp_coeff_1)

        units = b.params.get_metadata().derived_units
        return pyunits.convert(cp, units["heat_capacity_mole"])


class EnthMolCO2():
    # Only contribution to enthalpy is heat of absorption
    @staticmethod
    def build_parameters(cobj):
        cobj.dh_abs_co2 = Var(
                doc="Heat of absorption of CO2 @ Tref",
                units=pyunits.J/pyunits.mol)
        set_param_from_config(cobj, param="dh_abs_co2")

    @staticmethod
    def return_expression(b, cobj, T):
        return cobj.dh_abs_co2


class EnthMolSolvent():
    # Method to calculate specific enthalpy of solvents
    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mass_liq_comp_coeff_1"):
            CpMolSolvent.build_parameters(cobj)

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific enthalpy
        T = pyunits.convert(T, to_units=pyunits.K)
        Tr = pyunits.convert(b.params.temperature_ref, to_units=pyunits.K)

        units = b.params.get_metadata().derived_units

        h = (pyunits.convert(
            cobj.mw * (
                (cobj.cp_mass_liq_comp_coeff_5/5)*(T**5-Tr**5) +
                (cobj.cp_mass_liq_comp_coeff_4/4)*(T**4-Tr**4) +
                (cobj.cp_mass_liq_comp_coeff_3/3)*(T**3-Tr**3) +
                (cobj.cp_mass_liq_comp_coeff_2/2)*(T**2-Tr**2) +
                cobj.cp_mass_liq_comp_coeff_1*(T-Tr)),
            units["energy_mole"]))

        return h


class N2OAnalogy():
    # Henry's constant N2O Analogy Jiru et.al (2012)
    # Units of original expression are Pa*m^3/mol
    # TODO: Handle units of H
    def return_expression(b, p, j, T):
        t = T - 273.15*pyunits.K

        # Calculations require mass fraction of MEA and H2O on a CO2 free basis
        m_MEA = b.mole_frac_comp['MEA']*b.params.MEA.mw
        m_H2O = b.mole_frac_comp['H2O']*b.params.H2O.mw
        wt_MEA = m_MEA/(m_MEA+m_H2O)
        wt_H2O = m_H2O/(m_MEA+m_H2O)
        H_N2O_MEA = 2.448e5 * exp(-1348*pyunits.K / T)
        H_CO2_H2O = 3.52e6 * exp(-2113*pyunits.K / T)
        H_N2O_H2O = 8.449e6 * exp(-2283*pyunits.K / T)
        H_CO2_MEA = H_N2O_MEA * (H_CO2_H2O / H_N2O_H2O)
        lwm = (1.70981 + 0.03972*pyunits.K**-1 * t -
               4.3e-4*pyunits.K**-2 * t**2 -
               2.20377 * wt_H2O)

        return ((exp(wt_MEA * log(H_CO2_MEA) + wt_H2O * log(H_CO2_H2O) +
                     wt_MEA * wt_H2O * lwm)) *
                pyunits.Pa*pyunits.m**-3*pyunits.mol**-1)


class PressureSatSolvent():
    # Method for calculating saturation pressure ofsolvents
    @staticmethod
    def build_parameters(cobj):
        cobj.pressure_sat_comp_coeff_1 = Var(
                doc="Coefficient 1 for calculating Psat",
                units=pyunits.dimensionless)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="1")

        cobj.pressure_sat_comp_coeff_2 = Var(
                doc="Coefficient 2 for calculating Psat",
                units=pyunits.K)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="2")

        cobj.pressure_sat_comp_coeff_3 = Var(
                doc="Coefficient 3 for calculating Psat",
                units=pyunits.dimensionless)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="3")

        cobj.pressure_sat_comp_coeff_4 = Var(
                doc="Coefficient 4 for calculating Psat",
                units=pyunits.K**-2)
        set_param_from_config(cobj, param="pressure_sat_comp_coeff", index="4")

    @staticmethod
    def return_expression(b, cobj, T, dT=False):
        if dT:
            return PressureSatSolvent.dT_expression(b, cobj, T)

        return (exp(cobj.pressure_sat_comp_coeff_1 +
                    cobj.pressure_sat_comp_coeff_2/T +
                    cobj.pressure_sat_comp_coeff_3*log(T/pyunits.K) +
                    cobj.pressure_sat_comp_coeff_4*T**2))*pyunits.Pa


class VolMolSolvent():
    # Weiland Method for calculating molar volume of pure solvents [2]

    @staticmethod
    def build_parameters(cobj):
        cobj.dens_mol_liq_comp_coeff_1 = Var(
                doc="Parameter 1 for liquid phase molar density",
                units=pyunits.g/pyunits.mL/pyunits.K**2)
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="1")

        cobj.dens_mol_liq_comp_coeff_2 = Var(
                doc="Parameter 2 for liquid phase molar density",
                units=pyunits.g/pyunits.mL/pyunits.K)
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="2")

        cobj.dens_mol_liq_comp_coeff_3 = Var(
                doc="Parameter 3 for liquid phase molar density",
                units=pyunits.g/pyunits.mL)
        set_param_from_config(cobj, param="dens_mol_liq_comp_coeff", index="3")

    @staticmethod
    def return_expression(b, cobj, T):
        T = pyunits.convert(T, to_units=pyunits.K)

        rho = (cobj.dens_mol_liq_comp_coeff_1*T**2 +
               cobj.dens_mol_liq_comp_coeff_2*T +
               cobj.dens_mol_liq_comp_coeff_3)
        vol_mol = cobj.mw/rho

        units = b.params.get_metadata().derived_units

        return pyunits.convert(vol_mol, units["volume"]/units["amount"])


class VolMolCO2():
    # Weiland Method for calculating molar volume of disolved CO2 [2]

    @staticmethod
    def build_parameters(cobj):
        cobj.vol_mol_liq_comp_coeff_a = Var(
                doc="Parameter a for liquid phase molar volume",
                units=pyunits.mL/pyunits.mol)
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="a")

        cobj.vol_mol_liq_comp_coeff_b = Var(
                doc="Parameter b for liquid phase molar volume",
                units=pyunits.mL/pyunits.mol)
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="b")

        cobj.vol_mol_liq_comp_coeff_c = Var(
                doc="Parameter c for liquid phase molar volume",
                units=pyunits.mL/pyunits.mol)
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="c")

        cobj.vol_mol_liq_comp_coeff_d = Var(
                doc="Parameter d for liquid phase molar volume",
                units=pyunits.mL/pyunits.mol)
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="d")

        cobj.vol_mol_liq_comp_coeff_e = Var(
                doc="Parameter e for liquid phase molar volume",
                units=pyunits.mL/pyunits.mol)
        set_param_from_config(cobj, param="vol_mol_liq_comp_coeff", index="e")

    @staticmethod
    def return_expression(b, cobj, T):
        x = b.mole_frac_comp

        vol_mol = (cobj.vol_mol_liq_comp_coeff_a +
                   (cobj.vol_mol_liq_comp_coeff_b +
                    cobj.vol_mol_liq_comp_coeff_c * x["MEA"]) *
                   x["MEA"]*x["H2O"]/x["CO2"] +
                   (cobj.vol_mol_liq_comp_coeff_d +
                    cobj.vol_mol_liq_comp_coeff_e * x["MEA"]) *
                   x["MEA"])

        units = b.params.get_metadata().derived_units

        return pyunits.convert(vol_mol, units["volume"]/units["amount"])


# -----------------------------------------------------------------------------
# Transport property models
class Viscosity():
    def build_parameters(pobj):
        pobj.visc_d_coeff_a = Var(
                doc="Parameter a for liquid phase viscosity model",
                units=pyunits.K**-1)
        set_param_from_config(pobj, param="visc_d_coeff", index="a")

        pobj.visc_d_coeff_b = Var(
                doc="Parameter b for liquid phase viscosity model",
                units=pyunits.K**-1)
        set_param_from_config(pobj, param="visc_d_coeff", index="b")

        pobj.visc_d_coeff_c = Var(
                doc="Parameter c for liquid phase viscosity model",
                units=pyunits.dimensionless)
        set_param_from_config(pobj, param="visc_d_coeff", index="c")

        pobj.visc_d_coeff_d = Var(
                doc="Parameter d for liquid phase viscosity model",
                units=pyunits.dimensionless)
        set_param_from_config(pobj, param="visc_d_coeff", index="d")

        pobj.visc_d_coeff_e = Var(
                doc="Parameter e for liquid phase viscosity model",
                units=pyunits.dimensionless)
        set_param_from_config(pobj, param="visc_d_coeff", index="e")

        pobj.visc_d_coeff_f = Var(
                doc="Parameter f for liquid phase viscosity model",
                units=pyunits.K**-1)
        set_param_from_config(pobj, param="visc_d_coeff", index="f")

        pobj.visc_d_coeff_g = Var(
                doc="Parameter g for liquid phase viscosity model",
                units=pyunits.dimensionless)
        set_param_from_config(pobj, param="visc_d_coeff", index="g")

    def return_expression(blk, phase):
        pobj = blk.params.get_phase(phase)

        # Calculate mass fraction from mole fraction and molecular weights
        r = (blk.mole_frac_comp['MEA'] *
             blk.mw_phase["Liq"] / blk.mw_comp["MEA"] * 100)
        T = blk.temperature
        alpha = blk.mole_frac_comp['CO2'] / blk.mole_frac_comp['MEA']
        mu_H2O = (1.002e-3*pyunits.Pa/pyunits.s *
                  10**((1.3272 *
                        (293.15*pyunits.K - T -
                         0.001053*pyunits.K**-1 * (T - 293.15*pyunits.K)**2)) /
                       (T - 168.15*pyunits.K)))
        a = pobj.visc_d_coeff_a
        b = pobj.visc_d_coeff_b
        c = pobj.visc_d_coeff_c
        d = pobj.visc_d_coeff_d
        e = pobj.visc_d_coeff_e
        f = pobj.visc_d_coeff_f
        g = pobj.visc_d_coeff_g

        # Model appears to be entirely empirical, and units are not obvious
        # Assume each part of expression is unitless and use units to match
        return mu_H2O * exp(r *
                            (T * (a * r + b) + c * r + d) *
                            (alpha * (e * r + f * T + g) + 1) /
                            T**2 * pyunits.K**2)


# -----------------------------------------------------------------------------
# Configuration dictionary for aqueous MEA solvent
configuration = {
    # Specifying components
    "components": {
        'H2O': {"type": Solvent,
                "cp_mol_liq_comp": CpMolSolvent,
                "enth_mol_liq_comp": EnthMolSolvent,
                "pressure_sat_comp": PressureSatSolvent,
                "vol_mol_liq_comp": VolMolSolvent,
                "parameter_data": {
                    "mw": (0.01802, pyunits.kg/pyunits.mol),
                    "cp_mass_liq_comp_coeff": {
                        '1': 4.2107,
                        '2': -1.696e-3,
                        '3': 2.568e-5,
                        '4': -1.095e-7,
                        '5': 3.038e-10},
                    "dens_mol_liq_comp_coeff": {
                        '1': (-3.2484e-6, pyunits.g/pyunits.mL/pyunits.K**2),  # [2]
                        '2': (0.00165, pyunits.g/pyunits.mL/pyunits.K),
                        '3': (0.793, pyunits.g/pyunits.mL)},
                    "pressure_sat_comp_coeff": {
                        '1': 72.55,
                        '2': -7206.70,
                        '3': -7.1385,
                        '4': 4.05e-6}
                    }},
        'MEA': {"type": Solvent,
                "cp_mol_liq_comp": CpMolSolvent,
                "enth_mol_liq_comp": EnthMolSolvent,
                "pressure_sat_comp": PressureSatSolvent,
                "vol_mol_liq_comp": VolMolSolvent,
                "parameter_data": {
                    "mw": (0.06108, pyunits.kg/pyunits.mol),
                    "cp_mass_liq_comp_coeff": {
                        '1': 2.6161,
                        '2': 3.706e-3,
                        '3': 3.787e-6,
                        '4': 0.0,
                        '5': 0.0},
                    "dens_mol_liq_comp_coeff": {
                        '1': (-5.35162e-7, pyunits.g/pyunits.mL/pyunits.K**2),  # [2]
                        '2': (-4.51417e-4, pyunits.g/pyunits.mL/pyunits.K),
                        '3': (1.19451, pyunits.g/pyunits.mL)},
                    "pressure_sat_comp_coeff": {
                        '1': 172.78,
                        '2': -13492,
                        '3': -21.914,
                        '4': 1.38e-5}
                    }},
        'CO2': {"type": Solute,
                "cp_mol_liq_comp": CpMolCO2,
                "enth_mol_liq_comp": EnthMolCO2,
                "henry_component": {"Liq": N2OAnalogy},
                "vol_mol_liq_comp": VolMolCO2,
                "parameter_data": {
                    "mw": (0.04401, pyunits.kg/pyunits.mol),
                    "dh_abs_co2": -84000,
                    "vol_mol_liq_comp_coeff": {
                        'a': (10.2074, pyunits.mL/pyunits.mol),  # [2]
                        'b': (-2.2642, pyunits.mL/pyunits.mol),
                        'c': (3.0059, pyunits.mL/pyunits.mol),
                        'd': (207, pyunits.mL/pyunits.mol),
                        'e': (-563.3701, pyunits.mL/pyunits.mol)}}}},

    # Specifying phases
    "phases":  {'Liq': {"type": AqueousPhase,
                        "equation_of_state": Ideal,
                        "visc_d_phase": Viscosity,
                        "parameter_data": {
                            "visc_d_coeff": {
                                "a": -0.0838,
                                "b": 2.8817,
                                "c": 33.651,
                                "d": 1817.0,
                                "e": 0.00847,
                                "f": 0.0103,
                                "g": -2.3890}}}},

    # Set base units of measurement
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},

    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {"flow_mol": (0, 1, 1000, pyunits.mol/pyunits.s),
                     "temperature": (273.15, 298.15, 450, pyunits.K),
                     "pressure": (5e4, 101325, 1e6, pyunits.Pa)},
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K)}
