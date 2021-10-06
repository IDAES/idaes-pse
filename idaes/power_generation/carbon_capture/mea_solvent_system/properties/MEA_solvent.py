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
# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits, Var

# Import IDAES cores
from idaes.core import AqueousPhase, Solvent, Solute, Anion, Cation

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.eos.ideal import Ideal
from idaes.core.util.misc import set_param_from_config


# Set up logger
_log = logging.getLogger(__name__)


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
# Configuration dictionary for aqueous MEA solvent

configuration = {
    # Specifying components
    "components": {
        'H2O': {"type": Solvent,
                "cp_mol_liq_comp": CpMolSolvent,
                "enth_mol_liq_comp": EnthMolSolvent,
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
                        '3': (0.793, pyunits.g/pyunits.mL)}
                    }},
        'MEA': {"type": Solvent,
                "cp_mol_liq_comp": CpMolSolvent,
                "enth_mol_liq_comp": EnthMolSolvent,
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
                        '3': (1.19451, pyunits.g/pyunits.mL)}}},
        'CO2': {"type": Solute,
                "cp_mol_liq_comp": CpMolCO2,
                "enth_mol_liq_comp": EnthMolCO2,
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
                        "equation_of_state": Ideal}},

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
