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
the method used to calculate their vapor pressure:

    Carbon Dioxide (CO2) - Henry's Law
    Monoethanolamine (MEA) - non-volatile,
    Water (H2O) - Raoult's Law

Additionally, the following true ionic species are requried for calculating
transport properties:

    MEA+, MEACOO-, HCO3-

Assumptions:
    * No heats of formation
    * CO2 contribution to enthalpy is due to heat of absorption

References:
    [1] Hilliard thesis (1998)
    [2] Morgan et.al (2015)
    [3] NIST Webbook, https://webbook.nist.gov/
"""
# Import Pyomo units
from pyomo.environ import exp, log, units as pyunits, Var, Expression

# Import IDAES cores
from idaes.core import AqueousPhase, Solvent, Solute, Anion, Cation

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.base.generic_property import StateIndex
from idaes.models.properties.modular_properties.eos.ideal import Ideal

from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    log_power_law_equil,
)
from idaes.models.properties.modular_properties.base.utility import ConcentrationForm
from idaes.models.properties.modular_properties.phase_equil.henry import HenryType

from idaes.core.util.misc import set_param_from_config
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


# -----------------------------------------------------------------------------
# Pure Component Property methods for aqueous MEA
class CpMolCO2:
    # No contribution to enthalpy from CO2, but need to add a contribution
    # to cp to avoid artificial decrease in as CO2 is absorbed.
    # Assume cp contribution of CO2 is equal to weighted average of cp for
    # solvents
    @staticmethod
    def return_expression(b, cobj, T):
        return (
            (
                b.mole_frac_phase_comp["Liq", "H2O"] * b.cp_mol_phase_comp["Liq", "H2O"]
                + b.mole_frac_phase_comp["Liq", "MEA"]
                * b.cp_mol_phase_comp["Liq", "MEA"]
            )
            / b.mole_frac_phase_comp["Liq", "CO2"]
            * (1 / (1 - b.mass_frac_phase_comp["Liq", "CO2"]) - 1)
        )


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


class EnthMolCO2:
    # Only contribution to enthalpy is heat of absorption
    @staticmethod
    def build_parameters(cobj):
        cobj.dh_abs_co2 = Var(
            doc="Heat of absorption of CO2 @ Tref", units=pyunits.J / pyunits.mol
        )
        set_param_from_config(cobj, param="dh_abs_co2")

    @staticmethod
    def return_expression(b, cobj, T):
        return cobj.dh_abs_co2


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


class N2OAnalogy:
    # Henry's constant N2O Analogy Jiru et.al (2012)
    # Units of original expression are Pa*m^3/mol
    # TODO: Handle units of H
    @staticmethod
    def build_parameters(cobj, phase, h_type):
        cobj.lwm_coeff_1 = Var(
            doc="N2O Analogy Henry's constant coefficient 1",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="lwm_coeff", index="1")

        cobj.lwm_coeff_2 = Var(
            doc="N2O Analogy Henry's constant coefficient 2",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="lwm_coeff", index="2")

        cobj.lwm_coeff_3 = Var(
            doc="N2O Analogy Henry's constant coefficient 3",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="lwm_coeff", index="3")

        cobj.lwm_coeff_4 = Var(
            doc="N2O Analogy Henry's constant coefficient 4",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="lwm_coeff", index="4")

    @staticmethod
    def return_expression(b, p, j, T):
        cobj = b.params.get_component(j)
        t = T - 273.15 * pyunits.K

        # Calculations require mass fraction of MEA and H2O on a CO2 free basis
        m_MEA = b.mole_frac_comp["MEA"] * b.params.MEA.mw
        m_H2O = b.mole_frac_comp["H2O"] * b.params.H2O.mw
        wt_MEA = m_MEA / (m_MEA + m_H2O)
        wt_H2O = m_H2O / (m_MEA + m_H2O)
        H_N2O_MEA = 2.448e5 * exp(-1348 * pyunits.K / T)
        H_CO2_H2O = 3.52e6 * exp(-2113 * pyunits.K / T)
        H_N2O_H2O = 8.449e6 * exp(-2283 * pyunits.K / T)
        H_CO2_MEA = H_N2O_MEA * (H_CO2_H2O / H_N2O_H2O)
        lwm = (
            cobj.lwm_coeff_1
            + cobj.lwm_coeff_2 * pyunits.K**-1 * t
            + cobj.lwm_coeff_3 * pyunits.K**-2 * t**2
            + cobj.lwm_coeff_4 * wt_H2O
        )

        return (
            (
                exp(
                    wt_MEA * log(H_CO2_MEA)
                    + wt_H2O * log(H_CO2_H2O)
                    + wt_MEA * wt_H2O * lwm
                )
            )
            * pyunits.Pa
            * pyunits.m**3
            * pyunits.mol**-1
        )


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

        return (
            exp(
                cobj.pressure_sat_comp_coeff_1
                + cobj.pressure_sat_comp_coeff_2 / T
                + cobj.pressure_sat_comp_coeff_3 * log(T / pyunits.K)
                + cobj.pressure_sat_comp_coeff_4 * T**2
            )
        ) * pyunits.Pa


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
    # Weiland Method for calculating molar volume of disolved CO2 [2]

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


# -----------------------------------------------------------------------------
# Transport property models
class DiffusCO2:
    @staticmethod
    def build_parameters(cobj, phase):
        cobj.diffus_phase_comp_coeff_1 = Var(
            doc="Parameter 1 for liquid phase diffusivity model",
            units=pyunits.m**2 / pyunits.s,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="1")
        cobj.diffus_phase_comp_coeff_2 = Var(
            doc="Parameter 2 for liquid phase diffusivity model",
            units=pyunits.m**5 / pyunits.kmol / pyunits.s,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="2")
        cobj.diffus_phase_comp_coeff_3 = Var(
            doc="Parameter 3 for liquid phase diffusivity model",
            units=pyunits.m**8 / pyunits.kmol**2 / pyunits.s,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="3")
        cobj.diffus_phase_comp_coeff_4 = Var(
            doc="Parameter 4 for liquid phase diffusivity model", units=pyunits.K
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="4")
        cobj.diffus_phase_comp_coeff_5 = Var(
            doc="Parameter 5 for liquid phase diffusivity model",
            units=(pyunits.m**3) * pyunits.K / pyunits.kmol,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="5")

    @staticmethod
    def return_expression(blk, p, j, T):
        cobj = blk.params.get_component(j)
        C_MEA = blk.conc_mol_comp["MEA"] * 1e-3 * pyunits.kmol / pyunits.mol
        return (
            cobj.diffus_phase_comp_coeff_1
            + cobj.diffus_phase_comp_coeff_2 * C_MEA
            + cobj.diffus_phase_comp_coeff_3 * C_MEA**2
        ) * exp(
            (cobj.diffus_phase_comp_coeff_4 + (cobj.diffus_phase_comp_coeff_5 * C_MEA))
            / T
        )


class DiffusMEA:
    @staticmethod
    def build_parameters(cobj, phase):
        cobj.diffus_phase_comp_coeff_1 = Var(
            doc="Parameter 1 for liquid phase diffusivity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="1")
        cobj.diffus_phase_comp_coeff_2 = Var(
            doc="Parameter 2 for liquid phase diffusivity model", units=pyunits.K
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="2")
        cobj.diffus_phase_comp_coeff_3 = Var(
            doc="Parameter 3 for liquid phase diffusivity model",
            units=pyunits.m**3 / pyunits.mol,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="3")

    @staticmethod
    def return_expression(blk, p, j, T):
        cobj = blk.params.get_component(j)
        C_MEA = blk.conc_mol_comp["MEA"]
        return (
            exp(
                cobj.diffus_phase_comp_coeff_1
                + cobj.diffus_phase_comp_coeff_2 / T
                + cobj.diffus_phase_comp_coeff_3 * C_MEA
            )
            * pyunits.m**2
            / pyunits.s
        )


class DiffusIons:
    @staticmethod
    def build_parameters(cobj, phase):
        cobj.diffus_phase_comp_coeff_1 = Var(
            doc="Parameter 1 for liquid phase diffusivity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="1")
        cobj.diffus_phase_comp_coeff_2 = Var(
            doc="Parameter 2 for liquid phase diffusivity model", units=pyunits.K
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="2")
        cobj.diffus_phase_comp_coeff_3 = Var(
            doc="Parameter 3 for liquid phase diffusivity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="diffus_phase_comp_coeff", index="3")

    @staticmethod
    def return_expression(blk, p, j, T):
        cobj = blk.params.get_component(j)
        return (
            exp(
                cobj.diffus_phase_comp_coeff_1
                + cobj.diffus_phase_comp_coeff_2 / T
                + cobj.diffus_phase_comp_coeff_3
                * log(blk.visc_d_phase[p] / pyunits.get_units(blk.visc_d_phase[p]))
            )
            * pyunits.m**2
            / pyunits.s
        )


class DiffusNone:
    # placeholder method for components where diffusivity is not required
    @staticmethod
    def build_parameters(cobj, phase):
        pass

    @staticmethod
    def return_expression(blk, p, j, T):
        return Expression.Skip


class Viscosity:
    @staticmethod
    def build_parameters(pobj):
        pobj.visc_d_coeff_a = Var(
            doc="Parameter a for liquid phase viscosity model", units=pyunits.K**-1
        )
        set_param_from_config(pobj, param="visc_d_coeff", index="a")

        pobj.visc_d_coeff_b = Var(
            doc="Parameter b for liquid phase viscosity model", units=pyunits.K**-1
        )
        set_param_from_config(pobj, param="visc_d_coeff", index="b")

        pobj.visc_d_coeff_c = Var(
            doc="Parameter c for liquid phase viscosity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="visc_d_coeff", index="c")

        pobj.visc_d_coeff_d = Var(
            doc="Parameter d for liquid phase viscosity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="visc_d_coeff", index="d")

        pobj.visc_d_coeff_e = Var(
            doc="Parameter e for liquid phase viscosity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="visc_d_coeff", index="e")

        pobj.visc_d_coeff_f = Var(
            doc="Parameter f for liquid phase viscosity model", units=pyunits.K**-1
        )
        set_param_from_config(pobj, param="visc_d_coeff", index="f")

        pobj.visc_d_coeff_g = Var(
            doc="Parameter g for liquid phase viscosity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="visc_d_coeff", index="g")

    @staticmethod
    def return_expression(blk, phase):
        pobj = blk.params.get_phase(phase)

        r = (
            blk.mass_frac_phase_comp_apparent["Liq", "MEA"]
            / (
                blk.mass_frac_phase_comp_apparent["Liq", "MEA"]
                + blk.mass_frac_phase_comp_apparent["Liq", "H2O"]
            )
        ) * 100
        T = blk.temperature
        alpha = (
            blk.mole_frac_phase_comp_apparent["Liq", "CO2"]
            / blk.mole_frac_phase_comp_apparent["Liq", "MEA"]
        )
        mu_H2O = (
            1.002e-3
            * pyunits.Pa
            * pyunits.s
            * 10
            ** (
                (
                    1.3272
                    * (
                        293.15 * pyunits.K
                        - T
                        - 0.001053 * pyunits.K**-1 * (T - 293.15 * pyunits.K) ** 2
                    )
                )
                / (T - 168.15 * pyunits.K)
            )
        )
        a = pobj.visc_d_coeff_a
        b = pobj.visc_d_coeff_b
        c = pobj.visc_d_coeff_c
        d = pobj.visc_d_coeff_d
        e = pobj.visc_d_coeff_e
        f = pobj.visc_d_coeff_f
        g = pobj.visc_d_coeff_g

        # Model appears to be entirely empirical, and units are not obvious
        # Assume each part of expression is unitless and use units to match
        return mu_H2O * exp(
            r
            * (T * (a * r + b) + c * r + d)
            * (alpha * (e * r + f * T + g) + 1)
            / T**2
            * pyunits.K**2
        )


class ThermalCond:
    @staticmethod
    def build_parameters(pobj):
        pass

    @staticmethod
    def return_expression(blk, phase):
        Tb_MEA = 443 * pyunits.K
        Tc_MEA = 614.2 * pyunits.K
        T = blk.temperature
        Tr = T / Tc_MEA
        Tbr_MEA = Tb_MEA / Tc_MEA
        K_MEA = (
            1.1053152
            / (61.08**0.5)
            * (3 + 20 * (1 - Tr) ** (2 / 3))
            / (3 + 20 * (1 - Tbr_MEA) ** (2 / 3))
        )
        K_H2O = 0.6065 * (
            -1.48445
            + 4.12292 * T / (298.15 * pyunits.K)
            - 1.63866 * (T / (298.15 * pyunits.K)) ** 2
        )

        x = blk.mole_frac_phase_comp_apparent
        return (
            ((x["Liq", "H2O"] * K_H2O**-2 + x["Liq", "MEA"] * K_MEA**-2) ** -1)
            ** 0.5
            * pyunits.W
            / pyunits.m
            / pyunits.K
        )


class SurfTens:
    @staticmethod
    def build_parameters(pobj):
        pobj.surf_tens_H2O_coeff_1 = Var(
            doc="Parameter 1 for surface tension model of pure water",
            units=pyunits.N / pyunits.m,
        )
        set_param_from_config(pobj, param="surf_tens_H2O_coeff", index="1")
        pobj.surf_tens_H2O_coeff_2 = Var(
            doc="Parameter 2 for surface tension model of pure water",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_H2O_coeff", index="2")
        pobj.surf_tens_H2O_coeff_3 = Var(
            doc="Parameter 3 for surface tension model of pure water",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_H2O_coeff", index="3")
        pobj.surf_tens_H2O_coeff_4 = Var(
            doc="Parameter 4 for surface tension model of pure water",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_H2O_coeff", index="4")

        pobj.surf_tens_MEA_coeff_1 = Var(
            doc="Parameter 1 for surface tension model of pure MEA",
            units=pyunits.N / pyunits.m,
        )
        set_param_from_config(pobj, param="surf_tens_MEA_coeff", index="1")
        pobj.surf_tens_MEA_coeff_2 = Var(
            doc="Parameter 2 for surface tension model of pure MEA",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_MEA_coeff", index="2")
        pobj.surf_tens_MEA_coeff_3 = Var(
            doc="Parameter 3 for surface tension model of pure MEA",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_MEA_coeff", index="3")
        pobj.surf_tens_MEA_coeff_4 = Var(
            doc="Parameter 4 for surface tension model of pure MEA",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_MEA_coeff", index="4")

        pobj.surf_tens_CO2_coeff_1 = Var(
            doc="Parameter 1 for surface tension model of pure CO2",
            units=pyunits.N / pyunits.m,
        )
        set_param_from_config(pobj, param="surf_tens_CO2_coeff", index="1")
        pobj.surf_tens_CO2_coeff_2 = Var(
            doc="Parameter 2 for surface tension model of pure CO2",
            units=pyunits.N / pyunits.m,
        )
        set_param_from_config(pobj, param="surf_tens_CO2_coeff", index="2")
        pobj.surf_tens_CO2_coeff_3 = Var(
            doc="Parameter 3 for surface tension model of pure CO2",
            units=pyunits.N / pyunits.m,
        )
        set_param_from_config(pobj, param="surf_tens_CO2_coeff", index="3")
        pobj.surf_tens_CO2_coeff_4 = Var(
            doc="Parameter 4 for surface tension model of pure CO2",
            units=pyunits.N / pyunits.m / pyunits.K,
        )
        set_param_from_config(pobj, param="surf_tens_CO2_coeff", index="4")
        pobj.surf_tens_CO2_coeff_5 = Var(
            doc="Parameter 5 for surface tension model of pure CO2",
            units=pyunits.N / pyunits.m / pyunits.K,
        )
        set_param_from_config(pobj, param="surf_tens_CO2_coeff", index="5")
        pobj.surf_tens_CO2_coeff_6 = Var(
            doc="Parameter 6 for surface tension model of pure CO2",
            units=pyunits.N / pyunits.m / pyunits.K,
        )
        set_param_from_config(pobj, param="surf_tens_CO2_coeff", index="6")

        pobj.surf_tens_F_coeff_a = Var(
            doc="Parameter Fa for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="a")
        pobj.surf_tens_F_coeff_b = Var(
            doc="Parameter Fb for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="b")
        pobj.surf_tens_F_coeff_c = Var(
            doc="Parameter Fc for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="c")
        pobj.surf_tens_F_coeff_d = Var(
            doc="Parameter Fd for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="d")
        pobj.surf_tens_F_coeff_e = Var(
            doc="Parameter Fe for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="e")
        pobj.surf_tens_F_coeff_f = Var(
            doc="Parameter Ff for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="f")
        pobj.surf_tens_F_coeff_g = Var(
            doc="Parameter Fg for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="g")
        pobj.surf_tens_F_coeff_h = Var(
            doc="Parameter Fh for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="h")
        pobj.surf_tens_F_coeff_i = Var(
            doc="Parameter Fi for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="i")
        pobj.surf_tens_F_coeff_j = Var(
            doc="Parameter Fj for liquid phase surface tension model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="surf_tens_F_coeff", index="j")

    @staticmethod
    def return_expression(blk, phase):
        pobj = blk.params.get_phase(phase)

        T = blk.temperature
        r = blk.mass_frac_phase_comp["Liq", "MEA"] / (
            blk.mass_frac_phase_comp["Liq", "MEA"]
            + blk.mass_frac_phase_comp["Liq", "H2O"]
        )

        alpha = (
            blk.mole_frac_phase_comp_apparent["Liq", "CO2"]
            / blk.mole_frac_phase_comp_apparent["Liq", "MEA"]
        )

        Tc_h2o = blk.params.H2O.temperature_crit
        Tc_mea = blk.params.MEA.temperature_crit

        sigma_h2o = pobj.surf_tens_H2O_coeff_1 * (1 - T / Tc_h2o) ** (
            pobj.surf_tens_H2O_coeff_2
            + pobj.surf_tens_H2O_coeff_3 * T / Tc_h2o
            + pobj.surf_tens_H2O_coeff_4 * (T / Tc_h2o) ** 2
        )
        sigma_mea = pobj.surf_tens_MEA_coeff_1 * (1 - T / Tc_mea) ** (
            pobj.surf_tens_MEA_coeff_2
            + pobj.surf_tens_MEA_coeff_3 * T / Tc_mea
            + pobj.surf_tens_MEA_coeff_4 * (T / Tc_mea) ** 2
        )
        sigma_co2 = (
            pobj.surf_tens_CO2_coeff_1 * r**2
            + pobj.surf_tens_CO2_coeff_2 * r
            + pobj.surf_tens_CO2_coeff_3
            + T
            * (
                pobj.surf_tens_CO2_coeff_4 * r**2
                + pobj.surf_tens_CO2_coeff_5 * r
                + pobj.surf_tens_CO2_coeff_6
            )
        )

        Fa = pobj.surf_tens_F_coeff_a
        Fb = pobj.surf_tens_F_coeff_b
        Fc = pobj.surf_tens_F_coeff_c
        Fd = pobj.surf_tens_F_coeff_d
        Fe = pobj.surf_tens_F_coeff_e
        Ff = pobj.surf_tens_F_coeff_f
        Fg = pobj.surf_tens_F_coeff_g
        Fh = pobj.surf_tens_F_coeff_h
        Fi = pobj.surf_tens_F_coeff_i
        Fj = pobj.surf_tens_F_coeff_j

        return (
            sigma_h2o
            + (sigma_co2 - sigma_h2o)
            * blk.mole_frac_comp["CO2"]
            * (Fa + Fb * alpha + Fc * alpha**2 + Fd * r + Fe * r**2)
            + (sigma_mea - sigma_h2o)
            * (Ff + Fg * alpha + Fh * alpha**2 + Fi * r + Fj * r**2)
            * blk.mole_frac_comp["MEA"]
        )


# -----------------------------------------------------------------------------
# Equilibrium constant model
class k_eq:
    @staticmethod
    def build_parameters(rblock, config):
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

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        return exp(b.log_k_eq[r_idx]) * ((pyunits.m) ** 3 / pyunits.mol)

    @staticmethod
    def return_log_expression(b, rblock, r_idx, T):
        return b.log_k_eq[r_idx] == (
            rblock.k_eq_coeff_1
            + rblock.k_eq_coeff_2 / T
            + rblock.k_eq_coeff_3 * log(T / pyunits.K)
            + log(1e-3)
        )

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1


# -----------------------------------------------------------------------------
# Configuration dictionary for aqueous MEA solvent
configuration = {
    # Specifying components
    "components": {
        "H2O": {
            "type": Solvent,
            "cp_mol_liq_comp": CpMolSolvent,
            "diffus_phase_comp": {"Liq": DiffusNone},
            "enth_mol_liq_comp": EnthMolSolvent,
            "pressure_sat_comp": PressureSatSolvent,
            "vol_mol_liq_comp": VolMolSolvent,
            "parameter_data": {
                "mw": (0.01802, pyunits.kg / pyunits.mol),
                "cp_mass_liq_comp_coeff": {
                    "1": 4.2107,
                    "2": -1.696e-3,
                    "3": 2.568e-5,
                    "4": -1.095e-7,
                    "5": 3.038e-10,
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
                "temperature_crit": (647.13, pyunits.K),
            },
        },
        "MEA": {
            "type": Solvent,
            "cp_mol_liq_comp": CpMolSolvent,
            "diffus_phase_comp": {"Liq": DiffusMEA},
            "enth_mol_liq_comp": EnthMolSolvent,
            "pressure_sat_comp": PressureSatSolvent,
            "vol_mol_liq_comp": VolMolSolvent,
            "parameter_data": {
                "mw": (0.06108, pyunits.kg / pyunits.mol),
                "cp_mass_liq_comp_coeff": {
                    "1": 2.6161,
                    "2": 3.706e-3,
                    "3": 3.787e-6,
                    "4": 0.0,
                    "5": 0.0,
                },
                "dens_mol_liq_comp_coeff": {
                    "1": (-5.35162e-7, pyunits.g / pyunits.mL / pyunits.K**2),  # [2]
                    "2": (-4.51417e-4, pyunits.g / pyunits.mL / pyunits.K),
                    "3": (1.19451, pyunits.g / pyunits.mL),
                },
                "dh_vap": 58000,  # [3]
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
                "temperature_crit": (614.45, pyunits.K),
            },
        },
        "CO2": {
            "type": Solute,
            "cp_mol_liq_comp": CpMolCO2,
            "diffus_phase_comp": {"Liq": DiffusCO2},
            "enth_mol_liq_comp": EnthMolCO2,
            "henry_component": {
                "Liq": {
                    "method": N2OAnalogy,
                    "type": HenryType.Kpc,
                    "basis": StateIndex.true,
                }
            },
            "vol_mol_liq_comp": VolMolCO2,
            "parameter_data": {
                "mw": (0.04401, pyunits.kg / pyunits.mol),
                "dh_abs_co2": -84000,
                "lwm_coeff": {
                    "1": 1.70981,
                    "2": 0.03972,
                    "3": -4.3e-4,
                    "4": -2.20377,
                },
                "diffus_phase_comp_coeff": {
                    "1": 2.35e-6,
                    "2": 2.9837e-8,
                    "3": -9.7078e-9,
                    "4": -2119,
                    "5": -20.132,
                },
                "vol_mol_liq_comp_coeff": {
                    "a": (10.2074, pyunits.mL / pyunits.mol),  # [2]
                    "b": (-2.2642, pyunits.mL / pyunits.mol),
                    "c": (3.0059, pyunits.mL / pyunits.mol),
                    "d": (207, pyunits.mL / pyunits.mol),
                    "e": (-563.3701, pyunits.mL / pyunits.mol),
                },
            },
        },
        "MEA_+": {
            "type": Cation,
            "charge": +1,
            "diffus_phase_comp": {"Liq": DiffusIons},
            "parameter_data": {
                "diffus_phase_comp_coeff": {"1": -22.64, "2": -1000.0, "3": -0.7}
            },
        },
        "MEACOO_-": {
            "type": Anion,
            "charge": -1,
            "diffus_phase_comp": {"Liq": DiffusIons},
            "parameter_data": {
                "diffus_phase_comp_coeff": {"1": -22.64, "2": -1000.0, "3": -0.7}
            },
        },
        "HCO3_-": {
            "type": Anion,
            "charge": -1,
            "diffus_phase_comp": {"Liq": DiffusNone},
        },
    },
    # Specifying phases
    "phases": {
        "Liq": {
            "type": AqueousPhase,
            "equation_of_state": Ideal,
            "equation_of_state_options": {"property_basis": "apparent"},
            "surf_tens_phase": SurfTens,
            "therm_cond_phase": ThermalCond,
            "visc_d_phase": Viscosity,
            "parameter_data": {
                "surf_tens_H2O_coeff": {
                    "1": 0.18548,
                    "2": 2.717,
                    "3": -3.554,
                    "4": 2.047,
                },
                "surf_tens_MEA_coeff": {"1": 0.09945, "2": 1.067, "3": 0.0, "4": 0.0},
                "surf_tens_CO2_coeff": {
                    "1": -5.987,
                    "2": 3.7699,
                    "3": -0.43164,
                    "4": 0.018155,
                    "5": -0.01207,
                    "6": 0.002119,
                },
                "surf_tens_F_coeff": {
                    "a": 2.4558,
                    "b": -1.5311,
                    "c": 3.4994,
                    "d": -5.6398,
                    "e": 10.2109,
                    "f": 2.3122,
                    "g": 4.5608,
                    "h": -2.3924,
                    "i": 5.3324,
                    "j": -12.0494,
                },
                "visc_d_coeff": {
                    "a": -0.0838,
                    "b": 2.8817,
                    "c": 33.651,
                    "d": 1817.0,
                    "e": 0.00847,
                    "f": 0.0103,
                    "g": -2.3890,
                },
            },
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
        "temperature": (273.15, 298.15, 450, pyunits.K),
        "pressure": (5e4, 101325, 1e6, pyunits.Pa),
    },
    "state_components": StateIndex.apparent,
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    "inherent_reactions": {
        "carbamate": {
            "stoichiometry": {
                ("Liq", "MEA"): -2,
                ("Liq", "CO2"): -1,
                ("Liq", "MEA_+"): 1,
                ("Liq", "MEACOO_-"): 1,
            },
            "equilibrium_constant": k_eq,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {"k_eq_coeff": {"1": 233.4, "2": -3410, "3": -36.8}},
        },
        "bicarbonate": {
            "stoichiometry": {
                ("Liq", "MEA"): -1,
                ("Liq", "CO2"): -1,
                ("Liq", "H2O"): -1,
                ("Liq", "HCO3_-"): 1,
                ("Liq", "MEA_+"): 1,
            },
            "equilibrium_constant": k_eq,
            "equilibrium_form": log_power_law_equil,
            "concentration_form": ConcentrationForm.molarity,
            "parameter_data": {"k_eq_coeff": {"1": 176.72, "2": -2909, "3": -28.46}},
        },
    },
}
