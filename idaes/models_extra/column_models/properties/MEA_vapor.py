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
Ideal Vapor Phase properties for CO@ absorption via aqueous MEA.

The following apparent species are used to represent the mixture:

    Carbon Dioxide (CO2)
    Water (H2O)
    Oxygen (O2)
    Nitrogen (N2)

Assumptions:
    * No heats of formation


References:
    [1] Hilliard thesis (1998)
    [2] Morgan et.al (2015)
"""
# Import Pyomo units
from pyomo.environ import sqrt, units as pyunits, Var

# Import IDAES cores
from idaes.core import VaporPhase, Component
from idaes.core.util.constants import Constants

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.base.generic_property import StateIndex
from idaes.models.properties.modular_properties.eos.ideal import Ideal

from idaes.core.util.misc import set_param_from_config
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


# -----------------------------------------------------------------------------
# Pure Component Property methods for vapor mixture
class Cp:
    # Method to calculate cp for solvents [1]
    @staticmethod
    def build_parameters(cobj):
        cobj.cp_mol_vap_comp_coeff_1 = Var(
            doc="Parameter 1 for vapor phase molar specific heat capacity",
            units=pyunits.dimensionless,
        )
        set_param_from_config(cobj, param="cp_mol_vap_comp_coeff", index="1")

        cobj.cp_mol_vap_comp_coeff_2 = Var(
            doc="Parameter 2 for vapor phase molar specific heat capacity",
            units=pyunits.K**-1,
        )
        set_param_from_config(cobj, param="cp_mol_vap_comp_coeff", index="2")

        cobj.cp_mol_vap_comp_coeff_3 = Var(
            doc="Parameter 3 for vapor phase molar specific heat capacity",
            units=pyunits.K**2,
        )
        set_param_from_config(cobj, param="cp_mol_vap_comp_coeff", index="3")

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific heat capacity
        return (
            cobj.cp_mol_vap_comp_coeff_1
            + cobj.cp_mol_vap_comp_coeff_2 * T
            + cobj.cp_mol_vap_comp_coeff_3 * T**-2
        ) * Constants.gas_constant


class EnthMol:
    # Method to calculate specific enthalpy of solvents
    @staticmethod
    def build_parameters(cobj):
        if not hasattr(cobj, "cp_mol_vap_comp_coeff_1"):
            Cp.build_parameters(cobj)

    @staticmethod
    def return_expression(b, cobj, T):
        # Specific enthalpy
        Tref = b.params.temperature_ref
        return (
            cobj.cp_mol_vap_comp_coeff_1 * (T - Tref)
            + 0.5 * cobj.cp_mol_vap_comp_coeff_2 * (T**2 - Tref**2)
            - cobj.cp_mol_vap_comp_coeff_3 * (T**-1 - Tref**-1)
        ) * Constants.gas_constant


# -----------------------------------------------------------------------------
# Transport property models
class Diffus:
    @staticmethod
    def build_parameters(cobj, phase):
        pass

    @staticmethod
    def return_expression(blk, p, i, T):
        # Diffusion Coefficient(binary) constants -
        # Diffusion volumes in Fuller-Schettler-Giddings correlation
        # for estimating binary diffusivities of components in vapor phase
        # Reference: Table 3.1 pp 71 Seader Henley (2006)

        pobj = blk.params.get_phase(p)

        if not hasattr(pobj, "diffus_binary_param"):
            pobj.diffus_binary_param = Var(
                blk.component_list,
                initialize=pobj.config.parameter_data["diffus_binary_param"],
                units=pyunits.m**3,
                doc="Diffusion volume parameter for binary diffusivity",
            )
            pobj.diffus_binary_param.fix()

        binary_set = []
        for ii in blk.component_list:
            for jj in blk.component_list:
                if ii != jj and (jj, ii) not in binary_set:
                    binary_set.append((ii, jj))

        def diffus_binary(i, j):
            return (
                1.013e-2
                * (pyunits.kilogram / pyunits.m / pyunits.s**2)
                * pyunits.m**4
                * ((pyunits.kilogram / pyunits.mol) ** 0.5)
                / pyunits.s
                / pyunits.K**1.75
                * T**1.75
                / blk.pressure
                * sqrt(1e-3 * (1 / blk.mw_comp[i] + 1 / blk.mw_comp[j]))
                / (
                    pobj.diffus_binary_param[i] ** (1 / 3)
                    + pobj.diffus_binary_param[j] ** (1 / 3)
                )
                ** 2
            )

        return (1 - blk.mole_frac_comp[i]) / (
            sum(
                blk.mole_frac_comp[j] / diffus_binary(i, j)
                for j in blk.component_list
                if (i, j) in binary_set
            )
            + sum(
                blk.mole_frac_comp[j] / diffus_binary(j, i)
                for j in blk.component_list
                if (j, i) in binary_set
            )
        )


def visc_d_comp(blk, pobj, i):
    """
    Dynamic viscosity of vapor components
    Sutherland formula for N2 & O2
    DIPPR method for H2O & CO2
    """
    # DIPPR components
    if i == "H2O":
        return (
            (
                pobj.visc_d_h2o_coeff_1
                * blk.temperature**pobj.visc_d_h2o_coeff_2
                / (1 + pobj.visc_d_h2o_coeff_3 / blk.temperature)
            )
            * pyunits.Pa
            * pyunits.s
            * (1 / (pyunits.K) ** (pobj.visc_d_h2o_coeff_2))
        )
    elif i == "CO2":
        return (
            (
                pobj.visc_d_co2_coeff_1
                * blk.temperature**pobj.visc_d_co2_coeff_2
                / (1 + pobj.visc_d_co2_coeff_3 / blk.temperature)
            )
            * pyunits.Pa
            * pyunits.s
            * (1 / (pyunits.K) ** (pobj.visc_d_co2_coeff_2))
        )
    elif i == "N2":
        return (
            pobj.visc_d_n2_coeff_1
            * (pobj.visc_d_n2_coeff_2 + pobj.visc_d_n2_coeff_3)
            / (blk.temperature + pobj.visc_d_n2_coeff_3)
            * (blk.temperature / pobj.visc_d_n2_coeff_2) ** 1.5
        )
    elif i == "O2":
        return (
            pobj.visc_d_o2_coeff_1
            * (pobj.visc_d_o2_coeff_2 + pobj.visc_d_o2_coeff_3)
            / (blk.temperature + pobj.visc_d_o2_coeff_3)
            * (blk.temperature / pobj.visc_d_o2_coeff_2) ** 1.5
        )
    else:
        # Overrun error
        raise RuntimeError(f"Unrecognised component{i}")


class ThermalCond:
    @staticmethod
    def build_parameters(pobj):
        # Need viscosity parameters for thermal conductivity
        if not hasattr(pobj, "visc_d_h2o_coeff_1"):
            Viscosity.build_parameters(pobj)

        pobj.therm_cond_h2o_coeff_1 = Var(
            doc="Parameter 1 for H2O thermal conductivity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="therm_cond_h2o_coeff", index="1")
        pobj.therm_cond_h2o_coeff_2 = Var(
            doc="Parameter 2 for H2O thermal conductivity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="therm_cond_h2o_coeff", index="2")
        pobj.therm_cond_h2o_coeff_3 = Var(
            doc="Parameter 3 for H2O thermal conductivity model", units=pyunits.K
        )
        set_param_from_config(pobj, param="therm_cond_h2o_coeff", index="3")
        pobj.therm_cond_h2o_coeff_4 = Var(
            doc="Parameter 4 for H2O thermal conductivity model", units=pyunits.K**2
        )
        set_param_from_config(pobj, param="therm_cond_h2o_coeff", index="4")

        pobj.therm_cond_co2_coeff_1 = Var(
            doc="Parameter 1 for CO2 thermal conductivity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="therm_cond_co2_coeff", index="1")
        pobj.therm_cond_co2_coeff_2 = Var(
            doc="Parameter 2 for CO2 thermal conductivity model",
            units=pyunits.dimensionless,
        )
        set_param_from_config(pobj, param="therm_cond_co2_coeff", index="2")
        pobj.therm_cond_co2_coeff_3 = Var(
            doc="Parameter 3 for CO2 thermal conductivity model", units=pyunits.K
        )
        set_param_from_config(pobj, param="therm_cond_co2_coeff", index="3")
        pobj.therm_cond_co2_coeff_4 = Var(
            doc="Parameter 4 for CO2 thermal conductivity model", units=pyunits.K**2
        )
        set_param_from_config(pobj, param="therm_cond_co2_coeff", index="4")

        if "N2" in pobj.parent_block().component_list:
            pobj.therm_cond_n2_coeff_1 = Var(
                doc="Parameter 1 for N2 thermal conductivity model",
                units=pyunits.dimensionless,
            )
            set_param_from_config(pobj, param="therm_cond_n2_coeff", index="1")
            pobj.therm_cond_n2_coeff_2 = Var(
                doc="Parameter 2 for N2 thermal conductivity model",
                units=pyunits.dimensionless,
            )
            set_param_from_config(pobj, param="therm_cond_n2_coeff", index="2")
            pobj.therm_cond_n2_coeff_3 = Var(
                doc="Parameter 3 for N2 thermal conductivity model", units=pyunits.K
            )
            set_param_from_config(pobj, param="therm_cond_n2_coeff", index="3")
            pobj.therm_cond_n2_coeff_4 = Var(
                doc="Parameter 4 for N2 thermal conductivity model",
                units=pyunits.K**2,
            )
            set_param_from_config(pobj, param="therm_cond_n2_coeff", index="4")

            pobj.therm_cond_o2_coeff_1 = Var(
                doc="Parameter 1 for O2 thermal conductivity model",
                units=pyunits.dimensionless,
            )
            set_param_from_config(pobj, param="therm_cond_o2_coeff", index="1")
            pobj.therm_cond_o2_coeff_2 = Var(
                doc="Parameter 2 for O2 thermal conductivity model",
                units=pyunits.dimensionless,
            )
            set_param_from_config(pobj, param="therm_cond_o2_coeff", index="2")
            pobj.therm_cond_o2_coeff_3 = Var(
                doc="Parameter 3 for O2 thermal conductivity model", units=pyunits.K
            )
            set_param_from_config(pobj, param="therm_cond_o2_coeff", index="3")
            pobj.therm_cond_o2_coeff_4 = Var(
                doc="Parameter 4 for O2 thermal conductivity model",
                units=pyunits.K**2,
            )
            set_param_from_config(pobj, param="therm_cond_o2_coeff", index="4")

    @staticmethod
    def return_expression(blk, phase):
        pobj = blk.params.get_phase(phase)

        def therm_cond_comp(j):
            p1 = getattr(pobj, "therm_cond_" + j.lower() + "_coeff_1")
            p2 = getattr(pobj, "therm_cond_" + j.lower() + "_coeff_2")
            p3 = getattr(pobj, "therm_cond_" + j.lower() + "_coeff_3")
            p4 = getattr(pobj, "therm_cond_" + j.lower() + "_coeff_4")
            return (
                p1
                * (1 / (pyunits.K) ** (p2))
                * (pyunits.joule / (pyunits.m * pyunits.K * pyunits.s))
                * blk.temperature**p2
                / ((1 + p3 / blk.temperature) + (p4 / blk.temperature**2))
            )

        """
        Thermal conductivity of vapor phase
        Wassiljewa-Mason-Saxena mixing rule(low pressure)
        """
        k_vap = 0
        for i in blk.component_list:
            sumij = 0
            for j in blk.component_list:
                Aij = (
                    1
                    + (visc_d_comp(blk, pobj, i) / visc_d_comp(blk, pobj, j)) ** 0.5
                    * (blk.mw_comp[j] / blk.mw_comp[i]) ** 0.25
                ) ** 2 * (8 * (1 + blk.mw_comp[i] / blk.mw_comp[j])) ** -0.5
                sumij += blk.mole_frac_comp[j] * Aij
            k_vap += blk.mole_frac_comp[i] * therm_cond_comp(i) / sumij

        return k_vap


class Viscosity:
    @staticmethod
    def build_parameters(pobj):
        # Viscosity parameters are required for thermal conductivity, and
        # have likely already been built by the time this is triggered
        # for viscosity.
        # To avoid implict replacement, check to see if parameters already
        # exist
        if hasattr(pobj, "visc_d_h2o_coeff_1"):
            return None

        # Viscosity constants
        # CO2 & H2O calculated from:
        # Perry and Green Handbook; McGraw Hill, 8th edition 2008
        # C1*T^(C2)/(1+C3/T)
        # O2 & N2 Calculated from Sutherland Formula
        # C1*(C2 + C3)/(T+C3)*(T/C2)^1.5
        # constants C1, C2, C3 in Sutherlands' Formula are:
        # C1 = vis_d_ref
        # C2 = temperature_ref
        # C3 = sutherland constant

        pobj.visc_d_h2o_coeff_1 = Var(
            doc="Parameter 1 for H2O viscosity model", units=pyunits.dimensionless
        )
        set_param_from_config(pobj, param="visc_d_h2o_coeff", index="1")
        pobj.visc_d_h2o_coeff_2 = Var(
            doc="Parameter 2 for H2O viscosity model", units=pyunits.dimensionless
        )
        set_param_from_config(pobj, param="visc_d_h2o_coeff", index="2")
        pobj.visc_d_h2o_coeff_3 = Var(
            doc="Parameter 3 for H2O viscosity model", units=pyunits.K
        )
        set_param_from_config(pobj, param="visc_d_h2o_coeff", index="3")

        pobj.visc_d_co2_coeff_1 = Var(
            doc="Parameter 1 for CO2 viscosity model", units=pyunits.dimensionless
        )
        set_param_from_config(pobj, param="visc_d_co2_coeff", index="1")
        pobj.visc_d_co2_coeff_2 = Var(
            doc="Parameter 2 for CO2 viscosity model", units=pyunits.dimensionless
        )
        set_param_from_config(pobj, param="visc_d_co2_coeff", index="2")
        pobj.visc_d_co2_coeff_3 = Var(
            doc="Parameter 3 for CO2 viscosity model", units=pyunits.K
        )
        set_param_from_config(pobj, param="visc_d_co2_coeff", index="3")

        if "N2" in pobj.parent_block().component_list:
            pobj.visc_d_n2_coeff_1 = Var(
                doc="Parameter 1 for N2 viscosity model", units=pyunits.Pa * pyunits.s
            )
            set_param_from_config(pobj, param="visc_d_n2_coeff", index="1")
            pobj.visc_d_n2_coeff_2 = Var(
                doc="Parameter 2 for N2 viscosity model", units=pyunits.K
            )
            set_param_from_config(pobj, param="visc_d_n2_coeff", index="2")
            pobj.visc_d_n2_coeff_3 = Var(
                doc="Parameter 3 for N2 viscosity model", units=pyunits.K
            )
            set_param_from_config(pobj, param="visc_d_n2_coeff", index="3")

            pobj.visc_d_o2_coeff_1 = Var(
                doc="Parameter 1 for O2 viscosity model", units=pyunits.Pa * pyunits.s
            )
            set_param_from_config(pobj, param="visc_d_o2_coeff", index="1")
            pobj.visc_d_o2_coeff_2 = Var(
                doc="Parameter 2 for O2 viscosity model", units=pyunits.K
            )
            set_param_from_config(pobj, param="visc_d_o2_coeff", index="2")
            pobj.visc_d_o2_coeff_3 = Var(
                doc="Parameter 3 for O2 viscosity model", units=pyunits.K
            )
            set_param_from_config(pobj, param="visc_d_o2_coeff", index="3")

    @staticmethod
    def return_expression(blk, phase):
        pobj = blk.params.get_phase(phase)

        theta_ij = {}
        o = dict()
        for (i, j) in enumerate(blk.component_list, 1):
            o[i] = j

        for i in range(1, len(blk.component_list)):
            for j in range(i + 1, len(blk.component_list) + 1):
                theta_ij[o[i], o[j]] = (
                    1
                    + 2
                    * sqrt(visc_d_comp(blk, pobj, o[i]) / visc_d_comp(blk, pobj, o[j]))
                    * (blk.mw_comp[o[j]] / blk.mw_comp[o[i]]) ** 0.25
                    + visc_d_comp(blk, pobj, o[i])
                    / visc_d_comp(blk, pobj, o[j])
                    * (blk.mw_comp[o[j]] / blk.mw_comp[o[i]]) ** 0.5
                ) / (8 + 8 * blk.mw_comp[o[i]] / blk.mw_comp[o[j]]) ** 0.5

                theta_ij[o[j], o[i]] = (
                    visc_d_comp(blk, pobj, o[j])
                    / visc_d_comp(blk, pobj, o[i])
                    * blk.mw_comp[o[i]]
                    / blk.mw_comp[o[j]]
                    * theta_ij[o[i], o[j]]
                )

        for i in blk.component_list:
            for j in blk.component_list:
                if i == j:
                    theta_ij[i, j] = 1

        return sum(
            blk.mole_frac_comp[i]
            * visc_d_comp(blk, pobj, i)
            / sum(blk.mole_frac_comp[j] * theta_ij[i, j] for j in blk.component_list)
            for i in blk.component_list
        )


# -----------------------------------------------------------------------------
# Configuration dictionary for aqueous MEA solvent
flue_gas = {
    # Specifying components
    "components": {
        "CO2": {
            "type": Component,
            "cp_mol_ig_comp": Cp,
            "diffus_phase_comp": {"Vap": Diffus},
            "enth_mol_ig_comp": EnthMol,
            "parameter_data": {
                "mw": (0.04401, pyunits.kg / pyunits.mol),
                "cp_mol_vap_comp_coeff": {"1": 5.457, "2": 1.045e-3, "3": -1.157e5},
            },
        },
        "H2O": {
            "type": Component,
            "cp_mol_ig_comp": Cp,
            "diffus_phase_comp": {"Vap": Diffus},
            "enth_mol_ig_comp": EnthMol,
            "parameter_data": {
                "mw": (0.01802, pyunits.kg / pyunits.mol),
                "cp_mol_vap_comp_coeff": {"1": 3.47, "2": 1.45e-3, "3": 0.121e5},
            },
        },
        "N2": {
            "type": Component,
            "cp_mol_ig_comp": Cp,
            "diffus_phase_comp": {"Vap": Diffus},
            "enth_mol_ig_comp": EnthMol,
            "parameter_data": {
                "mw": (0.02801, pyunits.kg / pyunits.mol),
                "cp_mol_vap_comp_coeff": {"1": 3.28, "2": 0.593e-3, "3": 0.04e5},
            },
        },
        "O2": {
            "type": Component,
            "cp_mol_ig_comp": Cp,
            "diffus_phase_comp": {"Vap": Diffus},
            "enth_mol_ig_comp": EnthMol,
            "parameter_data": {
                "mw": (0.032, pyunits.kg / pyunits.mol),
                "cp_mol_vap_comp_coeff": {"1": 3.639, "2": 0.506e-3, "3": -0.227e5},
            },
        },
    },
    # Specifying phases
    "phases": {
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Ideal,
            "therm_cond_phase": ThermalCond,
            "visc_d_phase": Viscosity,
            "parameter_data": {
                "diffus_binary_param": {
                    "N2": 18.5,
                    "CO2": 26.7,
                    "H2O": 13.1,
                    "O2": 16.3,
                },
                "therm_cond_co2_coeff": {
                    "1": 3.69,
                    "2": -0.3838,
                    "3": 964,
                    "4": 1.86e6,
                },
                "therm_cond_h2o_coeff": {"1": 6.204e-6, "2": 1.3973, "3": 0.0, "4": 0},
                "therm_cond_n2_coeff": {
                    "1": 0.000331,
                    "2": 0.7722,
                    "3": 16.323,
                    "4": 373.72,
                },
                "therm_cond_o2_coeff": {"1": 0.00045, "2": 0.7456, "3": 56.699, "4": 0},
                "visc_d_co2_coeff": {"1": 2.148e-6, "2": 0.46, "3": 290},
                "visc_d_h2o_coeff": {"1": 1.7096e-8, "2": 1.1146, "3": 0.0},
                "visc_d_n2_coeff": {"1": 0.01781e-3, "2": 300.55, "3": 111},
                "visc_d_o2_coeff": {"1": 0.02018e-3, "2": 292.25, "3": 127},
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
}


wet_co2 = {
    # Specifying components
    "components": {
        "CO2": {
            "type": Component,
            "cp_mol_ig_comp": Cp,
            "diffus_phase_comp": {"Vap": Diffus},
            "enth_mol_ig_comp": EnthMol,
            "parameter_data": {
                "mw": (0.04401, pyunits.kg / pyunits.mol),
                "cp_mol_vap_comp_coeff": {"1": 5.457, "2": 1.045e-3, "3": -1.157e5},
            },
        },
        "H2O": {
            "type": Component,
            "cp_mol_ig_comp": Cp,
            "diffus_phase_comp": {"Vap": Diffus},
            "enth_mol_ig_comp": EnthMol,
            "parameter_data": {
                "mw": (0.01802, pyunits.kg / pyunits.mol),
                "cp_mol_vap_comp_coeff": {"1": 3.47, "2": 1.45e-3, "3": 0.121e5},
            },
        },
    },
    # Specifying phases
    "phases": {
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Ideal,
            "therm_cond_phase": ThermalCond,
            "visc_d_phase": Viscosity,
            "parameter_data": {
                "diffus_binary_param": {"CO2": 26.7, "H2O": 13.1},
                "therm_cond_co2_coeff": {
                    "1": 3.69,
                    "2": -0.3838,
                    "3": 964,
                    "4": 1.86e6,
                },
                "therm_cond_h2o_coeff": {"1": 6.204e-6, "2": 1.3973, "3": 0.0, "4": 0},
                "visc_d_co2_coeff": {"1": 2.148e-6, "2": 0.46, "3": 290},
                "visc_d_h2o_coeff": {"1": 1.7096e-8, "2": 1.1146, "3": 0.0},
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
}
