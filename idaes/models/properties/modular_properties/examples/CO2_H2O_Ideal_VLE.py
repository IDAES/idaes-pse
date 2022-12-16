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
Property package for MEA Carbon Capture System Condenser

CO2-H2O VLE system
Ideal VLE with non-condensables
"""
# Import Python libraries
import logging

from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component, PhaseType as PT

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity

from idaes.models.properties.modular_properties.pure.Perrys import Perrys
from idaes.models.properties.modular_properties.pure.NIST import NIST

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal CO2-H2O system
# Assumptions
# 1. Ideal gas-liquid system
# 2. CO2 conc. negligible in liquid phase, at equilibrium
# 3. Properties are applicable at standard pressure (1 atm)

# Reference Temperature : 298.15 K
# Reference Pressure : 101325 Pa

# Data Sources:

# [1] NIST Webbook, https://webbook.nist.gov/
#     Retrieved 27th November, 2020. Converted from bar to Pa
# [2] Perry's Chemical Engineers' Handbook 7th Ed.

# Initial temperature and pressure reference:
# Measurement and Modeling of the Phase Behavior of the
# (Carbon Dioxide + Water) Mixture at Temperatures From 298.15 K to 448.15
# KShu-Xin Hou, Geoffrey C. Maitland, and J. P. Martin Trusler*
# Qatar Carbonates and Carbon Storage Research Centre
# Department of Chemical Engineering Imperial College London,
# South Kensington Campus, London SW7 2AZ.  U.K

configuration = {
    # Specifying components
    "components": {
        "H2O": {
            "type": Component,
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": NIST,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (18.0153e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (220.64e5, pyunits.Pa),  # [1]
                "temperature_crit": (647, pyunits.K),  # [1]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (
                        5.459,
                        pyunits.kmol * pyunits.m**-3,
                    ),  # [2] pg. 2-98, temperature range 273.16 K - 333.15 K
                    "2": (0.30542, None),
                    "3": (647.13, pyunits.K),
                    "4": (0.081, None),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (
                        30.09200,
                        pyunits.J / pyunits.mol / pyunits.K,
                    ),  # [1] temperature range 500 K- 1700 K
                    "B": (
                        6.832514,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-1,
                    ),
                    "C": (
                        6.793435,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-2,
                    ),
                    "D": (
                        -2.534480,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-3,
                    ),
                    "E": (
                        0.082139,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**2,
                    ),
                    "F": (-250.8810, pyunits.kJ / pyunits.mol),
                    "G": (223.3967, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (
                        2.7637e5,
                        pyunits.J / pyunits.kmol / pyunits.K,
                    ),  # [2] pg 2-174, temperature range 273.16 K - 533.15 K
                    "2": (-2.0901e3, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (8.125, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (-1.4116e-2, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (9.3701e-6, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (
                    -285.83e3,
                    pyunits.J / pyunits.mol,
                ),  # [1]
                "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [1]
                "pressure_sat_comp_coeff": {
                    "A": (4.6543, None),  # [1], temperature range 255.9 K - 373 K
                    "B": (1435.264, pyunits.K),
                    "C": (-64.848, pyunits.K),
                },
            },
        },
        "CO2": {
            "type": Component,
            "valid_phase_types": PT.vaporPhase,
            "enth_mol_ig_comp": NIST,
            "parameter_data": {
                "mw": (44.0095e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (73.825e5, pyunits.Pa),  # [1]
                "temperature_crit": (304.23, pyunits.K),  # [1]
                "cp_mol_ig_comp_coeff": {  # [1], temperature range 298 K - 1200 K
                    "A": (24.99735, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (
                        55.18696,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-1,
                    ),
                    "C": (
                        -33.69137,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-2,
                    ),
                    "D": (
                        7.948387,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**-3,
                    ),
                    "E": (
                        -0.136638,
                        pyunits.J
                        * pyunits.mol**-1
                        * pyunits.K**-1
                        * pyunits.kiloK**2,
                    ),
                    "F": (-403.6075, pyunits.kJ / pyunits.mol),
                    "G": (228.2431, pyunits.J / pyunits.mol / pyunits.K),
                    "H": (0, pyunits.kJ / pyunits.mol),
                },
                "enth_mol_form_vap_comp_ref": (0, pyunits.J / pyunits.mol),  # [1]
            },
        },
    },
    # Specifying phases
    "phases": {
        "Liq": {"type": LiquidPhase, "equation_of_state": Ideal},
        "Vap": {"type": VaporPhase, "equation_of_state": Ideal},
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
        "flow_mol": (0, 10, 20, pyunits.mol / pyunits.s),
        "temperature": (273.15, 323.15, 1000, pyunits.K),
        "pressure": (5e4, 108900, 1e7, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}
