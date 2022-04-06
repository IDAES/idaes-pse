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
Benzene-Toluene phase equilibrium package for liquid and vapor using the
Peng-Robinson equation of state.

Example property package using the Generic Property Package Framework.
This exmample shows how to set up a property package to do benzene-toluene
phase equilibrium in the generic framework using the Peng-Robinson equation of
state for liquid and vapor phases along with methods drawn from the pre-built
IDAES property libraries.
"""
# Import Python libraries
import logging

# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
from idaes.models.properties.modular_properties.pure import RPP4


# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

configuration = {
    # Specifying components
    "components": {
        "benzene": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (78.1136e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                "temperature_crit": (562.2, pyunits.K),  # [1]
                "omega": 0.212,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (-3.392e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (4.739e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (-3.017e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (7.130e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (82.9e3, pyunits.J / pyunits.mol),  # [3]
                "entr_mol_form_vap_comp_ref": (
                    -269,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [3]
                "pressure_sat_comp_coeff": {
                    "A": (-6.98273, None),  # [1]
                    "B": (1.33213, None),
                    "C": (-2.62863, None),
                    "D": (-3.33399, None),
                },
            },
        },
        "toluene": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (92.1405e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (41e5, pyunits.Pa),  # [1]
                "temperature_crit": (591.8, pyunits.K),  # [1]
                "omega": 0.263,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (-2.435e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (5.125e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (-2.765e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (4.911e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (50.1e3, pyunits.J / pyunits.mol),  # [3]
                "entr_mol_form_vap_comp_ref": (
                    -321,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [3]
                "pressure_sat_comp_coeff": {
                    "A": (-7.28607, None),  # [1]
                    "B": (1.38091, None),
                    "C": (-2.83433, None),
                    "D": (-2.79168, None),
                },
            },
        },
    },
    # Specifying phases
    "phases": {
        "Liq": {
            "type": LiquidPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
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
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 500, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {
        "PR_kappa": {
            ("benzene", "benzene"): 0.000,
            ("benzene", "toluene"): 0.000,
            ("toluene", "benzene"): 0.000,
            ("toluene", "toluene"): 0.000,
        }
    },
}
