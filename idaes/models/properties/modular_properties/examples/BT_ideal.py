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
Benzene-Toluene phase equilibrium package using ideal liquid and vapor.

Example property package using the Generic Property Package Framework.
This exmample shows how to set up a property package to do benzene-toluene
phase equilibrium in the generic framework using ideal liquid and vapor
assumptions along with methods drawn from the pre-built IDAES property
libraries.
"""
# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import LiquidPhase, VaporPhase, Component
import idaes.logger as idaeslog

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from idaes.models.properties.modular_properties.pure import Perrys
from idaes.models.properties.modular_properties.pure import RPP4


# Set up logger
_log = idaeslog.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for an ideal Benzene-Toluene system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] Perry's Chemical Engineers' Handbook 7th Ed.
# [3] Engineering Toolbox, https://www.engineeringtoolbox.com
#     Retrieved 1st December, 2019

configuration = {
    # Specifying components
    "components": {
        "benzene": {
            "type": Component,
            "elemental_composition": {"C": 6, "H": 6},
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (78.1136e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                "temperature_crit": (562.2, pyunits.K),  # [1]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (1.0162, pyunits.kmol * pyunits.m**-3),  # [2] pg. 2-98
                    "2": (0.2655, None),
                    "3": (562.16, pyunits.K),
                    "4": (0.28212, None),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (-3.392e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (4.739e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (-3.017e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (7.130e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (1.29e2, pyunits.J / pyunits.kmol / pyunits.K),  # [2]
                    "2": (-1.7e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (6.48e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (49.0e3, pyunits.J / pyunits.mol),  # [3]
                "enth_mol_form_vap_comp_ref": (82.9e3, pyunits.J / pyunits.mol),  # [3]
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
            "elemental_composition": {"C": 7, "H": 8},
            "dens_mol_liq_comp": Perrys,
            "enth_mol_liq_comp": Perrys,
            "enth_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP4,
            "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            "parameter_data": {
                "mw": (92.1405e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (41e5, pyunits.Pa),  # [1]
                "temperature_crit": (591.8, pyunits.K),  # [1]
                "dens_mol_liq_comp_coeff": {
                    "eqn_type": 1,
                    "1": (0.8488, pyunits.kmol * pyunits.m**-3),  # [2] pg. 2-98
                    "2": (0.26655, None),
                    "3": (591.8, pyunits.K),
                    "4": (0.2878, None),
                },
                "cp_mol_ig_comp_coeff": {
                    "A": (-2.435e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (5.125e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (-2.765e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (4.911e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "cp_mol_liq_comp_coeff": {
                    "1": (1.40e2, pyunits.J / pyunits.kmol / pyunits.K),  # [2]
                    "2": (-1.52e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (6.95e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (12.0e3, pyunits.J / pyunits.mol),  # [3]
                "enth_mol_form_vap_comp_ref": (50.1e3, pyunits.J / pyunits.mol),  # [3]
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
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 450, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (1e5, pyunits.Pa),
    "temperature_ref": (300, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": IdealBubbleDew,
}
