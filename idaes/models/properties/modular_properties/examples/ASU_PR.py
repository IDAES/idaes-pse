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
Air separation phase equilibrium package using Peng-Robinson EoS.

Example property package using the Generic Property Package Framework.
This example shows how to set up a property package to do air separation
phase equilibrium in the generic framework using Peng-Robinson equation
along with methods drawn from the pre-built IDAES property libraries.

The example includes two dictionaries.

1. The dictionary named configuration contains parameters obtained from
The Properties of Gases and Liquids (1987) 4th edition and NIST.

2. The dictionary named configuration_Dowling_2015 contains parameters used in
A framework for efficient large scale equation-oriented flowsheet optimization
(2015) Dowling. The parameters are extracted from Properties of Gases and
Liquids (1977) 3rd edition for Antoine's vapor equation and acentric factors
and converted values from the Properties of Gases and Liquids (1977)
3rd edition to j.
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
from idaes.models.properties.modular_properties.pure import NIST
from idaes.models.properties.modular_properties.pure import RPP3

# Set up logger
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Configuration dictionary for a Peng-Robinson Oxygen-Argon-Nitrogen system

# Data Sources:
# [1] The Properties of Gases and Liquids (1987)
#     4th edition, Chemical Engineering Series - Robert C. Reid
# [2] NIST, https://webbook.nist.gov/
#     Retrieved 16th August, 2020
# [3] The Properties of Gases and Liquids (1987)
#     3rd edition, Chemical Engineering Series - Robert C. Reid
#     Cp parameters where converted to j in Dowling 2015
# [4] A framework for efficient large scale equation-oriented flowsheet optimization (2015)
#     Computers and Chemical Engineering - Alexander W. Dowling

configuration = {
    # Specifying components
    "components": {
        "nitrogen": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (28.0135e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (34e5, pyunits.Pa),  # [1]
                "temperature_crit": (126.2, pyunits.K),  # [1]
                "omega": 0.037,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (3.115e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (-1.357e-2, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (2.680e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (-1.168e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),  # [2]
                "entr_mol_form_vap_comp_ref": (
                    191.61,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [2]
                "pressure_sat_comp_coeff": {
                    "A": (3.7362, None),  # [2]
                    "B": (264.651, pyunits.K),
                    "C": (-6.788, pyunits.K),
                },
            },
        },
        "argon": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (39.948e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (48.98e5, pyunits.Pa),  # [1]
                "temperature_crit": (150.86, pyunits.K),  # [1]
                "omega": 0.001,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (2.050e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                    "B": (0.0, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (0.0, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (0.0, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),  # [2]
                "entr_mol_form_vap_comp_ref": (
                    154.8,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [2]
                "pressure_sat_comp_coeff": {
                    "A": (3.29555, None),  # [2]
                    "B": (215.24, pyunits.K),
                    "C": (-22.233, pyunits.K),
                },
            },
        },
        "oxygen": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": NIST,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (31.999e-3, pyunits.kg / pyunits.mol),  # [1]
                "pressure_crit": (50.43e5, pyunits.Pa),  # [1]
                "temperature_crit": (154.58, pyunits.K),  # [1]
                "omega": 0.025,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (2.811e1, pyunits.J / pyunits.mol / pyunits.K),
                    "B": (-3.680e-6, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (1.746e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (-1.065e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),  # [2]
                "entr_mol_form_vap_comp_ref": (
                    205.152,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [2]
                "pressure_sat_comp_coeff": {
                    "A": (3.85845, None),  # [2]
                    "B": (325.675, pyunits.K),
                    "C": (-5.667, pyunits.K),
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
        "temperature": (10, 300, 350, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {
        "PR_kappa": {
            ("nitrogen", "nitrogen"): 0.000,
            ("nitrogen", "argon"): -0.26e-2,
            ("nitrogen", "oxygen"): -0.119e-1,
            ("argon", "nitrogen"): -0.26e-2,
            ("argon", "argon"): 0.000,
            ("argon", "oxygen"): 0.104e-1,
            ("oxygen", "nitrogen"): -0.119e-1,
            ("oxygen", "argon"): 0.104e-1,
            ("oxygen", "oxygen"): 0.000,
        }
    },
}

configuration_Dowling_2015 = {
    # Specifying components
    "components": {
        "nitrogen": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP3,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (28.0135e-3, pyunits.kg / pyunits.mol),  # [3]
                "pressure_crit": (33.943875e5, pyunits.Pa),  # [4]
                "temperature_crit": (126.2, pyunits.K),  # [4]
                "omega": 0.04,  # [3]
                "cp_mol_ig_comp_coeff": {
                    "A": (3.112896e1, pyunits.J / pyunits.mol / pyunits.K),  # [3]
                    "B": (-1.356e-2, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (2.6878e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (-1.167e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),  # [2]
                "entr_mol_form_vap_comp_ref": (
                    191.61,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [2]
                "pressure_sat_comp_coeff": {
                    "A": (14.9342, None),  # [3]
                    "B": (588.72, pyunits.K),
                    "C": (-6.60, pyunits.K),
                },
            },
        },
        "argon": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP3,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (39.948e-3, pyunits.kg / pyunits.mol),  # [3]
                "pressure_crit": (48.737325e5, pyunits.Pa),  # [4]
                "temperature_crit": (150.86, pyunits.K),  # [4]
                "omega": -0.004,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (2.0790296e1, pyunits.J / pyunits.mol / pyunits.K),  # [3]
                    "B": (-3.209e-05, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (5.163e-08, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (0.0, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),  # [3]
                "entr_mol_form_vap_comp_ref": (
                    154.8,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [3]
                "pressure_sat_comp_coeff": {
                    "A": (15.2330, None),  # [3]
                    "B": (700.51, pyunits.K),
                    "C": (-5.84, pyunits.K),
                },
            },
        },
        "oxygen": {
            "type": Component,
            "enth_mol_ig_comp": RPP4,
            "entr_mol_ig_comp": RPP4,
            "pressure_sat_comp": RPP3,
            "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
            "parameter_data": {
                "mw": (31.999e-3, pyunits.kg / pyunits.mol),  # [3]
                "pressure_crit": (50.45985e5, pyunits.Pa),  # [4]
                "temperature_crit": (154.58, pyunits.K),  # [4]
                "omega": 0.021,  # [1]
                "cp_mol_ig_comp_coeff": {
                    "A": (2.8087192e1, pyunits.J / pyunits.mol / pyunits.K),  # [3]
                    "B": (-3.678e-6, pyunits.J / pyunits.mol / pyunits.K**2),
                    "C": (1.745e-5, pyunits.J / pyunits.mol / pyunits.K**3),
                    "D": (-1.064e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                },
                "enth_mol_form_vap_comp_ref": (0.0, pyunits.J / pyunits.mol),  # [2]
                "entr_mol_form_vap_comp_ref": (
                    205.152,
                    pyunits.J / pyunits.mol / pyunits.K,
                ),  # [2]
                "pressure_sat_comp_coeff": {
                    "A": (15.4075, None),  # [3]
                    "B": (734.55, pyunits.K),
                    "C": (-6.45, pyunits.K),
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
        "temperature": (10, 300, 350, pyunits.K),
        "pressure": (5e4, 1e5, 1e7, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    # Defining phase equilibria
    "phases_in_equilibrium": [("Vap", "Liq")],
    "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
    "bubble_dew_method": LogBubbleDew,
    "parameter_data": {
        "PR_kappa": {
            ("nitrogen", "nitrogen"): 0.000,
            ("nitrogen", "argon"): -0.26e-2,
            ("nitrogen", "oxygen"): -0.119e-1,
            ("argon", "nitrogen"): -0.26e-2,
            ("argon", "argon"): 0.000,
            ("argon", "oxygen"): 0.104e-1,
            ("oxygen", "nitrogen"): -0.119e-1,
            ("oxygen", "argon"): 0.104e-1,
            ("oxygen", "oxygen"): 0.000,
        }
    },
}
