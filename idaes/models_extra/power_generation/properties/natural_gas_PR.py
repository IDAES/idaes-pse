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
Natural gas property package for the vapor phase using Peng-Robinson equation
of state.
"""
# Import Python libraries
import logging
import copy
import enum

# Import Pyomo units
import pyomo.environ as pyo
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import VaporPhase, LiquidPhase, Component, PhaseType

from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.phase_equil.forms import (
    log_fugacity,
)
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.pure import NIST, RPP4, RPP5

from idaes.models.properties.modular_properties.base.generic_reaction import (
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_constant import arrhenius
from idaes.models.properties.modular_properties.reactions.rate_forms import (
    power_law_rate,
)

from idaes.core.util.exceptions import ConfigurationError

# Set up logger
_log = logging.getLogger(__name__)


class EosType(enum.Enum):
    PR = 1
    IDEAL = 2


# Property Sources

# Source: NIST webbook
# Properties: Heat capacity coefficients for all species except ethane,
# propane, and butane. Reference enthalpies and entropies for all species.

# Source: The Properties of Gases and Liquids (1987)
# 4th edition, Chemical Engineering Series - Robert C. Reid
# Properties: Critical temperatures and pressures. Omega.
# Heat capacity coefficients for ethane, propane, and butane.

_phase_dicts_pr = {
    "Vap": {
        "type": VaporPhase,
        "equation_of_state": Cubic,
        "equation_of_state_options": {"type": CubicType.PR},
    },
    "Liq": {
        "type": LiquidPhase,
        "equation_of_state": Cubic,
        "equation_of_state_options": {"type": CubicType.PR},
    },
}

_phase_dicts_ideal = {
    "Vap": {
        "type": VaporPhase,
        "equation_of_state": Ideal,
    },
}

_component_params = {
    "H2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"H": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0020159, pyunits.kg / pyunits.mol),
            "pressure_crit": (13e5, pyunits.Pa),
            "temperature_crit": (33.2, pyunits.K),
            "omega": -0.218,
            "cp_mol_ig_comp_coeff": {
                "A": 33.066178,
                "B": -11.363417,
                "C": 11.432816,
                "D": -2.772874,
                "E": -0.158558,
                "F": -9.980797,
                "G": 172.707974,
                "H": 0.0,
            },
        },
    },
    "CO": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 1, "O": 1},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0280101, pyunits.kg / pyunits.mol),
            "pressure_crit": (35e5, pyunits.Pa),
            "temperature_crit": (132.9, pyunits.K),
            "omega": 0.066,
            "cp_mol_ig_comp_coeff": {
                "A": 25.56759,
                "B": 6.09613,
                "C": 4.054656,
                "D": -2.671301,
                "E": 0.131021,
                "F": -118.0089,
                "G": 227.3665,
                "H": -110.5271,
            },
        },
    },
    "H2O": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase, PhaseType.liquidPhase],
        "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
        "elemental_composition": {"H": 2, "O": 1},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "pressure_sat_comp": NIST,
        "parameter_data": {
            "mw": (0.01801528, pyunits.kg / pyunits.mol),
            "pressure_crit": (221.2e5, pyunits.Pa),
            "temperature_crit": (647.3, pyunits.K),
            "omega": 0.344,
            "cp_mol_ig_comp_coeff": {
                "A": 30.092,
                "B": 6.832514,
                "C": 6.793435,
                "D": -2.53448,
                "E": 0.082139,
                "F": -250.881,
                "G": 223.3967,
                "H": -241.8264,
            },
            "pressure_sat_comp_coeff": {  # NIST <- Stull 1947
                "A": 4.6543,
                "B": 1435.264,
                "C": -64.848,
            },
        },
    },
    "CO2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 1, "O": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "cp_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.04401, pyunits.kg / pyunits.mol),
            "pressure_crit": (73.8e5, pyunits.Pa),
            "temperature_crit": (304.1, pyunits.K),
            "omega": 0.239,
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
        },
    },
    "O2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"O": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.031998, pyunits.kg / pyunits.mol),
            "pressure_crit": (50.4e5, pyunits.Pa),
            "temperature_crit": (154.6, pyunits.K),
            "omega": 0.025,
            "cp_mol_ig_comp_coeff": {
                "A": 30.03235,
                "B": 8.772972,
                "C": -3.988133,
                "D": 0.788313,
                "E": -0.741599,
                "F": -11.32468,
                "G": 236.1663,
                "H": 0.0,
            },
        },
    },
    "N2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"N": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0280134, pyunits.kg / pyunits.mol),
            "pressure_crit": (33.9e5, pyunits.Pa),
            "temperature_crit": (126.2, pyunits.K),
            "omega": 0.039,
            "cp_mol_ig_comp_coeff": {
                "A": 19.50583,
                "B": 19.88705,
                "C": -8.598535,
                "D": 1.369784,
                "E": 0.527601,
                "F": -4.935202,
                "G": 212.39,
                "H": 0.0,
            },
        },
    },
    "Ar": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"Ar": 1},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.039948, pyunits.kg / pyunits.mol),
            "pressure_crit": (48.7e5, pyunits.Pa),
            "temperature_crit": (150.8, pyunits.K),
            "omega": 0.001,
            "cp_mol_ig_comp_coeff": {
                "A": 20.786,
                "B": 2.825911e-07,
                "C": -1.464191e-07,
                "D": 1.092131e-08,
                "E": -3.661371e-08,
                "F": -6.19735,
                "G": 179.999,
                "H": 0.0,
            },
        },
    },
    "CH4": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 1, "H": 4},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0160425, pyunits.kg / pyunits.mol),
            "pressure_crit": (46e5, pyunits.Pa),
            "temperature_crit": (190.4, pyunits.K),
            "omega": 0.011,
            "cp_mol_ig_comp_coeff": {
                "A": -0.703029,
                "B": 108.4773,
                "C": -42.52157,
                "D": 5.862788,
                "E": 0.678565,
                "F": -76.84376,
                "G": 158.7163,
                "H": -74.8731,
            },
        },
    },
    "C2H6": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 2, "H": 6},
        "enth_mol_ig_comp": RPP4,
        "entr_mol_ig_comp": RPP4,
        "parameter_data": {
            "mw": (0.030069, pyunits.kg / pyunits.mol),
            "pressure_crit": (48.8e5, pyunits.Pa),
            "temperature_crit": (305.4, pyunits.K),
            "omega": 0.099,
            "cp_mol_ig_comp_coeff": {
                "A": 5.409,
                "B": 1.781e-1,
                "C": -6.938e-5,
                "D": 8.713e-9,
            },
            "enth_mol_form_vap_comp_ref": (-84000, pyunits.J / pyunits.mol),
            "entr_mol_form_vap_comp_ref": (229.2, pyunits.J / pyunits.mol / pyunits.K),
        },
    },
    "C3H8": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 3, "H": 8},
        "enth_mol_ig_comp": RPP4,
        "entr_mol_ig_comp": RPP4,
        "parameter_data": {
            "mw": (0.0320849, pyunits.kg / pyunits.mol),
            "pressure_crit": (42.5e5, pyunits.Pa),
            "temperature_crit": (369.8, pyunits.K),
            "omega": 0.153,
            "cp_mol_ig_comp_coeff": {
                "A": -4.224,
                "B": 3.063e-1,
                "C": -1.586e-4,
                "D": 3.215e-8,
            },
            "enth_mol_form_vap_comp_ref": (-104700, pyunits.J / pyunits.mol),
            "entr_mol_form_vap_comp_ref": (270.3, pyunits.J / pyunits.mol / pyunits.K),
        },
    },
    "C4H10": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 4, "H": 10},
        "enth_mol_ig_comp": RPP5,
        "entr_mol_ig_comp": RPP5,
        "parameter_data": {
            "mw": (0.058123, pyunits.kg / pyunits.mol),  # RPP5
            "pressure_crit": (37.96e5, pyunits.Pa),  # RPP5
            "temperature_crit": (425.12, pyunits.K),  # RPP5
            "omega": 0.200,  # RPP5
            "cp_mol_ig_comp_coeff": {  # RPP5
                "a0": 5.547,
                "a1": 5.536e-3,
                "a2": 8.057e-5,
                "a3": -10.571e-8,
                "a4": 4.134e-11,
            },
            "enth_mol_form_vap_comp_ref": (-125790, pyunits.J / pyunits.mol),  # RPP5
            "entr_mol_form_vap_comp_ref": (
                310.23,
                pyunits.J / pyunits.mol / pyunits.K,
            ),  # wikipedia data page
        },
    },
    "H2S": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"H": 2, "S": 1},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.034081, pyunits.kg / pyunits.mol),  # NIST
            "pressure_crit": (89.6291e5, pyunits.Pa),  # NIST <- Goodwin 1983
            "temperature_crit": (373.4, pyunits.K),  # NIST <- Goodwin 1983
            "omega": 0.090,  # RPP5
            "cp_mol_ig_comp_coeff": {  # NIST <- Chase 1998
                "A": 26.88412,
                "B": 18.67809,
                "C": 3.434203,
                "D": -3.378702,
                "E": 0.135882,
                "F": -28.91211,
                "G": 233.3747,
                "H": -20.50202,
            },
        },
    },
    "SO2": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"S": 1, "O": 2},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.064064, pyunits.kg / pyunits.mol),  # NIST
            "pressure_crit": (78.84e5, pyunits.Pa),  # RPP5
            "temperature_crit": (430.80, pyunits.K),  # RPP5
            "omega": 0.0,  # RPP5 (Not listed)
            "cp_mol_ig_comp_coeff": {  # NIST <- Chase 1998
                "A": 21.43049,
                "B": 74.35094,
                "C": -57.75217,
                "D": 16.35534,
                "E": 0.086731,
                "F": -305.7688,
                "G": 254.8872,
                "H": -296.8422,
            },
        },
    },
    "C2H4": {
        "type": Component,
        "valid_phase_types": [PhaseType.vaporPhase],
        "elemental_composition": {"C": 2, "H": 4},
        "enth_mol_ig_comp": NIST,
        "entr_mol_ig_comp": NIST,
        "parameter_data": {
            "mw": (0.0280532, pyunits.kg / pyunits.mol),  # NIST
            "pressure_crit": (50.6e5, pyunits.Pa),  # NIST
            "temperature_crit": (282.5, pyunits.K),  # NIST
            "omega": 0.087,  # RPP5
            "cp_mol_ig_comp_coeff": {  # NIST <- Chase 1998
                "A": -6.387880,
                "B": 184.4019,
                "C": -112.9718,
                "D": 28.49593,
                "E": 0.315540,
                "F": 48.17332,
                "G": 163.1568,
                "H": 52.46694,
            },
        },
    },
}


# returns a configuration dictionary for the list of specified components
def get_prop(components=None, phases="Vap", eos=EosType.PR, scaled=False):
    if components is None:
        components = list(_component_params.keys())
    configuration = {
        "components": {},  # fill in later based on selected components
        "parameter_data": {},
        "phases": {},
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
            "flow_mol": (0, 8000, 50000, pyunits.mol / pyunits.s),
            "temperature": (273.15, 500, 2500, pyunits.K),
            "pressure": (5e4, 1.3e5, 1e8, pyunits.Pa),
        },
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
    }

    c = configuration["components"]
    for comp in components:
        c[comp] = copy.deepcopy(_component_params[comp])
    if isinstance(phases, str):
        phases = [phases]
    for k in phases:
        if eos == EosType.PR:
            configuration["phases"][k] = copy.deepcopy(_phase_dicts_pr[k])
        elif eos == EosType.IDEAL:
            if k == "Liq":
                raise ConfigurationError(
                    "This parameter set does not support Ideal EOS with liquid"
                )
            configuration["phases"][k] = copy.deepcopy(_phase_dicts_ideal[k])
        else:
            raise ValueError("Invalid EoS.")
    if len(phases) > 1:
        p = tuple(phases)
        configuration["phases_in_equilibrium"] = [p]
        configuration["phase_equilibrium_state"] = {p: SmoothVLE}

    # Fill the binary parameters with zeros.
    d = configuration["parameter_data"]
    d["PR_kappa"] = {(a, b): 0 for a in c for b in c}

    # Change to scaled units if specified
    if scaled:
        configuration["base_units"]["mass"] = pyunits.Mg
        configuration["base_units"]["amount"] = pyunits.kmol

    return configuration


def get_rxn(property_package, reactions=None, scaled=False):
    rxns = {
        "property_package": property_package,
        "base_units": {
            "time": pyo.units.s,
            "length": pyo.units.m,
            "mass": pyo.units.kg,
            "amount": pyo.units.mol,
            "temperature": pyo.units.K,
        },
        "rate_reactions": {
            "ch4_cmb": {
                "stoichiometry": {
                    ("Vap", "CH4"): -1,
                    ("Vap", "O2"): -2,
                    ("Vap", "H2O"): 2,
                    ("Vap", "CO2"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
            "c2h6_cmb": {
                "stoichiometry": {
                    ("Vap", "C2H6"): -2,
                    ("Vap", "O2"): -7,
                    ("Vap", "H2O"): 6,
                    ("Vap", "CO2"): 4,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
            "c3h8_cmb": {
                "stoichiometry": {
                    ("Vap", "C3H8"): -1,
                    ("Vap", "O2"): -5,
                    ("Vap", "H2O"): 4,
                    ("Vap", "CO2"): 3,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
            "c4h10_cmb": {
                "stoichiometry": {
                    ("Vap", "C4H10"): -2,
                    ("Vap", "O2"): -13,
                    ("Vap", "H2O"): 10,
                    ("Vap", "CO2"): 8,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
            "co_cmb": {
                "stoichiometry": {
                    ("Vap", "CO"): -2,
                    ("Vap", "O2"): -1,
                    ("Vap", "CO2"): 2,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
            "h2_cmb": {
                "stoichiometry": {
                    ("Vap", "H2"): -2,
                    ("Vap", "O2"): -1,
                    ("Vap", "H2O"): 2,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
            "h2s_cmb": {
                "stoichiometry": {
                    ("Vap", "H2S"): -2,
                    ("Vap", "O2"): -3,
                    ("Vap", "H2O"): 2,
                    ("Vap", "SO2"): 2,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
            "c2h4_cmb": {
                "stoichiometry": {
                    ("Vap", "C2H4"): -1,
                    ("Vap", "O2"): -3,
                    ("Vap", "H2O"): 2,
                    ("Vap", "CO2"): 2,
                },
                "heat_of_reaction": constant_dh_rxn,
                "rate_constant": arrhenius,
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "dh_rxn_ref": 0,
                    "arrhenius_const": 0,
                    "energy_activation": 0,
                },
            },
        },
    }
    # Change to scaled units if specified
    if scaled:
        rxns["base_units"]["mass"] = pyunits.Mg
        rxns["base_units"]["amount"] = pyunits.kmol
    if reactions is None:
        return rxns
    else:
        for r in list(rxns["rate_reactions"].keys()):
            if r not in reactions:
                del rxns["rate_reactions"][r]
        return rxns
