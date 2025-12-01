#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Config dictionaries used in testing of electrolyte states with and without
inherent reactions.

Authors: Andrew Lee, Douglas Allan
"""
# Import Pyomo units
from pyomo.environ import units as pyunits

# Import IDAES cores
from idaes.core import AqueousPhase, VaporPhase
from idaes.core.base.components import (
    Solvent,
    Solute,
    Apparent,
    Cation,
    Anion,
    Component,
)

from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    power_law_equil,
)
from idaes.models.properties.modular_properties.base.utility import ConcentrationForm

from idaes.core.solvers import get_solver

solver = get_solver("ipopt_v2")


def dens_mol_H2O(*args, **kwargs):
    """
    Method returning density of water
    """
    return 55e3 * pyunits.mol / pyunits.m**3


def dummy_method(b, *args, **kwargs):
    """
    Dummy method for testing
    """
    return 42


def get_config_no_inherent_reactions():
    """
    Method to get a config file to test electrolyte states
    without inherent reactions.
    """
    config = {
        # Specifying components
        "components": {
            "H2O": {
                "type": Solvent,
                "enth_mol_ig_comp": dummy_method,
                "parameter_data": {"mw": (18e-3, pyunits.kg / pyunits.mol)},
            },
            "CO2": {
                "type": Solute,
                "enth_mol_ig_comp": dummy_method,
                "parameter_data": {"mw": (44e-3, pyunits.kg / pyunits.mol)},
            },
            "KHCO3": {
                "type": Apparent,
                "enth_mol_ig_comp": dummy_method,
                "dissociation_species": {"K+": 1, "HCO3-": 1},
                "parameter_data": {"mw": (100.1e-3, pyunits.kg / pyunits.mol)},
            },
            "K+": {
                "type": Cation,
                "charge": +1,
                "parameter_data": {"mw": (39.1e-3, pyunits.kg / pyunits.mol)},
            },
            "HCO3-": {
                "type": Anion,
                "charge": -1,
                "parameter_data": {"mw": (61e-3, pyunits.kg / pyunits.mol)},
            },
            "N2": {
                "type": Component,
                "enth_mol_ig_comp": dummy_method,
                "parameter_data": {"mw": (28e-3, pyunits.kg / pyunits.mol)},
            },
        },
        # Specifying phases
        "phases": {
            "Liq": {
                "type": AqueousPhase,
                "equation_of_state": DummyEoS,
                "equation_of_state_options": {"pH_range": "basic"},
            },
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
        # Leave the state definition and bounds to be
        # specified in each individual test file.
        # "state_definition": None,
        # "state_bounds": None,
        # "state_components": None,
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
    }
    return config


def get_config_with_inherent_reactions():
    """
    Method to get a config file to test electrolyte states
    with inherent reactions.
    """
    config = {
        # Specifying components
        "components": {
            "H2O": {
                "type": Solvent,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (18e-3, pyunits.kg / pyunits.mol)},
            },
            "KHCO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "HCO3-": 1},
                "parameter_data": {"mw": (100.1e-3, pyunits.kg / pyunits.mol)},
            },
            "K2CO3": {
                "type": Apparent,
                "dissociation_species": {"K+": 2, "CO3--": 1},
                "parameter_data": {"mw": (138.2e-3, pyunits.kg / pyunits.mol)},
            },
            "KOH": {
                "type": Apparent,
                "dissociation_species": {"K+": 1, "OH-": 1},
                "parameter_data": {"mw": (56.1e-3, pyunits.kg / pyunits.mol)},
            },
            "H+": {
                "type": Cation,
                "charge": +1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (1e-3, pyunits.kg / pyunits.mol)},
            },
            "K+": {
                "type": Cation,
                "charge": +1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (39.1e-3, pyunits.kg / pyunits.mol)},
            },
            "OH-": {
                "type": Anion,
                "charge": -1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (17e-3, pyunits.kg / pyunits.mol)},
            },
            "HCO3-": {
                "type": Anion,
                "charge": -1,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (61e-3, pyunits.kg / pyunits.mol)},
            },
            "CO3--": {
                "type": Anion,
                "charge": -2,
                "dens_mol_liq_comp": dens_mol_H2O,
                "parameter_data": {"mw": (60e-3, pyunits.kg / pyunits.mol)},
            },
        },
        # Specifying phases
        "phases": {
            "Liq": {
                "type": AqueousPhase,
                "equation_of_state": DummyEoS,
                "equation_of_state_options": {"pH_range": "basic"},
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
        # Leave the state definition and bounds to be
        # specified in each individual test file.
        # "state_definition": None,
        # "state_bounds": None,
        # "state_components": None,
        "pressure_ref": (101325, pyunits.Pa),
        "temperature_ref": (298.15, pyunits.K),
        "inherent_reactions": {
            "h2o_si": {
                "stoichiometry": {
                    ("Liq", "H2O"): -1,
                    ("Liq", "H+"): 1,
                    ("Liq", "OH-"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "equilibrium_constant": van_t_hoff,
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.molarity,
                "parameter_data": {
                    "reaction_order": {("Liq", "H+"): 1, ("Liq", "OH-"): 1},
                    "dh_rxn_ref": 1,
                    "k_eq_ref": 1e-14,
                    "T_eq_ref": 350,
                },
            },
            "co3_hco3": {
                "stoichiometry": {
                    ("Liq", "CO3--"): -1,
                    ("Liq", "H2O"): -1,
                    ("Liq", "HCO3-"): 1,
                    ("Liq", "OH-"): 1,
                },
                "heat_of_reaction": constant_dh_rxn,
                "equilibrium_constant": van_t_hoff,
                "equilibrium_form": power_law_equil,
                "concentration_form": ConcentrationForm.molarity,
                "parameter_data": {
                    "reaction_order": {
                        ("Liq", "CO3--"): -1,
                        ("Liq", "HCO3-"): 1,
                        ("Liq", "OH-"): 1,
                    },
                    "dh_rxn_ref": 1,
                    "k_eq_ref": 5e-11,
                    "T_eq_ref": 350,
                },
            },
        },
    }
    return config
