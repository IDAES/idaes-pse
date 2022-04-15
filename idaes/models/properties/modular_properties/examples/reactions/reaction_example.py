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
Example of defining a reaction package for the following reaction system:

R1)  A + B --> 2C  Power law rate reaction
R2)  B + C <-> D  Power law equilibrium reaction

Author: Andrew Lee
"""
from pyomo.environ import units as pyunits

from idaes.core import LiquidPhase, Component

from idaes.models.properties.modular_properties.state_definitions import FcTP
from idaes.models.properties.modular_properties.eos.ideal import Ideal

from idaes.models.properties.modular_properties.base.utility import ConcentrationForm
from idaes.models.properties.modular_properties.reactions.dh_rxn import constant_dh_rxn
from idaes.models.properties.modular_properties.reactions.rate_constant import arrhenius
from idaes.models.properties.modular_properties.reactions.rate_forms import (
    power_law_rate,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_constant import (
    van_t_hoff,
)
from idaes.models.properties.modular_properties.reactions.equilibrium_forms import (
    power_law_equil,
)


# First, create a thermophsyical property definition that will be used
# with the reactions
import idaes.models.properties.modular_properties.pure.Perrys as Perrys

# For this example, the thermophsycial properties will only define components
# and phases, but in practice users would also need to define some properties
thermo_configuration = {
    "components": {
        "A": {
            "type": Component,
            "enth_mol_liq_comp": Perrys,
            "parameter_data": {
                "cp_mol_liq_comp_coeff": {
                    "1": (1e2, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (6.5e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (5e4, pyunits.J / pyunits.mol),
            },
        },
        "B": {
            "type": Component,
            "enth_mol_liq_comp": Perrys,
            "parameter_data": {
                "cp_mol_liq_comp_coeff": {
                    "1": (1e2, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (6.5e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (5e4, pyunits.J / pyunits.mol),
            },
        },
        "C": {
            "type": Component,
            "enth_mol_liq_comp": Perrys,
            "parameter_data": {
                "cp_mol_liq_comp_coeff": {
                    "1": (1e2, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (6.5e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (5e4, pyunits.J / pyunits.mol),
            },
        },
        "D": {
            "type": Component,
            "enth_mol_liq_comp": Perrys,
            "parameter_data": {
                "cp_mol_liq_comp_coeff": {
                    "1": (1e2, pyunits.J / pyunits.kmol / pyunits.K),
                    "2": (-2e-1, pyunits.J / pyunits.kmol / pyunits.K**2),
                    "3": (6.5e-4, pyunits.J / pyunits.kmol / pyunits.K**3),
                    "4": (0, pyunits.J / pyunits.kmol / pyunits.K**4),
                    "5": (0, pyunits.J / pyunits.kmol / pyunits.K**5),
                },
                "enth_mol_form_liq_comp_ref": (5e4, pyunits.J / pyunits.mol),
            },
        },
    },
    "phases": {"Liq": {"type": LiquidPhase, "equation_of_state": Ideal}},
    "state_definition": FcTP,
    "state_bounds": {
        "flow_mol_comp": (0, 500, 1000),
        "temperature": (273.15, 300, 450),
        "pressure": (5e4, 1e5, 1e6),
    },
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
}


# Next, create the reaction property definition which describes the system on
# reactions to be modeled.
rxn_configuration = {
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "rate_reactions": {
        "R1": {
            "stoichiometry": {("Liq", "A"): -1, ("Liq", "B"): -1, ("Liq", "C"): 2},
            "heat_of_reaction": constant_dh_rxn,
            "rate_constant": arrhenius,
            "rate_form": power_law_rate,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (-10000, pyunits.J / pyunits.mol),
                "arrhenius_const": (1, pyunits.mol / pyunits.m**3 / pyunits.s),
                "energy_activation": (1000, pyunits.J / pyunits.mol),
            },
        }
    },
    "equilibrium_reactions": {
        "R2": {
            "stoichiometry": {("Liq", "B"): -1, ("Liq", "C"): -1, ("Liq", "D"): 1},
            "heat_of_reaction": constant_dh_rxn,
            "equilibrium_constant": van_t_hoff,
            "equilibrium_form": power_law_equil,
            "concentration_form": ConcentrationForm.moleFraction,
            "parameter_data": {
                "dh_rxn_ref": (-20000, pyunits.J / pyunits.mol),
                "k_eq_ref": (100, None),
                "T_eq_ref": (350, pyunits.K),
            },
        }
    },
}
