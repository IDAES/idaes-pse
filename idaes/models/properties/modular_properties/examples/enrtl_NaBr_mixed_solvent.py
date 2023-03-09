#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Example of using the eNRTL model to model solutions of NaBr in mixed solvents.

Reference:

[1] Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

Author: Andrew Lee
"""
import pyomo.environ as pyo

from idaes.core import AqueousPhase, Solvent, Apparent, Anion, Cation
from idaes.models.properties.modular_properties.eos.enrtl import ENRTL
from idaes.models.properties.modular_properties.base.generic_property import (
    StateIndex,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.pure.electrolyte import (
    relative_permittivity_constant,
)


def rho_H2O(b, *args, **kwargs):
    """Assume constant density of pure water"""
    return 1000 / 18e-3 * pyo.units.mol / pyo.units.m**3


def rho_MeOH(b, *args, **kwargs):
    """Assume constant density for ethanol"""
    return 792 / 32e-3 * pyo.units.mol / pyo.units.m**3


def rho_EtOH(b, *args, **kwargs):
    """Assume constant density for methanol"""
    return 789.45 / 46e-3 * pyo.units.mol / pyo.units.m**3


configuration = {
    "components": {
        "H2O": {
            "type": Solvent,
            "dens_mol_liq_comp": rho_H2O,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (18e-3, pyo.units.kg / pyo.units.mol),
                "relative_permittivity_liq_comp": 78.54,
            },
        },
        "MeOH": {
            "type": Solvent,
            "dens_mol_liq_comp": rho_MeOH,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (32e-3, pyo.units.kg / pyo.units.mol),
                "relative_permittivity_liq_comp": 32.6146,
            },
        },
        "EtOH": {
            "type": Solvent,
            "dens_mol_liq_comp": rho_EtOH,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (46e-3, pyo.units.kg / pyo.units.mol),
                "relative_permittivity_liq_comp": 24.113,
            },
        },
        "NaBr": {"type": Apparent, "dissociation_species": {"Na+": 1, "Br-": 1}},
        "Na+": {"type": Cation, "charge": +1},
        "Br-": {"type": Anion, "charge": -1},
    },
    "phases": {"Liq": {"type": AqueousPhase, "equation_of_state": ENRTL}},
    "base_units": {
        "time": pyo.units.s,
        "length": pyo.units.m,
        "mass": pyo.units.kg,
        "amount": pyo.units.mol,
        "temperature": pyo.units.K,
    },
    "state_definition": FTPx,
    "state_components": StateIndex.true,
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "parameter_data": {
        "Liq_alpha": {
            ("H2O", "EtOH"): 0.3031,
            ("H2O", "MeOH"): 0.2994,
            ("MeOH", "EtOH"): 0.3356,
            ("H2O", "Na+, Br-"): 0.2,
            ("MeOH", "Na+, Br-"): 0.2,
            ("EtOH", "Na+, Br-"): 0.1,
        },
        "Liq_tau": {
            ("H2O", "MeOH"): 1.4265,
            ("MeOH", "H2O"): -0.42864,
            ("H2O", "EtOH"): 2.2485,
            ("EtOH", "H2O"): -0.18514,
            ("MeOH", "EtOH"): -0.04394,
            ("EtOH", "MeOH"): 0.02147,
            ("H2O", "Na+, Br-"): 9.527,
            ("Na+, Br-", "H2O"): -4.790,
            ("MeOH", "Na+, Br-"): 5.910,
            ("Na+, Br-", "MeOH"): -3.863,
            ("EtOH", "Na+, Br-"): 6.118,
            ("Na+, Br-", "EtOH"): -4.450,
        },
    },
}
