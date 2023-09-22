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
Example for using the eNRTL model for aqueous solutions of NaCl and KCl.

Reference:

[1] Local Composition Model for Excess Gibbs Energy of Electrolyte Systems, Pt 1.
Chen, C.-C., Britt, H.I., Boston, J.F., Evans, L.B.,
AIChE Journal, 1982, Vol. 28(4), pgs. 588-596

[2] Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

[3] New Data on Activity Coefficients of Potassium, Nitrate, and Chloride Ions
in Aqueous Solutions of KNO3 and KCl by Ion Selective Electrodes
Dash, D., Kumar, S., Mallika, C., Kamachi Mudali, U.,
ISRN Chemical Engineering, 2012, doi:10.5402/2012/730154

Author: Andrew Lee
"""
import pyomo.environ as pyo

from idaes.core import AqueousPhase, Solvent, Apparent, Anion, Cation
from idaes.models.properties.modular_properties.eos.enrtl import ENRTL
from idaes.models.properties.modular_properties.eos.enrtl_reference_states import (
    Symmetric,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    StateIndex,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.pure.electrolyte import (
    relative_permittivity_constant,
)


def constant_density(b, *args, **kwargs):
    """Assume constant density of pure water"""
    return 1000 / 18e-3 * pyo.units.mol / pyo.units.m**3


configuration = {
    "components": {
        "H2O": {
            "type": Solvent,
            "dens_mol_liq_comp": constant_density,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (18e-3, pyo.units.kg / pyo.units.mol),
                "relative_permittivity_liq_comp": 78.54,
            },
        },
        "NaCl": {"type": Apparent, "dissociation_species": {"Na+": 1, "Cl-": 1}},
        "KCl": {"type": Apparent, "dissociation_species": {"K+": 1, "Cl-": 1}},
        "Na+": {"type": Cation, "charge": +1},
        "K+": {"type": Cation, "charge": +1},
        "Cl-": {"type": Anion, "charge": -1},
    },
    "phases": {
        "Liq": {
            "type": AqueousPhase,
            "equation_of_state": ENRTL,
            "equation_of_state_options": {"reference_state": Symmetric},
        }
    },
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
        "Liq_tau": {
            ("H2O", "Na+, Cl-"): 9.0234,
            ("Na+, Cl-", "H2O"): -4.5916,
            ("H2O", "K+, Cl-"): 8.1354,
            ("K+, Cl-", "H2O"): -4.1341,
        }
    },
}
