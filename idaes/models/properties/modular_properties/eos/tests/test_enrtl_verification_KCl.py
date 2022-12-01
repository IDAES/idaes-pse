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
Tests for eNRTL methods

References:

[1] Local Composition Model for Excess Gibbs Energy of Electrolyte Systems, Pt 1.
Chen, C.-C., Britt, H.I., Boston, J.F., Evans, L.B.,
AIChE Journal, 1982, Vol. 28(4), pgs. 588-596

[2] New Data on Activity Coefficients of Potassium, Nitrate, and Chloride Ions
in Aqueous Solutions of KNO3 and KCl by Ion Selective Electrodes
Dash, D., Kumar, S., Mallika, C., Kamachi Mudali, U.,
ISRN Chemical Engineering, 2012, doi:10.5402/2012/730154

Author: Andrew Lee
"""
import pytest
from math import log

from pyomo.environ import ConcreteModel, units as pyunits, value

from idaes.core import AqueousPhase, Solvent, Apparent, Anion, Cation
from idaes.models.properties.modular_properties.eos.enrtl import ENRTL
from idaes.models.properties.modular_properties.eos.enrtl_reference_states import (
    Unsymmetric,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    StateIndex,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.pure.electrolyte import (
    relative_permittivity_constant,
)


def dummy_method(b, *args, **kwargs):
    return 1000 / 18e-3 * pyunits.mol / pyunits.m**3


configuration = {
    "components": {
        "H2O": {
            "type": Solvent,
            "dens_mol_liq_comp": dummy_method,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (18e-3, pyunits.kg / pyunits.mol),
                "relative_permittivity_liq_comp": 78.54,
            },
        },
        "KCl": {"type": Apparent, "dissociation_species": {"K+": 1, "Cl-": 1}},
        "K+": {"type": Cation, "charge": +1},
        "Cl-": {"type": Anion, "charge": -1},
    },
    "phases": {
        "Liq": {
            "type": AqueousPhase,
            "equation_of_state": ENRTL,
            "equation_of_state_options": {"reference_state": Unsymmetric},
        }
    },
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    "state_definition": FTPx,
    "state_components": StateIndex.true,
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "parameter_data": {
        "Liq_tau": {
            ("H2O", "K+, Cl-"): 8.064,  # Table 1, [1]
            ("K+, Cl-", "H2O"): -4.107,
        }
    },
}  # Table 1, [1]


class TestStateBlockUnsymmetric(object):
    @pytest.fixture(scope="class")
    def model(self):
        config = dict(configuration)
        eos_opt = config["phases"]["Liq"]["equation_of_state_options"] = {}
        eos_opt["reference_state"] = Unsymmetric

        m = ConcreteModel()
        m.params = GenericParameterBlock(**config)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(298)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert model.params.Liq.tau["H2O", "K+, Cl-"].value == 8.064
        assert model.params.Liq.tau["K+, Cl-", "H2O"].value == -4.107

    @pytest.mark.unit
    def test_log_gamma(self, model):
        # start with pure water
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Data from [2] - Form {molality: gamma_KCl}
        data = {
            0.001: 0.965,
            0.002: 0.951,
            0.005: 0.927,
            0.01: 0.901,
            0.0201: 0.869,
            0.05: 0.816,
            0.0982: 0.768,
            0.1983: 0.717,
            0.3049: 0.687,
            0.3999: 0.665,
            0.5001: 0.649,
            0.6032: 0.636,
            0.7003: 0.626,
            0.7996: 0.617,
            0.8968: 0.61,
            0.9926: 0.604,
            1.1899: 0.594,
            1.3844: 0.586,
            1.5817: 0.58,
            1.7923: 0.576,
            1.9895: 0.573,
            2.5039: 0.568,
            2.9837: 0.568,
            3.4982: 0.571,
            3.994: 0.576,
            4.4897: 0.584,
            4.7909: 0.589,
            4.9908: 0.593,
        }

        for x, g in data.items():
            w = 1000 / 18

            model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(w / (w + 2 * x))
            model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(x / (w + 2 * x))
            model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(x / (w + 2 * x))

            conv = log(1 + 18 * 2 * x / 1000)  # Convert mole frac to molal basis
            assert pytest.approx(log(g), rel=3e-2, abs=6e-3) == value(
                model.state[1].Liq_log_gamma["K+"] - conv
            )
            assert pytest.approx(log(g), rel=3e-2, abs=6e-3) == value(
                model.state[1].Liq_log_gamma["Cl-"] - conv
            )
