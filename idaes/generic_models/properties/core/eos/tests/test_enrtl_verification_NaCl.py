##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for eNRTL methods

Parameter values from:

A Comprehensive Thermodynamic Model for High Salinity Produced Waters
Tanveer, S., and Chen, C.-C., AIChE Journal, 2019, Vol 66
Supplementary Information
 https://doi.org/10.1002/aic.16818

Author: Andrew Lee
"""
import pytest
from math import exp

from pyomo.environ import (ConcreteModel,
                           units as pyunits,
                           value)
from pyomo.util.check_units import assert_units_equivalent

from idaes.core import (AqueousPhase,
                        Solvent,
                        Apparent,
                        Anion,
                        Cation)
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock, StateIndex)
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.pure.electrolyte import \
    relative_permittivity_constant


def dummy_method(b, *args, **kwargs):
    return 42*pyunits.mol/pyunits.m**3


configuration = {
    "components": {
        "H2O": {"type": Solvent,
                "dens_mol_liq_comp": dummy_method,
                "relative_permittivity_liq_comp":
                    relative_permittivity_constant,
                "parameter_data": {"mw": (18E-3, pyunits.kg/pyunits.mol),
                                   "relative_permittivity_liq_comp": 101}},
        "NaCl": {"type": Apparent,
                 "dissociation_species": {"Na+": 1, "Cl-": 1}},
        "Na+": {"type": Cation,
                "charge": +1},
        "Cl-": {"type": Anion,
                "charge": -1}},
    "phases": {
        "Liq": {"type": AqueousPhase,
                "equation_of_state": ENRTL}},
    "base_units": {"time": pyunits.s,
                   "length": pyunits.m,
                   "mass": pyunits.kg,
                   "amount": pyunits.mol,
                   "temperature": pyunits.K},
    "state_definition": FTPx,
    "state_components": StateIndex.true,
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "parameter_data": {
        "Liq_tau": {
            ("H2O", "Na+, Cl-"): 8.865,
            ("Na+, Cl-", "H2O"): -4.541}}}


class TestStateBlockSymmetric(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(default=configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(300)

        return m

    @pytest.mark.unit
    def test_parameters(self, model):
        assert model.params.Liq.tau["H2O", "Na+, Cl-"].value == 8.865
        assert model.params.Liq.tau["Na+, Cl-", "H2O"].value == -4.541

    # TODO: Need to check effects of binary parameters
    @pytest.mark.unit
    def test_pure_water(self, model):
        # Start by setting all mole fractions to small number
        # Using 0 results in division by zero errors
        for k in model.state[1].mole_frac_phase_comp:
            model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

        # Test pure water
        model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(1)

        # Check mixing expressions
        assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(1, rel=1e-8)
        assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
            1e-12, rel=1e-8)
        assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
            -1e-12, rel=1e-8)

        for v in model.state[1].Liq_Y.values():
            assert value(v) == pytest.approx(1, rel=1e-8)

        for k, v in model.state[1].Liq_alpha.items():
            if k == ("H2O", "H2O"):
                assert value(v) == 0.3
            elif k in [("Na+", "Cl-"), ("Cl-", "Na+")]:
                assert value(v) == 0
            else:
                assert value(v) == 0.2

        assert value(model.state[1].Liq_G["H2O", "H2O"]) == 1
        assert value(model.state[1].Liq_G["H2O", "Na+"]) == exp(-0.2*8.865)
        assert value(model.state[1].Liq_G["Na+", "H2O"]) == exp(-0.2*-4.541)
        assert value(model.state[1].Liq_G["H2O", "Cl-"]) == exp(-0.2*8.865)
        assert value(model.state[1].Liq_G["Cl-", "H2O"]) == exp(-0.2*-4.541)
        assert value(model.state[1].Liq_G["Na+", "Cl-"]) == 1
        assert value(model.state[1].Liq_G["Cl-", "Na+"]) == 1

        assert value(model.state[1].Liq_tau["H2O", "H2O"]) == 0
        assert value(model.state[1].Liq_tau["H2O", "Na+"]) == pytest.approx(
            8.865, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "H2O"]) == pytest.approx(
            -4.541, rel=1e-8)
        assert value(model.state[1].Liq_tau["H2O", "Cl-"]) == pytest.approx(
            8.865, rel=1e-8)
        assert value(model.state[1].Liq_tau["Cl-", "H2O"]) == pytest.approx(
            -4.541, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "Cl-"]) == 0
        assert value(model.state[1].Liq_tau["Cl-", "Na+"]) == 0

        # Check activity coefficient contributions
        assert (value(model.state[1].Liq_log_gamma_pdh["H2O"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_lc["H2O"]) ==
                pytest.approx(0, abs=1e-10))
        assert (value(model.state[1].Liq_log_gamma_pdh["H2O"] +
                      model.state[1].Liq_log_gamma_lc["H2O"]) ==
                pytest.approx(0, abs=1e-10))

    @pytest.mark.unit
    def test_pure_NaCl(self, model):
        # Test pure NaCl
        model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(1e-12)
        model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(0.5)
        model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(0.5)

        # Check mixing expressions
        assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(
            1e-12, rel=1e-8)
        assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
            0.5, rel=1e-8)
        assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
            -0.5, rel=1e-8)

        for v in model.state[1].Liq_Y.values():
            assert value(v) == pytest.approx(1, rel=1e-8)

        for k, v in model.state[1].Liq_alpha.items():
            if k == ("H2O", "H2O"):
                assert value(v) == 0.3
            elif k in [("Na+", "Cl-"), ("Cl-", "Na+")]:
                assert value(v) == 0
            else:
                assert value(v) == 0.2

        assert value(model.state[1].Liq_G["H2O", "H2O"]) == 1
        assert value(model.state[1].Liq_G["H2O", "Na+"]) == exp(-0.2*8.865)
        assert value(model.state[1].Liq_G["Na+", "H2O"]) == exp(-0.2*-4.541)
        assert value(model.state[1].Liq_G["H2O", "Cl-"]) == exp(-0.2*8.865)
        assert value(model.state[1].Liq_G["Cl-", "H2O"]) == exp(-0.2*-4.541)
        assert value(model.state[1].Liq_G["Na+", "Cl-"]) == 1
        assert value(model.state[1].Liq_G["Cl-", "Na+"]) == 1

        assert value(model.state[1].Liq_tau["H2O", "H2O"]) == 0
        assert value(model.state[1].Liq_tau["H2O", "Na+"]) == pytest.approx(
            8.865, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "H2O"]) == pytest.approx(
            -4.541, rel=1e-8)
        assert value(model.state[1].Liq_tau["H2O", "Cl-"]) == pytest.approx(
            8.865, rel=1e-8)
        assert value(model.state[1].Liq_tau["Cl-", "H2O"]) == pytest.approx(
            -4.541, rel=1e-8)
        assert value(model.state[1].Liq_tau["Na+", "Cl-"]) == 0
        assert value(model.state[1].Liq_tau["Cl-", "Na+"]) == 0

        assert (value(model.state[1].Liq_log_gamma_pdh["Na+"]) ==
                pytest.approx(0, abs=1e-10))

        # assert (value(model.state[1].Liq_log_gamma_lc_I["Na+"]) ==
        #         value(model.state[1].Liq_log_gamma_lc_I0["Na+"]))

        # assert (value(model.state[1].Liq_log_gamma_lc["Na+"]) ==
        #         pytest.approx(0, abs=1e-10))
        # assert (value(model.state[1].Liq_log_gamma_pdh["Na+"] +
        #               model.state[1].Liq_log_gamma_lc["Na+"]) ==
        #         pytest.approx(0, abs=1e-10))

        assert (value(model.state[1].Liq_log_gamma_pdh["Cl-"]) ==
                pytest.approx(0, abs=1e-10))

    #     assert (value(model.state[1].Liq_log_gamma_lc_I["Cl-"]) ==
    #             value(model.state[1].Liq_log_gamma_lc_I0["Cl-"]))

    #     assert (value(model.state[1].Liq_log_gamma_lc["Cl-"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_pdh["Cl-"] +
    #                   model.state[1].Liq_log_gamma_lc["Cl-"]) ==
    #             pytest.approx(0, abs=1e-10))
