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

Reference:

Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797

Figures 1 and 2

Author: Andrew Lee
"""
import pytest
from math import exp, log

from pyomo.environ import (ConcreteModel,
                           units as pyunits,
                           value)

from idaes.core import (AqueousPhase,
                        Solvent,
                        Solute,
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
    return 53724.94095814129*pyunits.mol/pyunits.m**3


configuration = {
    "components": {
        "H2O": {"type": Solvent,
                "dens_mol_liq_comp": dummy_method,
                "relative_permittivity_liq_comp":
                    relative_permittivity_constant,
                "parameter_data": {
                    "mw": (18E-3, pyunits.kg/pyunits.mol),
                    "relative_permittivity_liq_comp": 73.41964622627609}},
        "MEA": {"type": Solvent,
                "dens_mol_liq_comp": dummy_method,
                "relative_permittivity_liq_comp":
                    relative_permittivity_constant,
                "parameter_data": {
                    "mw": (61.084E-3, pyunits.kg/pyunits.mol),
                    "relative_permittivity_liq_comp": 28.662707746759963}},
        "CO2": {"type": Solute,
                "dens_mol_liq_comp": dummy_method,
                "relative_permittivity_liq_comp":
                    relative_permittivity_constant,
                "parameter_data": {"mw": (44E-3, pyunits.kg/pyunits.mol),
                                   "relative_permittivity_liq_comp": 101}},
        "MEA+": {"type": Cation,
                 "charge": +1},
        "HCO3-": {"type": Anion,
                  "charge": -1},
        "MEACOO-": {"type": Anion,
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
        "Liq_alpha": {
            ('CO2', 'MEA+, HCO3-'): 0.1,
            ('CO2', 'MEA+, MEACOO-'): 0.1,
            ('H2O', 'MEA+, HCO3-'): 0.2,
            ('H2O', 'MEA+, MEACOO-'): 0.2,
            ('MEA', 'MEA+, HCO3-'): 0.1,
            ('MEA', 'MEA+, MEACOO-'): 0.1,
            ('MEA+, HCO3-', 'CO2'): 0.1,
            ('MEA+, HCO3-', 'H2O'): 0.2,
            ('MEA+, HCO3-', 'MEA'): 0.1,
            ('MEA+, MEACOO-', 'CO2'): 0.1,
            ('MEA+, MEACOO-', 'H2O'): 0.2,
            ('MEA+, MEACOO-', 'MEA'): 0.1,
            ('CO2', 'CO2'): 0.2,
            ('CO2', 'H2O'): 0.2,
            ('CO2', 'MEA'): 0.2,
            ('H2O', 'H2O'): 0.2,
            ('H2O', 'MEA'): 0.2,
            ('MEA', 'MEA'): 0.2},
        "Liq_tau": {
            ('H2O', 'MEA'): 0.509893610223642,
            ('MEA', 'H2O'): -1.38820670926518,
            ('CO2', 'MEA+, HCO3-'): 15.0,
            ('CO2', 'MEA+, MEACOO-'): 15.0,
            ('H2O', 'MEA+, HCO3-'): 8.5721,
            ('H2O', 'MEA+, MEACOO-'): 6.702476996805112,
            ('MEA', 'MEA+, HCO3-'): 8.0,
            ('MEA', 'MEA+, MEACOO-'): 8.0,
            ('MEA+, HCO3-', 'CO2'): -8.0,
            ('MEA+, HCO3-', 'H2O'): -4.0092,
            ('MEA+, HCO3-', 'MEA'): -4.0,
            ('MEA+, MEACOO-', 'CO2'): -8.0,
            ('MEA+, MEACOO-', 'H2O'): -3.036608306709265,
            ('MEA+, MEACOO-', 'MEA'): -4.0}}}


class Test_tau_0(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(default=configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(313)

        return m

    # # TODO: Need to check effects of binary parameters
    # @pytest.mark.unit
    # def test_pure_water(self, model):
    #     # Start by setting all mole fractions to small number
    #     # Using 0 results in division by zero errors
    #     for k in model.state[1].mole_frac_phase_comp:
    #         model.state[1].mole_frac_phase_comp[k].set_value(1e-12)

    #     # Test pure water
    #     model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(1)

    #     # Check mixing expressions
    #     assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(1, rel=1e-8)
    #     assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
    #         1e-12, rel=1e-8)
    #     assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
    #         1e-12, rel=1e-8)

    #     for k, v in model.state[1].Liq_Y.items():
    #         if k == "Cl-":
    #             assert value(v) == pytest.approx(1, rel=1e-8)
    #         else:
    #             assert value(v) == pytest.approx(0.5, rel=1e-8)

    #     for k, v in model.state[1].Liq_alpha.items():
    #         if k == ("H2O", "H2O"):
    #             assert value(v) == 0.3
    #         elif k in [("Na+", "Cl-"), ("K+", "Cl-")]:
    #             assert value(v) == 0.1
    #         elif k in [("Cl-", "Na+"), ("Cl-", "K+")]:
    #             assert value(v) == 0
    #         else:
    #             assert value(v) == 0.2

    #     assert value(model.state[1].Liq_G["Cl-", "H2O"]) == (
    #         0.5*exp(-0.2*-4.5916) + 0.5*exp(-0.2*-4.1341))
    #     assert value(model.state[1].Liq_G["Cl-", "K+"]) == 1
    #     assert value(model.state[1].Liq_G["Cl-", "Na+"]) == 1
    #     assert value(model.state[1].Liq_G["H2O", "Cl-"]) == (
    #         0.5*exp(-0.2*9.0234) + 0.5*exp(-0.2*8.1354))
    #     assert value(model.state[1].Liq_G["H2O", "H2O"]) == 1
    #     assert value(model.state[1].Liq_G["H2O", "K+"]) == exp(-0.2*8.1354)
    #     assert value(model.state[1].Liq_G["H2O", "Na+"]) == exp(-0.2*9.0234)
    #     assert value(model.state[1].Liq_G["K+", "Cl-"]) == 0.5
    #     assert value(model.state[1].Liq_G["K+", "H2O"]) == exp(-0.2*-4.1341)
    #     assert value(model.state[1].Liq_G["Na+", "Cl-"]) == 0.5
    #     assert value(model.state[1].Liq_G["Na+", "H2O"]) == exp(-0.2*-4.5916)

    #     assert value(model.state[1].Liq_tau["Cl-", "H2O"]) == (
    #         -log(0.5*exp(-0.2*-4.5916) + 0.5*exp(-0.2*-4.1341))/0.2)
    #     assert value(model.state[1].Liq_tau["Cl-", "K+"]) == 0
    #     assert value(model.state[1].Liq_tau["Cl-", "Na+"]) == 0
    #     assert value(model.state[1].Liq_tau["H2O", "Cl-"]) == (
    #         -log(0.5*exp(-0.2*9.0234) + 0.5*exp(-0.2*8.1354))/0.2)
    #     assert value(model.state[1].Liq_tau["H2O", "H2O"]) == 0
    #     assert value(model.state[1].Liq_tau["H2O", "K+"]) == 8.1354
    #     assert value(model.state[1].Liq_tau["H2O", "Na+"]) == 9.0234
    #     assert value(model.state[1].Liq_tau["K+", "Cl-"]) == -log(0.5)/0.1
    #     assert value(model.state[1].Liq_tau["K+", "H2O"]) == -4.1341
    #     assert value(model.state[1].Liq_tau["Na+", "Cl-"]) == -log(0.5)/0.1
    #     assert value(model.state[1].Liq_tau["Na+", "H2O"]) == -4.5916

    #     # Check activity coefficient contributions
    #     assert (value(model.state[1].Liq_log_gamma_pdh["H2O"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc["H2O"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_pdh["H2O"] +
    #                   model.state[1].Liq_log_gamma_lc["H2O"]) ==
    #             pytest.approx(0, abs=1e-10))

    # @pytest.mark.unit
    # def test_pure_salt(self, model):
    #     # Test pure NaCl
    #     model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(1e-12)
    #     model.state[1].mole_frac_phase_comp["Liq", "K+"].set_value(0.25)
    #     model.state[1].mole_frac_phase_comp["Liq", "Na+"].set_value(0.25)
    #     model.state[1].mole_frac_phase_comp["Liq", "Cl-"].set_value(0.5)

    #     # Check mixing expressions
    #     assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(
    #         1e-12, rel=1e-8)
    #     assert value(model.state[1].Liq_X["K+"]) == pytest.approx(
    #         0.25, rel=1e-8)
    #     assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
    #         0.25, rel=1e-8)
    #     assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
    #         0.5, rel=1e-8)
    #     assert value(model.state[1].Liq_X["H2O"]) == pytest.approx(
    #         value(model.state[1].Liq_X_ref["H2O"]), rel=1e-8)
    #     assert value(model.state[1].Liq_X["Na+"]) == pytest.approx(
    #         value(model.state[1].Liq_X_ref["Na+"]), rel=1e-8)
    #     assert value(model.state[1].Liq_X["Cl-"]) == pytest.approx(
    #         value(model.state[1].Liq_X_ref["Cl-"]), rel=1e-8)

    #     for k, v in model.state[1].Liq_Y.items():
    #         if k == "Cl-":
    #             assert value(v) == pytest.approx(1, rel=1e-8)
    #         else:
    #             assert value(v) == pytest.approx(0.5, rel=1e-8)

    #     for k, v in model.state[1].Liq_alpha.items():
    #         if k == ("H2O", "H2O"):
    #             assert value(v) == 0.3
    #         elif k in [("Na+", "Cl-"), ("K+", "Cl-")]:
    #             assert value(v) == 0.1
    #         elif k in [("Cl-", "Na+"), ("Cl-", "K+")]:
    #             assert value(v) == 0
    #         else:
    #             assert value(v) == 0.2

    #     assert value(model.state[1].Liq_G["Cl-", "H2O"]) == (
    #         0.5*exp(-0.2*-4.5916) + 0.5*exp(-0.2*-4.1341))
    #     assert value(model.state[1].Liq_G["Cl-", "K+"]) == 1
    #     assert value(model.state[1].Liq_G["Cl-", "Na+"]) == 1
    #     assert value(model.state[1].Liq_G["H2O", "Cl-"]) == (
    #         0.5*exp(-0.2*9.0234) + 0.5*exp(-0.2*8.1354))
    #     assert value(model.state[1].Liq_G["H2O", "H2O"]) == 1
    #     assert value(model.state[1].Liq_G["H2O", "K+"]) == exp(-0.2*8.1354)
    #     assert value(model.state[1].Liq_G["H2O", "Na+"]) == exp(-0.2*9.0234)
    #     assert value(model.state[1].Liq_G["K+", "Cl-"]) == 0.5
    #     assert value(model.state[1].Liq_G["K+", "H2O"]) == exp(-0.2*-4.1341)
    #     assert value(model.state[1].Liq_G["Na+", "Cl-"]) == 0.5
    #     assert value(model.state[1].Liq_G["Na+", "H2O"]) == exp(-0.2*-4.5916)

    #     assert value(model.state[1].Liq_tau["Cl-", "H2O"]) == (
    #         -log(0.5*exp(-0.2*-4.5916) + 0.5*exp(-0.2*-4.1341))/0.2)
    #     assert value(model.state[1].Liq_tau["Cl-", "K+"]) == 0
    #     assert value(model.state[1].Liq_tau["Cl-", "Na+"]) == 0
    #     assert value(model.state[1].Liq_tau["H2O", "Cl-"]) == (
    #         -log(0.5*exp(-0.2*9.0234) + 0.5*exp(-0.2*8.1354))/0.2)
    #     assert value(model.state[1].Liq_tau["H2O", "H2O"]) == 0
    #     assert value(model.state[1].Liq_tau["H2O", "K+"]) == 8.1354
    #     assert value(model.state[1].Liq_tau["H2O", "Na+"]) == 9.0234
    #     assert value(model.state[1].Liq_tau["K+", "Cl-"]) == -log(0.5)/0.1
    #     assert value(model.state[1].Liq_tau["K+", "H2O"]) == -4.1341
    #     assert value(model.state[1].Liq_tau["Na+", "Cl-"]) == -log(0.5)/0.1
    #     assert value(model.state[1].Liq_tau["Na+", "H2O"]) == -4.5916

    #     assert (value(model.state[1].Liq_log_gamma_pdh["Na+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I["Na+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I0["Na+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I["Na+"]) ==
    #             pytest.approx(value(model.state[1].Liq_log_gamma_lc_I0["Na+"]),
    #                           abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc["Na+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma["Na+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_pdh["Na+"] +
    #                   model.state[1].Liq_log_gamma_lc["Na+"]) ==
    #             pytest.approx(model.state[1].Liq_log_gamma["Na+"], abs=1e-10))

    #     assert (value(model.state[1].Liq_log_gamma_pdh["K+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I["K+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I0["K+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I["K+"]) ==
    #             pytest.approx(value(model.state[1].Liq_log_gamma_lc_I0["K+"]),
    #                           abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc["K+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma["K+"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_pdh["K+"] +
    #                   model.state[1].Liq_log_gamma_lc["K+"]) ==
    #             pytest.approx(model.state[1].Liq_log_gamma["K+"], abs=1e-10))

    #     assert (value(model.state[1].Liq_log_gamma_pdh["Cl-"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I["Cl-"]) ==
    #             pytest.approx(6.931471806, rel=1e-8))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I0["Cl-"]) ==
    #             pytest.approx(6.931471806, rel=1e-8))
    #     assert (value(model.state[1].Liq_log_gamma_lc_I["Cl-"]) ==
    #             pytest.approx(value(model.state[1].Liq_log_gamma_lc_I0["Cl-"]),
    #                           abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_lc["Cl-"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma["Cl-"]) ==
    #             pytest.approx(0, abs=1e-10))
    #     assert (value(model.state[1].Liq_log_gamma_pdh["Cl-"] +
    #                   model.state[1].Liq_log_gamma_lc["Cl-"]) ==
    #             pytest.approx(model.state[1].Liq_log_gamma["Cl-"], abs=1e-10))

    @pytest.mark.unit
    def test_1(self, model):
        # Set composition
        model.state[1].mole_frac_phase_comp["Liq", "CO2"].set_value(
            1.714017045e-5)
        model.state[1].mole_frac_phase_comp["Liq", "H2O"].set_value(
            0.83210883)
        model.state[1].mole_frac_phase_comp["Liq", "HCO3-"].set_value(
            0.0084691718)
        model.state[1].mole_frac_phase_comp["Liq", "MEA"].set_value(
            0.0085034521)
        model.state[1].mole_frac_phase_comp["Liq", "MEA+"].set_value(
            0.053123526)
        model.state[1].mole_frac_phase_comp["Liq", "MEACOO-"].set_value(
            0.044654355)

        # model.state[1].Liq_vol_mol_solvent.display()
        # model.state[1].Liq_relative_permittivity_solvent.display()
        # model.state[1].Liq_A_DH.display()
        # model.state[1].Liq_ionic_strength.display()
        # model.state[1].Liq_ionic_strength_ref.display()
        # model.state[1].Liq_log_gamma_pdh.display()

        # model.params.Liq.alpha.display()
        # model.params.Liq.tau.display()

        # model.state[1].Liq_X.display()
        # model.state[1].Liq_Y.display()
        model.state[1].Liq_alpha.display()
        model.state[1].Liq_tau.display()
        model.state[1].Liq_G.display()
        model.state[1].Liq_log_gamma_lc.display()

        assert False
