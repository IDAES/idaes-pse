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

Author: Andrew Lee
"""
import pytest

from pyomo.environ import (ConcreteModel,
                           Expression,
                           exp,
                           value,
                           Var,
                           units as pyunits)

from idaes.core import (AqueousPhase,
                        Solvent,
                        Solute,
                        Apparent,
                        Anion,
                        Cation)
from idaes.generic_models.properties.core.eos.enrtl import ENRTL
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)
from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.reactions.dh_rxn import \
    constant_dh_rxn
from idaes.generic_models.properties.core.reactions.equilibrium_constant import \
    van_t_hoff
from idaes.generic_models.properties.core.reactions.equilibrium_forms import \
    power_law_equil
from idaes.generic_models.properties.core.generic.generic_reaction import (
    ConcentrationForm)
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


configuration = {
    "components": {
        "H2O": {"type": Solvent},
        "C6H12": {"type": Solute},
        "NaCl": {"type": Apparent,
                 "dissociation_species": {"Na+": 1, "Cl-": 1}},
        "HCl": {"type": Apparent,
                "dissociation_species": {"H+": 1, "Cl-": 1}},
        "NaOH": {"type": Apparent,
                 "dissociation_species": {"Na+": 1, "OH-": 1}},
        "Na+": {"type": Cation,
                "charge": +1},
        "H+": {"type": Cation,
               "charge": +1},
        "Cl-": {"type": Anion,
                "charge": -1},
        "OH-": {"type": Anion,
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
    "pressure_ref": 1e5,
    "temperature_ref": 300,
    "inherent_reactions": {
        "h2o_si": {"stoichiometry": {("Liq", "H2O"): -1,
                                     ("Liq", "H+"): 1,
                                     ("Liq", "OH-"): 1},
                   "heat_of_reaction": constant_dh_rxn,
                   "equilibrium_constant": van_t_hoff,
                   "equilibrium_form": power_law_equil,
                   "concentration_form": ConcentrationForm.molarity,
                   "parameter_data": {
                       "reaction_order": {("Liq", "H+"): 1,
                                          ("Liq", "OH-"): 1},
                       "dh_rxn_ref": 1,
                       "k_eq_ref": 1e-14,
                       "T_eq_ref": 350}}}}


class TestParameters(object):
    @pytest.mark.unit
    def test_parameters_no_assignment(self):
        m = ConcreteModel()

        m.params = GenericParameterBlock(default=configuration)

        assert isinstance(m.params.Liq.alpha, Var)
        assert len(m.params.Liq.alpha) == 15
        for (i, j) in m.params.Liq.alpha:
            if i != j:
                assert (j, i) not in m.params.Liq.alpha
            if (i, j) in [
                    ("C6H12", "C6H12"), ("H2O", "H2O"), ("H2O", "C6H12")]:
                assert m.params.Liq.alpha[(i, j)].value == 0.3
                assert m.params.Liq.alpha[(i, j)].fixed
            else:
                assert m.params.Liq.alpha[(i, j)].value == 0.2
                assert m.params.Liq.alpha[(i, j)].fixed

        assert isinstance(m.params.Liq.tau, Var)
        assert len(m.params.Liq.tau) == 25
        for (i, j) in m.params.Liq.tau:
            assert m.params.Liq.tau[(i, j)].value == 0
            assert m.params.Liq.tau[(i, j)].fixed

    @pytest.mark.unit
    def test_parameters_assignment(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "NaCl")] = 0.6
        test_config["parameter_data"]["Liq_tau"] = {}
        test_config["parameter_data"]["Liq_tau"][("H2O", "NaCl")] = 0.1

        m = ConcreteModel()

        m.params = GenericParameterBlock(default=test_config)

        assert isinstance(m.params.Liq.alpha, Var)
        assert len(m.params.Liq.alpha) == 15
        for (i, j) in m.params.Liq.alpha:
            if i != j:
                assert (j, i) not in m.params.Liq.alpha
            if (i, j) == ("H2O", "NaCl"):
                assert m.params.Liq.alpha[(i, j)].value == 0.6
                assert m.params.Liq.alpha[(i, j)].fixed
            elif (i, j) in [
                    ("C6H12", "C6H12"), ("H2O", "H2O"), ("H2O", "C6H12")]:
                assert m.params.Liq.alpha[(i, j)].value == 0.3
                assert m.params.Liq.alpha[(i, j)].fixed
            else:
                assert m.params.Liq.alpha[(i, j)].value == 0.2
                assert m.params.Liq.alpha[(i, j)].fixed

        assert isinstance(m.params.Liq.tau, Var)
        assert len(m.params.Liq.tau) == 25
        for (i, j) in m.params.Liq.tau:
            print(i, j)
            if (i, j) == ("H2O", "NaCl"):
                assert m.params.Liq.tau[(i, j)].value == 0.1
                assert m.params.Liq.tau[(i, j)].fixed
            else:
                assert m.params.Liq.tau[(i, j)].value == 0
                assert m.params.Liq.tau[(i, j)].fixed

    @pytest.mark.unit
    def test_parameters_unsymmetric_alpha(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "NaCl")] = 0.6
        test_config["parameter_data"]["Liq_alpha"][("NaCl", "H2O")] = 0.8

        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params.Liq eNRTL alpha parameter assigned "
                           "non-symmetric value for pair \('H2O', 'NaCl'\). "
                           "Please assign only one value for component pair."):
            m.params = GenericParameterBlock(default=test_config)

    @pytest.mark.unit
    def test_parameters_alpha_symmetry_duplicate(self, caplog):
        caplog.set_level(
            idaeslog.INFO,
            logger=("idaes.generic_models.properties.core."
                    "generic.generic_property"))

        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "NaCl")] = 0.6
        test_config["parameter_data"]["Liq_alpha"][("NaCl", "H2O")] = 0.6

        m = ConcreteModel()

        m.params = GenericParameterBlock(default=test_config)

        assert ("eNRTL alpha value provided for both ('H2O', 'NaCl') and "
                "('NaCl', 'H2O'). It is only necessary to provide a "
                "value for one of these due to symmetry." in caplog.text)

    @pytest.mark.unit
    def test_parameters_alpha_unused_parameter(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "Na+")] = 0.6

        m = ConcreteModel()

        # TODO: Can't get regex to match tuple for some reason
        # For now, just chck the start of the expected string.
        with pytest.raises(ConfigurationError,
                           match="params.Liq eNRTL alpha parameter provided "
                           "for invalid component pair "):
            m.params = GenericParameterBlock(default=test_config)

    @pytest.mark.unit
    def test_parameters_tau_asymmetric(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_tau"] = {}
        test_config["parameter_data"]["Liq_tau"][("H2O", "NaCl")] = 0.1
        test_config["parameter_data"]["Liq_tau"][("NaCl", "H2O")] = -0.1

        m = ConcreteModel()

        m.params = GenericParameterBlock(default=test_config)

        assert isinstance(m.params.Liq.tau, Var)
        assert len(m.params.Liq.tau) == 25
        for (i, j) in m.params.Liq.tau:
            print(i, j)
            if (i, j) == ("H2O", "NaCl"):
                assert m.params.Liq.tau[(i, j)].value == 0.1
                assert m.params.Liq.tau[(i, j)].fixed
            elif (i, j) == ("NaCl", "H2O"):
                assert m.params.Liq.tau[(i, j)].value == -0.1
                assert m.params.Liq.tau[(i, j)].fixed
            else:
                assert m.params.Liq.tau[(i, j)].value == 0
                assert m.params.Liq.tau[(i, j)].fixed

    @pytest.mark.unit
    def test_parameters_tau_unused_parameter(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_tau"] = {}
        test_config["parameter_data"]["Liq_tau"][("H2O", "Na+")] = 0.6

        m = ConcreteModel()

        # TODO: Can't get regex to match tuple for some reason
        # For now, just chck the start of the expected string.
        with pytest.raises(ConfigurationError,
                           match="params.Liq eNRTL tau parameter provided "
                           "for invalid component pair "):
            m.params = GenericParameterBlock(default=test_config)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(default=configuration)

        m.state = m.params.build_state_block([1])

        return m

    @pytest.mark.unit
    def test_common(self, model):
        assert isinstance(model.state[1].Liq_X, Expression)
        assert len(model.state[1].Liq_X) == 6
        for j in model.state[1].Liq_X:
            if j in ["H2O", "C6H12"]:
                # _X should be mole_frac_phase_comp_true
                assert (
                    str(model.state[1].Liq_X[j]._expr) ==
                    str(model.state[1].mole_frac_phase_comp_true["Liq", j]))
            else:
                # _X should be mutiplied by charge
                assert (
                    str(model.state[1].Liq_X[j]._expr) ==
                    str(model.state[1].mole_frac_phase_comp_true["Liq", j] *
                        model.params.get_component(j).config.charge))

        assert isinstance(model.state[1].Liq_Y, Expression)
        assert len(model.state[1].Liq_Y) == 4
        for j in model.state[1].Liq_Y:
            if j in ["H+", "Na+"]:
                assert (str(model.state[1].Liq_Y[j]._expr) ==
                        str(model.state[1].Liq_X[j] /
                            (model.state[1].Liq_X["Na+"] +
                             model.state[1].Liq_X["H+"])))
            else:
                assert (str(model.state[1].Liq_Y[j]._expr) ==
                        str(model.state[1].Liq_X[j] /
                            (model.state[1].Liq_X["Cl-"] +
                             model.state[1].Liq_X["OH-"])))

    @pytest.mark.unit
    def test_alpha(self, model):
        assert isinstance(model.state[1].Liq_alpha, Expression)
        assert len(model.state[1].Liq_alpha) == 28

        # Molecule-molecule interactions
        assert (model.state[1].Liq_alpha["H2O", "H2O"].expr ==
                model.params.Liq.alpha["H2O", "H2O"])
        assert (model.state[1].Liq_alpha["H2O", "C6H12"].expr ==
                model.params.Liq.alpha["H2O", "C6H12"])
        assert (model.state[1].Liq_alpha["C6H12", "C6H12"].expr ==
                model.params.Liq.alpha["C6H12", "C6H12"])
        assert (model.state[1].Liq_alpha["C6H12", "H2O"].expr ==
                model.params.Liq.alpha["H2O", "C6H12"])

        # Molecule-ion interactions
        assert (model.state[1].Liq_alpha["H2O", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["H2O", "NaCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "NaOH"]))
        assert (model.state[1].Liq_alpha["H2O", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["H2O", "HCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "H2O"]))
        assert (model.state[1].Liq_alpha["Na+", "H2O"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["H2O", "NaCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "NaOH"]))
        assert (model.state[1].Liq_alpha["H+", "H2O"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["H2O", "HCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "H2O"]))
        assert (model.state[1].Liq_alpha["H2O", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["H2O", "NaCl"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "HCl"]))
        assert (model.state[1].Liq_alpha["H2O", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["H2O", "NaOH"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "H2O"]))
        assert (model.state[1].Liq_alpha["Cl-", "H2O"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["H2O", "NaCl"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "HCl"]))
        assert (model.state[1].Liq_alpha["OH-", "H2O"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["H2O", "NaOH"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "H2O"]))

        assert (model.state[1].Liq_alpha["C6H12", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["C6H12", "NaCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["C6H12", "NaOH"]))
        assert (model.state[1].Liq_alpha["C6H12", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["C6H12", "HCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "C6H12"]))
        assert (model.state[1].Liq_alpha["Na+", "C6H12"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["C6H12", "NaCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["C6H12", "NaOH"]))
        assert (model.state[1].Liq_alpha["H+", "C6H12"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["C6H12", "HCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "C6H12"]))
        assert (model.state[1].Liq_alpha["C6H12", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["C6H12", "NaCl"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["C6H12", "HCl"]))
        assert (model.state[1].Liq_alpha["C6H12", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["C6H12", "NaOH"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "C6H12"]))
        assert (model.state[1].Liq_alpha["Cl-", "C6H12"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["C6H12", "NaCl"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["C6H12", "HCl"]))
        assert (model.state[1].Liq_alpha["OH-", "C6H12"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["C6H12", "NaOH"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "C6H12"]))

        # Ion-ion interactions
        assert (model.state[1].Liq_alpha["Na+", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["NaCl", "NaCl"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["NaCl", "HCl"]))
        assert (model.state[1].Liq_alpha["Na+", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["NaOH", "NaOH"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "NaOH"]))
        assert (model.state[1].Liq_alpha["H+", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["NaCl", "HCl"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["HCl", "HCl"]))
        assert (model.state[1].Liq_alpha["H+", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 model.params.Liq.alpha["H2O", "NaOH"] +
                 model.state[1].Liq_Y["H+"] *
                 model.params.Liq.alpha["H2O", "H2O"]))
        assert (model.state[1].Liq_alpha["Cl-", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["NaCl", "NaCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["NaCl", "NaOH"]))
        assert (model.state[1].Liq_alpha["Cl-", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["HCl", "HCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "HCl"]))
        assert (model.state[1].Liq_alpha["OH-", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["NaCl", "NaOH"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["NaOH", "NaOH"]))
        assert (model.state[1].Liq_alpha["OH-", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 model.params.Liq.alpha["H2O", "HCl"] +
                 model.state[1].Liq_Y["OH-"] *
                 model.params.Liq.alpha["H2O", "H2O"]))

        # Like-ion interactions
        assert ("Na+", "Na+") not in model.state[1].Liq_alpha
        assert ("Na+", "H+") not in model.state[1].Liq_alpha
        assert ("H+", "Na+") not in model.state[1].Liq_alpha
        assert ("h+", "H+") not in model.state[1].Liq_alpha
        assert ("Cl-", "Cl-") not in model.state[1].Liq_alpha
        assert ("Cl-", "OH-") not in model.state[1].Liq_alpha
        assert ("OH-", "Cl-") not in model.state[1].Liq_alpha
        assert ("OH-", "OH-") not in model.state[1].Liq_alpha

    @pytest.mark.unit
    def test_G(self, model):
        assert isinstance(model.state[1].Liq_G, Expression)
        assert len(model.state[1].Liq_G) == 28

# exp(-alpha[i, j]*tau_rule(b, pobj, i, j, b.temperature))
        # Molecule-molecule interactions
        assert (model.state[1].Liq_G["H2O", "H2O"].expr ==
                exp(-model.params.Liq.alpha["H2O", "H2O"] *
                    model.params.Liq.tau["H2O", "H2O"]))
        assert (model.state[1].Liq_G["H2O", "C6H12"].expr ==
                exp(-model.params.Liq.alpha["H2O", "C6H12"] *
                    model.params.Liq.tau["H2O", "C6H12"]))
        assert (model.state[1].Liq_G["C6H12", "C6H12"].expr ==
                exp(-model.params.Liq.alpha["C6H12", "C6H12"] *
                    model.params.Liq.tau["C6H12", "C6H12"]))
        assert (model.state[1].Liq_G["C6H12", "H2O"].expr ==
                exp(-model.params.Liq.alpha["H2O", "C6H12"] *
                    model.params.Liq.tau["H2O", "C6H12"]))

        # Molecule-ion interactions
        assert (model.state[1].Liq_G["H2O", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["H2O", "NaCl"] *
                     model.params.Liq.tau["H2O", "NaCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "NaOH"] *
                     model.params.Liq.tau["H2O", "NaOH"])))
        assert (model.state[1].Liq_G["H2O", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["H2O", "HCl"] *
                     model.params.Liq.tau["H2O", "HCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "H2O"] *
                     model.params.Liq.tau["H2O", "H2O"])))
        assert (model.state[1].Liq_G["Na+", "H2O"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["H2O", "NaCl"] *
                     model.params.Liq.tau["H2O", "NaCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "NaOH"] *
                     model.params.Liq.tau["H2O", "NaOH"])))
        assert (model.state[1].Liq_G["H+", "H2O"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["H2O", "HCl"] *
                     model.params.Liq.tau["H2O", "HCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "H2O"] *
                     model.params.Liq.tau["H2O", "H2O"])))
        assert (model.state[1].Liq_G["H2O", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["H2O", "NaCl"] *
                     model.params.Liq.tau["H2O", "NaCl"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "HCl"] *
                     model.params.Liq.tau["H2O", "HCl"])))
        assert (model.state[1].Liq_G["H2O", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["H2O", "NaOH"] *
                     model.params.Liq.tau["H2O", "NaOH"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "H2O"] *
                     model.params.Liq.tau["H2O", "H2O"])))
        assert (model.state[1].Liq_G["Cl-", "H2O"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["H2O", "NaCl"] *
                     model.params.Liq.tau["H2O", "NaCl"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "HCl"] *
                     model.params.Liq.tau["H2O", "HCl"])))
        assert (model.state[1].Liq_G["OH-", "H2O"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["H2O", "NaOH"] *
                     model.params.Liq.tau["H2O", "NaOH"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "H2O"] *
                     model.params.Liq.tau["H2O", "H2O"])))

        assert (model.state[1].Liq_G["C6H12", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaCl"] *
                     model.params.Liq.tau["C6H12", "NaCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaOH"] *
                     model.params.Liq.tau["C6H12", "NaOH"])))
        assert (model.state[1].Liq_G["C6H12", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["C6H12", "HCl"] *
                     model.params.Liq.tau["C6H12", "HCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "C6H12"] *
                     model.params.Liq.tau["H2O", "C6H12"])))
        assert (model.state[1].Liq_G["Na+", "C6H12"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaCl"] *
                     model.params.Liq.tau["C6H12", "NaCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaOH"] *
                     model.params.Liq.tau["C6H12", "NaOH"])))
        assert (model.state[1].Liq_G["H+", "C6H12"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["C6H12", "HCl"] *
                     model.params.Liq.tau["C6H12", "HCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "C6H12"] *
                     model.params.Liq.tau["H2O", "C6H12"])))
        assert (model.state[1].Liq_G["C6H12", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaCl"] *
                     model.params.Liq.tau["C6H12", "NaCl"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["C6H12", "HCl"] *
                     model.params.Liq.tau["C6H12", "HCl"])))
        assert (model.state[1].Liq_G["C6H12", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaOH"] *
                     model.params.Liq.tau["C6H12", "NaOH"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "C6H12"] *
                     model.params.Liq.tau["H2O", "C6H12"])))
        assert (model.state[1].Liq_G["Cl-", "C6H12"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaCl"] *
                     model.params.Liq.tau["C6H12", "NaCl"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["C6H12", "HCl"] *
                     model.params.Liq.tau["C6H12", "HCl"])))
        assert (model.state[1].Liq_G["OH-", "C6H12"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["C6H12", "NaOH"] *
                     model.params.Liq.tau["C6H12", "NaOH"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "C6H12"] *
                     model.params.Liq.tau["H2O", "C6H12"])))

        # # Ion-ion interactions
        assert (model.state[1].Liq_G["Na+", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["NaCl", "NaCl"] *
                     model.params.Liq.tau["NaCl", "NaCl"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["NaCl", "HCl"] *
                     model.params.Liq.tau["NaCl", "HCl"])))
        assert (model.state[1].Liq_G["Na+", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["NaOH", "NaOH"] *
                     model.params.Liq.tau["NaOH", "NaOH"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "NaOH"] *
                     model.params.Liq.tau["H2O", "NaOH"])))
        assert (model.state[1].Liq_G["H+", "Cl-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["NaCl", "HCl"] *
                     model.params.Liq.tau["NaCl", "HCl"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["HCl", "HCl"] *
                     model.params.Liq.tau["HCl", "HCl"])))
        assert (model.state[1].Liq_G["H+", "OH-"].expr ==
                (model.state[1].Liq_Y["Na+"] *
                 exp(-model.params.Liq.alpha["H2O", "NaOH"] *
                     model.params.Liq.tau["H2O", "NaOH"]) +
                 model.state[1].Liq_Y["H+"] *
                 exp(-model.params.Liq.alpha["H2O", "H2O"] *
                     model.params.Liq.tau["H2O", "H2O"])))
        assert (model.state[1].Liq_G["Cl-", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["NaCl", "NaCl"] *
                     model.params.Liq.tau["NaCl", "NaCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["NaCl", "NaOH"] *
                     model.params.Liq.tau["NaCl", "NaOH"])))
        assert (model.state[1].Liq_G["Cl-", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["HCl", "HCl"] *
                     model.params.Liq.tau["HCl", "HCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "HCl"] *
                     model.params.Liq.tau["H2O", "HCl"])))
        assert (model.state[1].Liq_G["OH-", "Na+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["NaCl", "NaOH"] *
                     model.params.Liq.tau["NaCl", "NaOH"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["NaOH", "NaOH"] *
                     model.params.Liq.tau["NaOH", "NaOH"])))
        assert (model.state[1].Liq_G["OH-", "H+"].expr ==
                (model.state[1].Liq_Y["Cl-"] *
                 exp(-model.params.Liq.alpha["H2O", "HCl"] *
                     model.params.Liq.tau["H2O", "HCl"]) +
                 model.state[1].Liq_Y["OH-"] *
                 exp(-model.params.Liq.alpha["H2O", "H2O"] *
                     model.params.Liq.tau["H2O", "H2O"])))

        # Like-ion interactions
        assert ("Na+", "Na+") not in model.state[1].Liq_G
        assert ("Na+", "H+") not in model.state[1].Liq_G
        assert ("H+", "Na+") not in model.state[1].Liq_G
        assert ("h+", "H+") not in model.state[1].Liq_G
        assert ("Cl-", "Cl-") not in model.state[1].Liq_G
        assert ("Cl-", "OH-") not in model.state[1].Liq_G
        assert ("OH-", "Cl-") not in model.state[1].Liq_G
        assert ("OH-", "OH-") not in model.state[1].Liq_G
