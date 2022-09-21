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

Author: Andrew Lee
"""
import pytest

from pyomo.environ import (
    ConcreteModel,
    Expression,
    exp,
    log,
    Set,
    units as pyunits,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_equivalent

from idaes.core import AqueousPhase, Solvent, Solute, Apparent, Anion, Cation
from idaes.core.util.constants import Constants
from idaes.models.properties.modular_properties.eos.enrtl import ENRTL
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    StateIndex,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.pure.electrolyte import (
    relative_permittivity_constant,
)
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog


def dummy_method(b, *args, **kwargs):
    return 42 * pyunits.mol / pyunits.m**3


configuration = {
    "components": {
        "H2O": {
            "type": Solvent,
            "dens_mol_liq_comp": dummy_method,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (18e-3, pyunits.kg / pyunits.mol),
                "relative_permittivity_liq_comp": 101,
            },
        },
        "C6H12": {
            "type": Solute,
            "dens_mol_liq_comp": dummy_method,
            "relative_permittivity_liq_comp": relative_permittivity_constant,
            "parameter_data": {
                "mw": (84e-3, pyunits.kg / pyunits.mol),
                "relative_permittivity_liq_comp": 102,
            },
        },
        "NaCl": {"type": Apparent, "dissociation_species": {"Na+": 1, "Cl-": 1}},
        "HCl": {"type": Apparent, "dissociation_species": {"H+": 1, "Cl-": 1}},
        "NaOH": {"type": Apparent, "dissociation_species": {"Na+": 1, "OH-": 1}},
        "Na+": {"type": Cation, "charge": +1},
        "H+": {"type": Cation, "charge": +1},
        "Cl-": {"type": Anion, "charge": -1},
        "OH-": {"type": Anion, "charge": -1},
    },
    "phases": {"Liq": {"type": AqueousPhase, "equation_of_state": ENRTL}},
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
}


class TestParameters(object):
    @pytest.mark.unit
    def test_parameters_no_assignment(self):
        m = ConcreteModel()

        m.params = GenericParameterBlock(**configuration)

        assert isinstance(m.params.Liq.ion_pair_set, Set)
        assert len(m.params.Liq.ion_pair_set) == 4
        for p in m.params.Liq.ion_pair_set:
            assert p in [("Na+, Cl-"), ("Na+, OH-"), ("H+, Cl-"), ("H+, OH-")]

        assert isinstance(m.params.Liq.component_pair_set, Set)
        assert len(m.params.Liq.component_pair_set) == 32
        assert isinstance(m.params.Liq.component_pair_set_symmetric, Set)
        assert len(m.params.Liq.component_pair_set_symmetric) == 17

        assert isinstance(m.params.Liq.alpha, Var)
        assert len(m.params.Liq.alpha) == 17
        for (i, j) in m.params.Liq.alpha:
            if i != j:
                assert (j, i) not in m.params.Liq.alpha
            if (i, j) in [("C6H12", "C6H12"), ("H2O", "H2O"), ("H2O", "C6H12")]:
                assert m.params.Liq.alpha[(i, j)].value == 0.3
                assert m.params.Liq.alpha[(i, j)].fixed
            else:
                assert m.params.Liq.alpha[(i, j)].value == 0.2
                assert m.params.Liq.alpha[(i, j)].fixed

        assert isinstance(m.params.Liq.tau, Var)
        assert len(m.params.Liq.tau) == 32
        for (i, j) in m.params.Liq.tau:
            assert m.params.Liq.tau[(i, j)].value == 0
            assert m.params.Liq.tau[(i, j)].fixed

    @pytest.mark.unit
    def test_parameters_assignment(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "Na+, Cl-")] = 0.6
        test_config["parameter_data"]["Liq_tau"] = {}
        test_config["parameter_data"]["Liq_tau"][("H2O", "Na+, Cl-")] = 0.1

        m = ConcreteModel()

        m.params = GenericParameterBlock(**test_config)

        assert isinstance(m.params.Liq.alpha, Var)
        assert len(m.params.Liq.alpha) == 17
        for (i, j) in m.params.Liq.alpha:
            if i != j:
                assert (j, i) not in m.params.Liq.alpha
            if (i, j) == ("H2O", "Na+, Cl-"):
                assert m.params.Liq.alpha[(i, j)].value == 0.6
                assert m.params.Liq.alpha[(i, j)].fixed
            elif (i, j) in [("C6H12", "C6H12"), ("H2O", "H2O"), ("H2O", "C6H12")]:
                assert m.params.Liq.alpha[(i, j)].value == 0.3
                assert m.params.Liq.alpha[(i, j)].fixed
            else:
                assert m.params.Liq.alpha[(i, j)].value == 0.2
                assert m.params.Liq.alpha[(i, j)].fixed

        assert isinstance(m.params.Liq.tau, Var)
        assert len(m.params.Liq.tau) == 32
        for (i, j) in m.params.Liq.tau:
            print(i, j)
            if (i, j) == ("H2O", "Na+, Cl-"):
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
        test_config["parameter_data"]["Liq_alpha"][("H2O", "Na+, Cl-")] = 0.6
        test_config["parameter_data"]["Liq_alpha"][("Na+, Cl-", "H2O")] = 0.8

        m = ConcreteModel()

        # TODO: Having trouble getting regex to match component tuple
        # Using a wildcard for now
        with pytest.raises(
            ConfigurationError,
            match="params.Liq eNRTL alpha parameter assigned "
            "non-symmetric value for pair (.+?). Please assign "
            "only one value for component pair.",
        ):
            m.params = GenericParameterBlock(**test_config)

    @pytest.mark.unit
    def test_parameters_alpha_symmetry_duplicate(self, caplog):
        caplog.set_level(
            idaeslog.INFO,
            logger=(
                "idaes.models.properties.modular_properties." "eos.enrtl_parameters"
            ),
        )

        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "Na+, Cl-")] = 0.6
        test_config["parameter_data"]["Liq_alpha"][("Na+, Cl-", "H2O")] = 0.6

        m = ConcreteModel()

        m.params = GenericParameterBlock(**test_config)

        assert (
            "eNRTL alpha value provided for both ('H2O', 'Na+, Cl-') and "
            "('Na+, Cl-', 'H2O'). It is only necessary to provide a "
            "value for one of these due to symmetry." in caplog.text
        )

    @pytest.mark.unit
    def test_parameters_alpha_unused_parameter(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_alpha"] = {}
        test_config["parameter_data"]["Liq_alpha"][("H2O", "Na+")] = 0.6

        m = ConcreteModel()

        # TODO: Having trouble getting regex to match component tuple
        # Using a wildcard for now
        with pytest.raises(
            ConfigurationError,
            match="params.Liq eNRTL alpha parameter provided "
            "for invalid component pair (.+?). Please check "
            "typing and only provide parameters for valid "
            "species pairs.",
        ):
            m.params = GenericParameterBlock(**test_config)

    @pytest.mark.unit
    def test_parameters_tau_asymmetric(self):
        test_config = dict(configuration)
        test_config["parameter_data"] = {}
        test_config["parameter_data"]["Liq_tau"] = {}
        test_config["parameter_data"]["Liq_tau"][("H2O", "Na+, Cl-")] = 0.1
        test_config["parameter_data"]["Liq_tau"][("Na+, Cl-", "H2O")] = -0.1

        m = ConcreteModel()

        m.params = GenericParameterBlock(**test_config)

        assert isinstance(m.params.Liq.tau, Var)
        assert len(m.params.Liq.tau) == 32
        for (i, j) in m.params.Liq.tau:
            print(i, j)
            if (i, j) == ("H2O", "Na+, Cl-"):
                assert m.params.Liq.tau[(i, j)].value == 0.1
                assert m.params.Liq.tau[(i, j)].fixed
            elif (i, j) == ("Na+, Cl-", "H2O"):
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

        # TODO: Having trouble getting regex to match component tuple
        # Using a wildcard for now
        with pytest.raises(
            ConfigurationError,
            match="params.Liq eNRTL tau parameter provided for "
            "invalid component pair (.+?). Please check typing "
            "and only provide parameters for valid species "
            "pairs.",
        ):
            m.params = GenericParameterBlock(**test_config)


class TestStateBlockSymmetric(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(**configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(300)

        return m

    @pytest.mark.unit
    def test_common(self, model):
        # Reference state composition
        assert isinstance(model.state[1].Liq_x_ref, Expression)
        assert len(model.state[1].Liq_x_ref) == 6
        for k in model.state[1].Liq_x_ref:
            assert k in ["H2O", "C6H12", "Na+", "H+", "Cl-", "OH-"]
            if k in ["H2O", "C6H12"]:
                assert str(model.state[1].Liq_x_ref[k].expr) == str(0.0)
            else:
                assert str(model.state[1].Liq_x_ref[k].expr) == str(
                    model.state[1].mole_frac_phase_comp_true["Liq", k]
                    / (
                        model.state[1].mole_frac_phase_comp_true["Liq", "Cl-"]
                        + model.state[1].mole_frac_phase_comp_true["Liq", "OH-"]
                        + model.state[1].mole_frac_phase_comp_true["Liq", "Na+"]
                        + model.state[1].mole_frac_phase_comp_true["Liq", "H+"]
                    )
                )

        assert isinstance(model.state[1].Liq_X, Expression)
        assert len(model.state[1].Liq_X) == 6
        for j in model.state[1].Liq_X:
            if j in ["H2O", "C6H12"]:
                # _X should be mole_frac_phase_comp_true
                assert str(model.state[1].Liq_X[j]._expr) == str(
                    model.state[1].mole_frac_phase_comp_true["Liq", j]
                )
            else:
                # _X should be mutiplied by |charge|
                assert str(model.state[1].Liq_X[j]._expr) == str(
                    model.state[1].mole_frac_phase_comp_true["Liq", j]
                    * abs(model.params.get_component(j).config.charge)
                )

        assert isinstance(model.state[1].Liq_X_ref, Expression)
        assert len(model.state[1].Liq_X_ref) == 6
        for j in model.state[1].Liq_X_ref:
            if j in ["H2O", "C6H12"]:
                # _X should be mole_frac_phase_comp_true
                assert str(model.state[1].Liq_X_ref[j].expr) == str(
                    model.state[1].Liq_x_ref[j]
                )
            else:
                # _X should be mutiplied by |charge|
                assert str(model.state[1].Liq_X_ref[j]._expr) == str(
                    model.state[1].Liq_x_ref[j]
                    * abs(model.params.get_component(j).config.charge)
                )

        assert isinstance(model.state[1].Liq_Y, Expression)
        assert len(model.state[1].Liq_Y) == 4
        for j in model.state[1].Liq_Y:
            if j in ["H+", "Na+"]:
                assert str(model.state[1].Liq_Y[j]._expr) == str(
                    model.state[1].Liq_X[j]
                    / (model.state[1].Liq_X["Na+"] + model.state[1].Liq_X["H+"])
                )
            else:
                assert str(model.state[1].Liq_Y[j]._expr) == str(
                    model.state[1].Liq_X[j]
                    / (model.state[1].Liq_X["Cl-"] + model.state[1].Liq_X["OH-"])
                )

        assert isinstance(model.state[1].Liq_ionic_strength, Expression)
        assert len(model.state[1].Liq_ionic_strength) == 1
        assert str(model.state[1].Liq_ionic_strength.expr) == str(
            0.5
            * (
                model.params.get_component("Cl-").config.charge ** 2
                * model.state[1].mole_frac_phase_comp_true["Liq", "Cl-"]
                + model.params.get_component("OH-").config.charge ** 2
                * model.state[1].mole_frac_phase_comp_true["Liq", "OH-"]
                + model.params.get_component("Na+").config.charge ** 2
                * model.state[1].mole_frac_phase_comp_true["Liq", "Na+"]
                + model.params.get_component("H+").config.charge ** 2
                * model.state[1].mole_frac_phase_comp_true["Liq", "H+"]
            )
        )

        assert isinstance(model.state[1].Liq_ionic_strength_ref, Expression)
        assert len(model.state[1].Liq_ionic_strength_ref) == 1
        assert str(model.state[1].Liq_ionic_strength_ref.expr) == str(
            0.5
            * (
                model.params.get_component("Cl-").config.charge ** 2
                * model.state[1].Liq_x_ref["Cl-"]
                + model.params.get_component("OH-").config.charge ** 2
                * model.state[1].Liq_x_ref["OH-"]
                + model.params.get_component("Na+").config.charge ** 2
                * model.state[1].Liq_x_ref["Na+"]
                + model.params.get_component("H+").config.charge ** 2
                * model.state[1].Liq_x_ref["H+"]
            )
        )

        assert isinstance(model.state[1].Liq_vol_mol_solvent, Expression)
        assert len(model.state[1].Liq_vol_mol_solvent) == 1
        assert str(model.state[1].Liq_vol_mol_solvent.expr) == "1/(42*mol/m**3)"

        assert isinstance(model.state[1].Liq_relative_permittivity_solvent, Expression)
        assert len(model.state[1].Liq_relative_permittivity_solvent) == 1
        assert str(model.state[1].Liq_relative_permittivity_solvent.expr) == (
            str(model.params.get_component("H2O").relative_permittivity_liq_comp)
        )

        assert isinstance(model.state[1].Liq_A_DH, Expression)
        assert len(model.state[1].Liq_A_DH) == 1
        assert_units_equivalent(model.state[1].Liq_A_DH, pyunits.dimensionless)
        assert str(model.state[1].Liq_A_DH.expr) == str(
            (1 / 3)
            * (
                2
                * Constants.pi
                * Constants.avogadro_number
                / model.state[1].Liq_vol_mol_solvent
            )
            ** 0.5
            * (
                Constants.elemental_charge**2
                / (
                    4
                    * Constants.pi
                    * model.state[1].Liq_relative_permittivity_solvent
                    * Constants.vacuum_electric_permittivity
                    * Constants.boltzmann_constant
                    * model.state[1].temperature
                )
            )
            ** (3 / 2)
        )

        assert isinstance(model.state[1].Liq_log_gamma_pdh, Expression)
        assert len(model.state[1].Liq_log_gamma_pdh) == 6
        for j in model.state[1].Liq_log_gamma_pdh:
            assert j in ["H2O", "C6H12", "Na+", "H+", "Cl-", "OH-"]
            if j in ["H2O", "C6H12"]:
                assert str(model.state[1].Liq_log_gamma_pdh[j].expr) == str(
                    (
                        2
                        * model.state[1].Liq_A_DH
                        * model.state[1].Liq_ionic_strength ** (3 / 2)
                        / (1 + 14.9 * model.state[1].Liq_ionic_strength ** (1 / 2))
                    )
                )
            else:

                def ndxdn(j, k):
                    if j == k:
                        return (1 - model.state[1].Liq_x_ref[k]) / (
                            model.state[1].mole_frac_phase_comp_true["Liq", "Cl-"]
                            + model.state[1].mole_frac_phase_comp_true["Liq", "OH-"]
                            + model.state[1].mole_frac_phase_comp_true["Liq", "Na+"]
                            + model.state[1].mole_frac_phase_comp_true["Liq", "H+"]
                        )
                    else:
                        return -model.state[1].Liq_x_ref[k] / (
                            model.state[1].mole_frac_phase_comp_true["Liq", "Cl-"]
                            + model.state[1].mole_frac_phase_comp_true["Liq", "OH-"]
                            + model.state[1].mole_frac_phase_comp_true["Liq", "Na+"]
                            + model.state[1].mole_frac_phase_comp_true["Liq", "H+"]
                        )

                assert str(model.state[1].Liq_log_gamma_pdh[j].expr) == str(
                    -model.state[1].Liq_A_DH
                    * (
                        (2 * model.params.get_component(j).config.charge ** 2 / 14.9)
                        * log(
                            (1 + 14.9 * model.state[1].Liq_ionic_strength ** 0.5)
                            / (1 + 14.9 * model.state[1].Liq_ionic_strength_ref ** 0.5)
                        )
                        + (
                            model.params.get_component(j).config.charge ** 2
                            * model.state[1].Liq_ionic_strength ** 0.5
                            - 2 * model.state[1].Liq_ionic_strength ** 1.5
                        )
                        / (1 + 14.9 * model.state[1].Liq_ionic_strength ** 0.5)
                        - (
                            2
                            * model.state[1].Liq_ionic_strength
                            * model.state[1].Liq_ionic_strength_ref ** -0.5
                        )
                        / (1 + 14.9 * model.state[1].Liq_ionic_strength_ref ** 0.5)
                        * (
                            0.5
                            * (
                                model.params.get_component("Cl-").config.charge ** 2
                                * ndxdn(j, "Cl-")
                                + model.params.get_component("OH-").config.charge ** 2
                                * ndxdn(j, "OH-")
                                + model.params.get_component("Na+").config.charge ** 2
                                * ndxdn(j, "Na+")
                                + model.params.get_component("H+").config.charge ** 2
                                * ndxdn(j, "H+")
                            )
                        )
                    )
                )

        assert isinstance(model.state[1].Liq_log_gamma_lc_I, Expression)
        assert len(model.state[1].Liq_log_gamma_lc_I) == 6
        for k in model.state[1].Liq_log_gamma_lc_I:
            assert k in ["H2O", "C6H12", "Na+", "H+", "Cl-", "OH-"]

        assert isinstance(model.state[1].Liq_log_gamma_lc_I0, Expression)
        assert len(model.state[1].Liq_log_gamma_lc_I0) == 4
        for k in model.state[1].Liq_log_gamma_lc_I0:
            assert k in ["Na+", "H+", "Cl-", "OH-"]
            assert str(model.state[1].Liq_log_gamma_lc_I0[k].expr) != str(
                model.state[1].Liq_log_gamma_lc_I[k].expr
            )

        assert isinstance(model.state[1].Liq_log_gamma_lc, Expression)
        assert len(model.state[1].Liq_log_gamma_lc) == 6
        for k in model.state[1].Liq_log_gamma_lc:
            assert k in ["H2O", "C6H12", "Na+", "H+", "Cl-", "OH-"]
            if k in ["H2O", "C6H12"]:
                assert str(model.state[1].Liq_log_gamma_lc[k].expr) == str(
                    model.state[1].Liq_log_gamma_lc_I[k]
                )
            else:
                assert str(model.state[1].Liq_log_gamma_lc[k].expr) == str(
                    model.state[1].Liq_log_gamma_lc_I[k]
                    - model.state[1].Liq_log_gamma_lc_I0[k]
                )

        assert isinstance(model.state[1].Liq_log_gamma, Expression)
        assert len(model.state[1].Liq_log_gamma) == 6
        for k, v in model.state[1].Liq_log_gamma.items():
            assert str(model.state[1].Liq_log_gamma[k].expr) == str(
                model.state[1].Liq_log_gamma_pdh[k] + model.state[1].Liq_log_gamma_lc[k]
            )

    @pytest.mark.unit
    def test_alpha(self, model):
        assert isinstance(model.state[1].Liq_alpha, Expression)
        assert len(model.state[1].Liq_alpha) == 28

        # Molecule-molecule interactions
        assert str(model.state[1].Liq_alpha["H2O", "H2O"].expr) == str(
            model.params.Liq.alpha["H2O", "H2O"]
        )
        assert str(model.state[1].Liq_alpha["H2O", "C6H12"].expr) == str(
            model.params.Liq.alpha["H2O", "C6H12"]
        )
        assert str(model.state[1].Liq_alpha["C6H12", "H2O"].expr) == str(
            model.params.Liq.alpha["H2O", "C6H12"]
        )
        assert str(model.state[1].Liq_alpha["C6H12", "C6H12"].expr) == str(
            model.params.Liq.alpha["C6H12", "C6H12"]
        )

        # Molecule-ion interactions
        assert str(model.state[1].Liq_alpha["H2O", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * model.params.Liq.alpha["H2O", "Na+, Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["H2O", "Na+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["H2O", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * model.params.Liq.alpha["H2O", "H+, Cl-"]
                + model.state[1].Liq_Y["OH-"] * model.params.Liq.alpha["H2O", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["Na+", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * model.params.Liq.alpha["H2O", "Na+, Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["H2O", "Na+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["H+", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * model.params.Liq.alpha["H2O", "H+, Cl-"]
                + model.state[1].Liq_Y["OH-"] * model.params.Liq.alpha["H2O", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["H2O", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"] * model.params.Liq.alpha["H2O", "Na+, Cl-"]
                + model.state[1].Liq_Y["H+"] * model.params.Liq.alpha["H2O", "H+, Cl-"]
            )
        )
        assert str(model.state[1].Liq_alpha["H2O", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"] * model.params.Liq.alpha["H2O", "Na+, OH-"]
                + model.state[1].Liq_Y["H+"] * model.params.Liq.alpha["H2O", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["Cl-", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"] * model.params.Liq.alpha["H2O", "Na+, Cl-"]
                + model.state[1].Liq_Y["H+"] * model.params.Liq.alpha["H2O", "H+, Cl-"]
            )
        )
        assert str(model.state[1].Liq_alpha["OH-", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"] * model.params.Liq.alpha["H2O", "Na+, OH-"]
                + model.state[1].Liq_Y["H+"] * model.params.Liq.alpha["H2O", "H+, OH-"]
            )
        )

        assert str(model.state[1].Liq_alpha["C6H12", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["C6H12", "Na+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["C6H12", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * model.params.Liq.alpha["C6H12", "H+, Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["C6H12", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["Na+", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["C6H12", "Na+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["H+", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * model.params.Liq.alpha["C6H12", "H+, Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["C6H12", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["C6H12", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                + model.state[1].Liq_Y["H+"]
                * model.params.Liq.alpha["C6H12", "H+, Cl-"]
            )
        )
        assert str(model.state[1].Liq_alpha["C6H12", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * model.params.Liq.alpha["C6H12", "Na+, OH-"]
                + model.state[1].Liq_Y["H+"]
                * model.params.Liq.alpha["C6H12", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["Cl-", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                + model.state[1].Liq_Y["H+"]
                * model.params.Liq.alpha["C6H12", "H+, Cl-"]
            )
        )
        assert str(model.state[1].Liq_alpha["OH-", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * model.params.Liq.alpha["C6H12", "Na+, OH-"]
                + model.state[1].Liq_Y["H+"]
                * model.params.Liq.alpha["C6H12", "H+, OH-"]
            )
        )

        # Ion-ion interactions
        assert str(model.state[1].Liq_alpha["Na+", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"] * 0.2
                + model.state[1].Liq_Y["H+"]
                * model.params.Liq.alpha["Na+, Cl-", "H+, Cl-"]
            )
        )
        assert str(model.state[1].Liq_alpha["Na+", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"] * 0.2
                + model.state[1].Liq_Y["H+"]
                * model.params.Liq.alpha["Na+, OH-", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["H+", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * model.params.Liq.alpha["Na+, Cl-", "H+, Cl-"]
                + model.state[1].Liq_Y["H+"] * 0.2
            )
        )
        assert str(model.state[1].Liq_alpha["H+", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * model.params.Liq.alpha["Na+, OH-", "H+, OH-"]
                + model.state[1].Liq_Y["H+"] * 0.2
            )
        )
        assert str(model.state[1].Liq_alpha["Cl-", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * 0.2
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["Na+, Cl-", "Na+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["Cl-", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"] * 0.2
                + model.state[1].Liq_Y["OH-"]
                * model.params.Liq.alpha["H+, Cl-", "H+, OH-"]
            )
        )
        assert str(model.state[1].Liq_alpha["OH-", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * model.params.Liq.alpha["Na+, Cl-", "Na+, OH-"]
                + model.state[1].Liq_Y["OH-"] * 0.2
            )
        )
        assert str(model.state[1].Liq_alpha["OH-", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * model.params.Liq.alpha["H+, Cl-", "H+, OH-"]
                + model.state[1].Liq_Y["OH-"] * 0.2
            )
        )

        # Like species interactions
        assert ("Na+", "Na+") not in model.state[1].Liq_alpha
        assert ("Na+", "H+") not in model.state[1].Liq_alpha
        assert ("H+", "Na+") not in model.state[1].Liq_alpha
        assert ("H+", "H+") not in model.state[1].Liq_alpha
        assert ("Cl-", "Cl-") not in model.state[1].Liq_alpha
        assert ("Cl-", "OH-") not in model.state[1].Liq_alpha
        assert ("OH-", "Cl-") not in model.state[1].Liq_alpha
        assert ("OH-", "OH-") not in model.state[1].Liq_alpha

    @pytest.mark.unit
    def test_G(self, model):
        assert isinstance(model.state[1].Liq_G, Expression)
        assert len(model.state[1].Liq_G) == 28

        # Molecule-molecule interactions
        assert str(model.state[1].Liq_G["H2O", "H2O"].expr) == str(1.0)
        assert str(model.state[1].Liq_G["H2O", "C6H12"].expr) == str(
            exp(
                -model.params.Liq.alpha["H2O", "C6H12"]
                * model.params.Liq.tau["H2O", "C6H12"]
            )
        )
        assert str(model.state[1].Liq_G["C6H12", "H2O"].expr) == str(
            exp(
                -model.params.Liq.alpha["H2O", "C6H12"]
                * model.params.Liq.tau["C6H12", "H2O"]
            )
        )
        assert str(model.state[1].Liq_G["C6H12", "C6H12"].expr) == str(1.0)

        # Molecule-ion interactions
        assert str(model.state[1].Liq_G["H2O", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, Cl-"]
                    * model.params.Liq.tau["H2O", "Na+, Cl-"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, OH-"]
                    * model.params.Liq.tau["H2O", "Na+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["H2O", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, Cl-"]
                    * model.params.Liq.tau["H2O", "H+, Cl-"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, OH-"]
                    * model.params.Liq.tau["H2O", "H+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["Na+", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, Cl-"]
                    * model.params.Liq.tau["Na+, Cl-", "H2O"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, OH-"]
                    * model.params.Liq.tau["Na+, OH-", "H2O"]
                )
            )
        )
        assert str(model.state[1].Liq_G["H+", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, Cl-"]
                    * model.params.Liq.tau["H+, Cl-", "H2O"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, OH-"]
                    * model.params.Liq.tau["H+, OH-", "H2O"]
                )
            )
        )
        assert str(model.state[1].Liq_G["H2O", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, Cl-"]
                    * model.params.Liq.tau["H2O", "Na+, Cl-"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, Cl-"]
                    * model.params.Liq.tau["H2O", "H+, Cl-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["H2O", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, OH-"]
                    * model.params.Liq.tau["H2O", "Na+, OH-"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, OH-"]
                    * model.params.Liq.tau["H2O", "H+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["Cl-", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, Cl-"]
                    * model.params.Liq.tau["Na+, Cl-", "H2O"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, Cl-"]
                    * model.params.Liq.tau["H+, Cl-", "H2O"]
                )
            )
        )
        assert str(model.state[1].Liq_G["OH-", "H2O"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "Na+, OH-"]
                    * model.params.Liq.tau["Na+, OH-", "H2O"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["H2O", "H+, OH-"]
                    * model.params.Liq.tau["H+, OH-", "H2O"]
                )
            )
        )

        assert str(model.state[1].Liq_G["C6H12", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                    * model.params.Liq.tau["C6H12", "Na+, Cl-"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, OH-"]
                    * model.params.Liq.tau["C6H12", "Na+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["C6H12", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, Cl-"]
                    * model.params.Liq.tau["C6H12", "H+, Cl-"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, OH-"]
                    * model.params.Liq.tau["C6H12", "H+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["Na+", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                    * model.params.Liq.tau["Na+, Cl-", "C6H12"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, OH-"]
                    * model.params.Liq.tau["Na+, OH-", "C6H12"]
                )
            )
        )
        assert str(model.state[1].Liq_G["H+", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, Cl-"]
                    * model.params.Liq.tau["H+, Cl-", "C6H12"]
                )
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, OH-"]
                    * model.params.Liq.tau["H+, OH-", "C6H12"]
                )
            )
        )
        assert str(model.state[1].Liq_G["C6H12", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                    * model.params.Liq.tau["C6H12", "Na+, Cl-"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, Cl-"]
                    * model.params.Liq.tau["C6H12", "H+, Cl-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["C6H12", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, OH-"]
                    * model.params.Liq.tau["C6H12", "Na+, OH-"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, OH-"]
                    * model.params.Liq.tau["C6H12", "H+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["Cl-", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, Cl-"]
                    * model.params.Liq.tau["Na+, Cl-", "C6H12"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, Cl-"]
                    * model.params.Liq.tau["H+, Cl-", "C6H12"]
                )
            )
        )
        assert str(model.state[1].Liq_G["OH-", "C6H12"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "Na+, OH-"]
                    * model.params.Liq.tau["Na+, OH-", "C6H12"]
                )
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["C6H12", "H+, OH-"]
                    * model.params.Liq.tau["H+, OH-", "C6H12"]
                )
            )
        )

        # Ion-ion interactions
        assert str(model.state[1].Liq_G["Na+", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["Na+, Cl-", "H+, Cl-"]
                    * model.params.Liq.tau["Na+, Cl-", "H+, Cl-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["Na+", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                + model.state[1].Liq_Y["H+"]
                * exp(
                    -model.params.Liq.alpha["Na+, OH-", "H+, OH-"]
                    * model.params.Liq.tau["Na+, OH-", "H+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["H+", "Cl-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["Na+, Cl-", "H+, Cl-"]
                    * model.params.Liq.tau["H+, Cl-", "Na+, Cl-"]
                )
                + model.state[1].Liq_Y["H+"]
            )
        )
        assert str(model.state[1].Liq_G["H+", "OH-"].expr) == str(
            (
                model.state[1].Liq_Y["Na+"]
                * exp(
                    -model.params.Liq.alpha["Na+, OH-", "H+, OH-"]
                    * model.params.Liq.tau["H+, OH-", "Na+, OH-"]
                )
                + model.state[1].Liq_Y["H+"]
            )
        )
        assert str(model.state[1].Liq_G["Cl-", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["Na+, Cl-", "Na+, OH-"]
                    * model.params.Liq.tau["Na+, Cl-", "Na+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["Cl-", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                + model.state[1].Liq_Y["OH-"]
                * exp(
                    -model.params.Liq.alpha["H+, Cl-", "H+, OH-"]
                    * model.params.Liq.tau["H+, Cl-", "H+, OH-"]
                )
            )
        )
        assert str(model.state[1].Liq_G["OH-", "Na+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["Na+, Cl-", "Na+, OH-"]
                    * model.params.Liq.tau["Na+, OH-", "Na+, Cl-"]
                )
                + model.state[1].Liq_Y["OH-"]
            )
        )
        assert str(model.state[1].Liq_G["OH-", "H+"].expr) == str(
            (
                model.state[1].Liq_Y["Cl-"]
                * exp(
                    -model.params.Liq.alpha["H+, Cl-", "H+, OH-"]
                    * model.params.Liq.tau["H+, OH-", "H+, Cl-"]
                )
                + model.state[1].Liq_Y["OH-"]
            )
        )

        # Like species interactions
        assert ("Na+", "Na+") not in model.state[1].Liq_G
        assert ("Na+", "H+") not in model.state[1].Liq_G
        assert ("H+", "Na+") not in model.state[1].Liq_G
        assert ("H+", "H+") not in model.state[1].Liq_G
        assert ("Cl-", "Cl-") not in model.state[1].Liq_G
        assert ("Cl-", "OH-") not in model.state[1].Liq_G
        assert ("OH-", "Cl-") not in model.state[1].Liq_G
        assert ("OH-", "OH-") not in model.state[1].Liq_G

    @pytest.mark.unit
    def test_tau(self, model):
        assert isinstance(model.state[1].Liq_tau, Expression)
        assert len(model.state[1].Liq_tau) == 28

        # Molecule-molecule interactions
        assert str(model.state[1].Liq_tau["H2O", "H2O"].expr) == str(
            model.params.Liq.tau["H2O", "H2O"]
        )
        assert str(model.state[1].Liq_tau["H2O", "C6H12"].expr) == str(
            model.params.Liq.tau["H2O", "C6H12"]
        )
        assert str(model.state[1].Liq_tau["C6H12", "H2O"].expr) == str(
            model.params.Liq.tau["C6H12", "H2O"]
        )
        assert str(model.state[1].Liq_tau["C6H12", "C6H12"].expr) == str(
            model.params.Liq.tau["C6H12", "C6H12"]
        )

        for i, j in model.state[1].Liq_tau:
            if (i, j) not in [
                ("H2O", "H2O"),
                ("H2O", "C6H12"),
                ("C6H12", "H2O"),
                ("C6H12", "C6H12"),
            ]:
                assert str(model.state[1].Liq_tau[i, j].expr) == str(
                    -log(model.state[1].Liq_G[i, j]) / model.state[1].Liq_alpha[i, j]
                )

        # Like species interactions
        assert ("Na+", "Na+") not in model.state[1].Liq_tau
        assert ("Na+", "H+") not in model.state[1].Liq_tau
        assert ("H+", "Na+") not in model.state[1].Liq_tau
        assert ("H+", "H+") not in model.state[1].Liq_tau
        assert ("Cl-", "Cl-") not in model.state[1].Liq_tau
        assert ("Cl-", "OH-") not in model.state[1].Liq_tau
        assert ("OH-", "Cl-") not in model.state[1].Liq_tau
        assert ("OH-", "OH-") not in model.state[1].Liq_tau


class TestProperties(object):
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.params = GenericParameterBlock(**configuration)

        m.state = m.params.build_state_block([1])

        # Need to set a value of T for checking expressions later
        m.state[1].temperature.set_value(300)
        # Need to set a value of log_act_phase_solvents for checking osmotic pressure calculation
        m.state[1].log_act_phase_solvents["Liq"].set_value(-1.789578)

        return m

    @pytest.mark.unit
    def test_pressure_osm_phase(self, model):
        model.state[1].vol_mol_phase = Var(
            model.params.phase_list,
            initialize=18e-6,
            units=pyunits.m**3 / pyunits.mol,
        )

        assert_units_equivalent(model.state[1].pressure_osm_phase["Liq"], pyunits.Pa)
        assert len(model.state[1].pressure_osm_phase) == 1
        assert pytest.approx(
            value(-Constants.gas_constant * 300 * log(0.1670306) / 18e-6), rel=1e-6
        ) == value(model.state[1].pressure_osm_phase["Liq"])
