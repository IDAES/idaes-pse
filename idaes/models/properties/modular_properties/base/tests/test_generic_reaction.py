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
Tests for generic reaction package core code

Author: Andrew Lee
"""
import pytest
from sys import modules
from math import log

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    exp,
    Set,
    Var,
    value,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_equivalent

from idaes.core.base.property_meta import UnitSet
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS

from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)
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
    log_power_law_equil,
)

from idaes.core.util.constants import Constants as constants

from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
import idaes.logger as idaeslog


# -----------------------------------------------------------------------------
# Dummy methods to avoid errors
def set_metadata(b):
    pass


def define_state(b):
    b.temperature = Var(initialize=100)
    b.mole_frac_phase_comp = Var(
        b.params.phase_list, b.params.component_list, initialize=0.5
    )


# Declare a base units dict to save code later
base_units = {
    "time": pyunits.s,
    "length": pyunits.m,
    "mass": pyunits.kg,
    "amount": pyunits.mol,
    "temperature": pyunits.K,
}


# -----------------------------------------------------------------------------
class TestGenericReactionParameterBlock(object):
    @pytest.fixture
    def m(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = GenericParameterBlock(
            components={"c1": {}, "c2": {}},
            phases={
                "p1": {"equation_of_state": DummyEoS},
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        return m

    @pytest.mark.unit
    def test_rate_build(self, m):
        m.rxn_params = GenericReactionParameterBlock(
            property_package=m.params,
            rate_reactions={
                "r1": {
                    "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                    "heat_of_reaction": "foo",
                    "rate_form": "foo",
                }
            },
            base_units=base_units,
        )

        rxn_config = m.rxn_params.config.rate_reactions

        default_units = m.rxn_params.get_metadata().default_units
        assert isinstance(default_units, UnitSet)
        assert_units_equivalent(default_units.TIME, pyunits.s)
        assert_units_equivalent(default_units.LENGTH, pyunits.m)
        assert_units_equivalent(default_units.MASS, pyunits.kg)
        assert_units_equivalent(default_units.AMOUNT, pyunits.mol)
        assert_units_equivalent(default_units.TEMPERATURE, pyunits.K)
        assert_units_equivalent(default_units.CURRENT, pyunits.W)
        assert_units_equivalent(default_units.LUMINOUS_INTENSITY, pyunits.candela)

        assert isinstance(m.rxn_params.rate_reaction_idx, Set)
        assert len(m.rxn_params.rate_reaction_idx) == 1
        assert "r1" in m.rxn_params.rate_reaction_idx

        assert not hasattr(self, "equilibrium_reaction_idx")

        assert isinstance(m.rxn_params.rate_reaction_stoichiometry, dict)
        assert len(m.rxn_params.rate_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.rate_reaction_stoichiometry.items():
            if (k[1], k[2]) in rxn_config[k[0]].stoichiometry.keys():
                assert v == rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "equilibrium_reaction_stoichiometry")

        assert isinstance(m.rxn_params.reaction_idx, Set)
        assert m.rxn_params.reaction_idx == m.rxn_params.rate_reaction_idx

        assert isinstance(m.rxn_params.reaction_r1, Block)

    @pytest.mark.unit
    def test_invalid_unit(self, m):
        with pytest.raises(
            PropertyPackageError,
            match="Unrecognized units of measurement for quantity TIME " "\(foo\)",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                rate_reactions={
                    "r1": {
                        "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                        "heat_of_reaction": "foo",
                        "rate_form": "foo",
                    }
                },
                base_units={
                    "time": "foo",
                    "length": pyunits.m,
                    "mass": pyunits.kg,
                    "amount": pyunits.mol,
                    "temperature": pyunits.K,
                },
            )

    @pytest.mark.unit
    def test_rate_build_no_stoichiometry(self, m):
        with pytest.raises(
            ConfigurationError,
            match="rxn_params rate reaction r1 was not "
            "provided with a stoichiometry configuration "
            "argument.",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                base_units=base_units,
                rate_reactions={"r1": {"heat_of_reaction": "foo", "rate_form": "foo"}},
            )

    @pytest.mark.unit
    def test_rate_build_invalid_phase_stoichiometry(self, m):
        with pytest.raises(
            ConfigurationError,
            match="rxn_params stoichiometry for rate reaction "
            "r1 included unrecognised phase p7.",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                base_units=base_units,
                rate_reactions={
                    "r1": {
                        "stoichiometry": {("p7", "c1"): -1, ("p1", "c2"): 2},
                        "heat_of_reaction": "foo",
                        "rate_form": "foo",
                    }
                },
            )

    @pytest.mark.unit
    def test_rate_build_invalid_component_stoichiometry(self, m):
        with pytest.raises(
            ConfigurationError,
            match="rxn_params stoichiometry for rate reaction "
            "r1 included unrecognised component c7.",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                base_units=base_units,
                rate_reactions={
                    "r1": {
                        "stoichiometry": {("p1", "c7"): -1, ("p1", "c2"): 2},
                        "heat_of_reaction": "foo",
                        "rate_form": "foo",
                    }
                },
            )

    @pytest.mark.unit
    def test_rate_build_no_form(self, m, caplog):
        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.models.properties.modular_properties.base.generic_reaction"),
        )

        m.rxn_params = GenericReactionParameterBlock(
            property_package=m.params,
            base_units=base_units,
            rate_reactions={
                "r1": {
                    "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                    "heat_of_reaction": "foo",
                }
            },
        )

        assert (
            "rxn_params rate reaction r1 was not provided with a "
            "rate_form configuration argument. This is suitable for "
            "processes using stoichiometric reactors, but not for those "
            "using unit operations which rely on reaction rate." in caplog.text
        )

    @pytest.mark.unit
    def test_equil_build(self, m):
        m.rxn_params = GenericReactionParameterBlock(
            property_package=m.params,
            base_units=base_units,
            equilibrium_reactions={
                "e1": {
                    "stoichiometry": {("p2", "c1"): -3, ("p2", "c2"): 4},
                    "heat_of_reaction": "foo",
                    "equilibrium_form": "foo",
                }
            },
        )

        rxn_config = m.rxn_params.config.equilibrium_reactions

        assert isinstance(m.rxn_params.equilibrium_reaction_idx, Set)
        assert len(m.rxn_params.equilibrium_reaction_idx) == 1
        assert "e1" in m.rxn_params.equilibrium_reaction_idx

        assert not hasattr(self, "rate_reaction_idx")

        assert isinstance(m.rxn_params.equilibrium_reaction_stoichiometry, dict)
        assert len(m.rxn_params.equilibrium_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.equilibrium_reaction_stoichiometry.items():
            if (k[1], k[2]) in rxn_config[k[0]].stoichiometry.keys():
                assert v == rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "rate_reaction_stoichiometry")

        assert isinstance(m.rxn_params.reaction_idx, Set)
        assert m.rxn_params.reaction_idx == m.rxn_params.equilibrium_reaction_idx

        assert isinstance(m.rxn_params.reaction_e1, Block)

    @pytest.mark.unit
    def test_equil_build_no_stoichiometry(self, m):
        with pytest.raises(
            ConfigurationError,
            match="rxn_params equilibrium reaction e1 was not "
            "provided with a stoichiometry configuration "
            "argument.",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                base_units=base_units,
                equilibrium_reactions={
                    "e1": {"heat_of_reaction": "foo", "equilibrium_form": "foo"}
                },
            )

    @pytest.mark.unit
    def test_equil_build_invalid_phase_stoichiometry(self, m):
        with pytest.raises(
            ConfigurationError,
            match="rxn_params stoichiometry for equilibrium "
            "reaction e1 included unrecognised phase p7.",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                base_units=base_units,
                equilibrium_reactions={
                    "e1": {
                        "stoichiometry": {("p7", "c1"): -3, ("p2", "c2"): 4},
                        "heat_of_reaction": "foo",
                        "equilibrium_form": "foo",
                    }
                },
            )

    @pytest.mark.unit
    def test_equil_build_invalid_component_stoichiometry(self, m):
        with pytest.raises(
            ConfigurationError,
            match="rxn_params stoichiometry for equilibrium "
            "reaction e1 included unrecognised component c7.",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                base_units=base_units,
                equilibrium_reactions={
                    "e1": {
                        "stoichiometry": {("p2", "c7"): -3, ("p2", "c2"): 4},
                        "heat_of_reaction": "foo",
                        "equilibrium_form": "foo",
                    }
                },
            )

    @pytest.mark.unit
    def test_equil_build_no_form(self, m):
        with pytest.raises(
            ConfigurationError,
            match="rxn_params equilibrium reaction e1 was not "
            "provided with a equilibrium_form configuration "
            "argument.",
        ):
            m.rxn_params = GenericReactionParameterBlock(
                property_package=m.params,
                base_units=base_units,
                equilibrium_reactions={
                    "e1": {
                        "stoichiometry": {("p2", "c1"): -3, ("p2", "c2"): 4},
                        "heat_of_reaction": "foo",
                    }
                },
            )

    @pytest.mark.unit
    def test_rate_and_equil_build(self, m):
        m.rxn_params = GenericReactionParameterBlock(
            property_package=m.params,
            base_units=base_units,
            rate_reactions={
                "r1": {
                    "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                    "heat_of_reaction": "foo",
                    "rate_form": "foo",
                }
            },
            equilibrium_reactions={
                "e1": {
                    "stoichiometry": {("p2", "c1"): -3, ("p2", "c2"): 4},
                    "heat_of_reaction": "foo",
                    "equilibrium_form": "foo",
                }
            },
        )

        r_rxn_config = m.rxn_params.config.rate_reactions

        assert isinstance(m.rxn_params.rate_reaction_idx, Set)
        assert len(m.rxn_params.rate_reaction_idx) == 1
        assert "r1" in m.rxn_params.rate_reaction_idx

        assert not hasattr(self, "equilibrium_reaction_idx")

        assert isinstance(m.rxn_params.rate_reaction_stoichiometry, dict)
        assert len(m.rxn_params.rate_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.rate_reaction_stoichiometry.items():
            if (k[1], k[2]) in r_rxn_config[k[0]].stoichiometry.keys():
                assert v == r_rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "equilibrium_reaction_stoichiometry")

        e_rxn_config = m.rxn_params.config.equilibrium_reactions

        assert isinstance(m.rxn_params.equilibrium_reaction_idx, Set)
        assert len(m.rxn_params.equilibrium_reaction_idx) == 1
        assert "e1" in m.rxn_params.equilibrium_reaction_idx

        assert not hasattr(self, "rate_reaction_idx")

        assert isinstance(m.rxn_params.equilibrium_reaction_stoichiometry, dict)
        assert len(m.rxn_params.equilibrium_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.equilibrium_reaction_stoichiometry.items():
            if (k[1], k[2]) in e_rxn_config[k[0]].stoichiometry.keys():
                assert v == e_rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "rate_reaction_stoichiometry")

        assert isinstance(m.rxn_params.reaction_idx, Set)
        assert m.rxn_params.reaction_idx == (
            m.rxn_params.rate_reaction_idx | m.rxn_params.equilibrium_reaction_idx
        )
        assert len(m.rxn_params.reaction_idx) == 2

        assert isinstance(m.rxn_params.reaction_r1, Block)
        assert isinstance(m.rxn_params.reaction_e1, Block)

    @pytest.mark.unit
    def test_build_parameters(self, m):
        m.rxn_params = GenericReactionParameterBlock(
            property_package=m.params,
            base_units=base_units,
            rate_reactions={
                "r1": {
                    "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                    "heat_of_reaction": constant_dh_rxn,
                    "rate_form": "foo",
                    "parameter_data": {"dh_rxn_ref": -10000},
                }
            },
            equilibrium_reactions={
                "e1": {
                    "stoichiometry": {("p2", "c1"): -3, ("p2", "c2"): 4},
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_form": "foo",
                    "parameter_data": {"dh_rxn_ref": -20000},
                }
            },
        )

        assert isinstance(m.rxn_params.reaction_r1.dh_rxn_ref, Var)
        assert m.rxn_params.reaction_r1.dh_rxn_ref.fixed
        assert m.rxn_params.reaction_r1.dh_rxn_ref.value == -10000

        assert isinstance(m.rxn_params.reaction_e1.dh_rxn_ref, Var)
        assert m.rxn_params.reaction_e1.dh_rxn_ref.fixed
        assert m.rxn_params.reaction_e1.dh_rxn_ref.value == -20000


# -----------------------------------------------------------------------------
class TestGenericReactionBlock(object):
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = GenericParameterBlock(
            components={"c1": {}, "c2": {}},
            phases={
                "p1": {"equation_of_state": DummyEoS},
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.sblock = m.params.build_state_block([1])

        m.rxn_params = GenericReactionParameterBlock(
            property_package=m.params,
            base_units=base_units,
            rate_reactions={
                "r1": {
                    "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                    "heat_of_reaction": constant_dh_rxn,
                    "rate_constant": arrhenius,
                    "rate_form": power_law_rate,
                    "concentration_form": ConcentrationForm.moleFraction,
                    "parameter_data": {
                        "dh_rxn_ref": -10000,
                        "arrhenius_const": 1,
                        "energy_activation": 1000,
                    },
                }
            },
            equilibrium_reactions={
                "e1": {
                    "stoichiometry": {("p2", "c1"): -3, ("p2", "c2"): 4},
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_constant": van_t_hoff,
                    "equilibrium_form": power_law_equil,
                    "concentration_form": ConcentrationForm.moleFraction,
                    "parameter_data": {
                        "dh_rxn_ref": -20000,
                        "k_eq_ref": 100,
                        "T_eq_ref": 350,
                    },
                },
                "e2": {
                    "stoichiometry": {("p2", "c1"): -5, ("p2", "c2"): 6},
                    "heat_of_reaction": constant_dh_rxn,
                    "equilibrium_constant": van_t_hoff,
                    "equilibrium_form": log_power_law_equil,
                    "concentration_form": ConcentrationForm.moleFraction,
                    "parameter_data": {
                        "dh_rxn_ref": -20000,
                        "k_eq_ref": 100,
                        "T_eq_ref": 350,
                    },
                },
            },
        )

        m.rblock = m.rxn_params.build_reaction_block(
            [1], state_block=m.sblock, has_equilibrium=True
        )

        return m

    @pytest.mark.unit
    def test_component_phase_lists(self, model):
        assert model.rblock[1].component_list is model.params.component_list
        assert model.rblock[1].phase_list is model.params.phase_list
        assert model.rblock[1].phase_component_set is model.params._phase_component_set

    @pytest.mark.unit
    def test_dh_rxn(self, model):
        assert isinstance(model.rxn_params.reaction_r1.dh_rxn_ref, Var)
        assert isinstance(model.rxn_params.reaction_e1.dh_rxn_ref, Var)
        assert model.rxn_params.reaction_r1.dh_rxn_ref.value == -10000
        assert model.rxn_params.reaction_e1.dh_rxn_ref.value == -20000

        assert isinstance(model.rblock[1].dh_rxn, Expression)
        assert len(model.rblock[1].dh_rxn) == 3
        assert value(model.rblock[1].dh_rxn["r1"]) == -10000
        assert value(model.rblock[1].dh_rxn["e1"]) == -20000
        assert value(model.rblock[1].dh_rxn["e2"]) == -20000

    @pytest.mark.unit
    def test_rate_constant(self, model):
        assert isinstance(model.rxn_params.reaction_r1.arrhenius_const, Var)
        assert model.rxn_params.reaction_r1.arrhenius_const.value == 1
        assert isinstance(model.rxn_params.reaction_r1.energy_activation, Var)
        assert model.rxn_params.reaction_r1.energy_activation.value == 1000

        assert isinstance(model.rblock[1].k_rxn, Expression)
        assert len(model.rblock[1].k_rxn) == 1
        assert value(model.rblock[1].k_rxn["r1"]) == value(
            1 * exp(-1000 / (constants.gas_constant * model.sblock[1].temperature))
        )

    @pytest.mark.unit
    def test_reaction_rate(self, model):
        rblk = model.rxn_params.reaction_r1
        assert isinstance(rblk.reaction_order, Var)
        assert len(rblk.reaction_order) == 4
        assert rblk.reaction_order["p1", "c1"].value == 1
        assert rblk.reaction_order["p1", "c2"].value == 0
        assert rblk.reaction_order["p2", "c1"].value == 0
        assert rblk.reaction_order["p2", "c2"].value == 0

        assert isinstance(model.rblock[1].reaction_rate, Expression)
        assert len(model.rblock[1].reaction_rate) == 1
        assert value(model.rblock[1].reaction_rate["r1"]) == value(
            model.rblock[1].k_rxn["r1"]
            * model.sblock[1].mole_frac_phase_comp["p1", "c1"] ** 1
        )

    @pytest.mark.unit
    def test_reaction_rate_None(self, model):
        model.rxn_params.config.rate_reactions.r1.rate_form = None

        with pytest.raises(
            ConfigurationError,
            match="rblock\[1\] Generic Reaction r1 was not "
            "provided with a rate_form configuration "
            "argument.",
        ):
            model.rblock[1].reaction_rate

    @pytest.mark.unit
    def test_equilibrium_constant(self, model):
        assert isinstance(model.rxn_params.reaction_e1.k_eq_ref, Var)
        assert model.rxn_params.reaction_e1.k_eq_ref.value == 100
        assert isinstance(model.rxn_params.reaction_e1.T_eq_ref, Var)
        assert model.rxn_params.reaction_e1.T_eq_ref.value == 350

        assert isinstance(model.rblock[1].k_eq, Expression)
        assert len(model.rblock[1].k_eq) == 2
        assert value(model.rblock[1].k_eq["e1"]) == value(exp(1))
        assert value(model.rblock[1].k_eq["e2"]) == value(exp(1))

    @pytest.mark.unit
    def test_equilibrium_form(self, model):
        rblk = model.rxn_params.reaction_e1
        assert isinstance(rblk.reaction_order, Var)
        assert len(rblk.reaction_order) == 4
        assert rblk.reaction_order["p1", "c1"].value == 0
        assert rblk.reaction_order["p1", "c2"].value == 0
        assert rblk.reaction_order["p2", "c1"].value == -3
        assert rblk.reaction_order["p2", "c2"].value == 4

        rblk = model.rxn_params.reaction_e2
        assert isinstance(rblk.reaction_order, Var)
        assert len(rblk.reaction_order) == 4
        assert rblk.reaction_order["p1", "c1"].value == 0
        assert rblk.reaction_order["p1", "c2"].value == 0
        assert rblk.reaction_order["p2", "c1"].value == -5
        assert rblk.reaction_order["p2", "c2"].value == 6

        assert isinstance(model.rblock[1].equilibrium_constraint, Constraint)
        assert len(model.rblock[1].equilibrium_constraint) == 2
        assert value(model.rblock[1].equilibrium_constraint["e1"].body) == value(
            model.rblock[1].k_eq["e1"]
            - model.sblock[1].mole_frac_phase_comp["p2", "c1"] ** -3
            * model.sblock[1].mole_frac_phase_comp["p2", "c2"] ** 4
        )
        assert value(model.rblock[1].equilibrium_constraint["e2"].body) == value(
            model.rblock[1].log_k_eq["e2"]
            - model.sblock[1].log_mole_frac_phase_comp["p2", "c1"] * -5
            + model.sblock[1].log_mole_frac_phase_comp["p2", "c2"] * 6
        )

    @pytest.mark.unit
    def test_basic_scaling(self, model):
        model.rblock[1].calculate_scaling_factors()

        assert len(model.rblock[1].scaling_factor) == 7

        assert model.rblock[1].scaling_factor[model.rblock[1].dh_rxn["e1"]] == 5e-5
        assert model.rblock[1].scaling_factor[model.rblock[1].dh_rxn["r1"]] == 1e-4
        assert model.rblock[1].scaling_factor[model.rblock[1].k_eq["e1"]] == 1e-2
        assert model.rblock[1].scaling_factor[model.rblock[1].k_eq["e2"]] == 1e-2
        assert model.rblock[1].scaling_factor[model.rblock[1].log_k_eq["e1"]] == log(
            1e-2
        )
        assert model.rblock[1].scaling_factor[model.rblock[1].log_k_eq["e2"]] == log(
            1e-2
        )

        # Check that scaling factor was added to equilibrium constraint
        assert str(model.rblock[1].equilibrium_constraint["e1"].body) == str(
            (
                model.rblock[1].k_eq["e1"]
                - model.sblock[1].mole_frac_phase_comp["p2", "c1"]
                ** model.rxn_params.reaction_e1.reaction_order["p2", "c1"]
                * model.sblock[1].mole_frac_phase_comp["p2", "c2"]
                ** model.rxn_params.reaction_e1.reaction_order["p2", "c2"]
            )
            * 0.01
        )
        print(model.rblock[1].equilibrium_constraint["e2"].body)
        assert str(model.rblock[1].equilibrium_constraint["e2"].body) == str(
            (
                model.rblock[1].log_k_eq["e2"]
                - (
                    model.rxn_params.reaction_e2.reaction_order["p2", "c1"]
                    * model.sblock[1].log_mole_frac_phase_comp["p2", "c1"]
                    + model.rxn_params.reaction_e2.reaction_order["p2", "c2"]
                    * model.sblock[1].log_mole_frac_phase_comp["p2", "c2"]
                )
            )
            * log(0.01)
        )
