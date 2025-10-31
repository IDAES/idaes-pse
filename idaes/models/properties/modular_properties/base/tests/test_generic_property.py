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
Tests for generic property package core code

Author: Andrew Lee, Douglas Allan
"""
import functools
import pytest
import re
from sys import modules
from types import MethodType

from pyomo.environ import value, Block, ConcreteModel, Param, Set, Var, units as pyunits
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
    GenericParameterData,
    GenericStateBlock,
    _initialize_critical_props,
    ModularPropertiesInitializer,
    ModularPropertiesScaler,
)
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS

from idaes.core import (
    declare_process_block_class,
    Component,
    FlowsheetBlock,
    Phase,
    LiquidPhase,
    VaporPhase,
    MaterialBalanceType,
    MaterialFlowBasis,
    Solvent,
    PhaseType as PT,
)
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
from idaes.models.properties.modular_properties.phase_equil.henry import HenryType
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.core.base.property_meta import UnitSet

from idaes.models.properties.modular_properties.phase_equil.henry import HenryType
from idaes.models.properties.modular_properties.examples.BT_ideal import (
    configuration as BTconfig,
)
from idaes.models.unit_models.flash import Flash
from idaes.core.scaling import get_scaling_factor

import idaes.logger as idaeslog


# -----------------------------------------------------------------------------
# Dummy method to avoid errors when setting metadata dict
def set_metadata(b):
    pass


def calculate_scaling_factors(b):
    b.scaling_check = True


# Dummy build_parameter methods for tests
def build_parameters(cobj, p, *args, **kwargs):
    cobj.add_component("test_param_" + p, Var(initialize=42))


# Declare a base units dict to save code later
base_units = {
    "time": pyunits.s,
    "length": pyunits.m,
    "mass": pyunits.kg,
    "amount": pyunits.mol,
    "temperature": pyunits.K,
}


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True


class TestGenericParameterBlock(object):
    @pytest.mark.unit
    def test_build(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={"a": {}, "b": {}, "c": {}},
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "component_list": ["a", "b"],
                    "equation_of_state": DummyEoS,
                },
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        assert m.params.configured

        default_units = m.params.get_metadata().default_units
        assert isinstance(default_units, UnitSet)
        assert_units_equivalent(default_units.TIME, pyunits.s)
        assert_units_equivalent(default_units.LENGTH, pyunits.m)
        assert_units_equivalent(default_units.MASS, pyunits.kg)
        assert_units_equivalent(default_units.AMOUNT, pyunits.mol)
        assert_units_equivalent(default_units.TEMPERATURE, pyunits.K)
        assert_units_equivalent(default_units.CURRENT, pyunits.ampere)
        assert_units_equivalent(default_units.LUMINOUS_INTENSITY, pyunits.candela)

        assert isinstance(m.params.component_list, Set)
        assert len(m.params.component_list) == 3
        for j in m.params.component_list:
            assert j in ["a", "b", "c"]
            assert isinstance(m.params.get_component(j), Component)

        assert isinstance(m.params.phase_list, Set)
        assert len(m.params.phase_list) == 2
        for p in m.params.phase_list:
            assert p in ["p1", "p2"]
        assert isinstance(m.params.get_phase("p1"), LiquidPhase)
        assert isinstance(m.params.get_phase("p2"), Phase)
        assert m.params.p1.config.equation_of_state == DummyEoS
        assert m.params.p2.config.equation_of_state == DummyEoS

        assert isinstance(m.params._phase_component_set, Set)
        assert len(m.params._phase_component_set) == 5
        for i in m.params._phase_component_set:
            assert i in [
                ("p1", "a"),
                ("p1", "b"),
                ("p2", "a"),
                ("p2", "b"),
                ("p2", "c"),
            ]

        assert isinstance(m.params.pressure_ref, Param)
        assert m.params.pressure_ref.value == 1e5
        assert isinstance(m.params.temperature_ref, Param)
        assert m.params.temperature_ref.value == 300

        assert m.params.state_block_class is GenericStateBlock

        assert len(m.params.config.inherent_reactions) == 0
        assert m.params.config.reaction_basis == MaterialFlowBasis.molar

        assert not m.params.has_inherent_reactions

    @pytest.mark.unit
    def test_invalid_unit(self):
        m = ConcreteModel()

        with pytest.raises(
            PropertyPackageError,
            match=re.escape(
                "Unrecognized units of measurement for quantity TIME (foo)"
            ),
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units={
                    "time": "foo",
                    "length": pyunits.m,
                    "mass": pyunits.kg,
                    "amount": pyunits.mol,
                    "temperature": pyunits.K,
                },
            )

    @pytest.mark.unit
    def test_no_components(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params was not provided with a components argument.",
        ):
            m.params = DummyParameterBlock(
                phases={
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"},
                },
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_no_phases(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params was not provided with a phases argument. "
            "Did you forget to unpack the configurations dictionary?",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}}, base_units=base_units
            )

    @pytest.mark.unit
    def test_packed_dict(self):
        m = ConcreteModel()

        dummy_dict = {
            "phases": {
                "p1": {"equation_of_state": "foo"},
                "p2": {"equation_of_state": "bar"},
            },
        }

        with pytest.raises(
            ConfigurationError,
            match=re.escape(
                "params[phases] was not provided with a phases argument. "
                "Did you forget to unpack the configurations dictionary?"
            ),
        ):
            m.params = DummyParameterBlock(dummy_dict)

    @pytest.mark.unit
    def test_invalid_component_in_phase_component_list(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params phase-component list for phase p1 "
            "contained component d which is not in the master "
            "component list",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "d"],
                        "equation_of_state": "foo",
                    },
                    "p2": {"equation_of_state": "bar"},
                },
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_invalid_orphaned_component_in_phase_component_list(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Component c does not appear to be "
            "valid in any phase. Please check the component "
            "lists defined for each phase, and be sure you do "
            "not have generic Components in single-phase "
            "aqueous systems.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": "foo",
                    }
                },
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_invalid_component_in_phase_component_list_2(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params phase-component list for phase p1 "
            "contained component a, however this component is "
            "not valid for the given PhaseType",
        ):
            m.params = DummyParameterBlock(
                components={
                    "a": {"valid_phase_types": PT.solidPhase},
                    "b": {},
                    "c": {},
                },
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": "foo",
                    },
                    "p2": {"equation_of_state": "bar"},
                },
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_phase_component_set_from_valid_phases(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(
            components={
                "a": {"valid_phase_types": PT.liquidPhase},
                "b": {"valid_phase_types": PT.vaporPhase},
                "c": {"valid_phase_types": [PT.liquidPhase, PT.vaporPhase]},
            },
            phases={
                "p1": {"type": LiquidPhase, "equation_of_state": DummyEoS},
                "p2": {"type": VaporPhase, "equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        assert len(m.params._phase_component_set) == 4
        for i in m.params._phase_component_set:
            assert i in [("p1", "a"), ("p1", "c"), ("p2", "b"), ("p2", "c")]

    @pytest.mark.unit
    def test_no_state_definition(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Generic Property Package was not "
            "provided with a state_definition configuration "
            "argument. Please fix your property parameter "
            "definition to include this.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"},
                },
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_no_pressure_ref(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Generic Property Package was not "
            "provided with a pressure_ref configuration "
            "argument. Please fix your property parameter "
            "definition to include this.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={"p1": {"equation_of_state": "foo"}, "p2": {}},
                state_definition="baz",
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_convert_pressure_ref(self):
        m = ConcreteModel()

        # This will fail, but should set the reference pressure
        with pytest.raises(ConfigurationError):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={"p1": {"equation_of_state": "foo"}, "p2": {}},
                state_definition="baz",
                pressure_ref=(1, pyunits.bar),
                base_units=base_units,
            )
        assert m.params.pressure_ref.value == 1e5

    @pytest.mark.unit
    def test_log_no_units_pressure_ref(self, caplog):
        m = ConcreteModel()

        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.models.properties.modular_properties.base.generic_property"),
        )

        # This will fail, but should set the reference pressure
        with pytest.raises(ConfigurationError):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={"p1": {"equation_of_state": "foo"}, "p2": {}},
                state_definition="baz",
                pressure_ref=100000.0,
                base_units=base_units,
            )
        assert m.params.pressure_ref.value == 1e5
        assert (
            "params no units provided for parameter pressure_ref - "
            "assuming default units" in caplog.text
        )

    @pytest.mark.unit
    def test_no_temperature_ref(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Generic Property Package was not "
            "provided with a temperature_ref configuration "
            "argument. Please fix your property parameter "
            "definition to include this.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={"p1": {"equation_of_state": "foo"}, "p2": {}},
                state_definition="baz",
                pressure_ref=100000.0,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_log_no_units_temperature_ref(self, caplog):
        m = ConcreteModel()

        caplog.set_level(
            idaeslog.DEBUG,
            logger=("idaes.models.properties.modular_properties.base.generic_property"),
        )

        # This will fail, but should set the reference pressure
        with pytest.raises(ConfigurationError):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={"p1": {"equation_of_state": "foo"}, "p2": {}},
                state_definition="baz",
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )
        assert m.params.temperature_ref.value == 300
        assert (
            "params no units provided for parameter pressure_ref - "
            "assuming default units" in caplog.text
        )

    @pytest.mark.unit
    def test_convert_temperature_ref(self):
        m = ConcreteModel()

        # This will fail, but should set temerpature_ref
        with pytest.raises(ConfigurationError):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={"p1": {"equation_of_state": "foo"}, "p2": {}},
                state_definition="baz",
                pressure_ref=100000.0,
                temperature_ref=(540, pyunits.degR),
                base_units=base_units,
            )
        assert m.params.temperature_ref.value == 300

    @pytest.mark.unit
    def test_no_eos(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params phase p2 was not provided with an "
            "equation_of_state configuration argument. Please "
            "fix your property parameter definition to "
            "include this.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={"p1": {"equation_of_state": "foo"}, "p2": {}},
                state_definition="baz",
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_phases_in_equilibrium(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(
            components={
                "a": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                "b": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                "c": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
            },
            phases={
                "p1": {"equation_of_state": DummyEoS},
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            phases_in_equilibrium=[("p1", "p2")],
            phase_equilibrium_state={("p1", "p2"): "whoop"},
            base_units=base_units,
        )

        assert isinstance(m.params.phase_equilibrium_idx, Set)
        assert len(m.params.phase_equilibrium_idx) == 3
        for i in m.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert isinstance(m.params.phase_equilibrium_list, dict)
        assert m.params.phase_equilibrium_list == {
            "PE1": ["a", ("p1", "p2")],
            "PE2": ["b", ("p1", "p2")],
            "PE3": ["c", ("p1", "p2")],
        }

    @pytest.mark.unit
    def test_phases_in_equilibrium_no_form(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Generic Property Package component a "
            "is in equilibrium but phase_equilibrium_form"
            " was not specified.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"},
                },
                state_definition="baz",
                pressure_ref=100000.0,
                temperature_ref=300,
                phases_in_equilibrium=[("p1", "p2")],
                phase_equilibrium_state={("p1", "p2"): "whoop"},
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_phases_in_equilibrium_missing_pair_form(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Generic Property Package component b "
            "is in equilibrium but phase_equilibrium_form"
            " was not specified for all appropriate "
            "phase pairs.",
        ):
            # Also reverse order of phases for component a - this should pass
            # and component b should be flagged as missing
            m.params = DummyParameterBlock(
                components={
                    "a": {"phase_equilibrium_form": {("p2", "p1"): "foo"}},
                    "b": {"phase_equilibrium_form": {(1, 2): "foo"}},
                    "c": {},
                },
                phases={
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"},
                },
                state_definition="baz",
                pressure_ref=100000.0,
                temperature_ref=300,
                phases_in_equilibrium=[("p1", "p2")],
                phase_equilibrium_state={("p1", "p2"): "whoop"},
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_phases_in_equilibrium_no_formulation(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Generic Property Package provided "
            "with a phases_in_equilibrium argument but no "
            "method was specified for "
            "phase_equilibrium_state.",
        ):
            m.params = DummyParameterBlock(
                components={
                    "a": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "b": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "c": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                },
                phases={
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"},
                },
                state_definition="baz",
                pressure_ref=100000.0,
                temperature_ref=300,
                phases_in_equilibrium=[("p1", "p2")],
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_phases_in_equilibrium_missing_pair_formulation(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params Generic Property Package provided "
            "with a phases_in_equilibrium argument but "
            "phase_equilibrium_state was not specified "
            "for all phase pairs.",
        ):
            # Also reverse order of phases for component a - this should pass
            # and component b should be flagged as missing
            m.params = DummyParameterBlock(
                components={
                    "a": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "b": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "c": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                },
                phases={
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"},
                },
                state_definition="baz",
                pressure_ref=100000.0,
                temperature_ref=300,
                phases_in_equilibrium=[("p1", "p2")],
                phase_equilibrium_state={(1, 2): "whoop"},
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_parameter_construction_no_value(self):
        m = ConcreteModel()

        class test_class:
            # Mook up property method class for testing
            def build_parameters(c):
                c.test_var = Var()

        with pytest.raises(
            ConfigurationError,
            match="params parameter test_var was not assigned "
            "a value. Please check your configuration "
            "arguments.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {"dens_mol_liq_comp": test_class}, "b": {}, "c": {}},
                phases={
                    "p1": {"equation_of_state": DummyEoS},
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_parameter_construction_no_data(self):
        m = ConcreteModel()

        class test_class:
            # Mook up property method class for testing
            def build_parameters(c):
                c.config.parameter_data["test"]

        with pytest.raises(
            ConfigurationError,
            match="params values were not defined for "
            "parameter dens_mol_liq_comp in component a. "
            "Please check the parameter_data argument to "
            "ensure values are provided.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {"dens_mol_liq_comp": test_class}, "b": {}, "c": {}},
                phases={
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_no_elements(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={"a": {}, "b": {}, "c": {}},
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "component_list": ["a", "b"],
                    "equation_of_state": DummyEoS,
                },
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        assert not hasattr(m.params, "element_list")
        assert not hasattr(m.params, "element_comp")

    @pytest.mark.unit
    def test_partial_elements(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params not all Components declared an "
            "elemental_composition argument. Either all "
            "Components must declare this, or none.",
        ):
            m.params = DummyParameterBlock(
                components={
                    "a": {"elemental_composition": {"e1": 1}},
                    "b": {},
                    "c": {},
                },
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_elements_not_float(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match=re.escape(
                "params values in elemental_composition must "
                "be integers (not floats): e1: 2.0."
            ),
        ):
            m.params = DummyParameterBlock(
                components={
                    "a": {"elemental_composition": {"e1": 2.0}},
                    "b": {},
                    "c": {},
                },
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_elements(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(
            components={
                "a": {"elemental_composition": {"e1": 1, "e2": 2}},
                "b": {"elemental_composition": {"e3": 3, "e4": 4}},
                "c": {"elemental_composition": {"e1": 5, "e3": 6}},
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "component_list": ["a", "b"],
                    "equation_of_state": DummyEoS,
                },
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        assert isinstance(m.params.element_list, Set)
        assert len(m.params.element_list) == 4
        assert m.params.element_comp == {
            "a": {"e1": 1, "e2": 2, "e3": 0, "e4": 0},
            "b": {"e1": 0, "e2": 0, "e3": 3, "e4": 4},
            "c": {"e1": 5, "e2": 0, "e3": 6, "e4": 0},
        }

    @pytest.mark.unit
    def test_henry(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(
            components={
                "a": {
                    "henry_component": {
                        "p1": {"method": modules[__name__], "type": HenryType.Kpx}
                    }
                },
                "b": {},
                "c": {},
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "component_list": ["a", "b"],
                    "equation_of_state": DummyEoS,
                },
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        assert isinstance(m.params.a.test_param_p1, Var)
        assert m.params.a.test_param_p1.value == 42

    @pytest.mark.unit
    def test_henry_invalid_phase_name(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params component a was marked as a Henry's "
            "Law component in phase p3, but this is not a "
            "valid phase name.",
        ):
            m.params = DummyParameterBlock(
                components={
                    "a": {"henry_component": {"p3": modules[__name__]}},
                    "b": {},
                    "c": {},
                },
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_henry_invalid_phase_type(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params component a was marked as a Henry's "
            "Law component in phase p2, but this is not a "
            "Liquid phase.",
        ):
            m.params = DummyParameterBlock(
                components={
                    "a": {"henry_component": {"p2": modules[__name__]}},
                    "b": {},
                    "c": {},
                },
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
            )

    @pytest.mark.unit
    def test_inherent_reactions(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={"a": {}, "b": {}, "c": {}},
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "component_list": ["a", "b"],
                    "equation_of_state": DummyEoS,
                },
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
            inherent_reactions={
                "e1": {
                    "stoichiometry": {("p1", "a"): -3, ("p1", "b"): 4},
                    "heat_of_reaction": "foo",
                    "equilibrium_form": "foo",
                }
            },
        )

        assert m.params.has_inherent_reactions

        assert isinstance(m.params.inherent_reaction_idx, Set)
        assert len(m.params.inherent_reaction_idx) == 1
        assert m.params.inherent_reaction_idx == ["e1"]

        assert hasattr(m.params, "inherent_reaction_stoichiometry")
        assert len(m.params.inherent_reaction_stoichiometry) == 5
        assert m.params.inherent_reaction_stoichiometry == {
            ("e1", "p1", "a"): -3,
            ("e1", "p1", "b"): 4,
            ("e1", "p2", "a"): 0,
            ("e1", "p2", "b"): 0,
            ("e1", "p2", "c"): 0,
        }

        assert isinstance(m.params.reaction_e1, Block)
        assert isinstance(m.params.reaction_e1.reaction_order, Var)
        assert len(m.params.reaction_e1.reaction_order) == 5
        for i in m.params.reaction_e1.reaction_order:
            order = {
                ("p1", "a"): -3,
                ("p1", "b"): 4,
                ("p2", "a"): 0,
                ("p2", "b"): 0,
                ("p2", "c"): 0,
            }
            assert m.params.reaction_e1.reaction_order[i].value == order[i]

    @pytest.mark.unit
    def test_inherent_reactions_no_stoichiometry(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params inherent reaction e1 was not "
            "provided with a stoichiometry configuration "
            "argument.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
                inherent_reactions={
                    "e1": {"heat_of_reaction": "foo", "equilibrium_form": "foo"}
                },
            )

    @pytest.mark.unit
    def test_inherent_reactions_invalid_phase_stoichiometry(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params stoichiometry for inherent "
            "reaction e1 included unrecognised phase p7.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
                inherent_reactions={
                    "e1": {
                        "stoichiometry": {("p7", "a"): -3, ("p1", "b"): 4},
                        "heat_of_reaction": "foo",
                        "equilibrium_form": "foo",
                    }
                },
            )

    @pytest.mark.unit
    def test_inherent_reactions_invalid_component_stoichiometry(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params stoichiometry for inherent "
            "reaction e1 included unrecognised component c7.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
                inherent_reactions={
                    "e1": {
                        "stoichiometry": {("p1", "c7"): -3, ("p1", "b"): 4},
                        "heat_of_reaction": "foo",
                        "equilibrium_form": "foo",
                    }
                },
            )

    @pytest.mark.unit
    def test_inherent_reactions_no_form(self):
        m = ConcreteModel()

        with pytest.raises(
            ConfigurationError,
            match="params inherent reaction e1 was not "
            "provided with a equilibrium_form configuration "
            "argument.",
        ):
            m.params = DummyParameterBlock(
                components={"a": {}, "b": {}, "c": {}},
                phases={
                    "p1": {
                        "type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": DummyEoS,
                    },
                    "p2": {"equation_of_state": DummyEoS},
                },
                state_definition=modules[__name__],
                pressure_ref=100000.0,
                temperature_ref=300,
                base_units=base_units,
                inherent_reactions={
                    "e1": {"stoichiometry": {("p1", "a"): -3, ("p1", "b"): 4}}
                },
            )

    @pytest.mark.unit
    def test_default_scaling(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={"a": {}, "b": {}, "c": {}},
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "component_list": ["a", "b"],
                    "equation_of_state": DummyEoS,
                },
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        dsf = m.params.default_scaling_factor

        assert dsf[("temperature", None)] == 1e-2
        assert dsf[("pressure", None)] == 1e-5
        assert dsf[("dens_mol_phase", None)] == 1e-2
        assert dsf[("enth_mol", None)] == 1e-4
        assert dsf[("entr_mol", None)] == 1e-2
        assert dsf[("fug_phase_comp", None)] == 1e-4
        assert dsf[("fug_coeff_phase_comp", None)] == 1
        assert dsf[("gibbs_mol", None)] == 1e-4
        assert dsf[("mole_frac_comp", None)] == 1e3
        assert dsf[("mole_frac_phase_comp", None)] == 1e3
        assert dsf[("mw", None)] == 1e3
        assert dsf[("mw_comp", None)] == 1e3
        assert dsf[("mw_phase", None)] == 1e3

    @pytest.mark.unit
    def test_default_scaling_convert(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={"a": {}, "b": {}, "c": {}},
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "component_list": ["a", "b"],
                    "equation_of_state": DummyEoS,
                },
                "p2": {"equation_of_state": DummyEoS},
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.minutes,
                "length": pyunits.foot,
                "mass": pyunits.pound,
                "amount": pyunits.mol,
                "temperature": pyunits.degR,
            },
        )

        dsf = m.params.default_scaling_factor

        assert dsf[("temperature", None)] == 1e-2
        assert dsf[("pressure", None)] == 1e-8
        assert dsf[("dens_mol_phase", None)] == 1
        assert dsf[("enth_mol", None)] == 1e-9
        assert dsf[("entr_mol", None)] == 1e-7
        assert dsf[("fug_phase_comp", None)] == 1e-7
        assert dsf[("fug_coeff_phase_comp", None)] == 1
        assert dsf[("gibbs_mol", None)] == 1e-9
        assert dsf[("mole_frac_comp", None)] == 1e3
        assert dsf[("mole_frac_phase_comp", None)] == 1e3
        assert dsf[("mw", None)] == 1e3
        assert dsf[("mw_comp", None)] == 1e3
        assert dsf[("mw_phase", None)] == 1e3


# -----------------------------------------------------------------------------
# Dummy methods for testing build calls to sub-modules
def define_state(b):
    b.state_defined = True
    b.temperature = Var(initialize=100)


def get_material_flow_basis(self, *args, **kwargs):
    return MaterialFlowBasis.molar


def dummy_method(*args, **kwargs):
    return 42


class TestGenericStateBlock(object):
    @pytest.fixture()
    def frame(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {
                    "type": Solvent,
                    "diffus_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "visc_d_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "therm_cond_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "pressure_sat_comp": dummy_method,
                },
                "b": {
                    "diffus_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "visc_d_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "therm_cond_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "pressure_sat_comp": dummy_method,
                },
                "c": {
                    "diffus_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "visc_d_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "therm_cond_phase_comp": {"p1": dummy_method, "p2": dummy_method},
                    "pressure_sat_comp": dummy_method,
                },
            },
            phases={
                "p1": {
                    "equation_of_state": DummyEoS,
                    "therm_cond_phase": dummy_method,
                    "surf_tens_phase": dummy_method,
                    "visc_d_phase": dummy_method,
                },
                "p2": {
                    "equation_of_state": DummyEoS,
                    "therm_cond_phase": dummy_method,
                    "surf_tens_phase": dummy_method,
                    "visc_d_phase": dummy_method,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)
        m.props[1].get_material_flow_basis = MethodType(
            get_material_flow_basis, m.props[1].get_material_flow_basis
        )

        # Add necessary variables to state block
        m.props[1].flow_mol = Var(bounds=(0, 3000))
        m.props[1].pressure = Var(bounds=(1000, 3000))
        m.props[1].temperature = Var(bounds=(100, 200))
        m.props[1].mole_frac_phase_comp = Var(
            m.params.phase_list, m.params.component_list, initialize=0.3
        )
        m.props[1].phase_frac = Var(m.params.phase_list)

        return m

    @pytest.mark.unit
    def test_build(self, frame):
        assert isinstance(frame.props, Block)
        assert len(frame.props) == 1
        assert frame.props.default_initializer is ModularPropertiesInitializer

        # Check for expected behaviour for dummy methods
        assert frame.props[1].state_defined
        assert isinstance(frame.props[1].dummy_var, Var)
        assert frame.props[1].eos_common == 2

    @pytest.mark.unit
    def test_basic_scaling_legacy(self, frame):
        frame.props[1].calculate_scaling_factors()

        assert frame.props[1].scaling_check

        assert len(frame.props[1].scaling_factor) == 8
        assert frame.props[1].scaling_factor[frame.props[1].temperature] == 1e-2
        assert frame.props[1].scaling_factor[frame.props[1].pressure] == 1e-5
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["p1", "a"]
            ]
            == 1000
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["p1", "b"]
            ]
            == 1000
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["p1", "c"]
            ]
            == 1000
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["p2", "a"]
            ]
            == 1000
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["p2", "b"]
            ]
            == 1000
        )
        assert (
            frame.props[1].scaling_factor[
                frame.props[1].mole_frac_phase_comp["p2", "c"]
            ]
            == 1000
        )

    @pytest.mark.unit
    def test_basic_scaler_object(self, frame, caplog):
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 8
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10

        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_no_default_set(self, frame):
        scaler_obj = frame.props[1].default_scaler()
        with pytest.raises(
            ValueError,
            match=re.escape(
                "This scaler requires the user to provide a default "
                "scaling factor for props[1].flow_mol_phase[p1], but "
                "no default scaling factor was set."
            ),
        ):
            scaler_obj.scale_model(frame.props[1])

    @pytest.mark.unit
    def test_basic_scaler_object_enth_mol_phase_default_scaling(self, frame, caplog):
        frame.params.a.mw = 5e-3 * pyunits.kg / pyunits.mol
        frame.params.b.mw = 7e-3 * pyunits.kg / pyunits.mol
        frame.params.c.mw = 11e-3 * pyunits.kg / pyunits.mol
        frame.props[1].enth_mol_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.default_scaling_factors["enth_mol_phase"] = 1 / 19
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        for vardata in frame.props[1].enth_mol_phase.values():
            assert get_scaling_factor(vardata) == 1 / 19
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_enth_mol_phase_partial_scaling(self, frame, caplog):
        frame.params.a.mw = 5e-3 * pyunits.kg / pyunits.mol
        frame.params.b.mw = 7e-3 * pyunits.kg / pyunits.mol
        frame.params.c.mw = 11e-3 * pyunits.kg / pyunits.mol
        frame.props[1].enth_mol_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.default_scaling_factors["enth_mol_phase"] = 1 / 19
        scaler_obj.set_component_scaling_factor(
            frame.props[1].enth_mol_phase["p2"], 1 / 23
        )
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        get_scaling_factor(frame.props[1].enth_mol_phase["p1"]) == 1 / 19
        get_scaling_factor(frame.props[1].enth_mol_phase["p2"]) == 1 / 23
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_enth_mol_phase_manual_scaling(self, frame, caplog):
        frame.params.a.mw = 5e-3 * pyunits.kg / pyunits.mol
        frame.params.b.mw = 7e-3 * pyunits.kg / pyunits.mol
        frame.params.c.mw = 11e-3 * pyunits.kg / pyunits.mol
        frame.props[1].enth_mol_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.default_scaling_factors["enth_mol"] = 1 / 31  # Should be ignored
        scaler_obj.set_component_scaling_factor(
            frame.props[1].enth_mol_phase["p1"], 1 / 19
        )
        scaler_obj.set_component_scaling_factor(
            frame.props[1].enth_mol_phase["p2"], 1 / 23
        )
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        get_scaling_factor(frame.props[1].enth_mol_phase["p1"]) == 1 / 19
        get_scaling_factor(frame.props[1].enth_mol_phase["p2"]) == 1 / 23
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_enth_mol_phase_no_default_scaling_warning(
        self, frame, caplog
    ):
        frame.props[1].enth_mol_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.default_scaling_factors["enth_mol"] = 1 / 19  # Not enth_mol_phase
        with caplog.at_level(idaeslog.WARNING):
            scaler_obj.scale_model(frame.props[1])
        assert "Default scaling factor for molar enthalpy not set." in caplog.text

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        for vardata in frame.props[1].enth_mol_phase.values():
            assert get_scaling_factor(vardata) == 1
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_enth_mol_phase_partial_scaling_warning(
        self, frame, caplog
    ):
        frame.props[1].enth_mol_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.set_component_scaling_factor(
            frame.props[1].enth_mol_phase["p2"], 1 / 23
        )
        with caplog.at_level(idaeslog.WARNING):
            scaler_obj.scale_model(frame.props[1])
        assert "Default scaling factor for molar enthalpy not set." in caplog.text

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        assert get_scaling_factor(frame.props[1].enth_mol_phase["p1"]) == 1
        assert get_scaling_factor(frame.props[1].enth_mol_phase["p2"]) == 1 / 23
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_enth_mol_phase_mw_estimate(self, frame, caplog):
        frame.params.a.mw = 5e-3 * pyunits.kg / pyunits.mol
        frame.params.b.mw = 7e-3 * pyunits.kg / pyunits.mol
        frame.params.c.mw = 11e-3 * pyunits.kg / pyunits.mol
        frame.props[1].enth_mol_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        for vardata in frame.props[1].enth_mol_phase.values():
            assert get_scaling_factor(vardata) == pytest.approx(1 / (15 + 1 / 3) / 300)
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_enth_mol_phase_mw_estimate_partial_scaling(
        self, frame, caplog
    ):
        frame.params.a.mw = 5e-3 * pyunits.kg / pyunits.mol
        frame.params.b.mw = 7e-3 * pyunits.kg / pyunits.mol
        frame.params.c.mw = 11e-3 * pyunits.kg / pyunits.mol
        frame.props[1].enth_mol_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.set_component_scaling_factor(
            frame.props[1].enth_mol_phase["p2"], 1 / 23
        )
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        assert get_scaling_factor(frame.props[1].enth_mol_phase["p1"]) == pytest.approx(
            1 / (15 + 1 / 3) / 300
        )
        assert get_scaling_factor(frame.props[1].enth_mol_phase["p2"]) == 1 / 23
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_mw_phase(self, frame, caplog):
        frame.params.a.mw = 5e-3 * pyunits.kg / pyunits.mol
        frame.params.b.mw = 7e-3 * pyunits.kg / pyunits.mol
        frame.params.c.mw = 11e-3 * pyunits.kg / pyunits.mol
        frame.props[1].mw_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        assert get_scaling_factor(frame.props[1].mw_phase["p1"]) == pytest.approx(
            1e3 / (7 + 2 / 3)
        )
        assert get_scaling_factor(frame.props[1].mw_phase["p2"]) == pytest.approx(
            1e3 / (7 + 2 / 3)
        )
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1

    @pytest.mark.unit
    def test_basic_scaler_object_therm_cond_phase_default_scaling(self, frame, caplog):
        frame.props[1].therm_cond_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.default_scaling_factors["therm_cond_phase[p1]"] = 1 / 71
        scaler_obj.default_scaling_factors["therm_cond_phase[p2]"] = 1 / 479
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        for vardata in frame.props[1].enth_mol_phase.values():
            assert get_scaling_factor(vardata) == 1 / 19
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1
        assert get_scaling_factor(frame.props[1].therm_cond_phase["p1"]) == 1 / 71
        assert get_scaling_factor(frame.props[1].therm_cond_phase["p2"]) == 1 / 479

    @pytest.mark.unit
    def test_basic_scaler_object_therm_cond_phase_default_scaling(self, frame, caplog):
        frame.props[1].therm_cond_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        scaler_obj.default_scaling_factors["therm_cond_phase[p1]"] = 1 / 71
        scaler_obj.default_scaling_factors["therm_cond_phase[p2]"] = 1 / 479
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 10
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1
        assert get_scaling_factor(frame.props[1].therm_cond_phase["p1"]) == 1 / 71
        assert get_scaling_factor(frame.props[1].therm_cond_phase["p2"]) == 1 / 479

    @pytest.mark.unit
    def test_basic_scaler_object_therm_cond_phase_estimate_scaling(self, frame, caplog):
        frame.props[1].therm_cond_phase = Var(["p1", "p2"], initialize=0)
        scaler_class = frame.props[1].default_scaler
        assert scaler_class is ModularPropertiesScaler
        scaler_obj = scaler_class()
        scaler_obj.default_scaling_factors["flow_mol_phase"] = 1
        with caplog.at_level(idaeslog.INFO_HIGH):
            scaler_obj.scale_model(frame.props[1])
        assert len(caplog.text) == 0

        assert len(frame.props[1].scaling_factor) == 8
        assert len(frame.props[1].scaling_hint) == 2

        assert get_scaling_factor(frame.props[1].temperature) == 1 / 300
        assert get_scaling_factor(frame.props[1].pressure) == 1e-5
        for vardata in frame.props[1].mole_frac_phase_comp.values():
            assert get_scaling_factor(vardata) == 10
        for exprdata in frame.props[1].flow_mol_phase.values():
            assert get_scaling_factor(exprdata) == 1
        # p1 and p2 have an unknown phase type, so no scaling factor/hint is set.
        for vardata in frame.props[1].therm_cond_phase.values():
            assert get_scaling_factor(vardata) is None

    @pytest.mark.unit
    def test_build_all(self, frame):
        # Add necessary parameters
        frame.params.a.mw = 1
        frame.params.b.mw = 2
        frame.params.c.mw = 3

        # Add variables that would be built by state definition
        frame.props[1].flow_mol_phase_comp = Var(frame.props[1].phase_component_set)

        # Call all properties in metadata and assert they exist.
        for p in frame.params.get_metadata().properties.list_supported_properties():
            if p.name.endswith(("apparent", "true")):
                # True and apparent properties require electrolytes, which are
                # not tested here
                # Check that method exists and continue
                if frame.params.get_metadata().properties[p.name].method is not None:
                    assert hasattr(
                        frame.props[1],
                        frame.params.get_metadata().properties[p.name].method,
                    )
                continue
            elif p.name.endswith(("_crit")):
                # Critical properties will be tested elsewhere
                assert hasattr(
                    frame.props[1],
                    frame.params.get_metadata().properties[p.name].method,
                )
                continue
            elif p.name.endswith(("bubble", "bub", "dew")):
                # Bubble and dew properties require phase equilibria, which are
                # not tested here
                # Check that method exists and continue
                assert hasattr(
                    frame.props[1],
                    frame.params.get_metadata().properties[p.name].method,
                )
                continue
            elif p.name in ["dh_rxn", "log_k_eq"]:
                # Not testing inherent reactions here either
                # Check that method exists and continue
                assert hasattr(
                    frame.props[1],
                    frame.params.get_metadata().properties[p.name].method,
                )
                continue
            elif p.name in {
                "diffus_phase_comp",
            }:
                # phase indexed properties - these will be tested separately.
                continue
            elif p.supported:
                assert hasattr(frame.props[1], p.name)
            else:
                assert not hasattr(frame.props[1], p.name)

    @pytest.mark.unit
    def test_flows(self, frame):

        # Need to set the material flow basis for this test
        def set_flow_basis():
            return MaterialFlowBasis.molar

        frame.props[1].get_material_flow_basis = set_flow_basis

        frame.props[1].flow_mol.fix(100)
        frame.props[1].phase_frac["p1"].fix(0.75)
        frame.props[1].phase_frac["p2"].fix(0.25)
        frame.props[1].mole_frac_phase_comp["p1", "a"].fix(0.1)
        frame.props[1].mole_frac_phase_comp["p1", "b"].fix(0.2)
        frame.props[1].mole_frac_phase_comp["p1", "c"].fix(0.7)
        frame.props[1].mole_frac_phase_comp["p2", "c"].fix(0.1)
        frame.props[1].mole_frac_phase_comp["p2", "b"].fix(0.2)
        frame.props[1].mole_frac_phase_comp["p2", "a"].fix(0.7)
        frame.props[1].params.get_component("a").mw = 2
        frame.props[1].params.get_component("b").mw = 3
        frame.props[1].params.get_component("c").mw = 4

        assert value(frame.props[1].flow_vol) == pytest.approx(100 / 55e3, rel=1e-4)
        assert value(frame.props[1].flow_vol_phase["p1"]) == pytest.approx(
            100 / 55e3 * 0.75, rel=1e-4
        )
        assert value(frame.props[1].flow_vol_phase["p2"]) == pytest.approx(
            100 / 55e3 * 0.25, rel=1e-4
        )
        assert value(frame.props[1].flow_mol_comp["a"]) == pytest.approx(
            100 * 0.1 * 0.75 + 100 * 0.7 * 0.25, rel=1e-4
        )
        assert value(frame.props[1].flow_mol_comp["b"]) == pytest.approx(
            100 * 0.2 * 0.75 + 100 * 0.2 * 0.25, rel=1e-4
        )
        assert value(frame.props[1].flow_mol_comp["c"]) == pytest.approx(
            100 * 0.1 * 0.25 + 100 * 0.7 * 0.75, rel=1e-4
        )
        assert value(frame.props[1].mw_comp["a"]) == pytest.approx(2, rel=1e-4)
        assert value(frame.props[1].mw_phase["p1"]) == pytest.approx(
            2 * 0.1 + 3 * 0.2 + 4 * 0.7, rel=1e-4
        )
        assert value(frame.props[1].mw) == pytest.approx(
            0.75 * (2 * 0.1 + 3 * 0.2 + 4 * 0.7) + 0.25 * (2 * 0.7 + 3 * 0.2 + 4 * 0.1),
            rel=1e-4,
        )
        assert value(frame.props[1].flow_mass) == pytest.approx(
            100 * 0.75 * (2 * 0.1 + 3 * 0.2 + 4 * 0.7)
            + 100 * 0.25 * (2 * 0.7 + 3 * 0.2 + 4 * 0.1),
            rel=1e-4,
        )
        assert value(frame.props[1].flow_mass_phase["p1"]) == pytest.approx(
            100 * 0.75 * (2 * 0.1 + 3 * 0.2 + 4 * 0.7), rel=1e-4
        )
        assert value(frame.props[1].flow_mass_phase["p2"]) == pytest.approx(
            100 * 0.25 * (2 * 0.7 + 3 * 0.2 + 4 * 0.1), rel=1e-4
        )

    class dummy_prop:
        def return_expression(*args, **kwargs):
            return 4

        class visc_d_phase_comp(object):
            def return_expression(*args, **kwargs):
                return 5

        class therm_cond_phase_comp(object):
            def return_expression(*args, **kwargs):
                return 7

    @pytest.mark.unit
    def test_diffus_phase_comp(self, frame):
        frame.params.a.config.diffus_phase_comp = {
            "p1": TestGenericStateBlock.dummy_prop,
            "p2": TestGenericStateBlock.dummy_prop,
        }
        frame.params.b.config.diffus_phase_comp = {
            "p2": TestGenericStateBlock.dummy_prop
        }
        frame.params.c.config.diffus_phase_comp = {
            "p1": TestGenericStateBlock.dummy_prop
        }

        # There should be two skipped indices, so length should be 4
        assert len(frame.props[1].diffus_phase_comp) == 4
        assert value(frame.props[1].diffus_phase_comp["p1", "a"]) == 4
        assert value(frame.props[1].diffus_phase_comp["p2", "a"]) == 4
        assert value(frame.props[1].diffus_phase_comp["p2", "b"]) == 4
        assert value(frame.props[1].diffus_phase_comp["p1", "c"]) == 4

    @pytest.mark.unit
    def test_visc_d_phase_comp(self, frame):
        frame.params.a.config.visc_d_phase_comp = {
            "p1": TestGenericStateBlock.dummy_prop,
            "p2": TestGenericStateBlock.dummy_prop,
        }
        frame.params.b.config.visc_d_phase_comp = {
            "p2": TestGenericStateBlock.dummy_prop
        }
        frame.params.c.config.visc_d_phase_comp = {
            "p1": TestGenericStateBlock.dummy_prop
        }
        # Have to explicitly call constructor
        assert not hasattr(frame.props[1], "_visc_d_phase_comp")
        frame.props[1]._make_visc_d_phase_comp()

        # There should be two skipped indices, so length should be 4
        assert len(frame.props[1]._visc_d_phase_comp) == 4
        assert value(frame.props[1]._visc_d_phase_comp["p1", "a"]) == 5
        assert value(frame.props[1]._visc_d_phase_comp["p2", "a"]) == 5
        assert value(frame.props[1]._visc_d_phase_comp["p2", "b"]) == 5
        assert value(frame.props[1]._visc_d_phase_comp["p1", "c"]) == 5

    @pytest.mark.unit
    def test_therm_cond_phase_comp(self, frame):
        frame.params.a.config.therm_cond_phase_comp = {
            "p1": TestGenericStateBlock.dummy_prop,
            "p2": TestGenericStateBlock.dummy_prop,
        }
        frame.params.b.config.therm_cond_phase_comp = {
            "p2": TestGenericStateBlock.dummy_prop
        }
        frame.params.c.config.therm_cond_phase_comp = {
            "p1": TestGenericStateBlock.dummy_prop
        }
        # Have to explicitly call constructor
        assert not hasattr(frame.props[1], "_therm_cond_phase_comp")
        frame.props[1]._make_therm_cond_phase_comp()

        # There should be two skipped indices, so length should be 4
        assert len(frame.props[1]._therm_cond_phase_comp) == 4
        assert value(frame.props[1]._therm_cond_phase_comp["p1", "a"]) == 7
        assert value(frame.props[1]._therm_cond_phase_comp["p2", "a"]) == 7
        assert value(frame.props[1]._therm_cond_phase_comp["p2", "b"]) == 7
        assert value(frame.props[1]._therm_cond_phase_comp["p1", "c"]) == 7

    @pytest.mark.unit
    def test_surf_tens_phase(self, frame):
        # Surface tension is only built for Liquid Phases, and we have none.
        # All values should be None.
        for p in frame.props[1].phase_list:
            assert p not in frame.props[1].surf_tens_phase

    @pytest.mark.unit
    def test_therm_cond_phase(self, frame):
        frame.params.p1.config.therm_cond_phase = TestGenericStateBlock.dummy_prop
        frame.params.p2.config.therm_cond_phase = TestGenericStateBlock.dummy_prop

        for p in frame.props[1].phase_list:
            assert value(frame.props[1].therm_cond_phase[p]) == 4

    @pytest.mark.unit
    def test_visc_d_phase(self, frame):
        frame.params.p1.config.visc_d_phase = TestGenericStateBlock.dummy_prop
        frame.params.p2.config.visc_d_phase = TestGenericStateBlock.dummy_prop

        for p in frame.props[1].phase_list:
            assert value(frame.props[1].visc_d_phase[p]) == 4


class TestCriticalProps:
    @pytest.mark.unit
    def test_get_critical_ref_phase_exception(self):
        # No vapor or liquid phases
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "equation_of_state": DummyEoS,
                },
                "p2": {
                    "equation_of_state": DummyEoS,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        with pytest.raises(
            AttributeError,
            match="No liquid or vapor phase found to use as reference phase "
            "for calculating critical properties.",
        ):
            m.props[1]._get_critical_ref_phase()

    @pytest.mark.unit
    def test_get_critical_ref_phase_exception(self):
        # No vapor or liquid phases
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "equation_of_state": DummyEoS,
                },
                "p2": {
                    "equation_of_state": DummyEoS,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        with pytest.raises(
            PropertyPackageError,
            match="No liquid or vapor phase found to use as reference phase "
            "for calculating critical properties.",
        ):
            m.props[1]._get_critical_ref_phase()

    @pytest.mark.unit
    def test_get_critical_ref_phase_VL(self):
        # No vapor or liquid phases
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS,
                },
                "p2": {
                    "type": VaporPhase,
                    "equation_of_state": DummyEoS,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        assert m.props[1]._get_critical_ref_phase() == "p1"

    @pytest.mark.unit
    def test_get_critical_ref_phase_LL(self):
        # No vapor or liquid phases
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS,
                },
                "p2": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        assert m.props[1]._get_critical_ref_phase() == "p1"

    @pytest.mark.unit
    def test_get_critical_ref_phase_VLL(self):
        class DummyPE:
            def phase_equil(self, *args, **kwargs):
                pass

            @staticmethod
            def return_expression(b, *args, **kwargs):
                return b.temperature == 42

        # No vapor or liquid phases
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {
                    "phase_equilibrium_form": {("p2", "p3"): DummyPE},
                },
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS,
                },
                "p2": {
                    "type": VaporPhase,
                    "equation_of_state": DummyEoS,
                },
                "p3": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
            phases_in_equilibrium=[("p2", "p3")],
            phase_equilibrium_state={("p2", "p3"): DummyPE},
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        assert m.props[1]._get_critical_ref_phase() == "p3"

    @pytest.mark.unit
    def test_get_critical_ref_phase_VX(self):
        # No vapor or liquid phases
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "type": VaporPhase,
                    "equation_of_state": DummyEoS,
                },
                "p2": {
                    "equation_of_state": DummyEoS,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        assert m.props[1]._get_critical_ref_phase() == "p1"

    @pytest.mark.unit
    def test_critical_props_not_implemented(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        with pytest.raises(
            NotImplementedError,
            match=re.escape(
                "props[1] Equation of State module has not implemented a method for "
                "build_critical_properties. Please contact the EoS developer or use a "
                "different module."
            ),
        ):
            m.props[1]._critical_props()

    @pytest.mark.unit
    def test_critical_props_dummy_method(self):
        class DummyEoS2(DummyEoS):
            @staticmethod
            def build_critical_properties(b, *args, **kwargs):
                b._dummy_crit_executed = True

        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS2,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)
        m.props[1]._critical_props()

        assert m.props[1]._dummy_crit_executed

        assert isinstance(m.props[1].compress_fact_crit, Var)
        assert isinstance(m.props[1].dens_mol_crit, Var)
        assert isinstance(m.props[1].pressure_crit, Var)
        assert isinstance(m.props[1].temperature_crit, Var)

    @pytest.mark.unit
    def test_critical_props_attribute_error(self):
        class DummyEoS2(DummyEoS):
            @staticmethod
            def build_critical_properties(b, *args, **kwargs):
                raise AttributeError()

        m = ConcreteModel()
        m.params = DummyParameterBlock(
            components={
                "a": {},
            },
            phases={
                "p1": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS2,
                },
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units=base_units,
        )

        m.props = m.params.build_state_block([1], defined_state=False)
        try:
            m.props[1]._critical_props()
        except AttributeError:
            assert not hasattr(m.props[1], "compress_fact_crit")
            assert not hasattr(m.props[1], "dens_mol_crit")
            assert not hasattr(m.props[1], "pressure_crit")
            assert not hasattr(m.props[1], "temperature_crit")

    @pytest.mark.unit
    def test_initialize_critical_props(self):
        m = ConcreteModel()

        class DummyEoS2(DummyEoS):
            @staticmethod
            def build_critical_properties(b, *args, **kwargs):
                b._dummy_crit_executed = True

        # Dummy params block
        m.params = DummyParameterBlock(
            components={
                "a": {
                    "parameter_data": {
                        "compress_fact_crit": 0.1,
                        "dens_mol_crit": 1,
                        "pressure_crit": 1e5,
                        "temperature_crit": 100,
                    }
                },
                "b": {
                    "parameter_data": {
                        "compress_fact_crit": 0.2,
                        "dens_mol_crit": 2,
                        "pressure_crit": 2e5,
                        "temperature_crit": 200,
                    }
                },
                "c": {
                    "parameter_data": {
                        "compress_fact_crit": 0.3,
                        "dens_mol_crit": 3,
                        "pressure_crit": 3e5,
                        "temperature_crit": 300,
                    }
                },
            },
            phases={
                "Vap": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS2,
                },
                "Liq": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS2,
                },
            },
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            parameter_data={
                "PR_kappa": {
                    ("a", "a"): 0.0,
                    ("a", "b"): 0.0,
                    ("a", "c"): 0.0,
                    ("b", "a"): 0.0,
                    ("b", "b"): 0.0,
                    ("b", "c"): 0.0,
                    ("c", "a"): 0.0,
                    ("c", "b"): 0.0,
                    ("c", "c"): 0.0,
                },
            },
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        # Add common variables
        m.props[1].mole_frac_comp = Var(
            m.params.component_list, initialize=0.5, bounds=(1e-12, 1)
        )
        m.props[1].mole_frac_comp["a"].fix(0.4)
        m.props[1].mole_frac_comp["b"].fix(0.6)
        m.props[1].mole_frac_comp["c"].fix(1e-8)

        # Build critical props
        m.props[1]._critical_props()

        # Initialize critical props
        _initialize_critical_props(m.props[1])

        # these equalities are expected to be as exact as possible,
        # i.e. very small tolerances
        almost_exactly = functools.partial(pytest.approx, rel=1e-12)
        assert value(m.props[1].compress_fact_crit) == almost_exactly(
            0.4 * 0.1 + 0.6 * 0.2 + 1e-8 * 0.3
        )
        assert value(m.props[1].dens_mol_crit) == almost_exactly(
            0.4 * 1 + 0.6 * 2 + 1e-8 * 3
        )
        assert value(m.props[1].pressure_crit) == almost_exactly(
            0.4 * 1e5 + 0.6 * 2e5 + 1e-8 * 3e5
        )
        assert value(m.props[1].temperature_crit) == almost_exactly(
            0.4 * 1e2 + 0.6 * 2e2 + 1e-8 * 3e2
        )

    @pytest.mark.unit
    def test_initialize_critical_props_missing_value(self):
        m = ConcreteModel()

        class DummyEoS2(DummyEoS):
            @staticmethod
            def build_critical_properties(b, *args, **kwargs):
                b._dummy_crit_executed = True

        # Dummy params block
        m.params = DummyParameterBlock(
            components={
                "a": {
                    "parameter_data": {
                        "dens_mol_crit": 1,  # Missing Z_crit
                        "pressure_crit": 1e5,
                        "temperature_crit": 100,
                    }
                },
                "b": {
                    "parameter_data": {
                        "compress_fact_crit": 0.2,
                        "dens_mol_crit": 2,
                        "pressure_crit": 2e5,
                        "temperature_crit": 200,
                    }
                },
                "c": {
                    "parameter_data": {
                        "compress_fact_crit": 0.3,
                        "dens_mol_crit": 3,
                        "pressure_crit": 3e5,
                        "temperature_crit": 300,
                    }
                },
            },
            phases={
                "Vap": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS2,
                },
                "Liq": {
                    "type": LiquidPhase,
                    "equation_of_state": DummyEoS2,
                },
            },
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
            state_definition=modules[__name__],
            pressure_ref=100000.0,
            temperature_ref=300,
            parameter_data={
                "PR_kappa": {
                    ("a", "a"): 0.0,
                    ("a", "b"): 0.0,
                    ("a", "c"): 0.0,
                    ("b", "a"): 0.0,
                    ("b", "b"): 0.0,
                    ("b", "c"): 0.0,
                    ("c", "a"): 0.0,
                    ("c", "b"): 0.0,
                    ("c", "c"): 0.0,
                },
            },
        )

        m.props = m.params.build_state_block([1], defined_state=False)

        # Add common variables
        m.props[1].mole_frac_comp = Var(
            m.params.component_list, initialize=0.5, bounds=(1e-12, 1)
        )
        m.props[1].mole_frac_comp["a"].fix(0.4)
        m.props[1].mole_frac_comp["b"].fix(0.6)
        m.props[1].mole_frac_comp["c"].fix(1e-8)

        # Build critical props
        m.props[1]._critical_props()

        # Initialize critical props
        with pytest.raises(
            AttributeError,
            match="Missing attribute found when initializing compress_fact_crit. "
            "Make sure you have provided values for compress_fact_crit in all "
            "Component declarations.",
        ):
            _initialize_critical_props(m.props[1])


# Invalid property configuration to trigger configuration error
configuration = {
    # Specifying components
    "components": {
        "H2O": {
            "type": Component,
            "parameter_data": {
                "pressure_crit": (220.6e5, pyunits.Pa),
                "temperature_crit": (647, pyunits.K),
                "omega": 0.344,
            },
        },
    },
    # Specifying phases
    "phases": {
        "Liq": {
            "type": LiquidPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
        "Vap": {
            "type": VaporPhase,
            "equation_of_state": Cubic,
            "equation_of_state_options": {"type": CubicType.PR},
        },
    },
    # Set base units of measurement
    "base_units": {
        "time": pyunits.s,
        "length": pyunits.m,
        "mass": pyunits.kg,
        "amount": pyunits.mol,
        "temperature": pyunits.K,
    },
    # Specifying state definition
    "state_definition": FTPx,
    "state_bounds": {
        "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
        "temperature": (273.15, 300, 500, pyunits.K),
        "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
    },
    "pressure_ref": (101325, pyunits.Pa),
    "temperature_ref": (298.15, pyunits.K),
    "parameter_data": {
        "PR_kappa": {
            ("foo", "bar"): 0.000,
        }
    },
}


@pytest.mark.integration
def test_phase_component_flash():
    # Regression test for issue #1423
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.props = GenericParameterBlock(**BTconfig)

    m.fs.flash = Flash(
        property_package=m.fs.props,
        material_balance_type=MaterialBalanceType.componentPhase,
    )

    assert m.fs.props.phase_equilibrium_list == {
        "PE1": ["benzene", ("Vap", "Liq")],
        "PE2": ["toluene", ("Vap", "Liq")],
    }
    assert isinstance(m.fs.flash, Flash)
