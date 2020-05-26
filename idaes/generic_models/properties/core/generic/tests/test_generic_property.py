##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for generic property package core code

Author: Andrew Lee
"""
import pytest
from sys import modules

from pyomo.environ import Block, ConcreteModel, Param, Set, Var

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData,
        GenericStateBlock)
from idaes.generic_models.properties.core.generic.tests import dummy_eos

from idaes.core import (declare_process_block_class, Component,
                        Phase, LiquidPhase, VaporPhase)
from idaes.core.phases import PhaseType as PT
from idaes.core.util.exceptions import (ConfigurationError)


# -----------------------------------------------------------------------------
# Dummy method to avoid errors when setting metadata dict
def set_metadata(b):
    pass


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True


class TestGenericParameterBlock(object):
    def test_build(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"type": LiquidPhase,
                           "component_list": ["a", "b"],
                           "equation_of_state": dummy_eos},
                    "p2": {"equation_of_state": dummy_eos}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

        assert m.params.configured

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
        assert m.params.p1.config.equation_of_state == dummy_eos
        assert m.params.p2.config.equation_of_state == dummy_eos

        assert isinstance(m.params._phase_component_set, Set)
        assert len(m.params._phase_component_set) == 5
        for i in m.params._phase_component_set:
            assert i in [("p1", "a"), ("p1", "b"),
                         ("p2", "a"), ("p2", "b"), ("p2", "c")]

        assert isinstance(m.params.pressure_ref, Param)
        assert m.params.pressure_ref.value == 1e5
        assert isinstance(m.params.temperature_ref, Param)
        assert m.params.temperature_ref.value == 300

        assert m.params.state_block_class is GenericStateBlock

    def test_no_components(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params was not provided with a components "
                           "argument."):
            m.params = DummyParameterBlock(default={})

    def test_no_phases(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params was not provided with a phases "
                           "argument."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}}})

    def test_invalid_component_in_phase_component_list(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params phase-component list for phase p1 "
                           "contained component d which is not in the master "
                           "component list"):
            m.params = DummyParameterBlock(default={
                    "components": {"a": {}, "b": {}, "c": {}},
                    "phases": {
                        "p1": {"type": LiquidPhase,
                               "component_list": ["a", "d"],
                               "equation_of_state": "foo"},
                        "p2": {"equation_of_state": "bar"}}})

    def test_invalid_component_in_phase_component_list_2(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params phase-component list for phase p1 "
                           "contained component a, however this component is "
                           "not valid for the given PhaseType"):
            m.params = DummyParameterBlock(default={
                    "components": {
                        "a": {"valid_phase_types": PT.solidPhase},
                        "b": {},
                        "c": {}},
                    "phases": {
                        "p1": {"type": LiquidPhase,
                               "component_list": ["a", "b"],
                               "equation_of_state": "foo"},
                        "p2": {"equation_of_state": "bar"}}})

    def test_phase_component_set_from_valid_phases(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"valid_phase_types": PT.liquidPhase},
                    "b": {"valid_phase_types": PT.vaporPhase},
                    "c": {"valid_phase_types": [PT.liquidPhase,
                                                PT.vaporPhase]}},
                "phases": {
                    "p1": {"type": LiquidPhase,
                           "equation_of_state": dummy_eos},
                    "p2": {"type": VaporPhase,
                           "equation_of_state": dummy_eos}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

        assert len(m.params._phase_component_set) == 4
        for i in m.params._phase_component_set:
            assert i in [("p1", "a"), ("p1", "c"),
                         ("p2", "b"), ("p2", "c")]

    def test_no_state_definition(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was not "
                           "provided with a state_definition configuration "
                           "argument. Please fix your property parameter "
                           "definition to include this."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}}})

    def test_no_pressure_ref(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was not "
                           "provided with a pressure_ref configuration "
                           "argument. Please fix your property parameter "
                           "definition to include this."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {}},
                "state_definition": "baz"})

    def test_temperature_ref(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was not "
                           "provided with a temperature_ref configuration "
                           "argument. Please fix your property parameter "
                           "definition to include this."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {}},
                "state_definition": "baz",
                "pressure_ref": 1e5})

    def test_no_eos(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params phase p2 was not provided with an "
                           "equation_of_state configuration argument. Please "
                           "fix your property parameter definition to "
                           "include this."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300})

    def test_phases_in_equilibrium(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
            "components": {
                "a": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                "b": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                "c": {"phase_equilibrium_form": {("p1", "p2"): "foo"}}},
            "phases": {
                "p1": {"equation_of_state": dummy_eos},
                "p2": {"equation_of_state": dummy_eos}},
            "state_definition": modules[__name__],
            "pressure_ref": 1e5,
            "temperature_ref": 300,
            "phases_in_equilibrium": [("p1", "p2")],
            "phase_equilibrium_state": {("p1", "p2"): "whoop"}})

        assert isinstance(m.params.phase_equilibrium_idx, Set)
        assert len(m.params.phase_equilibrium_idx) == 3
        for i in m.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert isinstance(m.params.phase_equilibrium_list, dict)
        assert m.params.phase_equilibrium_list == {
            "PE1": {"a": ("p1", "p2")},
            "PE2": {"b": ("p1", "p2")},
            "PE3": {"c": ("p1", "p2")}}

    def test_phases_in_equilibrium_no_form(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package component a "
                                 "is in equilibrium but phase_equilibrium_form"
                                 " was not specified."):
            m.params = DummyParameterBlock(default={
                "components": {
                    "a": {},
                    "b": {},
                    "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("p1", "p2")],
                "phase_equilibrium_state": {("p1", "p2"): "whoop"}})

    def test_phases_in_equilibrium_missing_pair_form(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package component b "
                                 "is in equilibrium but phase_equilibrium_form"
                                 " was not specified for all appropriate "
                                 "phase pairs."):
            # Also reverse order of phases for component a - this should pass
            # and component b should be flagged as missing
            m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"phase_equilibrium_form": {("p2", "p1"): "foo"}},
                    "b": {"phase_equilibrium_form": {(1, 2): "foo"}},
                    "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("p1", "p2")],
                "phase_equilibrium_state": {("p1", "p2"): "whoop"}})

    def test_phases_in_equilibrium_no_formulation(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with a phases_in_equilibrium argument but no "
                           "method was specified for "
                           "phase_equilibrium_state."):
            m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "b": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "c": {"phase_equilibrium_form": {("p1", "p2"): "foo"}}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("p1", "p2")]})

    def test_phases_in_equilibrium_missing_pair_formulation(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with a phases_in_equilibrium argument but "
                           "phase_equilibrium_state was not specified "
                           "for all phase pairs."):
            # Also reverse order of phases for component a - this should pass
            # and component b should be flagged as missing
            m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "b": {"phase_equilibrium_form": {("p1", "p2"): "foo"}},
                    "c": {"phase_equilibrium_form": {("p1", "p2"): "foo"}}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("p1", "p2")],
                "phase_equilibrium_state": {(1, 2): "whoop"}})

    def test_parameter_construction_no_value(self):
        m = ConcreteModel()

        class test_class():
            # Mook up property method class for testing
            def build_parameters(c):
                c.test_var = Var()

        with pytest.raises(ConfigurationError,
                           match="params parameter test_var for component c "
                           "was not assigned a value. Please check your "
                           "configuration arguments."):
            m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"dens_mol_liq_comp": test_class},
                    "b": {},
                    "c": {}},
                "phases": {
                    "p1": {"equation_of_state": dummy_eos},
                    "p2": {"equation_of_state": dummy_eos}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

    def test_parameter_construction_no_data(self):
        m = ConcreteModel()

        class test_class():
            # Mook up property method class for testing
            def build_parameters(c):
                c.config.parameter_data["test"]

        with pytest.raises(ConfigurationError,
                           match="params values were not defined for "
                           "parameter dens_mol_liq_comp in component a. "
                           "Please check the parameter_data argument to "
                           "ensure values are provided."):
            m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"dens_mol_liq_comp": test_class},
                    "b": {},
                    "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

    def test_no_elements(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"type": LiquidPhase,
                           "component_list": ["a", "b"],
                           "equation_of_state": dummy_eos},
                    "p2": {"equation_of_state": dummy_eos}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

        assert not hasattr(m.params, "element_list")
        assert not hasattr(m.params, "element_comp")

    def test_partial_elements(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params not all Components declared an "
                           "elemental_composition argument. Either all "
                           "Components must declare this, or none."):
            m.params = DummyParameterBlock(default={
                    "components": {"a": {"elemental_composition": {"e1": 1}},
                                   "b": {},
                                   "c": {}},
                    "phases": {
                        "p1": {"type": LiquidPhase,
                               "component_list": ["a", "b"],
                               "equation_of_state": dummy_eos},
                        "p2": {"equation_of_state": dummy_eos}},
                    "state_definition": modules[__name__],
                    "pressure_ref": 1e5,
                    "temperature_ref": 300})

    def test_elements_not_float(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params values in elemental_composition must "
                           "be integers \(not floats\)\: e1\: 2.0."):
            m.params = DummyParameterBlock(default={
                    "components": {"a": {"elemental_composition": {"e1": 2.0}},
                                   "b": {},
                                   "c": {}},
                    "phases": {
                        "p1": {"type": LiquidPhase,
                               "component_list": ["a", "b"],
                               "equation_of_state": dummy_eos},
                        "p2": {"equation_of_state": dummy_eos}},
                    "state_definition": modules[__name__],
                    "pressure_ref": 1e5,
                    "temperature_ref": 300})

    def test_elements(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
                "components": {
                    "a": {"elemental_composition": {"e1": 1, "e2": 2}},
                    "b": {"elemental_composition": {"e3": 3, "e4": 4}},
                    "c": {"elemental_composition": {"e1": 5, "e3": 6}}},
                "phases": {
                    "p1": {"type": LiquidPhase,
                           "component_list": ["a", "b"],
                           "equation_of_state": dummy_eos},
                    "p2": {"equation_of_state": dummy_eos}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

        assert isinstance(m.params.element_list, Set)
        assert len(m.params.element_list) == 4
        assert m.params.element_comp == {
            "a": {"e1": 1, "e2": 2, "e3": 0, "e4": 0},
            "b": {"e1": 0, "e2": 0, "e3": 3, "e4": 4},
            "c": {"e1": 5, "e2": 0, "e3": 6, "e4": 0}}


# -----------------------------------------------------------------------------
# Dummy methods for testing build calls to sub-modules
def define_state(b):
    b.state_defined = True


class TestGenericStateBlock(object):
    @pytest.fixture()
    def frame(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": dummy_eos},
                    "p2": {"equation_of_state": dummy_eos}},
                "state_definition": modules[__name__],
                "pressure_ref": 1e5,
                "temperature_ref": 300})

        m.props = m.params.build_state_block([1],
                                             default={"defined_state": False})

        # Add necessary variables to state block
        m.props[1].pressure = Var(bounds=(1000, 3000))
        m.props[1].temperature = Var(bounds=(100, 200))
        m.props[1].mole_frac_phase_comp = Var(m.params.phase_list,
                                              m.params.component_list)
        m.props[1].phase_frac = Var(m.params.phase_list)

        return m

    def test_build(self, frame):
        assert isinstance(frame.props, Block)
        assert len(frame.props) == 1

        # Check for expected behaviour for dummy methods
        assert frame.props[1].state_defined
        assert isinstance(frame.props[1].dummy_var, Var)
        assert frame.props[1].eos_common == 2
