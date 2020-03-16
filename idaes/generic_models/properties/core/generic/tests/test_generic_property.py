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

from pyomo.environ import Block, ConcreteModel, Expression, Param, Set, Var
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericPropertyPackageError,
        get_method,
        GenericParameterData,
        GenericStateBlock)
from idaes.generic_models.properties.core.generic.tests import dummy_eos

from idaes.core import (declare_process_block_class, Component,
                        Phase, LiquidPhase)
from idaes.core.util.exceptions import (ConfigurationError,
                                        PropertyPackageError)
from idaes.core.util.misc import add_object_reference

# -----------------------------------------------------------------------------
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
                         "equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "foo",
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
        assert m.params.p1.config.equation_of_state == "foo"
        assert m.params.p2.config.equation_of_state == "bar"

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

    def test_phase_equilibrium_both_definitions(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was "
                           "provided with both a phases_in_equilibrium and a "
                           "phase_equilibrium_dict argument. Users should "
                           "provide only one of these."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("p1", "p2")],
                "phase_equilibrium_dict": {"r1": ("a", "p1", "p2")}})

    def test_phases_in_equilibrium(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
            "components": {"a": {}, "b": {}, "c": {}},
            "phases": {
                "p1": {"equation_of_state": "foo"},
                "p2": {"equation_of_state": "bar"}},
            "state_definition": "baz",
            "pressure_ref": 1e5,
            "temperature_ref": 300,
            "phases_in_equilibrium": [("p1", "p2")],
            "phase_equilibrium_formulation": "whoop"})

        assert isinstance(m.params.phase_equilibrium_idx, Set)
        assert len(m.params.phase_equilibrium_idx) == 3
        for i in m.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert isinstance(m.params.phase_equilibrium_list, dict)
        assert m.params.phase_equilibrium_list == {
            "PE1": {"a": ("p1", "p2")},
            "PE2": {"b": ("p1", "p2")},
            "PE3": {"c": ("p1", "p2")}}

    def test_phase_equilibrium_dict(self):
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
            "components": {"a": {}, "b": {}, "c": {}},
            "phases": {
                "p1": {"equation_of_state": "foo"},
                "p2": {"equation_of_state": "bar"}},
            "state_definition": "baz",
            "pressure_ref": 1e5,
            "temperature_ref": 300,
            "phase_equilibrium_dict": {"r1": ("a", "p1", "p2")},
            "phase_equilibrium_formulation": "whoop"})

        assert isinstance(m.params.phase_equilibrium_idx, Set)
        assert len(m.params.phase_equilibrium_idx) == 1
        for i in m.params.phase_equilibrium_idx:
            assert i in ["r1"]

        assert isinstance(m.params.phase_equilibrium_list, dict)
        assert m.params.phase_equilibrium_list == {
            "r1": {"a": ("p1", "p2")}}

    def test_phases_in_equilibrium_no_formulation(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with a phases_in_equilibrium or "
                           "phase_equilibrium_dict argument, but no method "
                           "was specified for phase_equilibrium_formulation."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phases_in_equilibrium": [("p1", "p2")]})

    def test_phases_equilibrium_dict_no_formulation(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with a phases_in_equilibrium or "
                           "phase_equilibrium_dict argument, but no method "
                           "was specified for phase_equilibrium_formulation."):
            m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    "p1": {"equation_of_state": "foo"},
                    "p2": {"equation_of_state": "bar"}},
                "state_definition": "baz",
                "pressure_ref": 1e5,
                "temperature_ref": 300,
                "phase_equilibrium_dict": {"r1": ("a", "p1", "p2")}})

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

        m.props = m.params.state_block_class([1],
                                             default={"defined_state": False,
                                                      "parameters": m.params})

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
        # assert frame.props[1].phase_equil_defined
