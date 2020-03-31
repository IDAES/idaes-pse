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
        GenericParameterData)
from idaes.generic_models.properties.core.generic.tests import dummy_eos

from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import (ConfigurationError,
                                        PropertyPackageError)
from idaes.core.util.misc import add_object_reference


# -----------------------------------------------------------------------------
supported_properties = ["phase_component_list",
                        "state_bounds",
                        "phase_equilibrium_formulation",
                        "phase_equilibrium_dict",
                        "temperature_bubble",
                        "temperature_dew",
                        "pressure_bubble",
                        "pressure_dew",
                        "dens_mol_liq_comp",
                        "enth_mol_liq_comp",
                        "enth_mol_ig_comp",
                        "entr_mol_liq_comp",
                        "entr_mol_ig_comp",
                        "pressure_sat_comp"]


# -----------------------------------------------------------------------------
def test_GenericPropertyPackageError():
    with pytest.raises(PropertyPackageError):
        raise GenericPropertyPackageError("block", "prop")


# Dummy method for testing get_method
def dummy_option(b):
    return b.dummy_response


class TestGetMethod(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = Block()

        # Add necessary parameters to parameter block
        m.params.config = ConfigBlock()
        m.params.config.declare("dummy_option", ConfigValue(default=None))

        # Create a dummy state block
        m.props = Block([1])
        m.props[1].config = ConfigBlock()
        add_object_reference(m.props[1], "params", m.params)

        m.props[1].dummy_response = "foo"

        return m

    def test_None(self, frame):
        with pytest.raises(GenericPropertyPackageError):
            get_method(frame.props[1], "dummy_option")

    def test_invlaid_option(self, frame):
        with pytest.raises(AttributeError):
            get_method(frame.props[1], "foo")

    def test_method(self, frame):
        # Test that get_method works when provided with the method directly
        frame.params.config.dummy_option = dummy_option

        assert get_method(frame.props[1], "dummy_option") is dummy_option
        assert get_method(frame.props[1], "dummy_option")(frame.props[1]) is \
            frame.props[1].dummy_response

    def test_module(self, frame):
        # Test that get_method works when pointed to a module with the method
        frame.params.config.dummy_option = modules[__name__]

        assert get_method(frame.props[1], "dummy_option") is dummy_option
        assert get_method(frame.props[1], "dummy_option")(frame.props[1]) is \
            frame.props[1].dummy_response

    def test_not_method_or_module(self, frame):
        # Test that get_method works when provided with the method directly
        frame.params.config.dummy_option = "foo"

        with pytest.raises(ConfigurationError):
            get_method(frame.props[1], "dummy_option")


# -----------------------------------------------------------------------------
@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True

    def parameters(self):
        self.parameters_set = True


class TestGenericParameterBlock(object):
    def test_configure(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"}})

        assert m.params.config.component_list == ["a", "b", "c"]
        assert m.params.config.phase_list == [1, 2]
        assert m.params.config.state_definition == "foo"
        assert m.params.config.equation_of_state == {1: "foo", 2: "bar"}

        # Tesrt number of config arguments. Note 1 inherited argument
        assert len(m.params.config) == len(supported_properties) + 4 + 1

        for i in supported_properties:
            assert m.params.config[i] is None

        assert m.params.configured

        assert isinstance(m.params.component_list, Set)
        assert len(m.params.component_list) == 3
        for j in m.params.component_list:
            assert j in ["a", "b", "c"]

        assert isinstance(m.params.phase_list, Set)
        assert len(m.params.phase_list) == 2
        for p in m.params.phase_list:
            assert p in [1, 2]

        assert isinstance(m.params._phase_component_set, Set)
        assert len(m.params._phase_component_set) == 6
        for v in m.params._phase_component_set:
            assert v[0] in [1, 2]
            assert v[1] in ["a", "b", "c"]

    def test_configure_phase_comp(self):
        # Phase-component list with invalid phase
        m = ConcreteModel()

        m.params = DummyParameterBlock(default={
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"],
                                         3: ["a", "b", "c"]},
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar", 3: "baz"}})

        assert isinstance(m.params.component_list, Set)
        assert len(m.params.component_list) == 3
        for j in m.params.component_list:
            assert j in ["a", "b", "c"]

        assert isinstance(m.params.phase_list, Set)
        assert len(m.params.phase_list) == 3
        for p in m.params.phase_list:
            assert p in [1, 2, 3]

        assert isinstance(m.params._phase_component_set, Set)
        assert len(m.params._phase_component_set) == 9
        for v in m.params._phase_component_set:
            assert v[0] in [1, 2, 3]
            assert v[1] in ["a", "b", "c"]

    def test_configure_no_components(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was not "
                           "provided with a component_list."):
            m.params = DummyParameterBlock()

    def test_configure_no_phases(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was not "
                           "provided with a phase_list."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"]})

    def test_configure_invalid_phase_comp_1(self):
        # Phase-component list with invalid phase
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params mismatch between phase_list and "
                            "phase_component_list. Phase 3 appears in "
                            "phase_component_list but not phase_list."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"],
                                         3: ["a", "b", "c"]}})

    def test_configure_invalid_phase_comp_2(self):
        # Phase-component list with invalid phase
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params mismatch between phase_list and "
                            "phase_component_list. Phase 3 appears in "
                            "phase_list but not phase_component_list."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2, 3],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"]}})

    def test_configure_invalid_phase_comp_3(self):
        # Phase-component list with invalid component
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params mismatch between component_list and "
                                "phase_component_list. Component d appears in "
                                "phase_component_list but not "
                                "component_list."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c", "d"],
                                         2: ["a", "b", "c", "d"]}})

    def test_configure_invalid_phase_comp_4(self):
        # Phase-component list with invalid component
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params mismatch between component_list and "
                                "phase_component_list. Component d appears in "
                                "component_list but not "
                                "phase_component_list."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c", "d"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"]}})

    def test_configure_no_state_definition(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was not "
                           "provided with a state_definition configuration "
                           "argument. Please fix your property parameter "
                           "definition to include this configuration "
                           "argument."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"]}})

    def test_configure_no_eos(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was not "
                           "provided with an equation_of_state configuration "
                           "argument. Please fix your property parameter "
                           "definition to include this configuration "
                           "argument."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"]},
                "state_definition": "foo"})

    def test_configure_eos_not_dict(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was provided"
                           " with an invalid equation_of_state configuration "
                           "argument. Argument must be a dict with phases as "
                           "keys."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"]},
                "state_definition": "foo",
                "equation_of_state": "bar"})

    def test_configure_eos_wrong_length(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package was provided"
                           " with an invalid equation_of_state configuration "
                           "argument. A value must be present for each "
                           "phase."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"]},
                "state_definition": "foo",
                "equation_of_state": {1: "bar"}})

    def test_configure_eos_invalid_phase(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package unrecognised"
                           " phase 3 in equation_of_state configuration "
                           "argument. Keys must be valid phases."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "phase_component_list": {1: ["a", "b", "c"],
                                         2: ["a", "b", "c"]},
                "state_definition": "foo",
                "equation_of_state": {1: "bar", 3: "baz"}})

    def test_configure_only_phase_equilibrium_formulation(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with only one of phase_equilibrium_formulation and"
                           " phase_equilibrium_dict. Either both of these "
                           "arguments need to be provided or neither."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo"})

    def test_configure_only_phase_equilibrium_dict(self):
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with only one of phase_equilibrium_formulation and"
                           " phase_equilibrium_dict. Either both of these "
                           "arguments need to be provided or neither."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_dict": {}})

    def test_configure_invalid_phase_equilibrium_dict_1(self):
        # Not dict
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with invalid phase_equilibrium_dict - value must "
                           "be a dict. Please see the documentation for the "
                           "correct form."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": "foo"})

    def test_configure_invalid_phase_equilibrium_dict_2(self):
        # Value not list
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with invalid phase_equilibrium_dict, foo. "
                           "Values in dict must be lists containing 2 "
                           "values."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": {1: "foo"}})

    def test_configure_invalid_phase_equilibrium_dict_3(self):
        # List with wrong number of values
        m = ConcreteModel()

        # Matching error text has been challenging for some reason
        with pytest.raises(ConfigurationError):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": {1: [1, 2, 3]}})

    def test_configure_invalid_phase_equilibrium_dict_4(self):
        # First value not component
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with invalid phase_equilibrium_dict. First value "
                           "in each list must be a valid component, received "
                           "foo."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": {1: ["foo", "bar"]}})

    def test_configure_invalid_phase_equilibrium_dict_5(self):
        # Second value tuple
        m = ConcreteModel()

        with pytest.raises(ConfigurationError,
                           match="params Generic Property Package provided "
                           "with invalid phase_equilibrium_dict. Second value "
                           "in each list must be a 2-tuple containing 2 valid "
                           "phases, received bar."):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": {1: ["a", "bar"]}})

    def test_configure_invalid_phase_equilibrium_dict_6(self):
        # Second value tuple, but wrong length
        m = ConcreteModel()

        # String matching difficult again
        with pytest.raises(ConfigurationError):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": {1: ["a", (1, 2, 3)]}})

    def test_configure_invalid_phase_equilibrium_dict_7(self):
        # Invalid phase in tuple
        m = ConcreteModel()

        # String matching difficult again
        with pytest.raises(ConfigurationError):
            m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": {1: ["a", (1, 3)]}})

    def test_configure_phase_equilibrium(self):
        # Invalid phase in tuple
        m = ConcreteModel()

        # String matching difficult again
        m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": "foo",
                "equation_of_state": {1: "foo", 2: "bar"},
                "phase_equilibrium_formulation": "foo",
                "phase_equilibrium_dict": {1: ["a", (1, 2)]}})

        assert m.params.phase_equilibrium_list == \
            m.params.config.phase_equilibrium_dict
        assert isinstance(m.params.phase_equilibrium_idx, Set)
        assert len(m.params.phase_equilibrium_idx) == \
            len(m.params.config.phase_equilibrium_dict)
        for k in m.params.phase_equilibrium_idx:
            assert k in m.params.config.phase_equilibrium_dict.keys()


# -----------------------------------------------------------------------------
# Dummy methods for testing build calls to sub-modules
def define_state(b):
    b.state_defined = True


def phase_equil(b):
    b.phase_equil_defined = True


def dummy_call(b):
    b.dummy_call = True


def pressure_sat_comp(b, j, T):
    return b.pressure


class TestGenericStateBlock(object):
    @pytest.fixture()
    def frame(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(default={
                "component_list": ["a", "b", "c"],
                "phase_list": [1, 2],
                "state_definition": modules[__name__],
                "equation_of_state": {1: dummy_eos,
                                      2: dummy_eos},
                "phase_equilibrium_formulation": modules[__name__],
                "phase_equilibrium_dict": {1: ["a", (1, 2)]}})

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
        assert frame.props[1].phase_equil_defined

    def test_temperture_bubble(self, frame):
        frame.params.config.temperature_bubble = dummy_call

        assert isinstance(frame.props[1].temperature_bubble, Var)
        assert frame.props[1].temperature_bubble.ub == 200
        assert frame.props[1].temperature_bubble.lb == 100

        assert isinstance(frame.props[1]._mole_frac_tbub, Var)
        assert len(frame.props[1]._mole_frac_tbub) == 3
        for i in frame.props[1]._mole_frac_tbub:
            assert i in frame.params.component_list

        assert frame.props[1].dummy_call

    def test_temperture_bubble_undefined(self, frame):
        with pytest.raises(GenericPropertyPackageError):
            frame.props[1].temperature_bubble()

    def test_temperture_dew(self, frame):
        frame.params.config.temperature_dew = dummy_call

        assert isinstance(frame.props[1].temperature_dew, Var)
        assert frame.props[1].temperature_dew.ub == 200
        assert frame.props[1].temperature_dew.lb == 100

        assert isinstance(frame.props[1]._mole_frac_tdew, Var)
        assert len(frame.props[1]._mole_frac_tdew) == 3
        for i in frame.props[1]._mole_frac_tdew:
            assert i in frame.params.component_list

        assert frame.props[1].dummy_call

    def test_temperture_dew_undefined(self, frame):
        with pytest.raises(GenericPropertyPackageError):
            frame.props[1].temperature_dew()

    def test_pressure_bubble(self, frame):
        frame.params.config.pressure_bubble = dummy_call

        assert isinstance(frame.props[1].pressure_bubble, Var)
        assert frame.props[1].pressure_bubble.ub == 3000
        assert frame.props[1].pressure_bubble.lb == 1000

        assert isinstance(frame.props[1]._mole_frac_pbub, Var)
        assert len(frame.props[1]._mole_frac_pbub) == 3
        for i in frame.props[1]._mole_frac_pbub:
            assert i in frame.params.component_list

        assert frame.props[1].dummy_call

    def test_pressure_bubble_undefined(self, frame):
        with pytest.raises(GenericPropertyPackageError):
            frame.props[1].pressure_bubble()

    def test_pressure_dew(self, frame):
        frame.params.config.pressure_dew = dummy_call

        assert isinstance(frame.props[1].pressure_dew, Var)
        assert frame.props[1].pressure_dew.ub == 3000
        assert frame.props[1].pressure_dew.lb == 1000

        assert isinstance(frame.props[1]._mole_frac_pdew, Var)
        assert len(frame.props[1]._mole_frac_pdew) == 3
        for i in frame.props[1]._mole_frac_pdew:
            assert i in frame.params.component_list

        assert frame.props[1].dummy_call

    def test_pressure_dew_undefined(self, frame):
        with pytest.raises(GenericPropertyPackageError):
            frame.props[1].pressure_dew()

    def test_dens_mass_phase(self, frame):
        assert isinstance(frame.props[1].dens_mass_phase, Expression)
        assert len(frame.props[1].dens_mass_phase) == 2
        for p in frame.props[1].dens_mass_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].dens_mass_phase[p].expr) == \
                str(frame.props[1].dummy_var)

    def test_dens_mass(self, frame):
        assert isinstance(frame.props[1].dens_mass, Expression)
        assert len(frame.props[1].dens_mass) == 1
        assert str(frame.props[1].dens_mass.expr) == \
            str(sum(frame.props[1].dens_mass_phase[p] *
                    frame.props[1].phase_frac[p]
                    for p in frame.props[1].params.phase_list))

        # Check that dependency variables also constructed properly
        assert isinstance(frame.props[1].dens_mass_phase, Expression)
        assert len(frame.props[1].dens_mass_phase) == 2
        for p in frame.props[1].dens_mass_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].dens_mass_phase[p].expr) == \
                str(frame.props[1].dummy_var)

    def test_dens_mol_phase(self, frame):
        assert isinstance(frame.props[1].dens_mol_phase, Expression)
        assert len(frame.props[1].dens_mol_phase) == 2
        for p in frame.props[1].dens_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].dens_mol_phase[p].expr) == \
                str(frame.props[1].dummy_var)

    def test_dens_mol(self, frame):
        assert isinstance(frame.props[1].dens_mol, Expression)
        assert len(frame.props[1].dens_mol) == 1
        assert str(frame.props[1].dens_mol.expr) == \
            str(sum(frame.props[1].dens_mol_phase[p] *
                    frame.props[1].phase_frac[p]
                    for p in frame.props[1].params.phase_list))

        # Check that dependency variables also constructed properly
        assert isinstance(frame.props[1].dens_mol_phase, Expression)
        assert len(frame.props[1].dens_mol_phase) == 2
        for p in frame.props[1].dens_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].dens_mol_phase[p].expr) == \
                str(frame.props[1].dummy_var)

    def test_enth_mol_phase_comp(self, frame):
        assert isinstance(frame.props[1].enth_mol_phase_comp, Expression)
        assert len(frame.props[1].enth_mol_phase_comp) == 6
        for k in frame.props[1].enth_mol_phase_comp:
            assert k[0] in frame.params.phase_list
            assert k[1] in frame.params.component_list
            assert str(frame.props[1].enth_mol_phase_comp[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_enth_mol_phase(self, frame):
        assert isinstance(frame.props[1].enth_mol_phase, Expression)
        assert len(frame.props[1].enth_mol_phase) == 2
        for k in frame.props[1].enth_mol_phase:
            assert k in frame.params.phase_list
            assert str(frame.props[1].enth_mol_phase[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_enth_mol(self, frame):
        assert isinstance(frame.props[1].enth_mol, Expression)
        assert len(frame.props[1].enth_mol) == 1
        assert str(frame.props[1].enth_mol.expr) == \
            str(sum(frame.props[1].enth_mol_phase[p] *
                    frame.props[1].phase_frac[p]
                    for p in frame.props[1].params.phase_list))

        # Check that dependency variables also constructed properly
        assert isinstance(frame.props[1].enth_mol_phase, Expression)
        assert len(frame.props[1].enth_mol_phase) == 2
        for p in frame.props[1].enth_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].enth_mol_phase[p].expr) == \
                str(frame.props[1].dummy_var)

    def test_entr_mol_phase_comp(self, frame):
        assert isinstance(frame.props[1].entr_mol_phase_comp, Expression)
        assert len(frame.props[1].entr_mol_phase_comp) == 6
        for k in frame.props[1].entr_mol_phase_comp:
            assert k[0] in frame.params.phase_list
            assert k[1] in frame.params.component_list
            assert str(frame.props[1].entr_mol_phase_comp[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_entr_mol_phase(self, frame):
        assert isinstance(frame.props[1].entr_mol_phase, Expression)
        assert len(frame.props[1].entr_mol_phase) == 2
        for k in frame.props[1].entr_mol_phase:
            assert k in frame.params.phase_list
            assert str(frame.props[1].entr_mol_phase[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_entr_mol(self, frame):
        assert isinstance(frame.props[1].entr_mol, Expression)
        assert len(frame.props[1].entr_mol) == 1
        assert str(frame.props[1].entr_mol.expr) == \
            str(sum(frame.props[1].entr_mol_phase[p] *
                    frame.props[1].phase_frac[p]
                    for p in frame.props[1].params.phase_list))

        # Check that dependency variables also constructed properly
        assert isinstance(frame.props[1].entr_mol_phase, Expression)
        assert len(frame.props[1].entr_mol_phase) == 2
        for p in frame.props[1].entr_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].entr_mol_phase[p].expr) == \
                str(frame.props[1].dummy_var)

    def test_fug_phase_comp(self, frame):
        assert isinstance(frame.props[1].fug_phase_comp, Expression)
        assert len(frame.props[1].fug_phase_comp) == 6
        for k in frame.props[1].fug_phase_comp:
            assert k[0] in frame.params.phase_list
            assert k[1] in frame.params.component_list
            assert str(frame.props[1].fug_phase_comp[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_fug_coeff_phase_comp(self, frame):
        assert isinstance(frame.props[1].fug_coeff_phase_comp, Expression)
        assert len(frame.props[1].fug_coeff_phase_comp) == 6
        for k in frame.props[1].fug_coeff_phase_comp:
            assert k[0] in frame.params.phase_list
            assert k[1] in frame.params.component_list
            assert str(frame.props[1].fug_coeff_phase_comp[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_gibbs_mol_phase_comp(self, frame):
        assert isinstance(frame.props[1].gibbs_mol_phase_comp, Expression)
        assert len(frame.props[1].gibbs_mol_phase_comp) == 6
        for k in frame.props[1].gibbs_mol_phase_comp:
            assert k[0] in frame.params.phase_list
            assert k[1] in frame.params.component_list
            assert str(frame.props[1].gibbs_mol_phase_comp[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_gibbs_mol_phase(self, frame):
        assert isinstance(frame.props[1].gibbs_mol_phase, Expression)
        assert len(frame.props[1].gibbs_mol_phase) == 2
        for k in frame.props[1].gibbs_mol_phase:
            assert k in frame.params.phase_list
            assert str(frame.props[1].gibbs_mol_phase[k].expr) == \
                str(frame.props[1].dummy_var)

    def test_gibbs_mol(self, frame):
        assert isinstance(frame.props[1].gibbs_mol, Expression)
        assert len(frame.props[1].gibbs_mol) == 1
        assert str(frame.props[1].gibbs_mol.expr) == \
            str(sum(frame.props[1].gibbs_mol_phase[p] *
                    frame.props[1].phase_frac[p]
                    for p in frame.props[1].params.phase_list))

        # Check that dependency variables also constructed properly
        assert isinstance(frame.props[1].gibbs_mol_phase, Expression)
        assert len(frame.props[1].gibbs_mol_phase) == 2
        for p in frame.props[1].gibbs_mol_phase:
            assert p in frame.params.phase_list
            assert str(frame.props[1].gibbs_mol_phase[p].expr) == \
                str(frame.props[1].dummy_var)

    def test_mw_phase(self, frame):
        frame.params.mw_comp = Param(frame.params.component_list,
                                     initialize={"a": 1, "b": 2, "c": 3})

        assert isinstance(frame.props[1].mw_phase, Expression)
        assert len(frame.props[1].mw_phase) == 2
        for k in frame.props[1].mw_phase:
            assert k in frame.params.phase_list
            assert str(frame.props[1].mw_phase[k].expr) == \
                str(sum(frame.props[1].mole_frac_phase_comp[k, j] *
                        frame.params.mw_comp[j]
                        for j in frame.params.component_list))

    def test_mw(self, frame):
        frame.params.mw_comp = Param(frame.params.component_list,
                                     initialize={"a": 1, "b": 2, "c": 3})

        assert isinstance(frame.props[1].mw, Expression)
        assert len(frame.props[1].mw) == 1
        assert str(frame.props[1].mw.expr) == \
            str(sum(frame.props[1].phase_frac[p] *
                    sum(frame.props[1].mole_frac_phase_comp[p, j] *
                        frame.params.mw_comp[j]
                        for j in frame.params.component_list)
                    for p in frame.params.phase_list))

    def test_pressure_sat_comp(self, frame):
        frame.params.config.pressure_sat_comp = modules[__name__]

        assert isinstance(frame.props[1].pressure_sat_comp, Expression)
        assert len(frame.props[1].pressure_sat_comp) == 3
        for j in frame.props[1].pressure_sat_comp:
            assert j in frame.params.component_list
            assert str(frame.props[1].pressure_sat_comp[j].expr) == \
                str(frame.props[1].pressure)
