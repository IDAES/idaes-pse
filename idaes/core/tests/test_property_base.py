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
Tests for flowsheet_model.

Author: Andrew Lee
"""
import pytest
import types

from pyomo.environ import ConcreteModel, Constraint, Set, Var
from pyomo.common.config import ConfigBlock

from idaes.core import (declare_process_block_class, PhysicalParameterBlock,
                        StateBlock, StateBlockData)
from idaes.core.phases import Phase
from idaes.core.components import Component
from idaes.core.util.exceptions import (PropertyPackageError,
                                        PropertyNotSupportedError)
from idaes.core.property_meta import PropertyClassMetadata


# -----------------------------------------------------------------------------
# Test ParameterBlock
@declare_process_block_class("ParameterBlock")
class _ParameterBlock(PhysicalParameterBlock):
    pass


@pytest.mark.unit
def test_config_block():
    # Test that PhysicalParameterBlock gets module information
    m = ConcreteModel()
    m.p = ParameterBlock()

    assert len(m.p.config) == 1
    assert isinstance(m.p.config.default_arguments, ConfigBlock)
    assert len(m.p.config.default_arguments) == 0


@pytest.mark.unit
def test_PhysicalParameterBlock():
    # Test that PhysicalParameterBlock builds correctly
    m = ConcreteModel()
    m.p = ParameterBlock()
    super(_ParameterBlock, m.p).build()


@pytest.mark.unit
def test_PhysicalParameter_NotImplementedErrors():
    # Test that class methods return NotImplementedError
    m = ConcreteModel()
    m.p = ParameterBlock()

    with pytest.raises(NotImplementedError):
        m.p.get_metadata()


@pytest.mark.unit
def test_get_phase_component_set():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.p.phase_list = Set(initialize=["1", "2", "3"])
    m.p.component_list = Set(initialize=["a", "b", "c"])

    # Need to call this to trigger build of Phase objects for now
    m.p._make_phase_objects()

    pc_set = m.p.get_phase_component_set()

    assert isinstance(m.p._phase_component_set, Set)
    assert len(m.p._phase_component_set) == 9
    for v in m.p._phase_component_set:
        assert v[0] in m.p.phase_list
        assert v[1] in m.p.component_list

    assert pc_set is m.p._phase_component_set

    # Check that method returns existing component
    # Delete phase list so that build will fail to make sure it isn't rebuilding
    m.p.del_component(m.p.phase_list)

    assert m.p.get_phase_component_set() is m.p._phase_component_set


@pytest.mark.unit
def test_get_phase_component_set_subset():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.p.phase_list = ["1", "2", "3"]
    m.p.phase_comp = {"1": ["a", "b", "c"],
                      "2": ["a", "b"],
                      "3": ["c"]}

    # Need to call this to trigger build of Phase objects for now
    m.p._make_phase_objects()

    pc_set = m.p.get_phase_component_set()

    assert isinstance(m.p._phase_component_set, Set)
    assert len(m.p._phase_component_set) == 6
    for v in m.p._phase_component_set:
        assert v[0] in m.p.phase_comp.keys()
        assert v[1] in m.p.phase_comp[v[0]]

    assert pc_set is m.p._phase_component_set


@pytest.mark.unit
def test_get_component():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.meta_object = PropertyClassMetadata()

    def get_metadata(self):
        return m.meta_object
    m.p.get_metadata = types.MethodType(get_metadata, m.p)

    with pytest.raises(AttributeError):
        m.p.get_component("foo")

    m.p.comp = Component()

    assert m.p.get_component("comp") is m.p.comp

    m.p.a = object()

    with pytest.raises(
            PropertyPackageError,
            match="p get_component found an attribute a, but it does not "
            "appear to be an instance of a Component object."):
        m.p.get_component("a")


@pytest.mark.unit
def test_get_phase():
    m = ConcreteModel()
    m.p = ParameterBlock()

    with pytest.raises(AttributeError):
        m.p.get_phase("foo")

    m.p.phase = Phase()

    assert m.p.get_phase("phase") is m.p.phase

    m.p.a = object()

    with pytest.raises(
            PropertyPackageError,
            match="p get_phase found an attribute a, but it does not "
            "appear to be an instance of a Phase object."):
        m.p.get_phase("a")


@pytest.mark.unit
def test_make_component_objects():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.meta_object = PropertyClassMetadata()

    def get_metadata(self):
        return m.meta_object
    m.p.get_metadata = types.MethodType(get_metadata, m.p)

    m.p.component_list = Set(initialize=["comp1", "comp2"])

    m.p._make_component_objects()

    assert isinstance(m.p.comp1, Component)
    assert isinstance(m.p.comp2, Component)


@pytest.mark.unit
def test_make_component_objects_already_exists():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.p.comp1 = object()

    m.p.component_list = Set(initialize=["comp1", "comp2"])

    with pytest.raises(PropertyPackageError,
                       match="p could not add Component object comp1 - an "
                       "object with that name already exists."):
        m.p._make_component_objects()


@pytest.mark.unit
def test_make_phase_objects():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.p.phase_list = Set(initialize=["phase1", "phase2"])
    m.p.phase_comp = {"phase1": [1, 2, 3],
                      "phase2": [4, 5, 6]}

    m.p._make_phase_objects()

    assert isinstance(m.p.phase1, Phase)
    assert isinstance(m.p.phase2, Phase)

    assert m.p.phase1.config.component_list == [1, 2, 3]
    assert m.p.phase2.config.component_list == [4, 5, 6]


@pytest.mark.unit
def test_make_phase_objects_already_exists():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.p.phase1 = object()

    m.p.phase_list = Set(initialize=["phase1", "phase2"])

    with pytest.raises(PropertyPackageError,
                       match="p could not add Phase object phase1 - an "
                       "object with that name already exists."):
        m.p._make_phase_objects()


@pytest.mark.unit
def test_validate_parameter_block_no_component_list():
    m = ConcreteModel()
    m.p = ParameterBlock()

    with pytest.raises(
            PropertyPackageError,
            match="Property package p has not defined a component list."):
        m.p._validate_parameter_block()


@pytest.mark.unit
def test_validate_parameter_block_no_phase_list():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.meta_object = PropertyClassMetadata()

    def get_metadata(self):
        return m.meta_object
    m.p.get_metadata = types.MethodType(get_metadata, m.p)

    m.p.component_list = Set(initialize=["c1", "c2"])

    with pytest.raises(
            PropertyPackageError,
            match="Property package p has not defined a phase list."):
        m.p._validate_parameter_block()


@pytest.mark.unit
def test_validate_parameter_block_invalid_component_object():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.p.component_list = Set(initialize=["foo"])

    m.p.phase_list = Set(initialize=["p1", "p2"])
    m.p.foo = object()

    with pytest.raises(TypeError):
        m.p._validate_parameter_block()


@pytest.mark.unit
def test_validate_parameter_block_invalid_phase_object():
    m = ConcreteModel()
    m.p = ParameterBlock()

    m.meta_object = PropertyClassMetadata()

    def get_metadata(self):
        return m.meta_object
    m.p.get_metadata = types.MethodType(get_metadata, m.p)

    m.p.component_list = Set(initialize=["c1", "c2"])

    m.p.phase_list = Set(initialize=["foo"])
    m.p.foo = object()

    with pytest.raises(TypeError):
        m.p._validate_parameter_block()


@pytest.mark.unit
def test_has_inherent_Reactions():
    m = ConcreteModel()
    m.p = ParameterBlock()

    assert not m.p.has_inherent_reactions

    m.p._has_inherent_reactions = True

    assert m.p.has_inherent_reactions


# -----------------------------------------------------------------------------
# Test StateBlock
@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class _StateBlockData(StateBlockData):
    def build(self):
        super(StateBlockData, self).build()


@declare_process_block_class("TestStateBlock2", block_class=StateBlock)
class _StateBlockData2(StateBlockData):
    def build(self):
        super(StateBlockData, self).build()

    def define_state_vars(self):
        return {}


@pytest.mark.unit
def test_StateBlockBase_initialize():
    # Test that StateBlock initialize method raises NotImplementedError
    m = ConcreteModel()
    # Need to index block so that it does not do multiple inheritance
    m.s = Set(initialize=[1, 2])
    m.p = TestStateBlock2(m.s)

    with pytest.raises(NotImplementedError):
        m.p.initialize()


@pytest.mark.unit
def test_StateBlockBase_report():
    # Test that StateBlock initialize method raises NotImplementedError
    m = ConcreteModel()
    # Need to index block so that it does not do multiple inheritance
    m.s = Set(initialize=[0, 1])
    m.p = TestStateBlock2(m.s)

    m.p.report(dof=True)


# -----------------------------------------------------------------------------
# Test StateBlockData
@pytest.mark.unit
def test_StateBlock_config():
    # Test that StateBlockData config has correct arguments
    m = ConcreteModel()
    m.p = TestStateBlock()

    assert len(m.p.config) == 3
    assert hasattr(m.p.config, "has_phase_equilibrium")
    assert hasattr(m.p.config, "defined_state")
    assert hasattr(m.p.config, "parameters")

    m.p.config.has_phase_equilibrium = True
    m.p.config.has_phase_equilibrium = False
    with pytest.raises(ValueError):
        m.p.config.has_phase_equilibrium = 'foo'
    with pytest.raises(ValueError):
        m.p.config.has_phase_equilibrium = 10

    m.p.config.defined_state = True
    m.p.config.defined_state = False
    with pytest.raises(ValueError):
        m.p.config.defined_state = 'foo'
    with pytest.raises(ValueError):
        m.p.config.defined_state = 10


@pytest.mark.unit
def test_StateBlock_NotImplementedErrors():
    # Test that placeholder methods return NotImplementedErrors
    m = ConcreteModel()
    m.p = TestStateBlock()

    with pytest.raises(NotImplementedError):
        m.p.define_state_vars()
    with pytest.raises(NotImplementedError):
        m.p.define_port_members()
    with pytest.raises(NotImplementedError):
        m.p.get_material_flow_terms()
    with pytest.raises(NotImplementedError):
        m.p.get_material_density_terms()
    with pytest.raises(NotImplementedError):
        m.p.get_material_diffusion_terms()
    with pytest.raises(NotImplementedError):
        m.p.get_enthalpy_flow_terms()
    with pytest.raises(NotImplementedError):
        m.p.get_energy_density_terms()
    with pytest.raises(NotImplementedError):
        m.p.get_energy_diffusion_terms()
    with pytest.raises(NotImplementedError):
        m.p.calculate_bubble_point_temperature()
    with pytest.raises(NotImplementedError):
        m.p.calculate_dew_point_temperature()
    with pytest.raises(NotImplementedError):
        m.p.calculate_bubble_point_pressure()
    with pytest.raises(NotImplementedError):
        m.p.calculate_dew_point_pressure()


# -----------------------------------------------------------------------------
# Test parameter block reference attribute
@declare_process_block_class("Parameters")
class _Parameters(PhysicalParameterBlock):
    def build(self):
        super(_Parameters, self).build()

        self.phase_list = ["p1"]
        self.component_list = ["c1"]

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({'a': {'method': 'a_method'},
                            'recursion1': {'method': '_recursion1'},
                            'recursion2': {'method': '_recursion2'},
                            'not_callable': {'method': 'test_obj'},
                            'raise_exception': {'method': '_raise_exception'},
                            'not_supported': {'method': False},
                            'does_not_create_component': {
                                'method': '_does_not_create_component'}})


@declare_process_block_class("StateTest", block_class=StateBlock)
class _StateTest(StateBlockData):
    def build(self):
        super(_StateTest, self).build()


@pytest.mark.unit
def test_param_ref():
    m = ConcreteModel()
    m.pb = Parameters()
    m.p = StateTest(default={"parameters": m.pb})

    assert m.p.params == m.p.config.parameters


@pytest.mark.unit
def test_validate_params():
    # Test that validate params has been triggered
    m = ConcreteModel()
    m.pb = Parameters()
    m.p = StateTest(default={"parameters": m.pb})

    # If validation has been triggered, Phase & Component objects should exist
    assert isinstance(m.pb.p1, Phase)
    assert isinstance(m.pb.c1, Component)


@pytest.mark.unit
def test_has_inherent_reactions_state_block():
    m = ConcreteModel()
    m.pb = Parameters()
    m.p = StateTest(default={"parameters": m.pb})

    assert not m.p.has_inherent_reactions

    m.pb._has_inherent_reactions = True

    assert m.p.has_inherent_reactions


# -----------------------------------------------------------------------------
# Test properties __getattr__ method
@declare_process_block_class("State", block_class=StateBlock)
class _State(StateBlockData):
    def build(self):
        super(StateBlockData, self).build()

        self.test_obj = 1

    def a_method(self):
        self.a = Var(initialize=1)

    def _recursion1(self):
        self.recursive_cons1 = Constraint(expr=self.recursion2 == 1)

    def _recursion2(self):
        self.recursive_cons2 = Constraint(expr=self.recursion1 == 1)

    def _raise_exception(self):
        raise Exception()

    def _does_not_create_component(self):
        pass


@pytest.fixture()
def m():
    m = ConcreteModel()
    m.pb = Parameters()
    m.p = State(default={"parameters": m.pb})

    return m


@pytest.mark.unit
def test_getattr_add_var(m):
    assert isinstance(m.p.a, Var)
    assert m.p.a.value == 1
    assert m.p.is_property_constructed("a")


@pytest.mark.unit
def test_lock_attributes(m):
    with pytest.raises(AttributeError):
        with m.p.lock_attribute_creation_context():
            m.p.a.value == 1
    # Make sure it unlocked
    assert m.p.a.value == 1


@pytest.mark.unit
def test_is_property_constructed(m):
    assert not m.p.is_property_constructed("a")
    assert m.p.a.value == 1
    assert m.p.is_property_constructed("a")


@pytest.mark.unit
def test_getattr_protected(m):
    with pytest.raises(PropertyPackageError):
        # Call a protected component that does not exist
        m.p.cons = Constraint(expr=m.p._foo == 1)


@pytest.mark.unit
def test_getattr_recursion(m):
    with pytest.raises(PropertyPackageError):
        # Call a component that triggers a recursive loop of calls
        m.p.cons = Constraint(expr=m.p.recursion1 == 1)


@pytest.mark.unit
def test_getattr_does_not_exist(m):
    with pytest.raises(PropertyNotSupportedError):
        m.p.cons = Constraint(expr=m.p.does_not_exist == 1)


@pytest.mark.unit
def test_getattr_not_callable(m):
    with pytest.raises(PropertyPackageError):
        m.p.cons = Constraint(expr=m.p.not_callable == 1)


@pytest.mark.unit
def test_getattr_not_supported(m):
    with pytest.raises(PropertyNotSupportedError):
        m.p.cons = Constraint(expr=m.p.not_supported == 1)


@pytest.mark.unit
def test_getattr_raise_exception(m):
    with pytest.raises(Exception):
        m.p.cons = Constraint(expr=m.p.raise_exception == 1)


# TODO : Need a test for cases where method does not create property
#@pytest.mark.unit
#def test_getattr_does_not_create_component(m):
#    with pytest.raises(PropertyPackageError):
#        m.p.cons = Constraint(expr=m.p.does_not_create_component == 1)
