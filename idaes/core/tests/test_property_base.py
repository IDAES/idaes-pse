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
Tests for flowsheet_model.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Constraint, Set, Var
from pyomo.common.config import ConfigBlock, ConfigValue
from idaes.core import (declare_process_block_class, PhysicalParameterBlock,
                        StateBlock, StateBlockData)
from idaes.core.util.exceptions import (PropertyPackageError,
                                        PropertyNotSupportedError)

# -----------------------------------------------------------------------------
# Test ParameterBlock
@declare_process_block_class("ParameterBlock")
class _ParameterBlock(PhysicalParameterBlock):
    pass


def test_config_block():
    # Test that PhysicalParameterBlock gets module information
    m = ConcreteModel()
    m.p = ParameterBlock()

    assert len(m.p.config) == 1
    assert isinstance(m.p.config.default_arguments, ConfigBlock)
    assert len(m.p.config.default_arguments) == 0


def test_PhysicalParameterBlock():
    # Test that PhysicalParameterBlock builds correctly
    m = ConcreteModel()
    m.p = ParameterBlock()
    super(_ParameterBlock, m.p).build()


def test_PhysicalParameter_NotImplementedErrors():
    # Test that class methods return NotImplementedError
    m = ConcreteModel()
    m.p = ParameterBlock()

    with pytest.raises(NotImplementedError):
        m.p.get_metadata()


def test_get_phase_component_set():
    m = ConcreteModel()
    m.p = ParameterBlock() 

    m.p.phase_list = Set(initialize=["1", "2", "3"])
    m.p.component_list = Set(initialize=["a", "b", "c"])

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


def test_get_phase_component_set_subset():
    m = ConcreteModel()
    m.p = ParameterBlock() 

    m.p.config = ConfigBlock()
    m.p.config.declare("phase_component_list",
                       ConfigValue(default={"1": ["a", "b", "c"],
                                            "2": ["a", "b"],
                                            "3": ["c"]}))

    pc_set = m.p.get_phase_component_set()

    assert isinstance(m.p._phase_component_set, Set)
    assert len(m.p._phase_component_set) == 6
    for v in m.p._phase_component_set:
        assert v[0] in m.p.config.phase_component_list.keys()
        assert v[1] in m.p.config.phase_component_list[v[0]]

    assert pc_set is m.p._phase_component_set


# -----------------------------------------------------------------------------
# Test StateBlock
@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class _StateBlockData(StateBlockData):
    pass


@declare_process_block_class("TestStateBlock2", block_class=StateBlock)
class _StateBlockData(StateBlockData):
    def define_state_vars(self):
        return {}


def test_StateBlockBase_initialize():
    # Test that StateBlock initialize method raises NotImplementedError
    m = ConcreteModel()
    # Need to index block so that it does not do multiple inheritance
    m.s = Set(initialize=[1, 2])
    m.p = TestStateBlock2(m.s)

    with pytest.raises(NotImplementedError):
        m.p.initialize()


def test_StateBlockBase_report():
    # Test that StateBlock initialize method raises NotImplementedError
    m = ConcreteModel()
    # Need to index block so that it does not do multiple inheritance
    m.s = Set(initialize=[0, 1])
    m.p = TestStateBlock2(m.s)

    m.p.report(dof=True)


# -----------------------------------------------------------------------------
# Test StateBlockData
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


def test_param_ref():
    m = ConcreteModel()
    m.pb = Parameters()
    m.p = StateTest(default={"parameters": m.pb})

    assert m.p.params == m.p.config.parameters


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


def test_getattr_add_var(m):
    assert isinstance(m.p.a, Var)
    assert m.p.a.value == 1


def test_getattr_protected(m):
    with pytest.raises(PropertyPackageError):
        # Call a protected component that does not exist
        m.p.cons = Constraint(expr=m.p._foo == 1)


def test_getattr_recursion(m):
    with pytest.raises(PropertyPackageError):
        # Call a component that triggers a recursive loop of calls
        m.p.cons = Constraint(expr=m.p.recursion1 == 1)


def test_getattr_does_not_exist(m):
    with pytest.raises(PropertyNotSupportedError):
        m.p.cons = Constraint(expr=m.p.does_not_exist == 1)


def test_getattr_not_callable(m):
    with pytest.raises(PropertyPackageError):
        m.p.cons = Constraint(expr=m.p.not_callable == 1)


def test_getattr_not_supported(m):
    with pytest.raises(PropertyNotSupportedError):
        m.p.cons = Constraint(expr=m.p.not_supported == 1)


def test_getattr_raise_exception(m):
    with pytest.raises(Exception):
        m.p.cons = Constraint(expr=m.p.raise_exception == 1)


# TODO : Need a test for cases where method does not create property
#def test_getattr_does_not_create_component(m):
#    with pytest.raises(PropertyPackageError):
#        m.p.cons = Constraint(expr=m.p.does_not_create_component == 1)
