##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for flowsheet_model.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Constraint, Var
from idaes.core import (declare_process_block_class, PropertyParameterBase,
                        StateBlockBase, StateBlockDataBase)


# -----------------------------------------------------------------------------
# Test ParameterBlock
@declare_process_block_class("ParameterBlock")
class _ParameterBlock(PropertyParameterBase):
    def build(self):
        pass


def test_PropertyParameterBase():
    # Test that PropertyParameterBase gets module information
    m = ConcreteModel()
    m.p = ParameterBlock()
    super(_ParameterBlock, m.p).build()

    assert hasattr(m.p, "property_module")


def test_NotImplementedErrors():
    # Test that class methods return NotImplementedError
    m = ConcreteModel()
    m.p = ParameterBlock()

    with pytest.raises(NotImplementedError):
        m.p.get_supported_properties()
    with pytest.raises(NotImplementedError):
        m.p.get_package_units()


# -----------------------------------------------------------------------------
# Test StateBlockBase
def test_StateBlockBase_initialize():
    # Test that StateBlockBase initialize method rasie NotImplementedError
    m = ConcreteModel()
    m.p = StateBlockData()

    with pytest.raises(NotImplementedError):
        m.p.initialize()


# -----------------------------------------------------------------------------
# Test StateBlockDataBase
@declare_process_block_class("BuildTest")
class _BuildTest(StateBlockDataBase):
    def build(self):
        super(_BuildTest, self).build()


def test_build_NotImplementedError():
    # Test that StateBlockDataBase build method returns NotImplementedErrors
    m = ConcreteModel()

    with pytest.raises(NotImplementedError):
        m.p = BuildTest()


@declare_process_block_class("StateBlockData", block_class=StateBlockBase)
class _StateBlockData(StateBlockDataBase):
    def build(self):
        pass

def test_StateBlock_config():
    # Test that StateBlockDataBase config has correct arguments
    m = ConcreteModel()
    m.p = StateBlockData()

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


def test_NotImplementedErrors():
    # Test that placeholder methods return NotImplementedErrors
    m = ConcreteModel()
    m.p = StateBlockData()

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
        m.p.get_enthlpy_flow_terms()
    with pytest.raises(NotImplementedError):
        m.p.get_enthlpy_density_terms()
    with pytest.raises(NotImplementedError):
        m.p.get_energy_diffusion_terms()


# -----------------------------------------------------------------------------
# Test properties __getattr__ method
@declare_process_block_class("Parameters")
class _Parameters(PropertyParameterBase):
    def build(self):
        super(_Parameters, self).build()

    @classmethod
    def get_supported_properties(self):
        return {'a': {'method': 'a_method', 'units': 'mol/s'},
                'recursion1': {'method': '_recursion1', 'units': 'Pa'},
                'recursion2': {'method': '_recursion2', 'units': 'K'},
                'raise_exception': {'method': '_raise_exception',
                                    'units': 'K'},
                'does_not_create_component': {
                        'method': '_does_not_create_component',
                        'units': 'K'}}


@declare_process_block_class("State", block_class=StateBlockBase)
class _State(StateBlockDataBase):
    def build(self):
        self.not_callable = 1

    def a_method(self):
        self.a = Var(initialize=1)


def test_getattr_add_var():
    m = ConcreteModel()
    m.pb = Parameters()
    m.p = State(parameters=m.pb)

    assert isinstance(m.p.a, Var)
    assert m.p.a.value == 1


#def test_getattr_protected():
#    m = ConcreteModel()
#    m.p = StateBlockData()
#
#    with pytest.raises(AttributeError):
#        # Call a protected component that does not exist
#        m.p.simple_cons = Constraint(expr=m.p._foo == 2)
#
#
#def test_getattr_recursion():
#    m = ConcreteModel()
#    m.p = StateBlockData()
#
#    with pytest.raises(Exception):
#        # Call a component that triggers a recursive loop of calls
#        m.p.simple_cons = Constraint(expr=m.p.recursion1 == 2)
#
#
#def test_getattr_does_not_exist():
#    m = ConcreteModel()
#    m.p = PropertyBlockData()
#    with pytest.raises(AttributeError):
#        m.p.simple_cons = Constraint(expr=m.p.does_not_exist == 2)
#
#
#def test_getattr_not_callable():
#    m = ConcreteModel()
#    m.p = PropertyBlockData()
#    with pytest.raises(AttributeError):
#        m.p.simple_cons = Constraint(expr=m.p.not_callable == 2)
#
#
#def test_getattr_raise_exception():
#    m = ConcreteModel()
#    m.pb = ParameterBlock()
#    super(_ParameterBlock, m.pb).build()
#    m.p = PropertyBlockData(parameters=m.pb)
#    with pytest.raises(KeyError):
#        m.p.simple_cons = Constraint(expr=m.p.raise_exception == 2)
#
#
#def test_getattr_does_not_create_component():
#    m = ConcreteModel()
#    m.p = PropertyBlockData()
#    with pytest.raises(Exception):
#        m.p.simple_cons = Constraint(expr=m.p.does_not_create_component == 2)

