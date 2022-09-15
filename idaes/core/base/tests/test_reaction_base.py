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
import inspect
from pyomo.environ import ConcreteModel, Constraint, Set, Var, units as pyunits
from pyomo.common.config import ConfigBlock
from idaes.core import (
    declare_process_block_class,
    ReactionParameterBlock,
    ReactionBlockBase,
    ReactionBlockDataBase,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
)
from idaes.core.util.exceptions import PropertyPackageError, PropertyNotSupportedError


# -----------------------------------------------------------------------------
# Test ParameterBlock
@declare_process_block_class("PropertyParameterBlock")
class _PropertyParameterBlock(PhysicalParameterBlock):
    def build(self):
        super(_PropertyParameterBlock, self).build()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {"prop1": {"method": None, "units": "m"}, "prop3": {"method": False}}
        )
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.g,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


@declare_process_block_class("ReactionParameterTestBlock")
class _ReactionParameterBlock(ReactionParameterBlock):
    def build(self):
        super(ReactionParameterBlock, self).build()


@pytest.mark.unit
def test_config_block():
    # Test that PhysicalParameterBlock gets module information
    m = ConcreteModel()
    m.r = ReactionParameterTestBlock()

    assert len(m.r.config) == 2
    assert isinstance(m.r.config.default_arguments, ConfigBlock)
    assert len(m.r.config.default_arguments) == 0


@pytest.mark.unit
def test_ReactionParameter_NotImplementedErrors():
    # Test that class methods return NotImplementedError
    m = ConcreteModel()
    m.r = ReactionParameterTestBlock()

    with pytest.raises(NotImplementedError):
        m.r.get_metadata()


@declare_process_block_class("ReactionParameterBlock2")
class _ReactionParameterBlock2(ReactionParameterBlock):
    def build(self):
        super(ReactionParameterBlock, self).build()

    @classmethod
    def get_required_properties(self):
        return ["prop1"]

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({"rxn1": {"method": None}})
        obj.add_default_units(
            {
                "time": pyunits.hr,
                "length": pyunits.m,
                "mass": pyunits.g,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )


@pytest.mark.unit
def test_validate_state_block_invalid_units():
    # Test validation of associated PropertyParameterBlock
    m = ConcreteModel()
    m.p = PropertyParameterBlock()
    m.r = ReactionParameterBlock2(property_package=m.p)

    with pytest.raises(PropertyPackageError):
        m.r._validate_property_parameter_units()


@declare_process_block_class("ReactionParameterBlock3")
class _ReactionParameterBlock3(ReactionParameterBlock):
    def build(self):
        super(ReactionParameterBlock, self).build()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({"rxn1": {"method": None}})
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.g,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )
        obj.add_required_properties({"prop2": "some"})


@pytest.mark.unit
def test_validate_state_block_unsupported_prop():
    # Test validation of associated PropertyParameterBlock
    m = ConcreteModel()
    m.p = PropertyParameterBlock()
    m.r = ReactionParameterBlock3(property_package=m.p)

    with pytest.raises(PropertyPackageError):
        m.r._validate_property_parameter_properties()


@declare_process_block_class("ReactionParameterBlock4")
class _ReactionParameterBlock4(ReactionParameterBlock):
    def build(self):
        super(ReactionParameterBlock, self).build()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({"rxn1": {"method": None}})
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.g,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )
        obj.add_required_properties({"prop3": "some"})


@pytest.mark.unit
def test_validate_state_block_unsupported_prop_False():
    # Test validation of associated PropertyParameterBlock
    m = ConcreteModel()
    m.p = PropertyParameterBlock()
    m.r = ReactionParameterBlock4(property_package=m.p)

    with pytest.raises(PropertyPackageError):
        m.r._validate_property_parameter_properties()


@declare_process_block_class("ReactionParameterBlock5")
class _ReactionParameterBlock5(ReactionParameterBlock):
    def build(self):
        super(ReactionParameterBlock, self).build()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({"rxn1": {"method": None}})
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.g,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )
        obj.add_required_properties({"prop1": "km"})


@pytest.mark.unit
def test_validate_state_block_req_prop_wrong_units():
    # Test validation of associated PropertyParameterBlock
    m = ConcreteModel()
    m.p = PropertyParameterBlock()
    m.r = ReactionParameterBlock5(property_package=m.p)

    with pytest.raises(PropertyPackageError):
        m.r._validate_property_parameter_properties()


@declare_process_block_class("ReactionParameterBlock6")
class _ReactionParameterBlock6(ReactionParameterBlock):
    def build(self):
        super(ReactionParameterBlock, self).build()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties({"rxn1": {"method": None}})
        obj.add_default_units(
            {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.g,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            }
        )
        obj.add_required_properties({"prop1": "m"})


@pytest.mark.unit
def test_ReactionParameterBase_build():
    # Test that ReactionParameterBlock builds correctly
    m = ConcreteModel()
    m.p = PropertyParameterBlock()
    m.r = ReactionParameterBlock6(property_package=m.p)
    super(_ReactionParameterBlock6, m.r).build()


# -----------------------------------------------------------------------------
# Test ReactionBlockBase
@declare_process_block_class("ReactionBlock", block_class=ReactionBlockBase)
class ReactionBlockData(ReactionBlockDataBase):
    def build(self):
        super(ReactionBlockDataBase, self).build()


@pytest.mark.unit
def test_ReactionBlockBase_initialize():
    # Test that ReactionBlockBase initialize method raises NotImplementedError
    m = ConcreteModel()
    # Need to index block so that it does not do multiple inheritance
    m.s = Set(initialize=[1, 2])
    m.r = ReactionBlock(m.s)

    with pytest.raises(NotImplementedError):
        m.r.initialize()


@pytest.mark.unit
def test_ReactionBlockBase_report():
    # Test that ReactionBlockBase report method raises NotImplementedError
    m = ConcreteModel()
    # Need to index block so that it does not do multiple inheritance
    m.s = Set(initialize=[1, 2])
    m.r = ReactionBlock(m.s)

    with pytest.raises(NotImplementedError):
        m.r.report()


# -----------------------------------------------------------------------------
# Test ReactionBlockDataBase
@pytest.mark.unit
def test_StateBlock_config():
    # Test that ReactionBlockDataBase config has correct arguments
    m = ConcreteModel()
    m.p = ReactionBlock()

    assert len(m.p.config) == 3
    assert hasattr(m.p.config, "parameters")
    assert hasattr(m.p.config, "state_block")
    assert hasattr(m.p.config, "has_equilibrium")

    m.p.config.has_equilibrium = True
    m.p.config.has_equilibrium = False
    with pytest.raises(ValueError):
        m.p.config.has_equilibrium = "foo"
    with pytest.raises(ValueError):
        m.p.config.has_equilibrium = 10


@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class StateTestBlockData(StateBlockData):
    def build(self):
        super(StateBlockData, self).build()


@pytest.mark.unit
def test_validate_state_block_fail():
    # Test that ReactionBlockDataBase validate_state_block returns error
    m = ConcreteModel()
    m.p = PropertyParameterBlock()
    m.p2 = PropertyParameterBlock()

    m.pb = TestStateBlock(parameters=m.p2)

    m.r = ReactionParameterBlock6(property_package=m.p)
    super(_ReactionParameterBlock6, m.r).build()

    m.rb = ReactionBlock(parameters=m.r, state_block=m.pb)

    with pytest.raises(PropertyPackageError):
        m.rb._validate_state_block()


@declare_process_block_class("ReactionBlock2", block_class=ReactionBlockBase)
class ReactionBlockData2(ReactionBlockDataBase):
    def build(self):
        super(ReactionBlockData2, self).build()


@pytest.mark.unit
def test_build():
    # Test that ReactionBlockDataBase builds correctly with good argumnets
    m = ConcreteModel()
    m.p = PropertyParameterBlock()

    m.pb = TestStateBlock(parameters=m.p)

    m.r = ReactionParameterBlock6(property_package=m.p)
    super(_ReactionParameterBlock6, m.r).build()

    m.rb = ReactionBlock2(parameters=m.r, state_block=m.pb)

    assert hasattr(m.rb, "state_ref")
    assert m.rb.params == m.rb.config.parameters


# -----------------------------------------------------------------------------
# Test reaction __getattr__ method
@declare_process_block_class("Parameters")
class _Parameters(ReactionParameterBlock):
    def build(self):
        super(ReactionParameterBlock, self).build()
        frm = inspect.stack()[1]
        self.__package_module = inspect.getmodule(frm[0])

    @classmethod
    def define_metadata(cls, obj):
        obj.add_properties(
            {
                "a": {"method": "a_method"},
                "recursion1": {"method": "_recursion1"},
                "recursion2": {"method": "_recursion2"},
                "not_callable": {"method": "test_obj"},
                "raise_exception": {"method": "_raise_exception"},
                "not_supported": {"method": False},
                "does_not_create_component": {"method": "_does_not_create_component"},
            }
        )


@declare_process_block_class("Reaction", block_class=ReactionBlockBase)
class _Reaction(ReactionBlockDataBase):
    def build(self):
        super(ReactionBlockDataBase, self).build()

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
    m.p = Reaction(parameters=m.pb)

    return m


@pytest.mark.unit
def test_getattr_add_var(m):
    assert isinstance(m.p.a, Var)
    assert m.p.a.value == 1
    assert m.p.is_property_constructed("a") == True


@pytest.mark.unit
def test_lock_attributes(m):
    with pytest.raises(AttributeError):
        with m.p.lock_attribute_creation_context():
            m.p.a.value == 1
    # Make sure it unlocked
    assert m.p.a.value == 1


@pytest.mark.unit
def test_is_property_constructed(m):
    assert m.p.is_property_constructed("a") == False
    assert m.p.a.value == 1
    assert m.p.is_property_constructed("a") == True


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
# @pytest.mark.unit
# def test_getattr_does_not_create_component(m):
#    with pytest.raises(PropertyPackageError):
#        m.p.cons = Constraint(expr=m.p.does_not_create_component == 1)
