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
Tests for ControlVolumeBase.

Author: Andrew Lee
"""
import inspect
import pytest
from pyomo.environ import ConcreteModel, Block, Set
from pyomo.common.config import ConfigBlock, ConfigValue
from idaes.core import (ControlVolumeBase, declare_process_block_class,
                        FlowsheetBlockData, UnitBlockData, useDefault,
                        PropertyParameterBase, ReactionParameterBase)
from idaes.core.util.exceptions import ConfigurationError, DynamicError

# TODO : Test CONFIG_Base, Enums


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("Unit")
class _UnitData(UnitBlockData):
    CONFIG = UnitBlockData.CONFIG()
    CONFIG.declare("property_package",
                   ConfigValue(default=None))
    CONFIG.declare("property_package_args",
                   ConfigValue(default={}))

    def build(self):
        super(_UnitData, self).build()


# -----------------------------------------------------------------------------
# Testing ControlVolumeBase
@declare_process_block_class("CVFrame")
class CVFrameData(ControlVolumeBase):
    def build(self):
        pass


def test_config_block():
    cv = CVFrame()

    assert len(cv.config) == 5
    assert cv.config.dynamic == useDefault
    assert cv.config.property_package == useDefault
    assert isinstance(cv.config.property_package_args, ConfigBlock)
    assert len(cv.config.property_package_args) == 0
    assert cv.config.reaction_package is None
    assert isinstance(cv.config.reaction_package_args, ConfigBlock)
    assert len(cv.config.reaction_package_args) == 0


# -----------------------------------------------------------------------------
# Test _setup_dynamics
def test_setup_dynamics_use_parent_value():
    # Test that dynamic = None works correctly
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.u = Unit(dynamic=False)
    m.fs.u.cv = CVFrame()

    m.fs.u.cv._setup_dynamics()

    assert m.fs.u.cv.config.dynamic is False
    assert m.fs.u.cv.time == [0]


def test_setup_dynamics_use_parent_value_fail_no_dynamic():
    # Test that dynamic = None works correctly
    fs = Flowsheet(dynamic=False, concrete=True)

    # Create a Block (with no dynamic attribute)
    fs.b = Block()
    # Add a time attribute to make sure the correct failure triggers
    fs.b.time = Set(initialize=[0])

    fs.b.cv = CVFrame()

    # _setup_dynamics should return DynamicError
    with pytest.raises(DynamicError):
        fs.b.cv._setup_dynamics()


def test_setup_dynamics_use_parent_value_fail_no_time():
    # Test that methods fails when no time domain in parent
    fs = Flowsheet(dynamic=False, concrete=True)

    # Mock-up an object with a dynamic flag but no time domain (using CVFrame)
    fs.u = CVFrame()
    # Set the config flag
    fs.u.config.dynamic = False

    # Add a second holdup block
    fs.u.cv = CVFrame()

    # _setup_dynamics should return DynamicError, as fs.u has no time
    with pytest.raises(DynamicError):
        fs.u.cv._setup_dynamics()


# -----------------------------------------------------------------------------
# Test _get_property_package
@declare_process_block_class("PropertyParameterBlock")
class _PropertyParameterBlock(PropertyParameterBase):
    def build(self):
        frm = inspect.stack()[1]
        self.property_module = inspect.getmodule(frm[0])

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])


def test_get_property_package_set():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock()
    m.cv = CVFrame(property_package=m.pp)
    m.cv._get_property_package()

    assert m.cv.property_module == m.pp.property_module


def test_get_property_package_default_args():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock(default_arguments={"test": "foo"})
    m.cv = CVFrame(property_package=m.pp)
    m.cv._get_property_package()

    assert m.cv.config.property_package_args["test"] == "foo"


def test_get_reaction_package_module_combine_args():
    # Test that local and default args combine correctly
    m = ConcreteModel()
    m.pp = PropertyParameterBlock(default_arguments={"test1": "foo",
                                                     "test2": "bar"})
    m.cv = CVFrame(property_package=m.pp,
                   property_package_args={"test2": "baz",
                                          "test3": "bar"})
    m.cv._get_property_package()

    assert m.cv.config.property_package_args["test1"] == "foo"
    assert m.cv.config.property_package_args["test2"] == "baz"
    assert m.cv.config.property_package_args["test3"] == "bar"


# -----------------------------------------------------------------------------
# Test _get_default_prop_pack
def test_get_default_prop_pack_works():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.config.default_property_package = m.fs.pp

    m.fs.cv = CVFrame()
    assert m.fs.cv._get_default_prop_pack() == m.fs.pp


# TODO : should test more failure modes
def test_get_default_prop_pack_no_default():
    m = ConcreteModel()
    m.fs = Flowsheet()

    m.fs.cv = CVFrame()
    with pytest.raises(ConfigurationError):
        m.fs.cv._get_default_prop_pack()


def test_get_property_package_call_to_get_default_prop_pack():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PropertyParameterBlock()
    m.fs.config.default_property_package = m.fs.pp

    m.fs.cv = CVFrame()
    m.fs.cv._get_property_package()
    assert m.fs.cv.property_module == m.fs.pp.property_module


# -----------------------------------------------------------------------------
# Test _get_indexing_sets
def test_get_indexing_sets():
    m = ConcreteModel()
    m.pp = PropertyParameterBlock()
    m.cv = CVFrame(property_package=m.pp)
    m.cv._get_property_package()
    m.cv._get_indexing_sets()

    assert hasattr(m.cv, "phase_list")
    assert hasattr(m.cv, "component_list")

# -----------------------------------------------------------------------------
# Test _get_reaction_package
def test_get_reaction_package_none():
    m = ConcreteModel()
    m.r = CVFrame()

    m.r._get_reaction_package()

    assert hasattr(m.r, "reaction_module") is False


@declare_process_block_class("ReactionParameterBlock")
class _ReactionParameterBlock(ReactionParameterBase):
    def build(self):
        frm = inspect.stack()[1]
        self.property_module = inspect.getmodule(frm[0])


def test_get_reaction_package_module():
    m = ConcreteModel()
    m.rp = ReactionParameterBlock(default_arguments={"test": "foo"})
    m.cv = CVFrame(reaction_package=m.rp)

    m.cv._get_reaction_package()

    assert m.cv.reaction_module == m.rp.property_module
    assert m.cv.config.reaction_package_args["test"] == "foo"


def test_get_reaction_package_module_default_args():
    # Test that local and default args combine correctly
    m = ConcreteModel()
    m.rp = ReactionParameterBlock(default_arguments={"test1": "foo",
                                                     "test2": "bar"})
    m.cv = CVFrame(reaction_package=m.rp,
                   reaction_package_args={"test2": "baz",
                                          "test3": "bar"})

    m.cv._get_reaction_package()

    assert m.cv.config.reaction_package_args["test1"] == "foo"
    assert m.cv.config.reaction_package_args["test2"] == "baz"
    assert m.cv.config.reaction_package_args["test3"] == "bar"
