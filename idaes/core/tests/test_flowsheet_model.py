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
from pyomo.environ import AbstractModel, Block, ConcreteModel, Set
from pyomo.dae import ContinuousSet
from idaes.core import FlowsheetBlockData, declare_process_block_class, \
                        PropertyParameterBase
from idaes.ui.report import degrees_of_freedom


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        pass


@declare_process_block_class("ParameterBlock")
class _ParameterBlock(PropertyParameterBase):
    def build(self):
        pass

def test_config_block():
    # Test that ConfigBlock has correct attributes
    fs = Flowsheet(dynamic=True, time_set=[2.0], concrete=True)

    assert fs.config.dynamic is True
    assert fs.config.time_set == [2.0]
    assert fs.config.default_property_package is None


def test_config_validation():
    # Test that ConfigBlock validation
    fs = Flowsheet(dynamic=True, time_set=[2.0], concrete=True)

    fs.p = ParameterBlock()

    assert len(fs.config) == 3

    # Test dynamic attribute - valid values
    fs.config.dynamic = False
    fs.config.dynamic = 'use_parent_value'
    # Test dynamic attribute - invalid values
    with pytest.raises(ValueError):
        fs.config.dynamic = "foo"  # invalid str
    with pytest.raises(ValueError):
        fs.config.dynamic = 5  # invalid int
    with pytest.raises(ValueError):
        fs.config.dynamic = 2.0  # invalid float
    with pytest.raises(ValueError):
        fs.config.dynamic = [2.0]  # invalid list
    with pytest.raises(ValueError):
        fs.config.dynamic = {'a': 2.0}  # invalid dict

    # Test time_set attribute - valid values
    fs.config.time_set = [1, 2, 3]
    fs.config.time_set = 5
    fs.config.time_set = 2.0
    # Test time_set attribute - invalid values
    with pytest.raises(ValueError):
        fs.config.time_set = "foo"  # invalid str
    with pytest.raises(ValueError):
        fs.config.time_set = {'a': 2.0}  # invalid dict

    # Test default_property_package attribute - valid values
    fs.config.default_property_package = fs.p
    # Test default_property_package - invalid values
    with pytest.raises(ValueError):
        fs.config.default_property_package = "foo"  # invalid str
    with pytest.raises(ValueError):
        fs.config.default_property_package = 5  # invalid int
    with pytest.raises(ValueError):
        fs.config.default_property_package = 2.0  # invalid float
    with pytest.raises(ValueError):
        fs.config.default_property_package = [2.0]  # invalid list
    with pytest.raises(ValueError):
        fs.config.default_property_package = {'a': 2.0}  # invalid dict


def test_setup_dynamics_top_level_concrete():
    # Test that method runs when flowsheet attached to ConcreteModel
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs._setup_dynamics()

    assert m.fs.time == [0]


def test_setup_dynamics_top_level_concrete2():
    # Test that method runs when flowsheet concrete flag set to True
    fs = Flowsheet(dynamic=False, concrete=True)
    fs._setup_dynamics()

    assert fs.time == [0]


def test_setup_dynamics_top_level_abstract():
    # Test that method returns AttributeError when flowsheet is not constructed
    fs = Flowsheet(dynamic=False)
    with pytest.raises(AttributeError):
        fs._setup_dynamics()


def test_dynamic_flowsheet_in_ss_flowsheet():
    # Test that declaring a dynamic flowsheet within a steady-state flowsheet
    # raises a ValueError
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs._setup_dynamics()
    m.fs.blk = Flowsheet(dynamic=True)
    with pytest.raises(ValueError):
        m.fs.blk._setup_dynamics()


def test_setup_dynamics_top_level_abstract2():
    # Test that method returns AttributeError when flowsheet attached to an
    # unconstructed AbstractModel
    m = AbstractModel()
    m.fs = Flowsheet(dynamic=False)
    with pytest.raises(AttributeError):
        m.fs._setup_dynamics()


def test_setup_dynamics_default():
    # Test that dynamic flag creates a ContinuousSet
    fs = Flowsheet(concrete=True)
    fs._setup_dynamics()

    assert not fs.config.dynamic
    assert isinstance(fs.time, Set)


def test_setup_dynamics_dynamic():
    # Test that dynamic flag creates a ContinuousSet
    fs = Flowsheet(dynamic=True, concrete=True)
    fs._setup_dynamics()

    assert isinstance(fs.time, ContinuousSet)
    assert fs.time == [0.0, 1.0]


def test_setup_dynamics_dynamic_invalid_time_set():
    # Test that dynamic flag creates a ContinuousSet
    fs = Flowsheet(dynamic=True, time_set=[1], concrete=True)
    with pytest.raises(ValueError):
        fs._setup_dynamics()


def test_setup_dynamics_steady_state():
    # Test that dynamic flag creates a Set
    fs = Flowsheet(dynamic=False, concrete=True)
    fs._setup_dynamics()

    assert isinstance(fs.time, Set)


def test_setup_dynamics_ConcreteModel_with_time():
    m = ConcreteModel()
    m.time = Set(initialize=[0])
    m.fs = Flowsheet(dynamic=False)
    m.fs._setup_dynamics()

    assert m.fs.time == m.time


def test_setup_dynamics_ConcreteModel_with_invalid_time():
    m = ConcreteModel()
    m.time = 1
    m.fs = Flowsheet(dynamic=False)
    with pytest.raises(TypeError):
        m.fs._setup_dynamics()


def test_setup_dynamics_subflowsheet():
    # Test that subflowsheets get parents time domain
    fs = Flowsheet(dynamic=True, concrete=True)
    fs._setup_dynamics()

    fs.sub = Flowsheet()
    fs.sub._setup_dynamics()

    assert isinstance(fs.time, ContinuousSet)
    assert fs.time == [0.0, 1.0]


def test_setup_dynamics_timeset_dynamic():
    # Test that timeset argument works for dynamic=True
    fs = Flowsheet(dynamic=True, time_set=[0.0, 5.0, 100.0], concrete=True)
    fs._setup_dynamics()

    assert isinstance(fs.time, ContinuousSet)
    assert fs.time == [0.0, 5.0, 100.0]


def test_setup_dynamics_timeset_steady_state():
    # Test that timeset argument works for dynamic=False
    fs = Flowsheet(dynamic=False, time_set=[0.0, 5.0, 100.0], concrete=True)
    fs._setup_dynamics()

    assert isinstance(fs.time, Set)
    assert fs.time == [0.0, 5.0, 100.0]


def test_setup_dynamics_subflowsheet_parent_with_no_time():
    # Test that subflowsheets get parents time domain
    fs = Flowsheet(dynamic=True, concrete=True)

    fs.sub = Flowsheet()
    with pytest.raises(AttributeError):
        fs.sub._setup_dynamics()


def test_build_method():
    # Test that build method automatically calls _setup_dynamics
    fs = Flowsheet(dynamic=False, concrete=True)
    super(_Flowsheet, fs).build()

    assert isinstance(fs.time, Set)
