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
import pytest
from pyomo.environ import ConcreteModel, Block, Set
from pyomo.common.config import ConfigValue
from idaes.core import (ControlVolumeBase, declare_process_block_class,
                        FlowsheetBlockData, UnitBlockData)
from idaes.core.util.exceptions import (ConfigurationError,
                                        DynamicError,
                                        PropertyPackageError)

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
    assert cv.config.dynamic is None
    assert cv.config.property_package is None
    assert cv.config.property_package_args is None
    assert cv.config.reaction_package is None
    assert cv.config.reaction_package_args is None


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

