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
Tests for unit_model.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Expression, Set, Var
from pyomo.common.config import ConfigValue

from idaes.core import (FlowsheetBlockData, declare_process_block_class,
                        UnitBlockData, PropertyParameterBase,
                        StateBlockDataBase, useDefault)
from idaes.core.util.exceptions import ConfigurationError, DynamicError


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("Unit")
class UnitData(UnitBlockData):
    def build(self):
        super(UnitBlockData, self).build()


def test_config_block():
    m = ConcreteModel()

    m.u = Unit()

    assert len(m.u. config) == 1
    assert m.u.config.dynamic == useDefault


def test_config_args():
    m = ConcreteModel()

    m.u = Unit(default={"dynamic": True})

    assert m.u.config.dynamic is True


def test_config_args_invalid():
    # Test validation of config arguments
    m = ConcreteModel()

    m.u = Unit()

    m.u.config.dynamic = True
    m.u.config.dynamic = False
    m.u.config.dynamic = None

    # Test that Value error raised when given invalid config arguments
    with pytest.raises(ValueError):
        m.u.config.dynamic = "foo"  # invalid str
    with pytest.raises(ValueError):
        m.u.config.dynamic = 5  # invalid int
    with pytest.raises(ValueError):
        m.u.config.dynamic = 2.0  # invalid float
    with pytest.raises(ValueError):
        m.u.config.dynamic = [2.0]  # invalid list
    with pytest.raises(ValueError):
        m.u.config.dynamic = {'a': 2.0}  # invalid dict


def test_setup_dynamics1():
    # Test that _setup_dynamics gets argument from parent
    m = ConcreteModel()

    m.fs = Flowsheet(default={"dynamic": False})

    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    assert m.fs.u.config.dynamic is False


def test_setup_dynamics2():
    # Test that _setup_dynamics returns an DynamicError when parent has no
    # dynamic config argument

    m = ConcreteModel()
    m.u = Unit()

    with pytest.raises(DynamicError):
        m.u._setup_dynamics()


def test_setup_dynamics_dynamic_in_steady_state():
    # Test that a DynamicError is raised when a dynamic models is placed in a
    # steady-state parent
    m = ConcreteModel()

    m.fs = Flowsheet(default={"dynamic": False})

    m.fs.u = Unit(default={"dynamic": True})
    with pytest.raises(DynamicError):
        m.fs.u._setup_dynamics()


def test_setup_dynamics_get_time():
    # Test that time domain is collected correctly
    m = ConcreteModel()

    m.fs = Flowsheet(default={"dynamic": False})

    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    assert m.fs.u.time == m.fs.time


def test_setup_dynamics_get_time_fails():
    # Test that DynamicError is raised when parent does not have time domain
    m = ConcreteModel()

    m.u = Unit()
    with pytest.raises(DynamicError):
        m.u._setup_dynamics()


def test_setup_dynamics_has_holdup():
    # Test that has_holdup argument is True when dynamic is True
    m = ConcreteModel()

    m.fs = Flowsheet(default={"dynamic": True})

    m.fs.u = Unit()
    m.fs.u.config.declare("has_holdup", ConfigValue(default=False))

    with pytest.raises(ConfigurationError):
        m.fs.u._setup_dynamics()
