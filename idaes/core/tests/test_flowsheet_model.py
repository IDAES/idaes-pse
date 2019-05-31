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
from pyomo.environ import ConcreteModel, Set
from pyomo.dae import ContinuousSet
from idaes.core import FlowsheetBlockData, declare_process_block_class, \
                        PhysicalParameterBlock, useDefault
from idaes.core.util.exceptions import DynamicError


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(FlowsheetBlockData, self).build()


@declare_process_block_class("ParameterBlock")
class _ParameterBlock(PhysicalParameterBlock):
    pass


def test_config_block():
    # Test that ConfigBlock has correct attributes
    fs = Flowsheet(default={"dynamic": True,
                            "time_set": [2.0]},
                   concrete=True)

    assert fs.config.dynamic is True
    assert fs.config.time_set == [2.0]
    assert fs.config.default_property_package is None


def test_is_flowsheet():
    # Test that flowsheet has is_flowsheet method and that it returns True
    fs = Flowsheet(default={"dynamic": True,
                            "time_set": [2.0]},
                   concrete=True)

    assert hasattr(fs, "is_flowsheet")
    assert fs.is_flowsheet()


def test_config_validation():
    # Test that ConfigBlock validation
    fs = Flowsheet(default={"dynamic": True,
                            "time_set": [2.0]},
                   concrete=True)

    fs.p = ParameterBlock()

    assert len(fs.config) == 4

    # Test dynamic attribute - valid values
    fs.config.dynamic = False
    fs.config.dynamic = None
    fs.config.dynamic = useDefault
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

    # Test time attribute
    assert fs.config.time is None

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


def test_setup_dynamics_defaults():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs._setup_dynamics()

    assert m.fs.config.dynamic is False
    assert isinstance(m.fs.time, Set)
    assert m.fs.time == [0]
    assert m.fs.config.time is m.fs.time


def test_setup_dynamics_ss_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs._setup_dynamics()

    assert m.fs.config.dynamic is False
    assert isinstance(m.fs.time, Set)
    assert m.fs.time == [0]
    assert m.fs.config.time is m.fs.time


def test_setup_dynamics_ss_time_set():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False, "time_set": [1, 2, 3]})
    m.fs._setup_dynamics()

    assert m.fs.config.dynamic is False
    assert isinstance(m.fs.time, Set)
    for t in m.fs.time:
        assert t in [1, 2, 3]
    assert len(m.fs.time) == 3
    assert m.fs.config.time is m.fs.time


def test_setup_dynamics_dynamic_default():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs._setup_dynamics()

    assert m.fs.config.dynamic is True
    assert isinstance(m.fs.time, ContinuousSet)
    for t in m.fs.time:
        assert t in [0, 1]
    assert m.fs.config.time is m.fs.time


def test_setup_dynamics_dynamic_time_set():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True, "time_set": [1, 2]})
    m.fs._setup_dynamics()

    assert m.fs.config.dynamic is True
    assert isinstance(m.fs.time, ContinuousSet)
    for t in m.fs.time:
        assert t in [1, 2]
    assert m.fs.config.time is m.fs.time


def test_setup_dynamics_ss_external():
    m = ConcreteModel()
    m.s = Set(initialize=[4, 5])
    m.fs = Flowsheet(default={"dynamic": False, "time": m.s})
    m.fs._setup_dynamics()

    assert m.fs.config.dynamic is False
    assert m.fs.config.time is m.s
    assert not hasattr(m.fs, "time")


def test_setup_dynamics_dynamic_external():
    m = ConcreteModel()
    m.s = ContinuousSet(initialize=[4, 5])
    m.fs = Flowsheet(default={"dynamic": True, "time": m.s})
    m.fs._setup_dynamics()

    assert m.fs.config.dynamic is True
    assert m.fs.config.time is m.s
    assert not hasattr(m.fs, "time")


def test_setup_dynamics_dynamic_invalid_external():
    m = ConcreteModel()
    m.s = Set(initialize=[4, 5])
    m.fs = Flowsheet(default={"dynamic": True, "time": m.s})
    with pytest.raises(DynamicError):
        m.fs._setup_dynamics()


def test_setup_dynamics_sub_defaults():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs._setup_dynamics()

    m.fs.sub = Flowsheet()
    m.fs.sub._setup_dynamics()

    assert m.fs.sub.config.dynamic is False
    assert m.fs.sub.config.time is m.fs.config.time


def test_setup_dynamics_sub_ss():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs._setup_dynamics()

    m.fs.sub = Flowsheet(default={"dynamic": False})
    m.fs.sub._setup_dynamics()

    assert m.fs.sub.config.dynamic is False
    assert m.fs.sub.config.time is m.fs.config.time


def test_setup_dynamics_sub_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs._setup_dynamics()

    m.fs.sub = Flowsheet(default={"dynamic": True})
    m.fs.sub._setup_dynamics()

    assert m.fs.sub.config.dynamic is True
    assert m.fs.sub.config.time is m.fs.config.time


def test_setup_dynamics_sub_dynamic_in_ss():
    m = ConcreteModel()
    m.fs = Flowsheet(default={"dynamic": False})
    m.fs._setup_dynamics()

    m.fs.sub = Flowsheet(default={"dynamic": True})
    with pytest.raises(DynamicError):
        m.fs.sub._setup_dynamics()


def test_setup_dynamics_sub_ss_external():
    m = ConcreteModel()
    m.s = Set(initialize=[4, 5])
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs._setup_dynamics()

    m.fs.sub = Flowsheet(default={"dynamic": False, "time": m.s})
    m.fs.sub._setup_dynamics()

    assert m.fs.sub.config.dynamic is False
    assert m.fs.sub.config.time is m.s


def test_setup_dynamics_sub_dynamic_external():
    m = ConcreteModel()
    m.s = ContinuousSet(initialize=[4, 5])
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs._setup_dynamics()

    m.fs.sub = Flowsheet(default={"dynamic": True, "time": m.s})
    m.fs.sub._setup_dynamics()

    assert m.fs.sub.config.dynamic is True
    assert m.fs.sub.config.time is m.s


def test_setup_dynamics_sub_dynamic_external_invalid():
    m = ConcreteModel()
    m.s = Set(initialize=[4, 5])
    m.fs = Flowsheet(default={"dynamic": True})
    m.fs._setup_dynamics()

    m.fs.sub = Flowsheet(default={"dynamic": True, "time": m.s})
    with pytest.raises(DynamicError):
        m.fs.sub._setup_dynamics()
