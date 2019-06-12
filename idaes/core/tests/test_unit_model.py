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
Tests for unit_model.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import Block, ConcreteModel, Set, Var
from pyomo.network import Port

from idaes.core import (FlowsheetBlockData, declare_process_block_class,
                        UnitModelBlockData, useDefault, PhysicalParameterBlock,
                        StateBlock, StateBlockData, ControlVolume0DBlock)
from idaes.core.util.exceptions import ConfigurationError, DynamicError


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("PhysicalParameterTestBlock")
class _PhysicalParameterBlock(PhysicalParameterBlock):
    def build(self):
        super(_PhysicalParameterBlock, self).build()

        self.phase_list = Set(initialize=["p1", "p2"])
        self.component_list = Set(initialize=["c1", "c2"])

        self.state_block_class = TestStateBlock


@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class StateTestBlockData(StateBlockData):
    def build(self):
        super(StateTestBlockData, self).build()

        self.a = Var(initialize=1)
        self.b = Var(initialize=2)
        self.c = Var(initialize=3)

    def define_state_vars(self):
        return {"a": self.a,
                "b": self.b,
                "c": self.c}


@declare_process_block_class("Unit")
class UnitData(UnitModelBlockData):
    def build(self):
        super(UnitModelBlockData, self).build()


def test_config_block():
    m = ConcreteModel()

    m.u = Unit()

    assert len(m.u. config) == 2
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
    # Intermediate Block is required, as ConcreteModels have a config attribute
    m.b = Block()
    m.b.u = Unit()

    with pytest.raises(DynamicError):
        m.b.u._setup_dynamics()


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


def test_setup_dynamics_has_holdup():
    # Test that has_holdup argument is True when dynamic is True
    m = ConcreteModel()

    m.fs = Flowsheet(default={"dynamic": True})

    m.fs.u = Unit()
    m.fs.u.config.has_holdup=False

    with pytest.raises(ConfigurationError):
        m.fs.u._setup_dynamics()


def test_add_port():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.prop = TestStateBlock(m.fs.time,
                             default={"parameters": m.fs.pp})

    p_obj = m.fs.u.add_port(name="test_port", block=m.fs.u.prop)

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "test_port")
    assert len(m.fs.u.test_port) == 1
    assert m.fs.u.test_port.a[0].value == m.fs.u.prop[0].a.value
    assert m.fs.u.test_port.b[0].value == m.fs.u.prop[0].b.value
    assert m.fs.u.test_port.c[0].value == m.fs.u.prop[0].c.value


def test_add_port_invalid_block():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.prop = TestStateBlock(m.fs.time,
                             default={"parameters": m.fs.pp})

    with pytest.raises(ConfigurationError):
        m.fs.u.add_port(name="test_port", block=m.fs.u)


def test_add_inlet_port_CV0D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_inlet_port()

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "inlet")
    assert len(m.fs.u.inlet) == 1

    # Set new inlet conditions to differentiate from outlet
    m.fs.u.control_volume.properties_in[0].a = 10
    m.fs.u.control_volume.properties_in[0].b = 20
    m.fs.u.control_volume.properties_in[0].c = 30

    assert m.fs.u.inlet.a[0].value == \
        m.fs.u.control_volume.properties_in[0].a.value
    assert m.fs.u.inlet.b[0].value == \
        m.fs.u.control_volume.properties_in[0].b.value
    assert m.fs.u.inlet.c[0].value == \
        m.fs.u.control_volume.properties_in[0].c.value


def test_add_inlet_port_CV0D_no_default_block():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    with pytest.raises(ConfigurationError):
        m.fs.u.add_inlet_port()


def test_add_inlet_port_CV0D_full_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.cv.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_inlet_port(name="test_port",
                                  block=m.fs.u.cv,
                                  doc="Test")

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "test_port")
    assert len(m.fs.u.test_port) == 1

    # Set new inlet conditions to differentiate from outlet
    m.fs.u.cv.properties_in[0].a = 10
    m.fs.u.cv.properties_in[0].b = 20
    m.fs.u.cv.properties_in[0].c = 30

    assert m.fs.u.test_port.a[0].value == \
        m.fs.u.cv.properties_in[0].a.value
    assert m.fs.u.test_port.b[0].value == \
        m.fs.u.cv.properties_in[0].b.value
    assert m.fs.u.test_port.c[0].value == \
        m.fs.u.cv.properties_in[0].c.value


def test_add_inlet_port_CV0D_part_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.cv.add_state_blocks()

    with pytest.raises(ConfigurationError):
        m.fs.u.add_inlet_port(block=m.fs.u.cv)

    with pytest.raises(ConfigurationError):
        m.fs.u.add_inlet_port(name="foo")


def test_add_outlet_port_CV0D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_outlet_port()

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "outlet")
    assert len(m.fs.u.outlet) == 1

    # Set new outlet conditions to differentiate from intlet
    m.fs.u.control_volume.properties_out[0].a = 10
    m.fs.u.control_volume.properties_out[0].b = 20
    m.fs.u.control_volume.properties_out[0].c = 30

    assert m.fs.u.outlet.a[0].value == \
        m.fs.u.control_volume.properties_out[0].a.value
    assert m.fs.u.outlet.b[0].value == \
        m.fs.u.control_volume.properties_out[0].b.value
    assert m.fs.u.outlet.c[0].value == \
        m.fs.u.control_volume.properties_out[0].c.value


def test_add_outlet_port_CV0D_no_default_block():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    with pytest.raises(ConfigurationError):
        m.fs.u.add_outlet_port()


def test_add_outlet_port_CV0D_full_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.cv.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_outlet_port(name="test_port",
                                  block=m.fs.u.cv,
                                  doc="Test")

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "test_port")
    assert len(m.fs.u.test_port) == 1

    # Set new outlet conditions to differentiate from inlet
    m.fs.u.cv.properties_out[0].a = 10
    m.fs.u.cv.properties_out[0].b = 20
    m.fs.u.cv.properties_out[0].c = 30

    assert m.fs.u.test_port.a[0].value == \
        m.fs.u.cv.properties_out[0].a.value
    assert m.fs.u.test_port.b[0].value == \
        m.fs.u.cv.properties_out[0].b.value
    assert m.fs.u.test_port.c[0].value == \
        m.fs.u.cv.properties_out[0].c.value


def test_add_inlet_port_CV0D_part_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.u.add_outlet_port(block=m.fs.u.cv)

    with pytest.raises(ConfigurationError):
        m.fs.u.add_outlet_port(name="foo")


def test_fix_unfix_initial_conditions():
    fs = Flowsheet(default={"dynamic": True, "time_set": [0, 1, 2]},
                   concrete=True)
    fs._setup_dynamics()

    fs.b = Unit()

    fs.b.material_accumulation = Var(fs.time, ["a", "b", "c"])
    fs.b.element_accumulation = Var(fs.time, ["a", "b", "c"])
    fs.b.enthalpy_accumulation = Var(fs.time)

    fs.fix_initial_conditions()

    for t in fs.time:
        for j in ["a", "b", "c"]:
            if t == 0:
                assert fs.b.material_accumulation[t, j].fixed
                assert fs.b.element_accumulation[t, j].fixed
                assert fs.b.enthalpy_accumulation[t].fixed
            else:
                assert fs.b.material_accumulation[t, j].fixed is False
                assert fs.b.element_accumulation[t, j].fixed is False
                assert fs.b.enthalpy_accumulation[t].fixed is False

    fs.unfix_initial_conditions()

    for t in fs.time:
        for j in ["a", "b", "c"]:
            assert fs.b.material_accumulation[t, j].fixed is False
            assert fs.b.element_accumulation[t, j].fixed is False
            assert fs.b.enthalpy_accumulation[t].fixed is False


def test_get_stream_table_contents_CV0D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    m.fs.u.add_inlet_port()
    m.fs.u.add_outlet_port()

    df = m.fs.u._get_stream_table_contents()

    assert df.loc["a"]["Inlet"] == 1
    assert df.loc["b"]["Inlet"] == 2
    assert df.loc["c"]["Inlet"] == 3

    assert df.loc["a"]["Outlet"] == 1
    assert df.loc["b"]["Outlet"] == 2
    assert df.loc["c"]["Outlet"] == 3


def test_get_stream_table_contents_CV0D_missing_default_port():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.u._get_stream_table_contents()


def test_report_CV0D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(
            default={"property_package": m.fs.pp})

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    m.fs.u.add_inlet_port()
    m.fs.u.add_outlet_port()

    m.fs.u.report()
