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
Tests for unit_model.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import Block, ConcreteModel, Constraint, Var, units
from pyomo.network import Port

from idaes.core import (
    FlowsheetBlockData,
    declare_process_block_class,
    UnitModelBlockData,
    useDefault,
    ControlVolume0DBlock,
    ControlVolume1DBlock,
    MaterialBalanceType,
    FlowDirection,
)
from idaes.core.util.exceptions import (
    BalanceTypeNotSupportedError,
    ConfigurationError,
    DynamicError,
)
from idaes.core.util.testing import PhysicalParameterTestBlock


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("Unit")
class UnitData(UnitModelBlockData):
    def build(self):
        super(UnitModelBlockData, self).build()


@pytest.mark.unit
def test_config_block():
    m = ConcreteModel()

    m.u = Unit()

    assert len(m.u.config) == 2
    assert m.u.config.dynamic == useDefault


@pytest.mark.unit
def test_config_args():
    m = ConcreteModel()

    m.u = Unit(dynamic=True)

    assert m.u.config.dynamic is True


@pytest.mark.unit
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
        m.u.config.dynamic = {"a": 2.0}  # invalid dict


@pytest.mark.unit
def test_setup_dynamics1():
    # Test that _setup_dynamics gets argument from parent
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=False)

    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    assert m.fs.u.config.dynamic is False


@pytest.mark.unit
def test_setup_dynamics2():
    # Test that _setup_dynamics returns an DynamicError when parent has no
    # dynamic config argument

    m = ConcreteModel()
    # Intermediate Block is required, as ConcreteModels have a config attribute
    m.b = Block()
    m.b.u = Unit()

    with pytest.raises(DynamicError):
        m.b.u._setup_dynamics()


@pytest.mark.unit
def test_setup_dynamics_dynamic_in_steady_state():
    # Test that a DynamicError is raised when a dynamic models is placed in a
    # steady-state parent
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=False)

    m.fs.u = Unit(dynamic=True)
    with pytest.raises(DynamicError):
        m.fs.u._setup_dynamics()


@pytest.mark.unit
def test_setup_dynamics_get_time():
    # Test that time domain is collected correctly
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=False)

    m.fs.u = Unit()
    m.fs.u._setup_dynamics()


@pytest.mark.unit
def test_setup_dynamics_has_holdup():
    # Test that has_holdup argument is True when dynamic is True
    m = ConcreteModel()

    m.fs = Flowsheet(dynamic=True, time_units=units.s)

    m.fs.u = Unit()
    m.fs.u.config.has_holdup = False

    with pytest.raises(ConfigurationError):
        m.fs.u._setup_dynamics()


@pytest.mark.unit
def test_add_port():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.prop = m.fs.pp.build_state_block(m.fs.time)

    p_obj = m.fs.u.add_port(name="test_port", block=m.fs.u.prop)

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "test_port")
    assert len(m.fs.u.test_port) == 1
    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert (
                m.fs.u.test_port.component_flow_phase[0, p, j].value
                == m.fs.u.prop[0].flow_mol_phase_comp[p, j].value
            )
    assert m.fs.u.test_port.pressure[0].value == m.fs.u.prop[0].pressure.value
    assert m.fs.u.test_port.temperature[0].value == m.fs.u.prop[0].temperature.value


@pytest.mark.unit
def test_add_port_invalid_block():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.prop = m.fs.pp.build_state_block(m.fs.time)

    with pytest.raises(
        ConfigurationError,
        match="fs.u block object provided to add_port method is not an "
        "instance of a StateBlock object \(does not have a build_port method\).",
    ):
        m.fs.u.add_port(name="test_port", block=m.fs.u)


@pytest.mark.unit
def test_add_inlet_port_CV0D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_inlet_port()

    assert isinstance(p_obj, Port)
    assert m.fs.u.inlet is p_obj
    assert len(m.fs.u.inlet) == 1

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert (
                m.fs.u.inlet.component_flow_phase[0, p, j]
                is m.fs.u.control_volume.properties_in[0].flow_mol_phase_comp[p, j]
            )
    assert m.fs.u.inlet.pressure[0] is m.fs.u.control_volume.properties_in[0].pressure
    assert (
        m.fs.u.inlet.temperature[0]
        is m.fs.u.control_volume.properties_in[0].temperature
    )


@pytest.mark.unit
def test_add_inlet_port_CV1D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.u.control_volume.add_geometry()
    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_inlet_port()

    assert isinstance(p_obj, Port)
    assert m.fs.u.inlet is p_obj
    assert len(m.fs.u.inlet) == 1

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert (
                m.fs.u.inlet.component_flow_phase[0, p, j]
                is m.fs.u.control_volume.properties[0, 0].flow_mol_phase_comp[p, j]
            )
    assert m.fs.u.inlet.pressure[0] is m.fs.u.control_volume.properties[0, 0].pressure
    assert (
        m.fs.u.inlet.temperature[0]
        is m.fs.u.control_volume.properties[0, 0].temperature
    )


@pytest.mark.unit
def test_add_inlet_port_CV1D_backward():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.u.control_volume.add_geometry(flow_direction=FlowDirection.backward)
    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_inlet_port()

    assert isinstance(p_obj, Port)
    assert m.fs.u.inlet is p_obj
    assert len(m.fs.u.inlet) == 1

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert (
                m.fs.u.inlet.component_flow_phase[0, p, j]
                is m.fs.u.control_volume.properties[0, 1].flow_mol_phase_comp[p, j]
            )
    assert m.fs.u.inlet.pressure[0] is m.fs.u.control_volume.properties[0, 1].pressure
    assert (
        m.fs.u.inlet.temperature[0]
        is m.fs.u.control_volume.properties[0, 1].temperature
    )


@pytest.mark.unit
def test_add_inlet_port_CV0D_no_default_block():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    with pytest.raises(
        ConfigurationError,
        match="fs.u add_inlet_port was called without a block argument"
        " but no default ControlVolume exists \(control_volume\). "
        "Please provide block to which the Port should be associated.",
    ):
        m.fs.u.add_inlet_port()


@pytest.mark.unit
def test_add_inlet_port_CV0D_full_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.cv.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_inlet_port(name="test_port", block=m.fs.u.cv, doc="Test")

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "test_port")
    assert len(m.fs.u.test_port) == 1

    # Set new inlet conditions to differentiate from outlet
    m.fs.u.cv.properties_in[0].a = 10
    m.fs.u.cv.properties_in[0].b = 20
    m.fs.u.cv.properties_in[0].c = 30

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert m.fs.u.test_port.component_flow_phase[0, p, j].value == (
                m.fs.u.cv.properties_in[0].flow_mol_phase_comp[p, j].value
            )
    assert (
        m.fs.u.test_port.pressure[0].value == m.fs.u.cv.properties_in[0].pressure.value
    )
    assert (
        m.fs.u.test_port.temperature[0].value
        == m.fs.u.cv.properties_in[0].temperature.value
    )


@pytest.mark.unit
def test_add_inlet_port_CV0D_part_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(
        ConfigurationError,
        match="fs.u add_inlet_port was called with a block argument, "
        "but a name argument was not provided. Either both "
        "a name and a block must be provided or neither.",
    ):
        m.fs.u.add_inlet_port(block=m.fs.u.cv)

    with pytest.raises(
        ConfigurationError,
        match="fs.u add_inlet_port was called without a block "
        "argument but a name argument was provided. Either "
        "both a name and a block must be provided or neither.",
    ):
        m.fs.u.add_inlet_port(name="foo")


@pytest.mark.unit
def test_add_inlet_port_CV0D_no_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(property_package=m.fs.pp)

    with pytest.raises(
        ConfigurationError,
        match="fs.u - control volume does not have expected "
        "names for StateBlocks. Please check that the control "
        "volume was constructed correctly.",
    ):
        m.fs.u.add_inlet_port()


@pytest.mark.unit
def test_add_outlet_port_CV0D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_outlet_port()

    assert isinstance(p_obj, Port)
    assert m.fs.u.outlet is p_obj
    assert len(m.fs.u.outlet) == 1

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert (
                m.fs.u.outlet.component_flow_phase[0, p, j]
                is m.fs.u.control_volume.properties_out[0].flow_mol_phase_comp[p, j]
            )
    assert m.fs.u.outlet.pressure[0] is m.fs.u.control_volume.properties_out[0].pressure
    assert (
        m.fs.u.outlet.temperature[0]
        is m.fs.u.control_volume.properties_out[0].temperature
    )


@pytest.mark.unit
def test_add_outlet_port_CV1D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.u.control_volume.add_geometry()
    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_outlet_port()

    assert isinstance(p_obj, Port)
    assert m.fs.u.outlet is p_obj
    assert len(m.fs.u.outlet) == 1

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert (
                m.fs.u.outlet.component_flow_phase[0, p, j]
                is m.fs.u.control_volume.properties[0, 1].flow_mol_phase_comp[p, j]
            )
    assert m.fs.u.outlet.pressure[0] is m.fs.u.control_volume.properties[0, 1].pressure
    assert (
        m.fs.u.outlet.temperature[0]
        is m.fs.u.control_volume.properties[0, 1].temperature
    )


@pytest.mark.unit
def test_add_outlet_port_CV1D_backward():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.u.control_volume.add_geometry(flow_direction=FlowDirection.backward)
    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_outlet_port()

    assert isinstance(p_obj, Port)
    assert m.fs.u.outlet is p_obj
    assert len(m.fs.u.outlet) == 1

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert (
                m.fs.u.outlet.component_flow_phase[0, p, j]
                is m.fs.u.control_volume.properties[0, 0].flow_mol_phase_comp[p, j]
            )
    assert m.fs.u.outlet.pressure[0] is m.fs.u.control_volume.properties[0, 0].pressure
    assert (
        m.fs.u.outlet.temperature[0]
        is m.fs.u.control_volume.properties[0, 0].temperature
    )


@pytest.mark.unit
def test_add_outlet_port_CV0D_no_default_block():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    with pytest.raises(
        ConfigurationError,
        match="fs.u add_outlet_port was called without a block "
        "argument but no default ControlVolume exists "
        "\(control_volume\). Please provide block to which the "
        "Port should be associated.",
    ):
        m.fs.u.add_outlet_port()


@pytest.mark.unit
def test_add_outlet_port_CV0D_full_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.cv.add_state_blocks(has_phase_equilibrium=False)

    p_obj = m.fs.u.add_outlet_port(name="test_port", block=m.fs.u.cv, doc="Test")

    assert isinstance(p_obj, Port)
    assert hasattr(m.fs.u, "test_port")
    assert len(m.fs.u.test_port) == 1

    # Set new outlet conditions to differentiate from inlet
    m.fs.u.cv.properties_out[0].a = 10
    m.fs.u.cv.properties_out[0].b = 20
    m.fs.u.cv.properties_out[0].c = 30

    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            assert m.fs.u.test_port.component_flow_phase[0, p, j].value == (
                m.fs.u.cv.properties_out[0].flow_mol_phase_comp[p, j].value
            )
    assert (
        m.fs.u.test_port.pressure[0].value == m.fs.u.cv.properties_out[0].pressure.value
    )
    assert (
        m.fs.u.test_port.temperature[0].value
        == m.fs.u.cv.properties_out[0].temperature.value
    )


@pytest.mark.unit
def test_add_outlet_port_CV0D_part_args():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(
        ConfigurationError,
        match="fs.u add_outlet_port was called with a block argument, "
        "but a name argument was not provided. Either both "
        "a name and a block must be provided or neither.",
    ):
        m.fs.u.add_outlet_port(block=m.fs.u.cv)

    with pytest.raises(
        ConfigurationError,
        match="fs.u add_outlet_port was called without a block "
        "argument but a name argument was provided. Either "
        "both a name and a block must be provided or neither.",
    ):
        m.fs.u.add_outlet_port(name="foo")


@pytest.mark.unit
def test_add_outlet_port_CV0D_no_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(property_package=m.fs.pp)

    with pytest.raises(
        ConfigurationError,
        match="fs.u - control volume does not have expected "
        "names for StateBlocks. Please check that the control "
        "volume was constructed correctly.",
    ):
        m.fs.u.add_outlet_port()


@pytest.mark.unit
def test_fix_unfix_initial_conditions():
    fs = Flowsheet(dynamic=True, time_set=[0, 1, 2], time_units=units.s, concrete=True)
    fs._setup_dynamics()

    fs.b = Unit()

    fs.b.material_accumulation = Var(fs.time, ["a", "b", "c"])
    fs.b.element_accumulation = Var(fs.time, ["a", "b", "c"])
    fs.b.energy_accumulation = Var(fs.time)

    fs.fix_initial_conditions()

    for t in fs.time:
        for j in ["a", "b", "c"]:
            if t == 0:
                assert fs.b.material_accumulation[t, j].fixed
                assert fs.b.element_accumulation[t, j].fixed
                assert fs.b.energy_accumulation[t].fixed
            else:
                assert fs.b.material_accumulation[t, j].fixed is False
                assert fs.b.element_accumulation[t, j].fixed is False
                assert fs.b.energy_accumulation[t].fixed is False

    fs.unfix_initial_conditions()

    for t in fs.time:
        for j in ["a", "b", "c"]:
            assert fs.b.material_accumulation[t, j].fixed is False
            assert fs.b.element_accumulation[t, j].fixed is False
            assert fs.b.energy_accumulation[t].fixed is False


@pytest.mark.unit
def test_get_stream_table_contents_CV0D():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    m.fs.u.add_inlet_port()
    m.fs.u.add_outlet_port()

    df = m.fs.u._get_stream_table_contents()

    assert df.loc["component_flow_phase ('p1', 'c1')"]["Inlet"] == 2
    assert df.loc["component_flow_phase ('p1', 'c2')"]["Inlet"] == 2
    assert df.loc["component_flow_phase ('p2', 'c1')"]["Inlet"] == 2
    assert df.loc["component_flow_phase ('p2', 'c2')"]["Inlet"] == 2
    assert df.loc["pressure"]["Inlet"] == 1e5
    assert df.loc["temperature"]["Inlet"] == 300

    assert df.loc["component_flow_phase ('p1', 'c1')"]["Outlet"] == 2
    assert df.loc["component_flow_phase ('p1', 'c2')"]["Outlet"] == 2
    assert df.loc["component_flow_phase ('p2', 'c1')"]["Outlet"] == 2
    assert df.loc["component_flow_phase ('p2', 'c2')"]["Outlet"] == 2
    assert df.loc["pressure"]["Outlet"] == 1e5
    assert df.loc["temperature"]["Outlet"] == 300


@pytest.mark.unit
def test_get_stream_table_contents_CV0D_missing_default_port():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()
    m.fs.u._setup_dynamics()

    m.fs.u.control_volume = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.u.control_volume.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.u._get_stream_table_contents()


@pytest.mark.unit
def test_add_state_material_balances_invalid_state_1():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb = object()

    with pytest.raises(
        ConfigurationError,
        match="fs.u state_1 argument to add_state_material_balances was "
        "not an instance of a State Block.",
    ):
        m.fs.u.add_state_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            state_1=m.fs.u.sb,
            state_2=m.fs.u.sb,
        )


@pytest.mark.unit
def test_add_state_material_balances_invalid_state_2():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp.build_state_block(m.fs.time)
    m.fs.u.sb2 = object()

    with pytest.raises(
        ConfigurationError,
        match="fs.u state_2 argument to add_state_material_balances was "
        "not an instance of a State Block.",
    ):
        m.fs.u.add_state_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            state_1=m.fs.u.sb1,
            state_2=m.fs.u.sb2,
        )


@pytest.mark.unit
def test_add_state_material_balances_mixed_states():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp1 = PhysicalParameterTestBlock()
    m.fs.pp2 = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp1.build_state_block(m.fs.time)
    m.fs.u.sb2 = m.fs.pp2.build_state_block(m.fs.time)

    with pytest.raises(
        ConfigurationError,
        match="fs.u add_state_material_balances method was provided with "
        "State Blocks are not linked to the same "
        "instance of a Physical Parameter Block. This method "
        "only supports linking State Blocks from the same "
        "Physical Parameter Block.",
    ):
        m.fs.u.add_state_material_balances(
            balance_type=MaterialBalanceType.componentTotal,
            state_1=m.fs.u.sb1,
            state_2=m.fs.u.sb2,
        )


@pytest.mark.unit
def test_add_state_material_balances_component_phase():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp.build_state_block(m.fs.time)
    m.fs.u.sb2 = m.fs.pp.build_state_block(m.fs.time)

    m.fs.u.add_state_material_balances(
        balance_type=MaterialBalanceType.componentPhase,
        state_1=m.fs.u.sb1,
        state_2=m.fs.u.sb2,
    )

    assert isinstance(m.fs.u.state_material_balances, Constraint)
    assert len(m.fs.u.state_material_balances) == 4
    for k in m.fs.u.state_material_balances:
        assert k in [
            (0.0, "p1", "c1"),
            (0.0, "p1", "c2"),
            (0.0, "p2", "c1"),
            (0.0, "p2", "c2"),
        ]


@pytest.mark.unit
def test_add_state_material_balances_component_total():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp.build_state_block(m.fs.time)
    m.fs.u.sb2 = m.fs.pp.build_state_block(m.fs.time)

    m.fs.u.add_state_material_balances(
        balance_type=MaterialBalanceType.componentTotal,
        state_1=m.fs.u.sb1,
        state_2=m.fs.u.sb2,
    )

    assert isinstance(m.fs.u.state_material_balances, Constraint)
    assert len(m.fs.u.state_material_balances) == 2
    for k in m.fs.u.state_material_balances:
        assert k in [(0.0, "c1"), (0.0, "c2")]


@pytest.mark.unit
def test_add_state_material_balances_total():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp.build_state_block(m.fs.time)
    m.fs.u.sb2 = m.fs.pp.build_state_block(m.fs.time)

    m.fs.u.add_state_material_balances(
        balance_type=MaterialBalanceType.total, state_1=m.fs.u.sb1, state_2=m.fs.u.sb2
    )

    assert isinstance(m.fs.u.state_material_balances, Constraint)
    assert len(m.fs.u.state_material_balances) == 1


@pytest.mark.unit
def test_add_state_material_balances_element_total():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp.build_state_block(m.fs.time)
    m.fs.u.sb2 = m.fs.pp.build_state_block(m.fs.time)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.u.add_state_material_balances(
            balance_type=MaterialBalanceType.elementTotal,
            state_1=m.fs.u.sb1,
            state_2=m.fs.u.sb2,
        )


@pytest.mark.unit
def test_add_state_material_balances_none():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp.build_state_block(m.fs.time)
    m.fs.u.sb2 = m.fs.pp.build_state_block(m.fs.time)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.u.add_state_material_balances(
            balance_type=MaterialBalanceType.none,
            state_1=m.fs.u.sb1,
            state_2=m.fs.u.sb2,
        )


@pytest.mark.unit
def test_add_state_material_balances_double_call():
    m = ConcreteModel()
    m.fs = Flowsheet()
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.u = Unit()

    m.fs.u.sb1 = m.fs.pp.build_state_block(m.fs.time)
    m.fs.u.sb2 = m.fs.pp.build_state_block(m.fs.time)

    m.fs.u.add_state_material_balances(
        balance_type=MaterialBalanceType.componentPhase,
        state_1=m.fs.u.sb1,
        state_2=m.fs.u.sb2,
    )

    with pytest.raises(AttributeError):
        m.fs.u.add_state_material_balances(
            balance_type=MaterialBalanceType.componentPhase,
            state_1=m.fs.u.sb1,
            state_2=m.fs.u.sb2,
        )
