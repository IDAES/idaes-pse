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
Tests for ControlVolumeBlockData.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Constraint, Expression, Set, units, Var
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from pyomo.common.config import ConfigBlock
from idaes.core import (
    ControlVolume0DBlock,
    ControlVolumeBlockData,
    FlowsheetBlockData,
    declare_process_block_class,
    FlowDirection,
    MaterialBalanceType,
    EnergyBalanceType,
)
from idaes.core.util.exceptions import (
    BalanceTypeNotSupportedError,
    ConfigurationError,
    PropertyNotSupportedError,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
)
import idaes.logger as idaeslog


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@declare_process_block_class("CVFrame")
class CVFrameData(ControlVolume0DBlock):
    def build(self):
        super(ControlVolumeBlockData, self).build()


# -----------------------------------------------------------------------------
# Basic tests
@pytest.mark.unit
def test_base_build():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    assert len(m.fs.cv.config) == 7
    assert m.fs.cv.config.dynamic is False
    assert m.fs.cv.config.has_holdup is False
    assert m.fs.cv.config.property_package == m.fs.pp
    assert isinstance(m.fs.cv.config.property_package_args, ConfigBlock)
    assert len(m.fs.cv.config.property_package_args) == 0
    assert m.fs.cv.config.reaction_package is None
    assert isinstance(m.fs.cv.config.reaction_package_args, ConfigBlock)
    assert len(m.fs.cv.config.reaction_package_args) == 0
    assert m.fs.cv.config.auto_construct is False

    assert hasattr(m.fs.config, "time")


# -----------------------------------------------------------------------------
# Test add_geometry
@pytest.mark.unit
def test_add_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_geometry()

    assert hasattr(m.fs.cv, "volume")
    assert len(m.fs.cv.volume) == 1.0
    assert m.fs.cv.volume[0].value == 1.0


# -----------------------------------------------------------------------------
# Test add_state_blocks
@pytest.mark.unit
def test_add_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    assert hasattr(m.fs.cv, "properties_in")
    assert len(m.fs.cv.properties_in[0].config) == 3
    assert m.fs.cv.properties_in[0].config.defined_state is True
    assert m.fs.cv.properties_in[0].config.has_phase_equilibrium is False
    assert m.fs.cv.properties_in[0].config.parameters == m.fs.pp

    assert hasattr(m.fs.cv, "properties_out")
    assert len(m.fs.cv.properties_out[0].config) == 3
    assert m.fs.cv.properties_out[0].config.defined_state is False
    assert m.fs.cv.properties_out[0].config.has_phase_equilibrium is False
    assert m.fs.cv.properties_out[0].config.parameters == m.fs.pp


@pytest.mark.unit
def test_add_state_block_forward_flow():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(
        information_flow=FlowDirection.forward, has_phase_equilibrium=False
    )

    assert m.fs.cv.properties_in[0].config.defined_state is True
    assert m.fs.cv.properties_out[0].config.defined_state is False


@pytest.mark.unit
def test_add_state_block_backward_flow():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(
        information_flow=FlowDirection.backward, has_phase_equilibrium=False
    )

    assert m.fs.cv.properties_in[0].config.defined_state is False
    assert m.fs.cv.properties_out[0].config.defined_state is True


@pytest.mark.unit
def test_add_state_blocks_has_phase_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)

    assert m.fs.cv.properties_in[0].config.has_phase_equilibrium is True
    assert m.fs.cv.properties_out[0].config.has_phase_equilibrium is True


@pytest.mark.unit
def test_add_state_blocks_no_has_phase_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_state_blocks()


@pytest.mark.unit
def test_add_state_blocks_custom_args():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, property_package_args={"test": "test"}
    )

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    assert len(m.fs.cv.properties_in[0].config) == 4
    assert m.fs.cv.properties_in[0].config.test == "test"

    assert len(m.fs.cv.properties_out[0].config) == 4
    assert m.fs.cv.properties_out[0].config.test == "test"


# -----------------------------------------------------------------------------
# Test add_reaction_blocks
@pytest.mark.unit
def test_add_reaction_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    assert hasattr(m.fs.cv, "reactions")
    assert len(m.fs.cv.reactions[0].config) == 3
    assert m.fs.cv.reactions[0].config.state_block == m.fs.cv.properties_out
    assert m.fs.cv.reactions[0].state_ref == m.fs.cv.properties_out[0]
    assert m.fs.cv.reactions[0].config.has_equilibrium is False
    assert m.fs.cv.reactions[0].config.parameters == m.fs.rp


@pytest.mark.unit
def test_add_reaction_blocks_has_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    assert m.fs.cv.reactions[0].config.has_equilibrium is True


@pytest.mark.unit
def test_add_reaction_blocks_no_has_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_reaction_blocks()


@pytest.mark.unit
def test_add_reaction_blocks_custom_args():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        reaction_package_args={"test1": 1},
    )

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    assert m.fs.cv.reactions[0].config.test1 == 1


# -----------------------------------------------------------------------------
# Test _add_phase_fractions
@pytest.mark.unit
def test_add_phase_fractions():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv._add_phase_fractions()

    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert len(m.fs.cv.phase_fraction) == 2
    assert hasattr(m.fs.cv, "sum_of_phase_fractions")


@pytest.mark.unit
def test_add_phase_fractions_single_phase():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_list)
    m.fs.pp.phase_list = Set(initialize=["p1"])

    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv._add_phase_fractions()

    assert isinstance(m.fs.cv.phase_fraction, Expression)
    assert len(m.fs.cv.phase_fraction) == 1
    assert not hasattr(m.fs.cv, "sum_of_phase_fractions")


# -----------------------------------------------------------------------------
# Test reaction rate conversion method
@pytest.mark.unit
def test_rxn_rate_conv_no_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 3
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for j in m.fs.pp.component_list:
            assert m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=False) == 1


@pytest.mark.unit
def test_rxn_rate_conv_property_basis_other():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 3
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for j in m.fs.pp.component_list:
            with pytest.raises(ConfigurationError):
                m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_reaction_basis_other():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.rp.basis_switch = 3

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for j in m.fs.pp.component_list:
            with pytest.raises(ConfigurationError):
                m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_both_molar():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for j in m.fs.pp.component_list:
            assert m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True) == 1


@pytest.mark.unit
def test_rxn_rate_conv_both_mass():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 2
    m.fs.rp.basis_switch = 2

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for j in m.fs.pp.component_list:
            assert m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True) == 1


@pytest.mark.unit
def test_rxn_rate_conv_mole_mass_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 1
    m.fs.rp.basis_switch = 2

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for j in m.fs.pp.component_list:
            with pytest.raises(PropertyNotSupportedError):
                m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_mass_mole_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 2
    m.fs.rp.basis_switch = 1

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for j in m.fs.pp.component_list:
            with pytest.raises(PropertyNotSupportedError):
                m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_mole_mass():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 1
    m.fs.rp.basis_switch = 2

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = {"c1": 2, "c2": 3}
        for j in m.fs.pp.component_list:
            assert (
                m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True)
                == 1 / m.fs.cv.properties_out[t].mw_comp[j]
            )


@pytest.mark.unit
def test_rxn_rate_conv_mass_mole():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 2
    m.fs.rp.basis_switch = 1

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = {"c1": 2, "c2": 3}
        for j in m.fs.pp.component_list:
            assert (
                m.fs.cv._rxn_rate_conv(t, j, has_rate_reactions=True)
                == m.fs.cv.properties_out[t].mw_comp[j]
            )


# -----------------------------------------------------------------------------
# Test add_material_balances default
@pytest.mark.unit
def test_add_material_balances_default_fail():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.pp.default_balance_switch = 2

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_material_balances(MaterialBalanceType.useDefault)


@pytest.mark.unit
def test_add_material_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_material_balances(MaterialBalanceType.useDefault)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_material_balances_rxn_molar():
    # use property package with mass basis to confirm correct rxn term units
    # add options so that all generation/extent terms exist
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    # Set property package to contain inherent reactions
    m.fs.pp._has_inherent_reactions = True

    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.rp.basis_switch = 1

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    units = m.fs.cv.config.property_package.get_metadata().get_derived_units
    pp_units = units("flow_mass")  # basis 2 is mass
    rp_units = units("flow_mole")  # basis 1 is molar

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    # add molecular weight variable to each time point, using correct units
    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = Var(
            m.fs.cv.properties_out[t].config.parameters.component_list,
            units=units("mass") / units("amount"),
        )

    # add material balances to control volume
    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentPhase,
        has_rate_reactions=True,
        has_equilibrium_reactions=True,
        has_phase_equilibrium=True,
    )

    assert_units_equivalent(m.fs.cv.rate_reaction_generation, rp_units)
    assert_units_equivalent(m.fs.cv.rate_reaction_extent, rp_units)
    assert_units_equivalent(m.fs.cv.equilibrium_reaction_generation, rp_units)
    assert_units_equivalent(m.fs.cv.equilibrium_reaction_extent, rp_units)
    assert_units_equivalent(m.fs.cv.inherent_reaction_generation, pp_units)
    assert_units_equivalent(m.fs.cv.inherent_reaction_extent, pp_units)
    assert_units_equivalent(m.fs.cv.phase_equilibrium_generation, pp_units)


@pytest.mark.unit
def test_add_material_balances_rxn_mass():
    # use property package with mass basis to confirm correct rxn term units
    # add options so that all generation/extent terms exist
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    # Set property package to contain inherent reactions
    m.fs.pp._has_inherent_reactions = True

    m.fs.pp.basis_switch = 1
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.rp.basis_switch = 2

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    units = m.fs.cv.config.property_package.get_metadata().get_derived_units
    pp_units = units("flow_mole")  # basis 2 is molar
    rp_units = units("flow_mass")  # basis 1 is mass

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    # add molecular weight variable to each time point, using correct units
    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = Var(
            m.fs.cv.properties_out[t].config.parameters.component_list,
            units=units("mass") / units("amount"),
        )

    # add material balances to control volume
    m.fs.cv.add_material_balances(
        balance_type=MaterialBalanceType.componentPhase,
        has_rate_reactions=True,
        has_equilibrium_reactions=True,
        has_phase_equilibrium=True,
    )

    assert_units_equivalent(m.fs.cv.rate_reaction_generation, rp_units)
    assert_units_equivalent(m.fs.cv.rate_reaction_extent, rp_units)
    assert_units_equivalent(m.fs.cv.equilibrium_reaction_generation, rp_units)
    assert_units_equivalent(m.fs.cv.equilibrium_reaction_extent, rp_units)
    assert_units_equivalent(m.fs.cv.inherent_reaction_generation, pp_units)
    assert_units_equivalent(m.fs.cv.inherent_reaction_extent, pp_units)
    assert_units_equivalent(m.fs.cv.phase_equilibrium_generation, pp_units)


# -----------------------------------------------------------------------------
# Test add_phase_component_balances
@pytest.mark.unit
def test_add_phase_component_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 4

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 8
    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert isinstance(m.fs.cv.material_holdup, Var)
    assert isinstance(m.fs.cv.material_accumulation, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_dynamic_no_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    # Do not add geometry
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances()


@pytest.mark.unit
def test_add_phase_component_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances(has_rate_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.rate_reaction_generation, Var)
    assert isinstance(m.fs.cv.rate_reaction_extent, Var)
    assert isinstance(m.fs.cv.rate_reaction_stoichiometry_constraint, Constraint)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_rate_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_rate_reactions=True)


@pytest.mark.unit
def test_add_phase_component_balances_rate_rxns_no_rxn_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.rp.del_component(m.fs.rp.rate_reaction_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_rate_reactions=True)


@pytest.mark.unit
def test_add_phase_component_balances_eq_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    mb = m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.equilibrium_reaction_generation, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_extent, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_stoichiometry_constraint, Constraint)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_eq_rxns_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_phase_component_balances_eq_rxns_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.rp.del_component(m.fs.rp.equilibrium_reaction_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_phase_component_balances_in_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    # Set property package to contain inherent reactions
    m.fs.pp._has_inherent_reactions = True

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.inherent_reaction_generation, Var)
    assert isinstance(m.fs.cv.inherent_reaction_extent, Var)
    assert isinstance(m.fs.cv.inherent_reaction_stoichiometry_constraint, Constraint)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_in_rxns_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    # Set property package to contain inherent reactions
    m.fs.pp._has_inherent_reactions = True
    # delete inherent_Reaction_dix to trigger exception
    m.fs.pp.del_component(m.fs.pp.inherent_reaction_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(
        PropertyNotSupportedError,
        match="fs.cv Property package does not contain a "
        "list of inherent reactions \(inherent_reaction_idx\), "
        "but include_inherent_reactions is True.",
    ):
        m.fs.cv.add_phase_component_balances()


@pytest.mark.unit
def test_add_phase_component_balances_eq_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_phase_component_balances_phase_eq():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances(has_phase_equilibrium=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.phase_equilibrium_generation, Var)
    assert isinstance(m.fs.cv.config.property_package.phase_equilibrium_idx, Set)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_phase_eq_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_phase_equilibrium=True)


@pytest.mark.unit
def test_add_phase_component_balances_phase_eq_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_phase_equilibrium=True)


@pytest.mark.unit
def test_add_phase_component_balances_mass_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances(has_mass_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.mass_transfer_term, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_custom_molar_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j] * units.mol / units.s

    mb = m.fs.cv.add_phase_component_balances(custom_molar_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_custom_molar_term_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j]

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(custom_molar_term=custom_method)


@pytest.mark.unit
def test_add_phase_component_balances_custom_molar_term_mass_flow_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j] * units.mol / units.s

    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = Var(
            m.fs.cv.properties_out[t].config.parameters.component_list,
            units=units.kg / units.mol,
        )

    mb = m.fs.cv.add_phase_component_balances(custom_molar_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_custom_molar_term_undefined_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 3
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j]

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(custom_molar_term=custom_method)


@pytest.mark.unit
def test_add_phase_component_balances_custom_mass_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j] * units.kg / units.s

    mb = m.fs.cv.add_phase_component_balances(custom_mass_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_custom_mass_term_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 1
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j]

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(custom_mass_term=custom_method)


@pytest.mark.unit
def test_add_phase_component_balances_custom_mass_term_mole_flow_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j] * units.kg / units.s

    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = Var(
            m.fs.cv.properties_out[t].config.parameters.component_list,
            units=units.kg / units.mol,
        )

    mb = m.fs.cv.add_phase_component_balances(custom_mass_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_custom_mass_term_undefined_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 3
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time,
        m.fs.cv.config.property_package.phase_list,
        m.fs.cv.config.property_package.component_list,
    )

    def custom_method(t, p, j):
        return m.fs.cv.test_var[t, p, j]

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(custom_mass_term=custom_method)


# -----------------------------------------------------------------------------
# Test add_total_component_balances
@pytest.mark.unit
def test_add_total_component_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 2

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert isinstance(m.fs.cv.material_holdup, Var)
    assert isinstance(m.fs.cv.material_accumulation, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_dynamic_no_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    # Do not add geometry
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances()


@pytest.mark.unit
def test_add_total_component_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_component_balances(has_rate_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.rate_reaction_generation, Var)
    assert isinstance(m.fs.cv.rate_reaction_extent, Var)
    assert isinstance(m.fs.cv.rate_reaction_stoichiometry_constraint, Constraint)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_rate_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_rate_reactions=True)


@pytest.mark.unit
def test_add_total_component_balances_rate_rxns_no_rxn_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.rp.del_component(m.fs.rp.rate_reaction_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(has_rate_reactions=True)


@pytest.mark.unit
def test_add_total_component_balances_eq_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    mb = m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.equilibrium_reaction_generation, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_extent, Var)
    assert isinstance(m.fs.cv.equilibrium_reaction_stoichiometry_constraint, Constraint)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_eq_rxns_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_total_component_balances_eq_rxns_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.rp.del_component(m.fs.rp.equilibrium_reaction_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_total_component_balances_eq_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_total_component_balances_in_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    # Set property package to contain inherent reactions
    m.fs.pp._has_inherent_reactions = True

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    mb = m.fs.cv.add_total_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.inherent_reaction_generation, Var)
    assert isinstance(m.fs.cv.inherent_reaction_extent, Var)
    assert isinstance(m.fs.cv.inherent_reaction_stoichiometry_constraint, Constraint)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_in_rxns_no_idx():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    # Set property package to contain inherent reactions
    m.fs.pp._has_inherent_reactions = True
    # delete inherent_Reaction_dix to trigger exception
    m.fs.pp.del_component(m.fs.pp.inherent_reaction_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(
        PropertyNotSupportedError,
        match="fs.cv Property package does not contain a "
        "list of inherent reactions \(inherent_reaction_idx\), "
        "but include_inherent_reactions is True.",
    ):
        m.fs.cv.add_total_component_balances()


@pytest.mark.unit
def test_add_total_component_balances_phase_eq_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(has_phase_equilibrium=True)


@pytest.mark.unit
def test_add_total_component_balances_mass_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_component_balances(has_mass_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.mass_transfer_term, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_custom_molar_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j] * units.mol / units.s

    mb = m.fs.cv.add_total_component_balances(custom_molar_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_custom_molar_term_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j]

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(custom_molar_term=custom_method)


@pytest.mark.unit
def test_add_total_component_balances_custom_molar_term_mass_flow_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j] * units.mol / units.s

    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = Var(
            m.fs.cv.properties_out[t].config.parameters.component_list,
            units=units.kg / units.mol,
        )

    mb = m.fs.cv.add_total_component_balances(custom_molar_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_custom_molar_term_undefined_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 3
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j]

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(custom_molar_term=custom_method)


@pytest.mark.unit
def test_add_total_component_balances_custom_mass_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j] * units.kg / units.s

    mb = m.fs.cv.add_total_component_balances(custom_mass_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_custom_mass_term_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 1
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j]

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_component_balances(custom_mass_term=custom_method)


@pytest.mark.unit
def test_add_total_component_balances_custom_mass_term_mole_flow_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 2
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j] * units.kg / units.s

    for t in m.fs.time:
        m.fs.cv.properties_out[t].mw_comp = Var(
            m.fs.cv.properties_out[t].config.parameters.component_list,
            units=units.kg / units.mol,
        )

    mb = m.fs.cv.add_total_component_balances(custom_mass_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 2

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_custom_mass_term_undefined_basis():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 3
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.cv.config.property_package.component_list
    )

    def custom_method(t, j):
        return m.fs.cv.test_var[t, j]

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_component_balances(custom_mass_term=custom_method)


# -----------------------------------------------------------------------------
# Test add_total_element_balances
@pytest.mark.unit
def test_add_total_element_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_element_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 3

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_element_balances_properties_not_supported():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.element_list)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_total_element_balances()


@pytest.mark.unit
def test_add_total_element_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_element_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 6
    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert isinstance(m.fs.cv.element_holdup, Var)
    assert isinstance(m.fs.cv.element_accumulation, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_element_balances_dynamic_no_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    # Do not add geometry
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances()


@pytest.mark.unit
def test_add_total_element_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_rate_reactions=True)


@pytest.mark.unit
def test_add_total_element_balances_eq_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_total_element_balances_phase_eq():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_phase_equilibrium=True)


@pytest.mark.unit
def test_add_total_element_balances_mass_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_element_balances(has_mass_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 3
    assert isinstance(m.fs.cv.elemental_mass_transfer_term, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_element_balances_custom_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.element_list)

    def custom_method(t, e):
        return m.fs.cv.test_var[t, e] * units.mol / units.s

    mb = m.fs.cv.add_total_element_balances(custom_elemental_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 3

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_element_balances_lineraly_dependent(caplog):
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    # Change elemental composition to introduce dependency
    m.fs.pp.element_comp = {
        "c1": {"H": 0, "He": 0, "Li": 1},
        "c2": {"H": 1, "He": 2, "Li": 0},
    }

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    mb = m.fs.cv.add_total_element_balances()
    # Check that logger message was recorded and has the right level
    msg = (
        "fs.cv detected linearly dependent element balance equations. "
        "Element balances will NOT be written for the following elements: "
        "['He']"
    )
    assert msg in caplog.text
    for record in caplog.records:
        if "['He']" in record.msg:
            assert record.levelno == idaeslog.INFO_LOW

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    for e in mb:
        # H and Li are not lineraly dependent and should have constraints
        assert e in [(0, "H"), (0, "Li")]
        # He is lineraly dependent on H and should be skipped

    assert_units_consistent(m)


# -----------------------------------------------------------------------------
# Test unsupported material balance types
@pytest.mark.unit
def test_add_total_material_balances():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.cv.add_total_material_balances()


# -----------------------------------------------------------------------------
# Test add_energy_balances default
@pytest.mark.unit
def test_add_energy_balances_default_fail():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.pp.default_balance_switch = 2

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_energy_balances(EnergyBalanceType.useDefault)


@pytest.mark.unit
def test_add_energy_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    eb = m.fs.cv.add_energy_balances(EnergyBalanceType.useDefault)

    assert isinstance(eb, Constraint)
    assert len(eb) == 1

    assert_units_consistent(m)


# -----------------------------------------------------------------------------
# Test phase enthalpy balances
@pytest.mark.unit
def test_add_total_enthalpy_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    eb = m.fs.cv.add_total_enthalpy_balances()

    assert isinstance(eb, Constraint)
    assert len(eb) == 1

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_enthalpy_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 2
    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert isinstance(m.fs.cv.energy_holdup, Var)
    assert isinstance(m.fs.cv.energy_accumulation, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_dynamic_no_geometry():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(
        property_package=m.fs.pp, reaction_package=m.fs.rp, dynamic=True
    )

    # Do not add geometry
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_enthalpy_balances()


@pytest.mark.unit
def test_add_total_enthalpy_balances_heat_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_enthalpy_balances(has_heat_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 1
    assert isinstance(m.fs.cv.heat, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_work_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_enthalpy_balances(has_work_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 1
    assert isinstance(m.fs.cv.work, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_enthalpy_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_enthalpy_balances(has_enthalpy_transfer=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 1
    assert isinstance(m.fs.cv.enthalpy_transfer, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_custom_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time)

    def custom_method(t):
        return m.fs.cv.test_var[t] * units.J / units.s

    mb = m.fs.cv.add_total_enthalpy_balances(custom_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 1

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_dh_rxn_no_extents():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_enthalpy_balances(has_heat_of_reaction=True)


@pytest.mark.unit
def test_add_total_enthalpy_balances_dh_rxn_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)
    m.fs.cv.add_phase_component_balances(has_rate_reactions=True)

    m.fs.cv.add_total_enthalpy_balances(has_heat_of_reaction=True)
    assert isinstance(m.fs.cv.heat_of_reaction, Expression)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_dh_rxn_equil_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)
    m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)

    m.fs.cv.add_total_enthalpy_balances(has_heat_of_reaction=True)
    assert isinstance(m.fs.cv.heat_of_reaction, Expression)

    assert_units_consistent(m)


# -----------------------------------------------------------------------------
# Test unsupported energy balance types
@pytest.mark.unit
def test_add_phase_enthalpy_balances():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.cv.add_phase_enthalpy_balances()


@pytest.mark.unit
def test_add_phase_energy_balances():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.cv.add_phase_energy_balances()


@pytest.mark.unit
def test_add_total_energy_balances():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.cv.add_total_energy_balances()


# -----------------------------------------------------------------------------
# Test add total pressure balances
@pytest.mark.unit
def test_add_total_pressure_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    eb = m.fs.cv.add_total_pressure_balances()

    assert isinstance(eb, Constraint)
    assert len(eb) == 1

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_pressure_balances_deltaP():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_pressure_balances(has_pressure_change=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 1
    assert isinstance(m.fs.cv.deltaP, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_pressure_balances_custom_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time)

    def custom_method(t):
        return m.fs.cv.test_var[t] * units.Pa

    mb = m.fs.cv.add_total_pressure_balances(custom_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 1

    assert_units_consistent(m)


# -----------------------------------------------------------------------------
# Test unsupported momentum balance types
@pytest.mark.unit
def test_add_phase_pressure_balances():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.cv.add_phase_pressure_balances()


@pytest.mark.unit
def test_add_phase_momentum_balances():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.cv.add_phase_momentum_balances()


@pytest.mark.unit
def test_add_total_momentum_balances():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(BalanceTypeNotSupportedError):
        m.fs.cv.add_total_momentum_balances()


# -----------------------------------------------------------------------------
# Test model checks, initialize and release_state
@pytest.mark.unit
def test_model_checks():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.model_check()

    for t in m.fs.time:
        assert m.fs.cv.properties_in[t].check is True
        assert m.fs.cv.properties_out[t].check is True
        assert m.fs.cv.reactions[t].check is True


@pytest.mark.unit
def test_initialize():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    f = m.fs.cv.initialize()

    for t in m.fs.time:
        assert m.fs.cv.properties_in[t].init_test is True
        assert m.fs.cv.properties_out[t].init_test is True
        assert m.fs.cv.properties_in[t].hold_state is True
        assert m.fs.cv.properties_out[t].hold_state is False
        assert m.fs.cv.reactions[t].init_test is True

    m.fs.cv.release_state(flags=f)

    for t in m.fs.time:
        assert m.fs.cv.properties_in[t].hold_state is False
        assert m.fs.cv.properties_out[t].hold_state is False


@pytest.mark.unit
def test_get_stream_table_contents():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)

    df = m.fs.cv._get_stream_table_contents()

    assert df.loc["component_flow_phase ('p1', 'c1')"]["In"] == 2
    assert df.loc["component_flow_phase ('p1', 'c2')"]["In"] == 2
    assert df.loc["component_flow_phase ('p2', 'c1')"]["In"] == 2
    assert df.loc["component_flow_phase ('p2', 'c2')"]["In"] == 2
    assert df.loc["pressure"]["In"] == 1e5
    assert df.loc["temperature"]["In"] == 300

    assert df.loc["component_flow_phase ('p1', 'c1')"]["Out"] == 2
    assert df.loc["component_flow_phase ('p1', 'c2')"]["Out"] == 2
    assert df.loc["component_flow_phase ('p2', 'c1')"]["Out"] == 2
    assert df.loc["component_flow_phase ('p2', 'c2')"]["Out"] == 2
    assert df.loc["pressure"]["Out"] == 1e5
    assert df.loc["temperature"]["Out"] == 300


@pytest.mark.unit
def test_get_performance_contents():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)
    m.fs.cv.add_material_balances(
        has_rate_reactions=True,
        has_equilibrium_reactions=True,
        has_phase_equilibrium=True,
        has_mass_transfer=True,
    )
    m.fs.cv.add_energy_balances(
        has_heat_of_reaction=True, has_work_transfer=True, has_heat_transfer=True
    )
    m.fs.cv.add_momentum_balances(has_pressure_change=True)

    dd = m.fs.cv._get_performance_contents()

    assert len(dd) == 3
    for k in dd.keys():
        assert k in ("vars", "exprs", "params")
    assert len(dd["vars"]) == 36
    for k in dd["vars"].keys():
        assert k in [
            "Volume",
            "Heat Transfer",
            "Work Transfer",
            "Pressure Change",
            "Phase Fraction [p1]",
            "Phase Fraction [p2]",
            "Energy Holdup [p1]",
            "Energy Holdup [p2]",
            "Energy Accumulation [p1]",
            "Energy Accumulation [p2]",
            "Material Holdup [p1, c1]",
            "Material Holdup [p1, c2]",
            "Material Holdup [p2, c1]",
            "Material Holdup [p2, c2]",
            "Material Accumulation [p1, c1]",
            "Material Accumulation [p1, c2]",
            "Material Accumulation [p2, c1]",
            "Material Accumulation [p2, c2]",
            "Rate Reaction Generation [p1, c1]",
            "Rate Reaction Generation [p1, c2]",
            "Rate Reaction Generation [p2, c1]",
            "Rate Reaction Generation [p2, c2]",
            "Equilibrium Reaction Generation [p1, c1]",
            "Equilibrium Reaction Generation [p1, c2]",
            "Equilibrium Reaction Generation [p2, c1]",
            "Equilibrium Reaction Generation [p2, c2]",
            "Mass Transfer Term [p1, c1]",
            "Mass Transfer Term [p1, c2]",
            "Mass Transfer Term [p2, c1]",
            "Mass Transfer Term [p2, c2]",
            "Rate Reaction Extent [r1]",
            "Rate Reaction Extent [r2]",
            "Equilibrium Reaction Extent [e1]",
            "Equilibrium Reaction Extent [e2]",
            "Phase Equilibrium Generation [e1]",
            "Phase Equilibrium Generation [e2]",
        ]

    assert len(dd["exprs"]) == 1
    for k in dd["exprs"].keys():
        assert k in ["Heat of Reaction Term"]

    assert len(dd["params"]) == 0


@pytest.mark.unit
def test_get_performance_contents_elemental():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)
    m.fs.cv.add_total_element_balances(has_mass_transfer=True)
    m.fs.cv.add_energy_balances(
        has_heat_of_reaction=False, has_work_transfer=True, has_heat_transfer=True
    )
    m.fs.cv.add_momentum_balances(has_pressure_change=True)

    dd = m.fs.cv._get_performance_contents()

    assert len(dd) == 3
    for k in dd.keys():
        assert k in ("vars", "exprs", "params")
    assert len(dd["vars"]) == 19
    for k in dd["vars"].keys():
        assert k in [
            "Volume",
            "Heat Transfer",
            "Work Transfer",
            "Pressure Change",
            "Phase Fraction [p1]",
            "Phase Fraction [p2]",
            "Energy Holdup [p1]",
            "Energy Holdup [p2]",
            "Energy Accumulation [p1]",
            "Energy Accumulation [p2]",
            "Elemental Holdup [H]",
            "Elemental Holdup [He]",
            "Elemental Holdup [Li]",
            "Elemental Accumulation [H]",
            "Elemental Accumulation [He]",
            "Elemental Accumulation [Li]",
            "Elemental Transfer Term [H]",
            "Elemental Transfer Term [He]",
            "Elemental Transfer Term [Li]",
        ]

    assert len(dd["exprs"]) == 12
    for k in dd["exprs"].keys():
        assert k in [
            "Element Flow In [p1, H]",
            "Element Flow In [p1, He]",
            "Element Flow In [p1, Li]",
            "Element Flow In [p2, H]",
            "Element Flow In [p2, He]",
            "Element Flow In [p2, Li]",
            "Element Flow Out [p1, H]",
            "Element Flow Out [p1, He]",
            "Element Flow Out [p1, Li]",
            "Element Flow Out [p2, H]",
            "Element Flow Out [p2, He]",
            "Element Flow Out [p2, Li]",
        ]

    assert len(dd["params"]) == 0


@pytest.mark.unit
def test_reports():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume0DBlock(property_package=m.fs.pp, reaction_package=m.fs.rp)

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)
    m.fs.cv.add_material_balances(
        has_rate_reactions=True,
        has_equilibrium_reactions=True,
        has_phase_equilibrium=True,
    )
    m.fs.cv.add_energy_balances(has_heat_of_reaction=True, has_heat_transfer=True)
    m.fs.cv.add_momentum_balances(has_pressure_change=True)

    m.fs.cv.report()
