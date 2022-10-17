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
from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Param,
    Set,
    units,
    value,
    Var,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.common.config import ConfigBlock
from pyomo.core.base.constraint import _GeneralConstraintData
from idaes.core import (
    ControlVolume1DBlock,
    FlowsheetBlockData,
    declare_process_block_class,
    FlowDirection,
    MaterialBalanceType,
    EnergyBalanceType,
    DistributedVars,
)
from idaes.core.base.control_volume1d import ControlVolume1DBlockData
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
class CVFrameData(ControlVolume1DBlockData):
    def build(self):
        super(ControlVolume1DBlockData, self).build()


# -----------------------------------------------------------------------------
# Test DistributedVars Enum
@pytest.mark.unit
def test_DistributedVars():
    assert len(DistributedVars) == 2


# -----------------------------------------------------------------------------
# Basic tests
@pytest.mark.unit
def test_base_build():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = CVFrame(property_package=m.fs.pp)

    assert len(m.fs.cv.config) == 12
    assert m.fs.cv.config.dynamic is False
    assert m.fs.cv.config.has_holdup is False
    assert m.fs.cv.config.property_package == m.fs.pp
    assert isinstance(m.fs.cv.config.property_package_args, ConfigBlock)
    assert len(m.fs.cv.config.property_package_args) == 0
    assert m.fs.cv.config.reaction_package is None
    assert isinstance(m.fs.cv.config.reaction_package_args, ConfigBlock)
    assert len(m.fs.cv.config.reaction_package_args) == 0
    assert m.fs.cv.config.auto_construct is False
    assert m.fs.cv.config.area_definition == DistributedVars.uniform
    assert m.fs.cv.config.transformation_method is None
    assert m.fs.cv.config.transformation_scheme is None
    assert m.fs.cv.config.finite_elements is None
    assert m.fs.cv.config.collocation_points is None


@pytest.mark.unit
def test_validate_config_args_transformation_method_none():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    with pytest.raises(ConfigurationError):
        m.fs.cv = ControlVolume1DBlock(property_package=m.fs.pp)


@pytest.mark.unit
def test_validate_config_args_transformation_scheme_none():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    with pytest.raises(ConfigurationError):
        m.fs.cv = ControlVolume1DBlock(
            property_package=m.fs.pp, transformation_method="dae.finite_difference"
        )


@pytest.mark.unit
def test_validate_config_args_transformation_scheme_invalid_1():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    with pytest.raises(ConfigurationError):
        m.fs.cv = ControlVolume1DBlock(
            property_package=m.fs.pp,
            transformation_method="dae.finite_difference",
            transformation_scheme="LAGRANGE-RADAU",
        )


@pytest.mark.unit
def test_validate_config_args_transformation_scheme_invalid_2():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    with pytest.raises(ConfigurationError):
        m.fs.cv = ControlVolume1DBlock(
            property_package=m.fs.pp,
            transformation_method="dae.finite_difference",
            transformation_scheme="LAGRANGE-LEGENDRE",
        )


@pytest.mark.unit
def test_validate_config_args_transformation_scheme_invalid_3():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    with pytest.raises(ConfigurationError):
        m.fs.cv = ControlVolume1DBlock(
            property_package=m.fs.pp,
            transformation_method="dae.collocation",
            transformation_scheme="BACKWARD",
        )


@pytest.mark.unit
def test_validate_config_args_transformation_scheme_invalid_4():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    with pytest.raises(ConfigurationError):
        m.fs.cv = ControlVolume1DBlock(
            property_package=m.fs.pp,
            transformation_method="dae.collocation",
            transformation_scheme="FORWARD",
        )


# -----------------------------------------------------------------------------
# Test add_geometry
@pytest.mark.unit
def test_add_geometry_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()

    assert isinstance(m.fs.cv.length_domain, ContinuousSet)
    assert len(m.fs.cv.length_domain) == 2
    assert isinstance(m.fs.cv.area, Var)
    assert len(m.fs.cv.area) == 1.0
    assert m.fs.cv.area.value == 1.0
    assert isinstance(m.fs.cv.length, Var)
    assert len(m.fs.cv.length) == 1.0
    assert m.fs.cv.length.value == 1.0
    assert m.fs.cv._flow_direction == FlowDirection.forward


@pytest.mark.unit
def test_add_geometry_inherited_domain():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.domain = ContinuousSet(bounds=(0, 1))

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.cv.add_geometry(length_domain=m.fs.domain)

    assert m.fs.cv.length_domain == m.fs.domain


@pytest.mark.unit
def test_add_geometry_length_domain_set():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.cv.add_geometry(length_domain_set=[0.0, 0.2, 0.7, 1.0])

    assert len(m.fs.cv.length_domain) == 4
    for p in m.fs.cv.length_domain:
        assert p in [0.0, 0.2, 0.7, 1.0]


@pytest.mark.unit
def test_add_geometry_flow_direction():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.cv.add_geometry(flow_direction=FlowDirection.backward)

    assert m.fs.cv._flow_direction == FlowDirection.backward


@pytest.mark.unit
def test_add_geometry_flow_direction_invalid():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_geometry(flow_direction="foo")


@pytest.mark.unit
def test_add_geometry_discretized_area():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        area_definition=DistributedVars.variant,
    )

    m.fs.cv.add_geometry()

    assert len(m.fs.cv.area) == 2


@pytest.mark.unit
def test_add_geometry_length_var_Var():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.length = Var(initialize=4)

    m.fs.cv.add_geometry(length_var=m.fs.length)

    assert isinstance(m.fs.cv.length_domain, ContinuousSet)
    assert len(m.fs.cv.length_domain) == 2
    assert isinstance(m.fs.cv.area, Var)
    assert len(m.fs.cv.area) == 1.0
    assert m.fs.cv.area.value == 1.0
    assert m.fs.cv.length is m.fs.length
    assert isinstance(m.fs.cv.length, Var)
    assert len(m.fs.cv.length) == 1.0
    assert m.fs.cv.length.value == 4
    assert m.fs.cv._flow_direction == FlowDirection.forward


@pytest.mark.unit
def test_add_geometry_length_var_Param():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.length = Param(initialize=4)

    m.fs.cv.add_geometry(length_var=m.fs.length)

    assert isinstance(m.fs.cv.length_domain, ContinuousSet)
    assert len(m.fs.cv.length_domain) == 2
    assert isinstance(m.fs.cv.area, Var)
    assert len(m.fs.cv.area) == 1.0
    assert m.fs.cv.area.value == 1.0
    assert m.fs.cv.length is m.fs.length
    assert isinstance(m.fs.cv.length, Param)
    assert len(m.fs.cv.length) == 1.0
    assert m.fs.cv.length.value == 4
    assert m.fs.cv._flow_direction == FlowDirection.forward


@pytest.mark.unit
def test_add_geometry_length_var_Expression():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.length = Expression(expr=4)

    m.fs.cv.add_geometry(length_var=m.fs.length)

    assert isinstance(m.fs.cv.length_domain, ContinuousSet)
    assert len(m.fs.cv.length_domain) == 2
    assert isinstance(m.fs.cv.area, Var)
    assert len(m.fs.cv.area) == 1.0
    assert m.fs.cv.area.value == 1.0
    assert m.fs.cv.length is m.fs.length
    assert isinstance(m.fs.cv.length, Expression)
    assert len(m.fs.cv.length) == 1.0
    assert value(m.fs.cv.length) == 4
    assert m.fs.cv._flow_direction == FlowDirection.forward


@pytest.mark.unit
def test_add_geometry_length_var_invalid_type():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    with pytest.raises(
        ConfigurationError,
        match="fs.cv length_var must be a Pyomo Var, Param or " "Expression.",
    ):
        m.fs.cv.add_geometry(length_var="foo")


@pytest.mark.unit
def test_add_geometry_length_var_indexed():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.length = Var([1, 2, 3, 4])

    with pytest.raises(
        ConfigurationError,
        match="fs.cv length_var must be a scalar \(unindexed\) " "component.",
    ):
        m.fs.cv.add_geometry(length_var=m.fs.length)


# -----------------------------------------------------------------------------
# Test apply_transformation
@pytest.mark.unit
def test_apply_transformation_finite_elements_none():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
    )

    m.fs.cv.add_geometry()
    with pytest.raises(ConfigurationError):
        m.fs.cv.apply_transformation()


@pytest.mark.unit
def test_apply_transformation_collocation_points_none():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.collocation",
        transformation_scheme="LAGRANGE-RADAU",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    with pytest.raises(ConfigurationError):
        m.fs.cv.apply_transformation()


@pytest.mark.unit
def test_apply_transformation_BFD_10():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.apply_transformation()

    assert len(m.fs.cv.length_domain) == 11


@pytest.mark.unit
def test_apply_transformation_FFD_12():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="FORWARD",
        finite_elements=12,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.apply_transformation()

    assert len(m.fs.cv.length_domain) == 13


@pytest.mark.unit
def test_apply_transformation_Lagrange_Radau_8_3():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.collocation",
        transformation_scheme="LAGRANGE-RADAU",
        finite_elements=8,
        collocation_points=3,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.apply_transformation()

    assert len(m.fs.cv.length_domain) == 25


@pytest.mark.unit
def test_apply_transformation_Lagrange_Legendre_3_7():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.collocation",
        transformation_scheme="LAGRANGE-LEGENDRE",
        finite_elements=9,
        collocation_points=4,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.apply_transformation()

    assert len(m.fs.cv.length_domain) == 46


@pytest.mark.unit
def test_apply_transformation_external_domain():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)

    m.fs.cset = ContinuousSet(bounds=(0, 1))

    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry(length_domain=m.fs.cset)
    with pytest.raises(ConfigurationError):
        m.fs.cv.apply_transformation()


# -----------------------------------------------------------------------------
# Test add_state_blocks
@pytest.mark.unit
def test_add_state_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    assert hasattr(m.fs.cv, "properties")
    assert len(m.fs.cv.properties) == 2

    for x in m.fs.cv.length_domain:
        assert len(m.fs.cv.properties[0, x].config) == 3
        if x == 0:
            assert m.fs.cv.properties[0, x].config.defined_state is True
        else:
            assert m.fs.cv.properties[0, x].config.defined_state is False
        assert m.fs.cv.properties[0, x].config.has_phase_equilibrium is False
        assert m.fs.cv.properties[0, x].config.parameters == m.fs.pp


@pytest.mark.unit
def test_add_state_block_forward_flow():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(
        information_flow=FlowDirection.forward, has_phase_equilibrium=False
    )

    assert m.fs.cv.properties[0, 0].config.defined_state is True
    assert m.fs.cv.properties[0, 1].config.defined_state is False


@pytest.mark.unit
def test_add_state_block_backward_flow():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(
        information_flow=FlowDirection.backward, has_phase_equilibrium=False
    )

    assert m.fs.cv.properties[0, 0].config.defined_state is False
    assert m.fs.cv.properties[0, 1].config.defined_state is True


@pytest.mark.unit
def test_add_state_blocks_has_phase_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)

    for x in m.fs.cv.length_domain:
        assert m.fs.cv.properties[0, x].config.has_phase_equilibrium is True


@pytest.mark.unit
def test_add_state_blocks_no_has_phase_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_state_blocks()


@pytest.mark.unit
def test_add_state_blocks_custom_args():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        property_package_args={"test": "test"},
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    for x in m.fs.cv.length_domain:
        assert len(m.fs.cv.properties[0, x].config) == 4
        assert m.fs.cv.properties[0, x].config.test == "test"


# -----------------------------------------------------------------------------
# Test add_reaction_blocks
@pytest.mark.unit
def test_add_reaction_blocks():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    assert hasattr(m.fs.cv, "reactions")
    assert len(m.fs.cv.reactions) == 2
    assert len(m.fs.cv.reactions[0, 0].config) == 3
    assert m.fs.cv.reactions[0, 0].config.state_block == m.fs.cv.properties
    assert m.fs.cv.reactions[0, 0].state_ref == m.fs.cv.properties[0, 0]
    assert m.fs.cv.reactions[0, 0].config.has_equilibrium is False
    assert m.fs.cv.reactions[0, 0].config.parameters == m.fs.rp


@pytest.mark.unit
def test_add_reaction_blocks_has_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    assert m.fs.cv.reactions[0, 0].config.has_equilibrium is True


@pytest.mark.unit
def test_add_reaction_blocks_no_has_equilibrium():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_reaction_blocks()


@pytest.mark.unit
def test_add_reaction_blocks_custom_args():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        reaction_package_args={"test1": 1},
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    assert m.fs.cv.reactions[0, 0].config.test1 == 1


# -----------------------------------------------------------------------------
# Test _add_phase_fractions
@pytest.mark.unit
def test_add_phase_fractions():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv._add_phase_fractions()

    assert isinstance(m.fs.cv.phase_fraction, Var)
    assert len(m.fs.cv.phase_fraction) == 4
    assert isinstance(m.fs.cv.sum_of_phase_fractions, Constraint)


@pytest.mark.unit
def test_add_phase_fractions_single_phase():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.del_component(m.fs.pp.phase_list)
    m.fs.pp.phase_list = Set(initialize=["p1"])

    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv._add_phase_fractions()

    assert isinstance(m.fs.cv.phase_fraction, Expression)
    assert len(m.fs.cv.phase_fraction) == 2
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            for j in m.fs.pp.component_list:
                assert m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=False) == 1


@pytest.mark.unit
def test_rxn_rate_conv_property_basis_other():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 3
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            for j in m.fs.pp.component_list:
                with pytest.raises(ConfigurationError):
                    m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_reaction_basis_other():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.rp.basis_switch = 3

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            for j in m.fs.pp.component_list:
                with pytest.raises(ConfigurationError):
                    m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_both_molar():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            for j in m.fs.pp.component_list:
                assert m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True) == 1


@pytest.mark.unit
def test_rxn_rate_conv_both_mass():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 2
    m.fs.rp.basis_switch = 2

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            for j in m.fs.pp.component_list:
                assert m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True) == 1


@pytest.mark.unit
def test_rxn_rate_conv_mole_mass_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 1
    m.fs.rp.basis_switch = 2

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            for j in m.fs.pp.component_list:
                with pytest.raises(PropertyNotSupportedError):
                    m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_mass_mole_no_mw():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 2
    m.fs.rp.basis_switch = 1

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            for j in m.fs.pp.component_list:
                with pytest.raises(PropertyNotSupportedError):
                    m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True)


@pytest.mark.unit
def test_rxn_rate_conv_mole_mass():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 1
    m.fs.rp.basis_switch = 2

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = {"c1": 2, "c2": 3}
            for j in m.fs.pp.component_list:
                assert (
                    m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True)
                    == 1 / m.fs.cv.properties[t, x].mw_comp[j]
                )


@pytest.mark.unit
def test_rxn_rate_conv_mass_mole():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.basis_switch = 2
    m.fs.rp.basis_switch = 1

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = {"c1": 2, "c2": 3}
            for j in m.fs.pp.component_list:
                assert (
                    m.fs.cv._rxn_rate_conv(t, x, j, has_rate_reactions=True)
                    == m.fs.cv.properties[t, x].mw_comp[j]
                )


# -----------------------------------------------------------------------------
# Test add_material_balances default
@pytest.mark.unit
def test_add_material_balances_default_fail():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_material_balances(MaterialBalanceType.useDefault)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            with pytest.raises(KeyError):
                assert m.fs.cv.material_balances[0, 0, p, j]
            assert type(m.fs.cv.material_balances[0, 1, p, j]) is _GeneralConstraintData

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    units = m.fs.cv.config.property_package.get_metadata().get_derived_units
    pp_units = units("flow_mass") / units("length")  # basis 2 is mass
    rp_units = units("flow_mole") / units("length")  # basis 1 is molar

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    # add molecular weight variable to each time point, using correct units
    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = Var(
                m.fs.cv.properties[t, x].config.parameters.component_list,
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    units = m.fs.cv.config.property_package.get_metadata().get_derived_units
    pp_units = units("flow_mole") / units("length")  # basis 2 is molar
    rp_units = units("flow_mass") / units("length")  # basis 1 is mass

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    # add molecular weight variable to each time point, using correct units
    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = Var(
                m.fs.cv.properties[t, x].config.parameters.component_list,
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            with pytest.raises(KeyError):
                assert m.fs.cv.material_balances[0, 0, p, j]
            assert type(m.fs.cv.material_balances[0, 1, p, j]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_default_FFD():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="FORWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    for p in m.fs.pp.phase_list:
        for j in m.fs.pp.component_list:
            with pytest.raises(KeyError):
                assert m.fs.cv.material_balances[0, 1, p, j]
            assert type(m.fs.cv.material_balances[0, 0, p, j]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_distrubuted_area():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        area_definition=DistributedVars.variant,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        dynamic=True,
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
def test_add_phase_component_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=True)

    with pytest.raises(PropertyNotSupportedError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_phase_component_balances_eq_rxns_no_ReactionBlock():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_phase_component_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_phase_component_balances_phase_eq():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_phase_component_balances(has_phase_equilibrium=True)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4
    assert isinstance(m.fs.cv.phase_equilibrium_generation, Var)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_phase_eq_not_active():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
        return m.fs.cv.test_var[t, p, j] * units.mol / units.s / units.m

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
        return m.fs.cv.test_var[t, p, j] * units.mol / units.s / units.m

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = Var(
                m.fs.cv.properties[t, x].config.parameters.component_list,
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
        return m.fs.cv.test_var[t, p, j] * units.kg / units.s / units.m

    mb = m.fs.cv.add_phase_component_balances(custom_mass_term=custom_method)

    assert isinstance(mb, Constraint)
    assert len(mb) == 4

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_phase_component_balances_custom_mass_term_no_mw_comp():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.pp.basis_switch = 1
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
        return m.fs.cv.test_var[t, p, j] * units.kg / units.s / units.m

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = Var(
                m.fs.cv.properties[t, x].config.parameters.component_list,
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(
        m.fs.cv.flowsheet().time, m.fs.pp.phase_list, m.fs.pp.component_list
    )

    def custom_method(t, x, p, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 2

    for j in m.fs.pp.component_list:
        with pytest.raises(KeyError):
            assert m.fs.cv.material_balances[0, 0, j]
        assert type(m.fs.cv.material_balances[0, 1, j]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_default_FFD():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="FORWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_component_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 2

    for j in m.fs.pp.component_list:
        with pytest.raises(KeyError):
            assert m.fs.cv.material_balances[0, 1, j]
        assert type(m.fs.cv.material_balances[0, 0, j]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_component_balances_distrubuted_area():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        area_definition=DistributedVars.variant,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        dynamic=True,
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
def test_add_total_component_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
        return m.fs.cv.test_var[t, j] * units.mol / units.s / units.m

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
        return m.fs.cv.test_var[t, j] * units.mol / units.s / units.m

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = Var(
                m.fs.cv.properties[t, x].config.parameters.component_list,
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
        return m.fs.cv.test_var[t, j] * units.kg / units.s / units.m

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
        return m.fs.cv.test_var[t, j] * units.kg / units.s / units.m

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            m.fs.cv.properties[t, x].mw_comp = Var(
                m.fs.cv.properties[t, x].config.parameters.component_list,
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.component_list)

    def custom_method(t, x, j):
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_element_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 3

    for j in m.fs.pp.element_list:
        with pytest.raises(KeyError):
            assert m.fs.cv.element_balances[0, 0, j]
        assert type(m.fs.cv.element_balances[0, 1, j]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_element_balances_default_FFD():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="FORWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_element_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 3

    for j in m.fs.pp.element_list:
        with pytest.raises(KeyError):
            assert m.fs.cv.element_balances[0, 1, j]
        assert type(m.fs.cv.element_balances[0, 0, j]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_element_balances_distrubuted_area():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        area_definition=DistributedVars.variant,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        dynamic=True,
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
def test_add_total_element_balances_rate_rxns():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_equilibrium_reactions=True)


@pytest.mark.unit
def test_add_total_element_balances_phase_eq():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    with pytest.raises(ConfigurationError):
        m.fs.cv.add_total_element_balances(has_phase_equilibrium=True)


@pytest.mark.unit
def test_add_total_element_balances_mass_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time, m.fs.pp.element_list)

    def custom_method(t, x, e):
        return m.fs.cv.test_var[t, e] * units.mol / units.s / units.m

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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
    for i in mb:
        # Should be no constraints at x = 0
        # H and Li are not lineraly dependent and should have constraints
        # He is lineraly dependent on H and should be skipped
        assert i in [(0, 1, "H"), (0, 1, "Li")]

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    eb = m.fs.cv.add_energy_balances(EnergyBalanceType.useDefault)

    assert isinstance(eb, Constraint)
    assert len(eb) == 1
    assert isinstance(m.fs.cv._enthalpy_flow, Var)
    assert isinstance(m.fs.cv.enthalpy_flow_linking_constraint, Constraint)
    assert isinstance(m.fs.cv.enthalpy_flow_dx, DerivativeVar)

    with pytest.raises(KeyError):
        assert m.fs.cv.enthalpy_balances[0, 0]
    assert type(m.fs.cv.enthalpy_balances[0, 1]) is _GeneralConstraintData

    assert_units_consistent(m)


# -----------------------------------------------------------------------------
# Test phase enthalpy balances
@pytest.mark.unit
def test_add_total_enthalpy_balances_default():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    eb = m.fs.cv.add_total_enthalpy_balances()

    assert isinstance(eb, Constraint)
    assert len(eb) == 1
    assert isinstance(m.fs.cv._enthalpy_flow, Var)
    assert isinstance(m.fs.cv.enthalpy_flow_linking_constraint, Constraint)
    assert isinstance(m.fs.cv.enthalpy_flow_dx, DerivativeVar)

    with pytest.raises(KeyError):
        assert m.fs.cv.enthalpy_balances[0, 0]
    assert type(m.fs.cv.enthalpy_balances[0, 1]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_default_FFD():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="FORWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    eb = m.fs.cv.add_total_enthalpy_balances()

    assert isinstance(eb, Constraint)
    assert len(eb) == 1
    assert isinstance(m.fs.cv._enthalpy_flow, Var)
    assert isinstance(m.fs.cv.enthalpy_flow_linking_constraint, Constraint)
    assert isinstance(m.fs.cv.enthalpy_flow_dx, DerivativeVar)

    with pytest.raises(KeyError):
        assert m.fs.cv.enthalpy_balances[0, 1]
    assert type(m.fs.cv.enthalpy_balances[0, 0]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_distributed_area():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        area_definition=DistributedVars.variant,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    eb = m.fs.cv.add_total_enthalpy_balances()

    assert isinstance(eb, Constraint)
    assert len(eb) == 1
    assert isinstance(m.fs.cv._enthalpy_flow, Var)
    assert isinstance(m.fs.cv.enthalpy_flow_linking_constraint, Constraint)
    assert isinstance(m.fs.cv.enthalpy_flow_dx, DerivativeVar)

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_enthalpy_balances_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
        dynamic=True,
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
def test_add_total_enthalpy_balances_heat_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time)

    def custom_method(t, x):
        return m.fs.cv.test_var[t] * units.J / units.s / units.m

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_pressure_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 1
    assert isinstance(m.fs.cv.pressure, Var)  # Reference to state block pressure
    assert isinstance(m.fs.cv.pressure_dx, DerivativeVar)

    with pytest.raises(KeyError):
        assert m.fs.cv.pressure_balance[0, 0]
    assert type(m.fs.cv.pressure_balance[0, 1]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_pressure_balances_default_FFD():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="FORWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    mb = m.fs.cv.add_total_pressure_balances()

    assert isinstance(mb, Constraint)
    assert len(mb) == 1
    assert isinstance(m.fs.cv.pressure, Var)  # Reference to state block pressure
    assert isinstance(m.fs.cv.pressure_dx, DerivativeVar)

    with pytest.raises(KeyError):
        assert m.fs.cv.pressure_balance[0, 1]
    assert type(m.fs.cv.pressure_balance[0, 0]) is _GeneralConstraintData

    assert_units_consistent(m)


@pytest.mark.unit
def test_add_total_pressure_balances_deltaP():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.test_var = Var(m.fs.cv.flowsheet().time)

    def custom_method(t, x):
        return m.fs.cv.test_var[t] * units.Pa / units.m

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

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

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.model_check()

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            assert m.fs.cv.properties[t, x].check is True
            assert m.fs.cv.reactions[t, x].check is True


@pytest.mark.unit
def test_initialize():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=True)
    m.fs.cv.add_reaction_blocks(has_equilibrium=False)

    m.fs.cv.initialize()

    for t in m.fs.time:
        for x in m.fs.cv.length_domain:
            assert m.fs.cv.properties[t, x].init_test is True
            assert m.fs.cv.reactions[t, x].init_test is True


@pytest.mark.unit
def test_report():
    # Test that calling report method on a 1D control volume returns a
    # NotImplementedError to inform the user that we don't support reports
    # on 1D models yet. This is because it is difficult to concisely report
    # distributed data.
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()
    m.fs.rp = ReactionParameterTestBlock(property_package=m.fs.pp)
    m.fs.pp.del_component(m.fs.pp.phase_equilibrium_idx)

    m.fs.cv = ControlVolume1DBlock(
        property_package=m.fs.pp,
        reaction_package=m.fs.rp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    with pytest.raises(NotImplementedError):
        m.fs.cv.report()
