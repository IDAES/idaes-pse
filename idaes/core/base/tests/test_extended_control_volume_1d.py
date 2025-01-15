#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for ExtendedControlVolumeBlockData.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Constraint, units
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    ExtendedControlVolume1DBlock,
    FlowsheetBlockData,
    declare_process_block_class,
)
from idaes.core.util.exceptions import (
    ConfigurationError,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
)


# -----------------------------------------------------------------------------
# Mockup classes for testing
@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(_Flowsheet, self).build()


@pytest.mark.unit
def test_add_isothermal_constraint():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ExtendedControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    m.fs.cv.add_isothermal_constraint()

    assert isinstance(m.fs.cv.enthalpy_balances, Constraint)
    assert len(m.fs.cv.enthalpy_balances) == 1  # x==0 is skipped
    assert str(m.fs.cv.enthalpy_balances[0, 1].expr) == str(
        m.fs.cv.properties[0, 0].temperature == m.fs.cv.properties[0, 1].temperature
    )
    assert (0, 0) not in m.fs.cv.enthalpy_balances

    assert_units_consistent(m.fs.cv)


@pytest.mark.unit
def test_add_isothermal_constraint_dynamic():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=True, time_set=[0, 1, 2, 3], time_units=units.s)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ExtendedControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    m.fs.cv.add_geometry()
    m.fs.cv.add_state_blocks(has_phase_equilibrium=False)

    m.fs.cv.add_isothermal_constraint()

    assert isinstance(m.fs.cv.enthalpy_balances, Constraint)
    assert len(m.fs.cv.enthalpy_balances) == 4  # x==0 is skipped
    for t in m.fs.time:
        assert str(m.fs.cv.enthalpy_balances[t, 1].expr) == str(
            m.fs.cv.properties[t, 0].temperature == m.fs.cv.properties[t, 1].temperature
        )
        assert (t, 0) not in m.fs.cv.enthalpy_balances

    assert_units_consistent(m.fs.cv)


@pytest.mark.unit
def test_add_isothermal_constraint_heat_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ExtendedControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    with pytest.raises(
        ConfigurationError,
        match="fs.cv: isothermal energy balance option requires that has_heat_transfer is False. "
        "If you are trying to solve for heat duty to achieve isothermal operation, please use "
        "a full energy balance and add a constraint to equate inlet and outlet temperatures.",
    ):
        m.fs.cv.add_isothermal_constraint(has_heat_transfer=True)


@pytest.mark.unit
def test_add_isothermal_constraint_work_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ExtendedControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    with pytest.raises(
        ConfigurationError,
        match="fs.cv: isothermal energy balance option requires that has_work_transfer is False. "
        "If you are trying to solve for work under isothermal operation, please use "
        "a full energy balance and add a constraint to equate inlet and outlet temperatures.",
    ):
        m.fs.cv.add_isothermal_constraint(has_work_transfer=True)


@pytest.mark.unit
def test_add_isothermal_constraint_enthalpy_transfer():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ExtendedControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    with pytest.raises(
        ConfigurationError,
        match="fs.cv: isothermal energy balance option does not support enthalpy transfer.",
    ):
        m.fs.cv.add_isothermal_constraint(has_enthalpy_transfer=True)


@pytest.mark.unit
def test_add_isothermal_constraint_heat_of_rxn():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ExtendedControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    with pytest.raises(
        ConfigurationError,
        match="fs.cv: isothermal energy balance option requires that has_heat_of_reaction is False. "
        "If you are trying to solve for heat duty to achieve isothermal operation, please use "
        "a full energy balance and add a constraint to equate inlet and outlet temperatures.",
    ):
        m.fs.cv.add_isothermal_constraint(has_heat_of_reaction=True)


@pytest.mark.unit
def test_add_isothermal_constraint_custom_term():
    m = ConcreteModel()
    m.fs = Flowsheet(dynamic=False)
    m.fs.pp = PhysicalParameterTestBlock()

    m.fs.cv = ExtendedControlVolume1DBlock(
        property_package=m.fs.pp,
        transformation_method="dae.finite_difference",
        transformation_scheme="BACKWARD",
        finite_elements=10,
    )

    with pytest.raises(
        ConfigurationError,
        match="fs.cv: isothermal energy balance option does not support custom terms.",
    ):
        m.fs.cv.add_isothermal_constraint(custom_term="foo")
