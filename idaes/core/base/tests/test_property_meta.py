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
Tests for process_meta

Author: Andrew Lee
"""
import pytest

from pyomo.environ import units

from idaes.core.base.property_meta import UnitSet
from idaes.core.util.exceptions import PropertyPackageError


@pytest.mark.unit
def test_invalid_base_quantity():
    with pytest.raises(
        PropertyPackageError, match="Unrecognized base quantities: \['foo'\]"
    ):
        UnitSet(foo=units.s)


@pytest.mark.unit
def test_invalid_require_base_quantity():
    with pytest.raises(
        PropertyPackageError,
        match="Unrecognized units of measurement for quantity time \(foo\)",
    ):
        UnitSet(time="foo")


@pytest.mark.unit
def test_missing_require_base_quantity():
    with pytest.raises(
        PropertyPackageError,
        match="Units of measurement not provided for base quantity time. "
        "Units must be provided for all base quantities except for current "
        "and luminous intensity.",
    ):
        UnitSet()


@pytest.fixture(scope="module")
def unit_set():
    return UnitSet(
        time=units.s,
        length=units.m,
        mass=units.kg,
        amount=units.mol,
        temperature=units.K,
        current=units.W,
        luminous_intensity=units.candela,
    )


@pytest.mark.unit
def test_time(unit_set):
    assert unit_set.TIME == units.s
    assert unit_set["time"] == units.s


@pytest.mark.unit
def test_length(unit_set):
    assert unit_set.LENGTH == units.m
    assert unit_set["length"] == units.m


@pytest.mark.unit
def test_mass(unit_set):
    assert unit_set.MASS == units.kg
    assert unit_set["mass"] == units.kg


@pytest.mark.unit
def test_amount(unit_set):
    assert unit_set.AMOUNT == units.mol
    assert unit_set["amount"] == units.mol


@pytest.mark.unit
def test_temperature(unit_set):
    assert unit_set.TEMPERATURE == units.K
    assert unit_set["temperature"] == units.K


@pytest.mark.unit
def test_current(unit_set):
    assert unit_set.CURRENT == units.W
    assert unit_set["current"] == units.W


@pytest.mark.unit
def test_luminous_intensity(unit_set):
    assert unit_set.LUMINOUS_INTENSITY == units.candela
    assert unit_set["luminous_intensity"] == units.candela
    # Also check with space
    assert unit_set["luminous intensity"] == units.candela
