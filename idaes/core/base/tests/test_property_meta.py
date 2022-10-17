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
from pyomo.util.check_units import assert_units_equivalent

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


# ---------------------------------------------------------------------------------------------
# Test derived units
derived_quantities = {
    "area": units.m**2,
    "volume": units.m**3,
    "flow_mass": units.kg * units.s**-1,
    "flow_mole": units.mol * units.s**-1,
    "flow_vol": units.m**3 * units.s**-1,
    "flux_mass": (units.kg * units.s**-1 * units.m**-2),
    "flux_mole": (units.mol * units.s**-1 * units.m**-2),
    "flux_energy": (units.kg * units.s**-3),
    "velocity": (units.m * units.s**-1),
    "acceleration": (units.m * units.s**-2),
    "density_mass": (units.kg * units.m**-3),
    "density_mole": (units.mol * units.m**-3),
    "molecular_weight": (units.kg / units.mol),
    "energy": (units.kg * units.m**2 * units.s**-2),
    "energy_mass": (units.m**2 * units.s**-2),
    "energy_mole": (units.kg * units.m**2 * units.s**-2 * units.mol**-1),
    "dynamic_viscosity": (units.kg * units.m**-1 * units.s**-1),
    "entropy": (units.kg * units.m**2 * units.s**-2 * units.K**-1),
    "entropy_mass": (units.m**2 * units.s**-2 * units.K**-1),
    "entropy_mole": (
        units.kg * units.m**2 * units.s**-2 * units.K**-1 * units.mol**-1
    ),
    "power": (units.kg * units.m**2 * units.s**-3),
    "pressure": (units.kg * units.m**-1 * units.s**-2),
    "heat_capacity_mass": (units.m**2 * units.s**-2 * units.K**-1),
    "heat_capacity_mole": (
        units.kg * units.m**2 * units.s**-2 * units.K**-1 * units.mol**-1
    ),
    "heat_transfer_coefficient": (units.kg * units.s**-3 * units.K**-1),
    "thermal_conductivity": (units.kg * units.m * units.s**-3 * units.K**-1),
    "gas_constant": (
        units.kg * units.m**2 * units.s**-2 * units.K**-1 * units.mol**-1
    ),
}


@pytest.mark.unit
@pytest.mark.parametrize("quantity", derived_quantities.keys())
def test_derived_units(unit_set, quantity):
    assert_units_equivalent(unit_set[quantity], derived_quantities[quantity])
