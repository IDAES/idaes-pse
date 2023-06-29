#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for process_meta

Author: Andrew Lee
"""
import pytest

from pyomo.environ import ConcreteModel, units
from pyomo.util.check_units import assert_units_equivalent

from idaes.core.base.property_meta import PropertyClassMetadata, UnitSet
from idaes.core.util.exceptions import PropertyPackageError


@pytest.mark.unit
def test_invalid_require_base_quantity():
    with pytest.raises(
        PropertyPackageError,
        match="Unrecognized units of measurement for quantity TIME \(foo\)",
    ):
        us = UnitSet()
        us.set_units(time="foo")


@pytest.mark.unit
def test_mismatched_length_units():
    with pytest.raises(
        PropertyPackageError,
        match="Invalid units of measurement for quantity LENGTH \(s\). "
        "Please ensure units provided are valid for this quantity.",
    ):
        us = UnitSet()
        us.set_units(length=units.s)


@pytest.mark.unit
def test_mismatched_mass_units():
    with pytest.raises(
        PropertyPackageError,
        match="Invalid units of measurement for quantity MASS \(m\). "
        "Please ensure units provided are valid for this quantity.",
    ):
        us = UnitSet()
        us.set_units(mass=units.m)


@pytest.mark.unit
def test_mismatched_amount_units():
    with pytest.raises(
        PropertyPackageError,
        match="Invalid units of measurement for quantity AMOUNT \(m\). "
        "Please ensure units provided are valid for this quantity.",
    ):
        us = UnitSet()
        us.set_units(amount=units.m)


@pytest.mark.unit
def test_mismatched_temperature_units():
    with pytest.raises(
        PropertyPackageError,
        match="Invalid units of measurement for quantity TEMPERATURE \(m\). "
        "Please ensure units provided are valid for this quantity.",
    ):
        us = UnitSet()
        us.set_units(temperature=units.m)


@pytest.mark.unit
def test_mismatched_current_units():
    with pytest.raises(
        PropertyPackageError,
        match="Invalid units of measurement for quantity CURRENT \(m\). "
        "Please ensure units provided are valid for this quantity.",
    ):
        us = UnitSet()
        us.set_units(current=units.m)


@pytest.mark.unit
def test_mismatched_luminous_intensity_units():
    with pytest.raises(
        PropertyPackageError,
        match="Invalid units of measurement for quantity LUMINOUS_INTENSITY \(m\). "
        "Please ensure units provided are valid for this quantity.",
    ):
        us = UnitSet()
        us.set_units(luminous_intensity=units.m)


@pytest.mark.unit
def test_mismatched_time_units():
    with pytest.raises(
        PropertyPackageError,
        match="Invalid units of measurement for quantity TIME \(m\). "
        "Please ensure units provided are valid for this quantity.",
    ):
        us = UnitSet()
        us.set_units(time=units.m)


@pytest.fixture(scope="module")
def unit_set():
    us = UnitSet()
    us.set_units(
        time=units.s,
        length=units.m,
        mass=units.kg,
        amount=units.mol,
        temperature=units.K,
        current=units.ampere,
        luminous_intensity=units.candela,
    )
    return us


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
    assert unit_set.CURRENT == units.ampere
    assert unit_set["current"] == units.ampere


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
    "voltage": (units.kg * units.m**2 * units.s**-3 * units.ampere**-1),
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


@pytest.mark.unit
def test_add_default_units_extra_arg():
    m = ConcreteModel()

    m.meta_object = PropertyClassMetadata()

    with pytest.raises(
        TypeError,
        match="Unexpected argument for base quantities found when creating "
        "UnitSet. Please ensure that units are only defined for the seven "
        "base quantities.",
    ):
        m.meta_object.add_default_units(
            {"foo": "bar"},
        )


@pytest.mark.unit
def test_add_properties():
    m = ConcreteModel()
    m.meta_object = PropertyClassMetadata()

    m.meta_object.add_properties(
        {
            "flow_vol": {"method": "foo", "supported": True},
            "conc_mol_comp": {"method": "bar", "supported": True},
            "mw_phase": {"required": True, "supported": False},
            "dens_mol_phase_comp": {
                "method": "baz",
                "supported": True,
                "valid_range": (1, 10),
            },
        }
    )

    assert m.meta_object._properties.flow_vol._none.method == "foo"
    assert m.meta_object._properties.flow_vol._none.supported
    assert not m.meta_object._properties.flow_vol._none.required
    assert m.meta_object._properties.flow_vol._none.valid_range is None

    assert m.meta_object._properties.flow_vol._comp.method is None
    assert not m.meta_object._properties.flow_vol._comp.supported
    assert not m.meta_object._properties.flow_vol._comp.required
    assert m.meta_object._properties.flow_vol._comp.valid_range is None

    assert m.meta_object._properties.flow_vol._phase.method is None
    assert not m.meta_object._properties.flow_vol._phase.supported
    assert not m.meta_object._properties.flow_vol._phase.required
    assert m.meta_object._properties.flow_vol._phase.valid_range is None

    assert m.meta_object._properties.flow_vol._phase_comp.method is None
    assert not m.meta_object._properties.flow_vol._phase_comp.supported
    assert not m.meta_object._properties.flow_vol._phase_comp.required
    assert m.meta_object._properties.flow_vol._phase_comp.valid_range is None

    assert m.meta_object._properties.conc_mol._none.method is None
    assert not m.meta_object._properties.conc_mol._none.supported
    assert not m.meta_object._properties.conc_mol._none.required
    assert m.meta_object._properties.conc_mol._none.valid_range is None

    assert m.meta_object._properties.conc_mol._comp.method == "bar"
    assert m.meta_object._properties.conc_mol._comp.supported
    assert not m.meta_object._properties.conc_mol._comp.required
    assert m.meta_object._properties.conc_mol._comp.valid_range is None

    assert m.meta_object._properties.conc_mol._phase.method is None
    assert not m.meta_object._properties.conc_mol._phase.supported
    assert not m.meta_object._properties.conc_mol._phase.required
    assert m.meta_object._properties.conc_mol._phase.valid_range is None

    assert m.meta_object._properties.conc_mol._phase_comp.method is None
    assert not m.meta_object._properties.conc_mol._phase_comp.supported
    assert not m.meta_object._properties.conc_mol._phase_comp.required
    assert m.meta_object._properties.conc_mol._phase_comp.valid_range is None

    assert m.meta_object._properties.mw._none.method is None
    assert not m.meta_object._properties.mw._none.supported
    assert not m.meta_object._properties.mw._none.required
    assert m.meta_object._properties.mw._none.valid_range is None

    assert m.meta_object._properties.mw._comp.method is None
    assert not m.meta_object._properties.mw._comp.supported
    assert not m.meta_object._properties.mw._comp.required
    assert m.meta_object._properties.mw._comp.valid_range is None

    assert m.meta_object._properties.mw._phase.method is None
    assert not m.meta_object._properties.mw._phase.supported
    assert m.meta_object._properties.mw._phase.required
    assert m.meta_object._properties.mw._phase.valid_range is None

    assert m.meta_object._properties.mw._phase_comp.method is None
    assert not m.meta_object._properties.mw._phase_comp.supported
    assert not m.meta_object._properties.mw._phase_comp.required
    assert m.meta_object._properties.mw._phase_comp.valid_range is None

    assert m.meta_object._properties.dens_mol._none.method is None
    assert not m.meta_object._properties.dens_mol._none.supported
    assert not m.meta_object._properties.dens_mol._none.required
    assert m.meta_object._properties.dens_mol._none.valid_range is None

    assert m.meta_object._properties.dens_mol._comp.method is None
    assert not m.meta_object._properties.dens_mol._comp.supported
    assert not m.meta_object._properties.dens_mol._comp.required
    assert m.meta_object._properties.dens_mol._comp.valid_range is None

    assert m.meta_object._properties.dens_mol._phase.method is None
    assert not m.meta_object._properties.dens_mol._phase.supported
    assert not m.meta_object._properties.dens_mol._phase.required
    assert m.meta_object._properties.dens_mol._phase.valid_range is None

    assert m.meta_object._properties.dens_mol._phase_comp.method == "baz"
    assert m.meta_object._properties.dens_mol._phase_comp.supported
    assert not m.meta_object._properties.dens_mol._phase_comp.required
    assert m.meta_object._properties.dens_mol._phase_comp.valid_range == (1, 10)
