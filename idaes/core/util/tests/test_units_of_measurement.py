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
Tests for units of measurement utility functions
"""
import pytest

from pyomo.environ import as_quantity, units

from idaes.core.util.units_of_measurement import convert_quantity_to_reporting_units


@pytest.mark.unit
def test_convert_quantity_pressure():
    q = as_quantity(1 * units.atm)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == 101325
    assert str(q2.u) == "pascal"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_energy():
    q = as_quantity(1 * units.BTU)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == 1055.056
    assert str(q2.u) == "joule"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_energy_mass():
    q = as_quantity(1 * units.BTU / units.pound)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == pytest.approx(1055.056 / 0.453592, rel=1e-6)
    assert str(q2.u) == "joule / kilogram"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_energy_mole():
    q = as_quantity(1 * units.BTU / units.mol)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == 1055.056
    assert str(q2.u) == "joule / mole"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_entropy():
    q = as_quantity(1 * units.BTU / units.degR)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == pytest.approx(1055.056 * 1.8, rel=1e-6)
    assert str(q2.u) == "joule / kelvin"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_entropy_mass():
    q = as_quantity(1 * units.BTU / units.pound / units.degR)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == pytest.approx(1055.056 / 0.453592 * 1.8, rel=1e-6)
    assert str(q2.u) == "joule / kelvin / kilogram"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_entropy_mole():
    q = as_quantity(1 * units.BTU / units.mol / units.degR)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == pytest.approx(1055.056 * 1.8, rel=1e-6)
    assert str(q2.u) == "joule / kelvin / mole"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_power():
    q = as_quantity(1 * units.BTU / units.hour)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == pytest.approx(0.2930711, rel=1e-6)
    assert str(q2.u) == "watt"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_other():
    q = as_quantity(1 * units.foot**0.5)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == 0.3048**0.5
    assert str(q2.u) == "meter ** 0.5"
    assert q is not q2


@pytest.mark.unit
def test_convert_quantity_dimensionless():
    q = as_quantity(1 * units.dimensionless)

    q2 = convert_quantity_to_reporting_units(q)

    assert q2.m == 1
    assert str(q2.u) == "dimensionless"
    assert q is q2
