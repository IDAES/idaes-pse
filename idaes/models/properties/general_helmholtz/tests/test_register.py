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
"""Test the parts of component registry that aren't usually hit."""

import pytest

from idaes.models.properties.general_helmholtz.components import (
    register_helmholtz_component,
    remove_component,
    registered_components,
    component_registered,
    viscosity_available,
    thermal_conductivity_available,
    surface_tension_available,
    component_registered,
    eos_reference,
    viscosity_reference,
    thermal_conductivity_reference,
    surface_tension_reference,
)


@pytest.mark.unit
def test_not_registered():
    assert not viscosity_available("not real")
    assert not thermal_conductivity_available("not real")
    assert not surface_tension_available("not real")
    assert not component_registered("not real")


@pytest.mark.unit
def test_remove_not_registered():
    with pytest.raises(KeyError):
        remove_component("not real")


@pytest.mark.unit
def test_remove_registered():
    register_helmholtz_component("fooponent")
    assert component_registered("fooponent")
    remove_component("fooponent")
    assert not component_registered("fooponent")


@pytest.mark.unit
def test_registered_components():
    """Use some standard components to test list"""
    assert "h2o" in registered_components()
    assert "co2" in registered_components()
    assert "r134a" in registered_components()


@pytest.mark.unit
def test_register_and_get():
    register_helmholtz_component(
        "fooponent",
        viscosity=True,
        thermal_conductivity=True,
        surface_tension=True,
        eos_ref="cite eos",
        viscosity_ref="cite viscosity",
        thermal_conductivity_ref="cite thermal conductivity",
        surface_tension_ref="cite surface tension",
    )
    assert component_registered("fooponent")
    assert viscosity_available("fooponent")
    assert thermal_conductivity_available("fooponent")
    assert surface_tension_available("fooponent")
    assert eos_reference("fooponent") == "cite eos"
    assert viscosity_reference("fooponent") == "cite viscosity"
    assert thermal_conductivity_reference("fooponent") == "cite thermal conductivity"
    assert surface_tension_reference("fooponent") == "cite surface tension"
    remove_component("fooponent")
    assert not component_registered("fooponent")


@pytest.mark.unit
def test_register_and_get_list_refs():
    register_helmholtz_component(
        "fooponent",
        viscosity=True,
        thermal_conductivity=True,
        surface_tension=True,
        eos_ref=["cite", "eos"],
        viscosity_ref=["cite", "viscosity"],
        thermal_conductivity_ref=["cite", "thermal conductivity"],
        surface_tension_ref=["cite", "surface tension"],
    )
    assert component_registered("fooponent")
    assert viscosity_available("fooponent")
    assert thermal_conductivity_available("fooponent")
    assert surface_tension_available("fooponent")
    assert eos_reference("fooponent") == "cite\neos"
    assert viscosity_reference("fooponent") == "cite\nviscosity"
    assert thermal_conductivity_reference("fooponent") == "cite\nthermal conductivity"
    assert surface_tension_reference("fooponent") == "cite\nsurface tension"
    remove_component("fooponent")
    assert not component_registered("fooponent")


@pytest.mark.unit
def test_not_registered_get():
    assert not component_registered("fooponent")
    assert not viscosity_available("fooponent")
    assert not thermal_conductivity_available("fooponent")
    assert not surface_tension_available("fooponent")
    assert eos_reference("fooponent") is None
    assert viscosity_reference("fooponent") is None
    assert thermal_conductivity_reference("fooponent") is None
    assert surface_tension_reference("fooponent") is None
