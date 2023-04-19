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
"""This module provides functions to register and retrieve component information
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

import os
import idaes

_default_data_dir = os.path.join(idaes.bin_directory, "helm_data")
_default_data_dir = os.path.join(_default_data_dir, "")

_components = {}


class _ComponentStruct(object):
    def __init__(
        self,
        transport_module=None,
        viscosity=False,
        thermal_conductivity=False,
        surface_tension=False,
    ):
        self.transport_module = transport_module
        # these are true if external functions to calculate are available
        self.viscosity = viscosity
        self.thermal_conductivity = thermal_conductivity
        self.surface_tension = surface_tension


def register_helmholtz_component(
    comp_str,
    viscosity=False,
    thermal_conductivity=False,
    surface_tension=False,
):
    _components[comp_str] = _ComponentStruct(
        viscosity=viscosity,
        thermal_conductivity=thermal_conductivity,
        surface_tension=surface_tension,
    )


def viscosity_available(comp_str):
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].viscosity
    return False


def thermal_conductivity_available(comp_str):
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].thermal_conductivity
    return False


def surface_tension_available(comp_str):
    comp_str = comp_str.lower()
    if component_registered(comp_str):
        return _components[comp_str].surface_tension
    return False


def component_registered(comp_str):
    comp_str = comp_str.lower()
    return comp_str in _components
