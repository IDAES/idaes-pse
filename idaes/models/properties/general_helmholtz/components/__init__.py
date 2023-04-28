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
"""This module provides a list of supported component and Pyomo expressions for
some properties not implemented as external functions.
"""

__author__ = "John Eslick"

from idaes.models.properties.general_helmholtz.components.registry import (
    register_helmholtz_component,
    viscosity_available,
    thermal_conductivity_available,
    surface_tension_available,
    component_registered,
)

register_helmholtz_component(
    "h2o", viscosity=True, thermal_conductivity=True, surface_tension=True
)
register_helmholtz_component(
    "co2", viscosity=True, thermal_conductivity=True, surface_tension=True
)
register_helmholtz_component(
    "r134a", viscosity=True, thermal_conductivity=True, surface_tension=True
)
register_helmholtz_component(
    "r1234ze", viscosity=True, thermal_conductivity=True, surface_tension=False
)
register_helmholtz_component(
    "r227ea", viscosity=False, thermal_conductivity=False, surface_tension=False
)
