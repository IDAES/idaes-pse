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
    viscosity_available,
    thermal_conductivity_available,
    surface_tension_available,
    component_registered,
)


@pytest.mark.unit
def test_not_registered():
    assert not viscosity_available("not real")
    assert not thermal_conductivity_available("not real")
    assert not surface_tension_available("not real")
    assert not component_registered("not real")
