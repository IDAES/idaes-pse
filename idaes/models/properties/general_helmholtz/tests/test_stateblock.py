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
"""Test the parts of helmholtz_state that aren't hit by other tests in this."""

import pytest

import pyomo.environ as pyo

import idaes.models.properties.general_helmholtz.helmholtz_state as hstate
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    AmountBasis,
)
from idaes.core import FlowsheetBlock


@pytest.mark.unit
def test_set_not_fixed():
    m = pyo.ConcreteModel()
    m.x = pyo.Var()
    m.y = pyo.Var(initialize=1)
    state = {"x": 5}
    hstate._StateBlock._set_not_fixed(m.x, state, key="x", hold=True)
    assert pyo.value(m.x) == 5
    assert m.x.fixed
    hstate._StateBlock._set_not_fixed(m.y, state, key="y", hold=True)
    assert pyo.value(m.y) == 1


@pytest.mark.unit
def test_misc_methods():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.prop_water = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE
    )
    m.fs.prop_water_mass = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MASS
    )

    m.fs.sb = m.fs.prop_water.build_state_block()
    m.fs.sb_mass = m.fs.prop_water_mass.build_state_block()

    assert m.fs.sb.model_check() is None
    assert "Molar Flow" in m.fs.sb.define_display_vars()
    assert "Mass Flow" in m.fs.sb.define_display_vars()
    assert "T" in m.fs.sb.define_display_vars()
    assert "P" in m.fs.sb.define_display_vars()
    assert "Vapor Fraction" in m.fs.sb.define_display_vars()
    assert "Molar Enthalpy" in m.fs.sb.define_display_vars()

    assert m.fs.sb.flow_mol in m.fs.sb.extensive_state_vars()
    assert m.fs.sb.pressure in m.fs.sb.intensive_state_vars()

    assert "Molar Flow" in m.fs.sb_mass.define_display_vars()
    assert "Mass Flow" in m.fs.sb_mass.define_display_vars()
    assert "T" in m.fs.sb_mass.define_display_vars()
    assert "P" in m.fs.sb_mass.define_display_vars()
    assert "Vapor Fraction" in m.fs.sb_mass.define_display_vars()
    assert "Mass Enthalpy" in m.fs.sb_mass.define_display_vars()
