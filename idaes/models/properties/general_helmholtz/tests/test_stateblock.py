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
"""Test the parts of helmholtz_state that aren't hit by other tests in this."""

import re

import pytest

import pyomo.environ as pyo

import idaes.models.properties.general_helmholtz.helmholtz_state as hstate
from idaes.models.properties.general_helmholtz import (
    HelmholtzParameterBlock,
    AmountBasis,
    StateVars,
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
def test_scaler_object_mole_basis():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.prop_water = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE
    )

    m.fs.sb = m.fs.prop_water.build_state_block([0])

    assert m.fs.sb[0].default_scaler is hstate.HelmholtzEoSScaler
    scaler_obj = m.fs.sb[0].default_scaler()

    with pytest.raises(
        ValueError,
        match=re.escape(
            "This scaler requires the user to provide a default scaling factor for fs.sb[0].flow_mol, but no default scaling factor was set."
        ),
    ):
        scaler_obj.scale_model(m.fs.sb[0])

    scaler_obj.default_scaling_factors["flow_mol"] = 1 / 137
    scaler_obj.scale_model(m.fs.sb[0])

    sb = m.fs.sb[0]
    assert len(sb.scaling_factor) == 3

    # Variables
    assert sb.scaling_factor[sb.flow_mol] == 1 / 137
    assert sb.scaling_factor[sb.pressure] == 1e-6
    assert sb.scaling_factor[sb.enth_mol] == 1e-3

    # No constraints

    # Expressions
    assert len(sb.scaling_hint) == 49
    # Test a subset of values
    assert sb.scaling_hint[sb.temperature] == 1e-1
    assert sb.scaling_hint[sb.vapor_frac] == 10
    assert sb.scaling_hint[sb.enth_mol_sat_phase["Liq"]] == 1e-2
    assert sb.scaling_hint[sb.entr_mol_phase["Vap"]] == 1e-1
    assert sb.scaling_hint[sb.cp_mol_phase["Vap"]] == 1e-2
    assert sb.scaling_hint[sb.cv_mass_phase["Liq"]] == 1e-3
    assert sb.scaling_hint[sb.dens_mol_phase["Vap"]] == 1
    assert sb.scaling_hint[sb.energy_internal_mass] == 1e-3
    assert sb.scaling_hint[sb.flow_vol] == pytest.approx(
        1 / (137 * 18.015 / 1000 * pyo.sqrt(5 * 1000)), rel=1e-3
    )

    # No constraints, so we cannot test condition number


@pytest.mark.unit
def test_scaler_object_mass_basis():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.prop_water = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MASS
    )

    m.fs.sb = m.fs.prop_water.build_state_block([0])

    assert m.fs.sb[0].default_scaler is hstate.HelmholtzEoSScaler
    scaler_obj = m.fs.sb[0].default_scaler()

    with pytest.raises(
        ValueError,
        match=re.escape(
            "This scaler requires the user to provide a default scaling factor for fs.sb[0].flow_mass, but no default scaling factor was set."
        ),
    ):
        scaler_obj.scale_model(m.fs.sb[0])

    scaler_obj.default_scaling_factors["flow_mass"] = 1 / 137
    scaler_obj.scale_model(m.fs.sb[0])

    from idaes.core.scaling import report_scaling_factors

    report_scaling_factors(m.fs.sb)
    sb = m.fs.sb[0]
    assert len(sb.scaling_factor) == 3

    # Variables
    assert sb.scaling_factor[sb.flow_mass] == 1 / 137
    assert sb.scaling_factor[sb.pressure] == 1e-6
    assert sb.scaling_factor[sb.enth_mass] == 1e-3

    # No constraints

    # Expressions
    assert len(sb.scaling_hint) == 49
    # Test a subset of values
    assert sb.scaling_hint[sb.temperature] == 1e-1
    assert sb.scaling_hint[sb.vapor_frac] == 10
    assert sb.scaling_hint[sb.enth_mol_sat_phase["Liq"]] == 1e-2
    assert sb.scaling_hint[sb.entr_mol_phase["Vap"]] == 1e-1
    assert sb.scaling_hint[sb.cp_mol_phase["Vap"]] == 1e-2
    assert sb.scaling_hint[sb.cv_mass_phase["Liq"]] == 1e-3
    assert sb.scaling_hint[sb.dens_mol_phase["Vap"]] == 1
    assert sb.scaling_hint[sb.energy_internal_mass] == 1e-3
    assert sb.scaling_hint[sb.flow_vol] == pytest.approx(
        1 / (137 * pyo.sqrt(5 * 1000)), rel=1e-3
    )

    # No constraints, so we cannot test condition number


@pytest.mark.unit
def test_scaler_object_mole_basis_TPX():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.prop_water = HelmholtzParameterBlock(
        pure_component="h2o", amount_basis=AmountBasis.MOLE, state_vars=StateVars.TPX
    )

    m.fs.sb = m.fs.prop_water.build_state_block([0])

    assert m.fs.sb[0].default_scaler is hstate.HelmholtzEoSScaler
    scaler_obj = m.fs.sb[0].default_scaler()

    scaler_obj.default_scaling_factors["flow_mol"] = 1 / 137
    scaler_obj.scale_model(m.fs.sb[0])

    sb = m.fs.sb[0]
    assert len(sb.scaling_factor) == 6

    # Variables
    assert sb.scaling_factor[sb.flow_mol] == 1 / 137
    assert sb.scaling_factor[sb.pressure] == 1e-6
    assert sb.scaling_factor[sb.temperature] == 1e-1
    assert sb.scaling_factor[sb.vapor_frac] == 10

    # Constraints
    assert sb.scaling_factor[sb.eq_complementarity] == pytest.approx(1e-7)
    assert sb.scaling_factor[sb.eq_sat] == pytest.approx(1e-9)

    # Expressions
    assert len(sb.scaling_hint) == 48
    # Test a subset of values
    assert sb.scaling_hint[sb.enth_mol] == 1e-3
    assert sb.scaling_hint[sb.enth_mol_sat_phase["Liq"]] == 1e-2
    assert sb.scaling_hint[sb.entr_mol_phase["Vap"]] == 1e-1
    assert sb.scaling_hint[sb.cp_mol_phase["Vap"]] == 1e-2
    assert sb.scaling_hint[sb.cv_mass_phase["Liq"]] == 1e-3
    assert sb.scaling_hint[sb.dens_mol_phase["Vap"]] == 1
    assert sb.scaling_hint[sb.energy_internal_mass] == 1e-3
    assert sb.scaling_hint[sb.flow_vol] == pytest.approx(
        1 / (137 * 18.015 / 1000 * pyo.sqrt(5 * 1000)), rel=1e-3
    )

    # Only a single active constraint, the condition number
    # is trivially 1


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
