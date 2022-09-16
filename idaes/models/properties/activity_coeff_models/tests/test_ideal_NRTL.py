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
Tests for Ideal + NRTL Liquid activity coefficient state block;
only tests for construction as parameters need to be provided or estimated
from VLE data to compute the activity coefficients.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import ConcreteModel
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom

# -----------------------------------------------------------------------------
# Create a flowsheet for test
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)

# vapor-liquid (NRTL)
m.fs.properties_NRTL_vl = BTXParameterBlock(
    valid_phase=("Liq", "Vap"), activity_coeff_model="NRTL"
)
m.fs.state_block_NRTL_vl = m.fs.properties_NRTL_vl.build_state_block(defined_state=True)

# liquid only (NRTL)
m.fs.properties_NRTL_l = BTXParameterBlock(
    valid_phase="Liq", activity_coeff_model="NRTL"
)
m.fs.state_block_NRTL_l = m.fs.properties_NRTL_l.build_state_block(
    has_phase_equilibrium=False, defined_state=True
)

# vapour only (NRTL)
m.fs.properties_NRTL_v = BTXParameterBlock(
    valid_phase="Vap", activity_coeff_model="NRTL"
)
m.fs.state_block_NRTL_v = m.fs.properties_NRTL_v.build_state_block(
    has_phase_equilibrium=False, defined_state=True
)


@pytest.mark.unit
def test_build_inlet_state_block():
    assert len(m.fs.properties_NRTL_vl.config) == 4

    # vapor-liquid (NRTL)
    assert m.fs.properties_NRTL_vl.config.valid_phase == (
        "Vap",
        "Liq",
    ) or m.fs.properties_NRTL_vl.config.valid_phase == ("Liq", "Vap")
    assert len(m.fs.properties_NRTL_vl.phase_list) == 2
    assert m.fs.properties_NRTL_vl.phase_list == ["Liq", "Vap"]
    assert m.fs.state_block_NRTL_vl.config.defined_state
    assert hasattr(m.fs.state_block_NRTL_vl, "eq_phase_equilibrium")
    assert hasattr(m.fs.state_block_NRTL_vl, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_NRTL_vl, "eq_mol_frac_out")

    # liquid only (NRTL)
    assert len(m.fs.properties_NRTL_l.config) == 4

    assert m.fs.properties_NRTL_l.config.valid_phase == "Liq"
    assert len(m.fs.properties_NRTL_l.phase_list) == 1
    assert m.fs.properties_NRTL_l.phase_list == ["Liq"]
    assert m.fs.state_block_NRTL_l.config.defined_state
    assert not hasattr(m.fs.state_block_NRTL_l, "eq_phase_equilibrium")
    assert not hasattr(m.fs.state_block_NRTL_l, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_NRTL_l, "eq_mol_frac_out")

    # vapor only (NRTL)
    assert len(m.fs.properties_NRTL_v.config) == 4

    assert m.fs.properties_NRTL_v.config.valid_phase == "Vap"
    assert len(m.fs.properties_NRTL_v.phase_list) == 1
    assert m.fs.properties_NRTL_v.phase_list == ["Vap"]
    assert m.fs.state_block_NRTL_v.config.defined_state
    assert not hasattr(m.fs.state_block_NRTL_v, "eq_phase_equilibrium")
    assert not hasattr(m.fs.state_block_NRTL_v, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_NRTL_v, "eq_mol_frac_out")


@pytest.mark.unit
def test_setInputs_inlet_state_block():

    # vapor-liquid (NRTL)
    m.fs.state_block_NRTL_vl.flow_mol.fix(1)
    m.fs.state_block_NRTL_vl.temperature.fix(368)
    m.fs.state_block_NRTL_vl.pressure.fix(101325)
    m.fs.state_block_NRTL_vl.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block_NRTL_vl.mole_frac_comp["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_NRTL_vl) == 6

    # Fix NRTL specific variables and check if DOF is 0
    m.fs.properties_NRTL_vl.alpha.fix()
    m.fs.properties_NRTL_vl.tau.fix()

    assert degrees_of_freedom(m.fs.state_block_NRTL_vl) == 0

    # liquid only (NRTL)
    m.fs.state_block_NRTL_l.flow_mol.fix(1)
    m.fs.state_block_NRTL_l.temperature.fix(368)
    m.fs.state_block_NRTL_l.pressure.fix(101325)
    m.fs.state_block_NRTL_l.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block_NRTL_l.mole_frac_comp["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_NRTL_l) == 0

    # vapour only (NRTL)
    m.fs.state_block_NRTL_v.flow_mol.fix(1)
    m.fs.state_block_NRTL_v.temperature.fix(368)
    m.fs.state_block_NRTL_v.pressure.fix(101325)
    m.fs.state_block_NRTL_v.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block_NRTL_v.mole_frac_comp["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_NRTL_v) == 0


# Create a flowsheet object to test outlet state blocks
m.fs1 = FlowsheetBlock(dynamic=False)

# vapor-liquid (NRTL)
m.fs1.properties_NRTL_vl = BTXParameterBlock(
    valid_phase=("Liq", "Vap"), activity_coeff_model="NRTL"
)
m.fs1.state_block_NRTL_vl = m.fs1.properties_NRTL_vl.build_state_block(
    defined_state=False
)

# liquid only (NRTL)
m.fs1.properties_NRTL_l = BTXParameterBlock(
    valid_phase="Liq", activity_coeff_model="NRTL"
)
m.fs1.state_block_NRTL_l = m.fs1.properties_NRTL_l.build_state_block(
    has_phase_equilibrium=False, defined_state=False
)

# vapour only (NRTL)
m.fs1.properties_NRTL_v = BTXParameterBlock(
    valid_phase="Vap", activity_coeff_model="NRTL"
)
m.fs1.state_block_NRTL_v = m.fs1.properties_NRTL_v.build_state_block(
    has_phase_equilibrium=False, defined_state=False
)


@pytest.mark.unit
def test_build_outlet_state_block():
    assert len(m.fs.properties_NRTL_vl.config) == 4

    # vapor-liquid (NRTL)
    assert m.fs1.properties_NRTL_vl.config.valid_phase == (
        "Vap",
        "Liq",
    ) or m.fs1.properties_NRTL_vl.config.valid_phase == ("Liq", "Vap")
    assert len(m.fs1.properties_NRTL_vl.phase_list) == 2
    assert m.fs1.properties_NRTL_vl.phase_list == ["Liq", "Vap"]
    assert not m.fs1.state_block_NRTL_vl.config.defined_state
    assert hasattr(m.fs1.state_block_NRTL_vl, "eq_phase_equilibrium")
    assert hasattr(m.fs1.state_block_NRTL_vl, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_NRTL_vl, "eq_mol_frac_out")

    # liquid only (NRTL)
    assert len(m.fs1.properties_NRTL_l.config) == 4

    assert m.fs1.properties_NRTL_l.config.valid_phase == "Liq"
    assert len(m.fs1.properties_NRTL_l.phase_list) == 1
    assert m.fs1.properties_NRTL_l.phase_list == ["Liq"]
    assert not m.fs1.state_block_NRTL_l.config.defined_state
    assert not hasattr(m.fs1.state_block_NRTL_l, "eq_phase_equilibrium")
    assert not hasattr(m.fs1.state_block_NRTL_l, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_NRTL_l, "eq_mol_frac_out")

    # vapour only (NRTL)
    assert len(m.fs1.properties_NRTL_v.config) == 4

    assert m.fs1.properties_NRTL_v.config.valid_phase == "Vap"
    assert len(m.fs1.properties_NRTL_v.phase_list) == 1
    assert m.fs1.properties_NRTL_v.phase_list == ["Vap"]
    assert not m.fs1.state_block_NRTL_v.config.defined_state
    assert not hasattr(m.fs1.state_block_NRTL_v, "eq_phase_equilibrium")
    assert not hasattr(m.fs1.state_block_NRTL_v, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_NRTL_v, "eq_mol_frac_out")


@pytest.mark.unit
def test_setInputs_outlet_state_block():

    # vapor-liquid (NRTL)
    m.fs1.state_block_NRTL_vl.flow_mol.fix(1)
    m.fs1.state_block_NRTL_vl.temperature.fix(368)
    m.fs1.state_block_NRTL_vl.pressure.fix(101325)
    m.fs1.state_block_NRTL_vl.mole_frac_comp["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs1.state_block_NRTL_vl) == 6

    # liquid only (NRTL)
    m.fs1.state_block_NRTL_l.flow_mol.fix(1)
    m.fs1.state_block_NRTL_l.temperature.fix(368)
    m.fs1.state_block_NRTL_l.pressure.fix(101325)
    m.fs1.state_block_NRTL_l.mole_frac_comp["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs1.state_block_NRTL_l) == 0

    # vapour only (NRTL)
    m.fs1.state_block_NRTL_v.flow_mol.fix(1)
    m.fs1.state_block_NRTL_v.temperature.fix(368)
    m.fs1.state_block_NRTL_v.pressure.fix(101325)
    m.fs1.state_block_NRTL_v.mole_frac_comp["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs1.state_block_NRTL_v) == 0


@pytest.mark.integration
def test_units_consistent():
    assert_units_consistent(m.fs)
    assert_units_consistent(m.fs1)
