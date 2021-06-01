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
Tests for Ideal + Ideal Liquid (i.e. no activity coefficient) state block;
only tests for construction as parameters need to be provided or estimated
from VLE data to compute the activity coefficients.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import ConcreteModel, TerminationCondition, \
    SolverStatus, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import get_solver

solver = get_solver()

# -----------------------------------------------------------------------------
# Create a flowsheet for test
m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

# vapor-liquid (ideal) - FTPz
m.fs.properties_ideal_vl_FTPz = BTXParameterBlock(
    default={"valid_phase": ('Liq', 'Vap'),
             "activity_coeff_model": "Ideal",
             "state_vars": "FTPz"})
m.fs.state_block_ideal_vl_FTPz = m.fs.properties_ideal_vl_FTPz.build_state_block(
    default={"defined_state": True})

# liquid only (ideal)
m.fs.properties_ideal_l = BTXParameterBlock(default={"valid_phase":
                                                     'Liq',
                                                     "activity_coeff_model":
                                                     "Ideal"})
m.fs.state_block_ideal_l = m.fs.properties_ideal_l.build_state_block(
    default={"has_phase_equilibrium": False,
             "defined_state": True})

# vapour only (ideal)
m.fs.properties_ideal_v = BTXParameterBlock(default={"valid_phase":
                                                     'Vap',
                                                     "activity_coeff_model":
                                                     "Ideal"})
m.fs.state_block_ideal_v = m.fs.properties_ideal_v.build_state_block(
    default={"has_phase_equilibrium": False,
             "defined_state": True})


@pytest.mark.unit
def test_build_inlet_state_block():
    assert len(m.fs.properties_ideal_vl_FTPz.config) == 4

    # vapor-liquid (ideal)
    assert m.fs.properties_ideal_vl_FTPz.config.valid_phase == ('Vap', 'Liq') or \
        m.fs.properties_ideal_vl_FTPz.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs.properties_ideal_vl_FTPz.phase_list) == 2
    assert m.fs.properties_ideal_vl_FTPz.phase_list == ["Liq", "Vap"]
    assert m.fs.state_block_ideal_vl_FTPz.config.defined_state
    assert hasattr(m.fs.state_block_ideal_vl_FTPz, "eq_phase_equilibrium")
    assert not hasattr(m.fs.state_block_ideal_vl_FTPz, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_ideal_vl_FTPz, "eq_mol_frac_out")

    # liquid only (ideal)
    assert len(m.fs.properties_ideal_l.config) == 4

    assert m.fs.properties_ideal_l.config.valid_phase == 'Liq'
    assert len(m.fs.properties_ideal_l.phase_list) == 1
    assert m.fs.properties_ideal_l.phase_list == ["Liq"]
    assert m.fs.state_block_ideal_l.config.defined_state
    assert not hasattr(m.fs.state_block_ideal_l, "eq_phase_equilibrium")
    assert not hasattr(m.fs.state_block_ideal_l, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_ideal_l, "eq_mol_frac_out")

    # vapor only (ideal)
    assert len(m.fs.properties_ideal_v.config) == 4

    assert m.fs.properties_ideal_v.config.valid_phase == 'Vap'
    assert len(m.fs.properties_ideal_v.phase_list) == 1
    assert m.fs.properties_ideal_v.phase_list == ["Vap"]
    assert m.fs.state_block_ideal_v.config.defined_state
    assert not hasattr(m.fs.state_block_ideal_v, "eq_phase_equilibrium")
    assert not hasattr(m.fs.state_block_ideal_v, "eq_activity_coeff")
    assert not hasattr(m.fs.state_block_ideal_v, "eq_mol_frac_out")


@pytest.mark.unit
def test_setInputs_inlet_state_block():

    # vapor-liquid (ideal)
    m.fs.state_block_ideal_vl_FTPz.flow_mol.fix(1)
    m.fs.state_block_ideal_vl_FTPz.temperature.fix(368)
    m.fs.state_block_ideal_vl_FTPz.pressure.fix(101325)
    m.fs.state_block_ideal_vl_FTPz.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block_ideal_vl_FTPz.mole_frac_comp["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_ideal_vl_FTPz) == 0

    # liquid only (ideal)
    m.fs.state_block_ideal_l.flow_mol.fix(1)
    m.fs.state_block_ideal_l.temperature.fix(368)
    m.fs.state_block_ideal_l.pressure.fix(101325)
    m.fs.state_block_ideal_l.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block_ideal_l.mole_frac_comp["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_ideal_l) == 0

    # vapour only (ideal)
    m.fs.state_block_ideal_v.flow_mol.fix(1)
    m.fs.state_block_ideal_v.temperature.fix(368)
    m.fs.state_block_ideal_v.pressure.fix(101325)
    m.fs.state_block_ideal_v.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block_ideal_v.mole_frac_comp["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block_ideal_v) == 0


@pytest.mark.unit
def test_solve():
    # vapor-liquid
    m.fs.state_block_ideal_vl_FTPz.initialize()
    results = solver.solve(m.fs.state_block_ideal_vl_FTPz, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    # Check for VLE results
    assert value(m.fs.state_block_ideal_vl_FTPz.mole_frac_phase_comp['Liq',
                                                                'benzene']) == \
        pytest.approx(0.4121, abs=1e-3)
    assert value(m.fs.state_block_ideal_vl_FTPz.mole_frac_phase_comp['Vap',
                                                                'benzene']) == \
        pytest.approx(0.6339, abs=1e-3)

    # liquid only
    m.fs.state_block_ideal_l.initialize()
    results = solver.solve(m.fs.state_block_ideal_l, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    # Check for results
    assert value(m.fs.state_block_ideal_l.mole_frac_phase_comp['Liq',
                                                          'benzene']) == \
        pytest.approx(0.5, abs=1e-3)
    assert value(m.fs.state_block_ideal_l.mole_frac_phase_comp['Liq',
                                                          'toluene']) == \
        pytest.approx(0.5, abs=1e-3)

    # vapor only
    m.fs.state_block_ideal_v.initialize()
    results = solver.solve(m.fs.state_block_ideal_v, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    # Check for results
    assert value(m.fs.state_block_ideal_v.mole_frac_phase_comp['Vap',
                                                          'benzene']) == \
        pytest.approx(0.5, abs=1e-3)
    assert value(m.fs.state_block_ideal_v.mole_frac_phase_comp['Vap',
                                                          'toluene']) == \
        pytest.approx(0.5, abs=1e-3)


# Create a flowsheet object to test outlet state blocks
m.fs1 = FlowsheetBlock(default={"dynamic": False})

# vapor-liquid (ideal)
m.fs1.properties_ideal_vl = BTXParameterBlock(default={"valid_phase":
                                                       ('Liq', 'Vap'),
                                                       "activity_coeff_model":
                                                       "Ideal"})
m.fs1.state_block_ideal_vl = m.fs1.properties_ideal_vl.build_state_block(
    default={"defined_state": False})

# liquid only (ideal)
m.fs1.properties_ideal_l = BTXParameterBlock(default={"valid_phase":
                                                      "Liq",
                                                      "activity_coeff_model":
                                                      "Ideal"})
m.fs1.state_block_ideal_l = m.fs1.properties_ideal_l.build_state_block(
    default={"has_phase_equilibrium": False,
             "defined_state": False})

# vapour only (ideal)
m.fs1.properties_ideal_v = BTXParameterBlock(default={"valid_phase":
                                                      "Vap",
                                                      "activity_coeff_model":
                                                      "Ideal"})
m.fs1.state_block_ideal_v = m.fs1.properties_ideal_v.build_state_block(
    default={"has_phase_equilibrium": False,
             "defined_state": False})


@pytest.mark.unit
def test_build_outlet_state_block():
    assert len(m.fs1.properties_ideal_vl.config) == 4

    # vapor-liquid (ideal)
    assert m.fs1.properties_ideal_vl.config.valid_phase == ('Vap', 'Liq') or \
        m.fs1.properties_ideal_vl.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs1.properties_ideal_vl.phase_list) == 2
    assert m.fs1.properties_ideal_vl.phase_list == ["Liq", "Vap"]
    assert not m.fs1.state_block_ideal_vl.config.defined_state
    assert hasattr(m.fs1.state_block_ideal_vl, "eq_phase_equilibrium")
    assert not hasattr(m.fs1.state_block_ideal_vl, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_ideal_vl, "eq_mol_frac_out")

    # liquid only (ideal)
    assert len(m.fs1.properties_ideal_l.config) == 4

    assert m.fs1.properties_ideal_l.config.valid_phase == 'Liq'
    assert len(m.fs1.properties_ideal_l.phase_list) == 1
    assert m.fs1.properties_ideal_l.phase_list == ["Liq"]
    assert not m.fs1.state_block_ideal_l.config.defined_state
    assert not hasattr(m.fs1.state_block_ideal_l, "eq_phase_equilibrium")
    assert not hasattr(m.fs1.state_block_ideal_l, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_ideal_l, "eq_mol_frac_out")

    # vapour only (ideal)
    assert len(m.fs1.properties_ideal_v.config) == 4

    assert m.fs1.properties_ideal_v.config.valid_phase == 'Vap'
    assert len(m.fs1.properties_ideal_v.phase_list) == 1
    assert m.fs1.properties_ideal_v.phase_list == ["Vap"]
    assert not m.fs1.state_block_ideal_v.config.defined_state
    assert not hasattr(m.fs1.state_block_ideal_v, "eq_phase_equilibrium")
    assert not hasattr(m.fs1.state_block_ideal_v, "eq_activity_coeff")
    assert hasattr(m.fs1.state_block_ideal_v, "eq_mol_frac_out")


@pytest.mark.unit
def test_setInputs_outlet_state_block():

    # vapor-liquid (ideal)
    m.fs1.state_block_ideal_vl.flow_mol.fix(1)
    m.fs1.state_block_ideal_vl.temperature.fix(368)
    m.fs1.state_block_ideal_vl.pressure.fix(101325)
    m.fs1.state_block_ideal_vl.mole_frac_comp["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs1.state_block_ideal_vl) == 0

    # liquid only (ideal)
    m.fs1.state_block_ideal_l.flow_mol.fix(1)
    m.fs1.state_block_ideal_l.temperature.fix(368)
    m.fs1.state_block_ideal_l.pressure.fix(101325)
    m.fs1.state_block_ideal_l.mole_frac_comp["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs1.state_block_ideal_l) == 0

    # vapour only (ideal)
    m.fs1.state_block_ideal_v.flow_mol.fix(1)
    m.fs1.state_block_ideal_v.temperature.fix(368)
    m.fs1.state_block_ideal_v.pressure.fix(101325)
    m.fs1.state_block_ideal_v.mole_frac_comp["benzene"].fix(0.5)

    assert degrees_of_freedom(m.fs1.state_block_ideal_v) == 0


@pytest.mark.integration
def test_units_consistent():
    assert_units_consistent(m.fs)
    assert_units_consistent(m.fs1)
