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
Tests for CLC gas phase thermo state block; tests for construction and solve
Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.testing import initialization_tester
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import GasPhaseParameterBlock

# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.fixture(scope="class")
def gas_prop():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # gas properties and state inlet block
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.unit = m.fs.properties.build_state_block(
        default={"parameters": m.fs.properties,
                 "defined_state": True})

    m.fs.unit.flow_mol.fix(1)
    m.fs.unit.temperature.fix(450)
    m.fs.unit.pressure.fix(1.60E5)
    m.fs.unit.mole_frac_comp["CO2"].fix(0.4772)
    m.fs.unit.mole_frac_comp["H2O"].fix(0.0646)
    m.fs.unit.mole_frac_comp["CH4"].fix(0.4582)

    return m


@pytest.mark.unit
def test_build_inlet_state_block(gas_prop):
    assert isinstance(gas_prop.fs.unit.mw, Var)
    assert isinstance(gas_prop.fs.unit.dens_mol, Var)
    assert isinstance(gas_prop.fs.unit.dens_mol_comp, Var)
    assert isinstance(gas_prop.fs.unit.dens_mass, Var)
    assert isinstance(gas_prop.fs.unit.cp_mol_comp, Var)
    assert isinstance(gas_prop.fs.unit.cp_mol, Var)
    assert isinstance(gas_prop.fs.unit.cp_mass, Var)
    assert isinstance(gas_prop.fs.unit.visc_d, Var)
    assert isinstance(gas_prop.fs.unit.enth_mol, Var)
    assert isinstance(gas_prop.fs.unit.entr_mol, Var)


@pytest.mark.unit
def test_setInputs_state_block(gas_prop):
    assert degrees_of_freedom(gas_prop.fs.unit) == 0


@pytest.fixture(scope="class")
def gas_prop_unscaled():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # gas properties and state inlet block
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.unit = m.fs.properties.build_state_block(
        default={"parameters": m.fs.properties,
                 "defined_state": True})

    m.fs.unit.flow_mol.fix(1)
    m.fs.unit.temperature.fix(450)
    m.fs.unit.pressure.fix(1.60E5)
    m.fs.unit.mole_frac_comp["CO2"].fix(0.4772)
    m.fs.unit.mole_frac_comp["H2O"].fix(0.0646)
    m.fs.unit.mole_frac_comp["CH4"].fix(0.4582)

    return m


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unscaled(gas_prop_unscaled):
    initialization_tester(
            gas_prop_unscaled)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve_unscaled(gas_prop_unscaled):

    assert hasattr(gas_prop_unscaled.fs.unit, "mw")
    assert hasattr(gas_prop_unscaled.fs.unit, "dens_mol")
    assert hasattr(gas_prop_unscaled.fs.unit, "dens_mol_comp")
    assert hasattr(gas_prop_unscaled.fs.unit, "dens_mass")
    assert hasattr(gas_prop_unscaled.fs.unit, "visc_d")
    assert hasattr(gas_prop_unscaled.fs.unit, "therm_cond")
    assert hasattr(gas_prop_unscaled.fs.unit, "diffusion_comp")
    assert hasattr(gas_prop_unscaled.fs.unit, "cp_mol_comp")
    assert hasattr(gas_prop_unscaled.fs.unit, "cp_mol")
    assert hasattr(gas_prop_unscaled.fs.unit, "cp_mass")
    assert hasattr(gas_prop_unscaled.fs.unit, "enth_mol")
    assert hasattr(gas_prop_unscaled.fs.unit, "enth_mol_comp")
    assert hasattr(gas_prop_unscaled.fs.unit, "entr_mol")

    assert hasattr(gas_prop_unscaled.fs.unit, "get_material_flow_terms")
    assert hasattr(gas_prop_unscaled.fs.unit, "get_enthalpy_flow_terms")
    assert hasattr(gas_prop_unscaled.fs.unit, "get_material_density_terms")
    assert hasattr(gas_prop_unscaled.fs.unit, "get_energy_density_terms")

    results = solver.solve(gas_prop_unscaled)

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.component
def test_scaling(gas_prop):
    # Calculate scaling factors

    # Construct property methods to build the constraints
    assert hasattr(gas_prop.fs.unit, "mw")
    assert hasattr(gas_prop.fs.unit, "dens_mol")
    assert hasattr(gas_prop.fs.unit, "dens_mol_comp")
    assert hasattr(gas_prop.fs.unit, "dens_mass")
    assert hasattr(gas_prop.fs.unit, "visc_d")
    assert hasattr(gas_prop.fs.unit, "therm_cond")
    assert hasattr(gas_prop.fs.unit, "diffusion_comp")
    assert hasattr(gas_prop.fs.unit, "cp_mol_comp")
    assert hasattr(gas_prop.fs.unit, "cp_mol")
    assert hasattr(gas_prop.fs.unit, "cp_mass")
    assert hasattr(gas_prop.fs.unit, "enth_mol")
    assert hasattr(gas_prop.fs.unit, "enth_mol_comp")
    assert hasattr(gas_prop.fs.unit, "entr_mol")

    # Call flow and density methods to construct flow and density expressions
    for i in gas_prop.fs.unit._params.component_list:
        gas_prop.fs.unit.get_material_flow_terms('Vap', i)
        gas_prop.fs.unit.get_material_density_terms('Vap', i)
    gas_prop.fs.unit.get_enthalpy_flow_terms('Vap')
    gas_prop.fs.unit.get_energy_density_terms('Vap')

    # Calculate scaling factors now that constraints/expressions are built
    iscale.calculate_scaling_factors(gas_prop)

    # Test scaling
    assert (pytest.approx(1e-3, abs=1e-2) ==
            iscale.get_scaling_factor(gas_prop.fs.unit.dens_mol))

    for i, c in gas_prop.fs.unit.material_flow_terms.items():
        assert (pytest.approx(1e-2, abs=1e-2) ==
                iscale.get_scaling_factor(c))
    for i, c in gas_prop.fs.unit.material_density_terms.items():
        assert (pytest.approx(1e-2, abs=1e-3) ==
                iscale.get_scaling_factor(c))
    for i, c in gas_prop.fs.unit.energy_density_terms.items():
        assert (pytest.approx(1e-6, abs=1e-2) ==
                iscale.get_scaling_factor(c))
    for i, c in gas_prop.fs.unit.enthalpy_flow_terms.items():
        assert (pytest.approx(1e-9, abs=1e-2) ==
                iscale.get_scaling_factor(c))

    assert (pytest.approx(1e2, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                gas_prop.fs.unit.mw_eqn))
    assert (pytest.approx(1e-5, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                gas_prop.fs.unit.ideal_gas))
    for i, c in gas_prop.fs.unit.comp_conc_eqn.items():
        assert (pytest.approx(1e-2, abs=1e-2) ==
                iscale.get_constraint_transform_applied_scaling_factor(c))
    for i, c in gas_prop.fs.unit.visc_d_constraint.items():
        assert (pytest.approx(1e5, abs=1e-2) ==
                iscale.get_constraint_transform_applied_scaling_factor(c))
    for i, c in gas_prop.fs.unit.diffusion_comp_constraint.items():
        assert (pytest.approx(1e5, abs=1e-2) ==
                iscale.get_constraint_transform_applied_scaling_factor(c))
    assert (pytest.approx(1, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                gas_prop.fs.unit.therm_cond_constraint))
    for i, c in gas_prop.fs.unit.cp_shomate_eqn.items():
        assert (pytest.approx(1e-6, abs=1e-2) ==
                iscale.get_constraint_transform_applied_scaling_factor(c))
    assert (pytest.approx(1e-6, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                gas_prop.fs.unit.mixture_heat_capacity_eqn))
    assert (pytest.approx(1e-6, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                gas_prop.fs.unit.cp_mass_basis))
    for i, c in gas_prop.fs.unit.enthalpy_shomate_eqn.items():
        assert (pytest.approx(1e-6, abs=1e-2) ==
                iscale.get_constraint_transform_applied_scaling_factor(c))
    assert (pytest.approx(1e-6, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                gas_prop.fs.unit.mixture_enthalpy_eqn))
    assert (pytest.approx(1e-2, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                gas_prop.fs.unit.entropy_correlation))


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(gas_prop):
    initialization_tester(
            gas_prop)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(gas_prop):

    assert hasattr(gas_prop.fs.unit, "mw")
    assert hasattr(gas_prop.fs.unit, "dens_mol")
    assert hasattr(gas_prop.fs.unit, "cp_mol")
    assert hasattr(gas_prop.fs.unit, "visc_d")
    assert hasattr(gas_prop.fs.unit, "enth_mol")
    assert hasattr(gas_prop.fs.unit, "entr_mol")

    results = solver.solve(gas_prop)

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(gas_prop):
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.mw.value)
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.dens_mol.value)
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.cp_mol.value)
    assert (pytest.approx(1e-5, abs=1e-2) ==
            gas_prop.fs.unit.visc_d.value)
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.enth_mol.value)
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.entr_mol.value)


@pytest.mark.component
def test_units_consistent(gas_prop):

    # Construct property methods to build the constraints
    assert hasattr(gas_prop.fs.unit, "mw")
    assert hasattr(gas_prop.fs.unit, "dens_mol")
    assert hasattr(gas_prop.fs.unit, "dens_mol_comp")
    assert hasattr(gas_prop.fs.unit, "dens_mass")
    assert hasattr(gas_prop.fs.unit, "visc_d")
    assert hasattr(gas_prop.fs.unit, "therm_cond")
    assert hasattr(gas_prop.fs.unit, "diffusion_comp")
    assert hasattr(gas_prop.fs.unit, "cp_mol_comp")
    assert hasattr(gas_prop.fs.unit, "cp_mol")
    assert hasattr(gas_prop.fs.unit, "cp_mass")
    assert hasattr(gas_prop.fs.unit, "enth_mol")
    assert hasattr(gas_prop.fs.unit, "enth_mol_comp")
    assert hasattr(gas_prop.fs.unit, "entr_mol")

    # Call flow and density methods to construct flow and density expressions
    for i in gas_prop.fs.unit._params.component_list:
        gas_prop.fs.unit.get_material_flow_terms('Vap', i)
        gas_prop.fs.unit.get_material_density_terms('Vap', i)
    gas_prop.fs.unit.get_enthalpy_flow_terms('Vap')
    gas_prop.fs.unit.get_energy_density_terms('Vap')

    assert_units_consistent(gas_prop)
