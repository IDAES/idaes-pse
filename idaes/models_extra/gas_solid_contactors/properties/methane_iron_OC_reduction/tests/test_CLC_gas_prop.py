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

from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)

from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_phase_thermo import (
    GasPhaseParameterBlock,
)

# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.fixture(scope="class")
def gas_prop():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # gas properties and state inlet block
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.unit = m.fs.properties.build_state_block(
        [0], parameters=m.fs.properties, defined_state=True
    )

    m.fs.unit[0].flow_mol.fix(1)
    m.fs.unit[0].temperature.fix(450)
    m.fs.unit[0].pressure.fix(1.60e5)
    m.fs.unit[0].mole_frac_comp["CO2"].fix(0.4772)
    m.fs.unit[0].mole_frac_comp["H2O"].fix(0.0646)
    m.fs.unit[0].mole_frac_comp["CH4"].fix(0.4582)

    return m


@pytest.mark.unit
def test_build_inlet_state_block(gas_prop):
    assert isinstance(gas_prop.fs.unit[0].mw, Var)
    assert isinstance(gas_prop.fs.unit[0].dens_mol, Var)
    assert isinstance(gas_prop.fs.unit[0].dens_mol_comp, Var)
    assert isinstance(gas_prop.fs.unit[0].dens_mass, Var)
    assert isinstance(gas_prop.fs.unit[0].cp_mol_comp, Var)
    assert isinstance(gas_prop.fs.unit[0].cp_mol, Var)
    assert isinstance(gas_prop.fs.unit[0].cp_mass, Var)
    assert isinstance(gas_prop.fs.unit[0].visc_d, Var)
    assert isinstance(gas_prop.fs.unit[0].enth_mol, Var)
    assert isinstance(gas_prop.fs.unit[0].entr_mol, Var)


@pytest.mark.unit
def test_setInputs_state_block(gas_prop):
    assert degrees_of_freedom(gas_prop.fs.unit[0]) == 0


@pytest.fixture(scope="class")
def gas_prop_unscaled():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # gas properties and state inlet block
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.unit = m.fs.properties.build_state_block(
        [0], parameters=m.fs.properties, defined_state=True
    )

    m.fs.unit[0].flow_mol.fix(1)
    m.fs.unit[0].temperature.fix(450)
    m.fs.unit[0].pressure.fix(1.60e5)
    m.fs.unit[0].mole_frac_comp["CO2"].fix(0.4772)
    m.fs.unit[0].mole_frac_comp["H2O"].fix(0.0646)
    m.fs.unit[0].mole_frac_comp["CH4"].fix(0.4582)

    return m


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unscaled(gas_prop_unscaled):
    orig_fixed_vars = fixed_variables_set(gas_prop_unscaled)
    orig_act_consts = activated_constraints_set(gas_prop_unscaled)

    gas_prop_unscaled.fs.unit.initialize()

    assert degrees_of_freedom(gas_prop_unscaled) == 0

    fin_fixed_vars = fixed_variables_set(gas_prop_unscaled)
    fin_act_consts = activated_constraints_set(gas_prop_unscaled)

    assert len(fin_act_consts) == len(orig_act_consts)
    assert len(fin_fixed_vars) == len(orig_fixed_vars)

    for c in fin_act_consts:
        assert c in orig_act_consts
    for v in fin_fixed_vars:
        assert v in orig_fixed_vars


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve_unscaled(gas_prop_unscaled):

    assert hasattr(gas_prop_unscaled.fs.unit[0], "mw")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "dens_mol")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "dens_mol_comp")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "dens_mass")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "visc_d")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "therm_cond")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "diffusion_comp")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "cp_mol_comp")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "cp_mol")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "cp_mass")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "enth_mol")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "enth_mol_comp")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "entr_mol")

    assert hasattr(gas_prop_unscaled.fs.unit[0], "get_material_flow_terms")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "get_enthalpy_flow_terms")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "get_material_density_terms")
    assert hasattr(gas_prop_unscaled.fs.unit[0], "get_energy_density_terms")

    results = solver.solve(gas_prop_unscaled)

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.component
def test_scaling(gas_prop):
    # Calculate scaling factors

    # Construct property methods to build the constraints
    assert hasattr(gas_prop.fs.unit[0], "mw")
    assert hasattr(gas_prop.fs.unit[0], "dens_mol")
    assert hasattr(gas_prop.fs.unit[0], "dens_mol_comp")
    assert hasattr(gas_prop.fs.unit[0], "dens_mass")
    assert hasattr(gas_prop.fs.unit[0], "visc_d")
    assert hasattr(gas_prop.fs.unit[0], "therm_cond")
    assert hasattr(gas_prop.fs.unit[0], "diffusion_comp")
    assert hasattr(gas_prop.fs.unit[0], "cp_mol_comp")
    assert hasattr(gas_prop.fs.unit[0], "cp_mol")
    assert hasattr(gas_prop.fs.unit[0], "cp_mass")
    assert hasattr(gas_prop.fs.unit[0], "enth_mol")
    assert hasattr(gas_prop.fs.unit[0], "enth_mol_comp")
    assert hasattr(gas_prop.fs.unit[0], "entr_mol")

    # Call flow and density methods to construct flow and density expressions
    for i in gas_prop.fs.unit[0]._params.component_list:
        gas_prop.fs.unit[0].get_material_flow_terms("Vap", i)
        gas_prop.fs.unit[0].get_material_density_terms("Vap", i)
    gas_prop.fs.unit[0].get_enthalpy_flow_terms("Vap")
    gas_prop.fs.unit[0].get_energy_density_terms("Vap")

    # Calculate scaling factors now that constraints/expressions are built
    iscale.calculate_scaling_factors(gas_prop)

    # Test scaling
    assert pytest.approx(1e-3, abs=1e-2) == iscale.get_scaling_factor(
        gas_prop.fs.unit[0].dens_mol
    )

    for i, c in gas_prop.fs.unit[0].material_flow_terms.items():
        assert pytest.approx(1e-2, abs=1e-2) == iscale.get_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].material_density_terms.items():
        assert pytest.approx(1e-2, abs=1e-3) == iscale.get_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].energy_density_terms.items():
        assert pytest.approx(1e-6, abs=1e-2) == iscale.get_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].enthalpy_flow_terms.items():
        assert pytest.approx(1e-9, abs=1e-2) == iscale.get_scaling_factor(c)

    assert pytest.approx(
        1e2, abs=1e-2
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].mw_eqn
    )
    assert pytest.approx(
        1e-5, abs=1e-2
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].ideal_gas
    )
    for i, c in gas_prop.fs.unit[0].comp_conc_eqn.items():
        assert pytest.approx(
            1e-2, abs=1e-2
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].visc_d_constraint.items():
        assert pytest.approx(
            1e5, abs=1e-2
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].diffusion_comp_constraint.items():
        assert pytest.approx(
            1e5, abs=1e-2
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    assert pytest.approx(
        1, abs=1e-2
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].therm_cond_constraint
    )
    for i, c in gas_prop.fs.unit[0].cp_shomate_eqn.items():
        assert pytest.approx(
            1e-6, abs=1e-2
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    assert pytest.approx(
        1e-6, abs=1e-2
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].mixture_heat_capacity_eqn
    )
    assert pytest.approx(
        1e-6, abs=1e-2
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].cp_mass_basis
    )
    for i, c in gas_prop.fs.unit[0].enthalpy_shomate_eqn.items():
        assert pytest.approx(
            1e-6, abs=1e-2
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    assert pytest.approx(
        1e-6, abs=1e-2
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].mixture_enthalpy_eqn
    )
    assert pytest.approx(
        1e-2, abs=1e-2
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].entropy_correlation
    )


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(gas_prop):
    orig_fixed_vars = fixed_variables_set(gas_prop)
    orig_act_consts = activated_constraints_set(gas_prop)

    gas_prop.fs.unit.initialize()

    assert degrees_of_freedom(gas_prop) == 0

    fin_fixed_vars = fixed_variables_set(gas_prop)
    fin_act_consts = activated_constraints_set(gas_prop)

    assert len(fin_act_consts) == len(orig_act_consts)
    assert len(fin_fixed_vars) == len(orig_fixed_vars)

    for c in fin_act_consts:
        assert c in orig_act_consts
    for v in fin_fixed_vars:
        assert v in orig_fixed_vars


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(gas_prop):

    assert hasattr(gas_prop.fs.unit[0], "mw")
    assert hasattr(gas_prop.fs.unit[0], "dens_mol")
    assert hasattr(gas_prop.fs.unit[0], "cp_mol")
    assert hasattr(gas_prop.fs.unit[0], "visc_d")
    assert hasattr(gas_prop.fs.unit[0], "enth_mol")
    assert hasattr(gas_prop.fs.unit[0], "entr_mol")

    results = solver.solve(gas_prop)

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(gas_prop):
    assert pytest.approx(1, abs=1e-2) == gas_prop.fs.unit[0].mw.value
    assert pytest.approx(1, abs=1e-2) == gas_prop.fs.unit[0].dens_mol.value
    assert pytest.approx(1, abs=1e-2) == gas_prop.fs.unit[0].cp_mol.value
    assert pytest.approx(1e-5, abs=1e-2) == gas_prop.fs.unit[0].visc_d.value
    assert pytest.approx(1, abs=1e-2) == gas_prop.fs.unit[0].enth_mol.value
    assert pytest.approx(1, abs=1e-2) == gas_prop.fs.unit[0].entr_mol.value


@pytest.mark.component
def test_units_consistent(gas_prop):

    # Construct property methods to build the constraints
    assert hasattr(gas_prop.fs.unit[0], "mw")
    assert hasattr(gas_prop.fs.unit[0], "dens_mol")
    assert hasattr(gas_prop.fs.unit[0], "dens_mol_comp")
    assert hasattr(gas_prop.fs.unit[0], "dens_mass")
    assert hasattr(gas_prop.fs.unit[0], "visc_d")
    assert hasattr(gas_prop.fs.unit[0], "therm_cond")
    assert hasattr(gas_prop.fs.unit[0], "diffusion_comp")
    assert hasattr(gas_prop.fs.unit[0], "cp_mol_comp")
    assert hasattr(gas_prop.fs.unit[0], "cp_mol")
    assert hasattr(gas_prop.fs.unit[0], "cp_mass")
    assert hasattr(gas_prop.fs.unit[0], "enth_mol")
    assert hasattr(gas_prop.fs.unit[0], "enth_mol_comp")
    assert hasattr(gas_prop.fs.unit[0], "entr_mol")

    # Call flow and density methods to construct flow and density expressions
    for i in gas_prop.fs.unit[0]._params.component_list:
        gas_prop.fs.unit[0].get_material_flow_terms("Vap", i)
        gas_prop.fs.unit[0].get_material_density_terms("Vap", i)
    gas_prop.fs.unit[0].get_enthalpy_flow_terms("Vap")
    gas_prop.fs.unit[0].get_energy_density_terms("Vap")

    assert_units_consistent(gas_prop)
