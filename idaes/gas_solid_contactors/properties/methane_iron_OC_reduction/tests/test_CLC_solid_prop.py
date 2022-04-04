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
Tests for CLC solid phase thermo state block; tests for construction and solve
Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.testing import initialization_tester
from idaes.core.util import get_solver, scaling as iscale

from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    solid_phase_thermo import SolidPhaseParameterBlock

# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.fixture(scope="class")
def solid_prop():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # solid properties and state inlet block
    m.fs.properties = SolidPhaseParameterBlock()
    m.fs.unit = m.fs.properties.build_state_block(
        default={"parameters": m.fs.properties,
                 "defined_state": True})

    m.fs.unit.flow_mass.fix(1)
    m.fs.unit.particle_porosity.fix(0.27)
    m.fs.unit.temperature.fix(1183.15)
    m.fs.unit.mass_frac_comp["Fe2O3"].fix(0.45)
    m.fs.unit.mass_frac_comp["Fe3O4"].fix(1e-9)
    m.fs.unit.mass_frac_comp["Al2O3"].fix(0.55)

    return m


@pytest.mark.unit
def test_build_inlet_state_block(solid_prop):
    assert isinstance(solid_prop.fs.unit.dens_mass_skeletal, Var)
    assert isinstance(solid_prop.fs.unit.enth_mol_comp, Var)
    assert isinstance(solid_prop.fs.unit.enth_mass, Var)
    assert isinstance(solid_prop.fs.unit.cp_mol_comp, Var)
    assert isinstance(solid_prop.fs.unit.cp_mass, Var)


@pytest.mark.unit
def test_setInputs_state_block(solid_prop):
    assert degrees_of_freedom(solid_prop.fs.unit) == 0


@pytest.fixture(scope="class")
def solid_prop_unscaled():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # solid properties and state inlet block
    m.fs.properties = SolidPhaseParameterBlock()
    m.fs.unit = m.fs.properties.build_state_block(
        default={"parameters": m.fs.properties,
                 "defined_state": True})

    m.fs.unit.flow_mass.fix(1)
    m.fs.unit.particle_porosity.fix(0.27)
    m.fs.unit.temperature.fix(1183.15)
    m.fs.unit.mass_frac_comp["Fe2O3"].fix(0.45)
    m.fs.unit.mass_frac_comp["Fe3O4"].fix(1e-9)
    m.fs.unit.mass_frac_comp["Al2O3"].fix(0.55)

    return m


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unscaled(solid_prop_unscaled):
    initialization_tester(
            solid_prop_unscaled)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve_unscaled(solid_prop_unscaled):

    assert hasattr(solid_prop_unscaled.fs.unit, "dens_mass_skeletal")
    assert hasattr(solid_prop_unscaled.fs.unit, "cp_mass")
    assert hasattr(solid_prop_unscaled.fs.unit, "enth_mass")

    results = solver.solve(solid_prop_unscaled)

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.component
def test_scaling(solid_prop):
    # Calculate scaling factors

    # Construct property methods to build the constraints

    assert hasattr(solid_prop.fs.unit, "flow_mass")
    assert hasattr(solid_prop.fs.unit, "particle_porosity")
    assert hasattr(solid_prop.fs.unit, "temperature")
    assert hasattr(solid_prop.fs.unit, "mass_frac_comp")
    assert hasattr(solid_prop.fs.unit, "dens_mass_skeletal")
    assert hasattr(solid_prop.fs.unit, "dens_mass_particle")
    assert hasattr(solid_prop.fs.unit, "cp_mol_comp")
    assert hasattr(solid_prop.fs.unit, "cp_mass")
    assert hasattr(solid_prop.fs.unit, "enth_mass")
    assert hasattr(solid_prop.fs.unit, "enth_mol_comp")

    # Call flow and density methods to construct flow and density expressions
    for i in solid_prop.fs.unit._params.component_list:
        solid_prop.fs.unit.get_material_flow_terms('Sol', i)
        solid_prop.fs.unit.get_material_density_terms('Sol', i)
    solid_prop.fs.unit.get_enthalpy_flow_terms('Sol')
    solid_prop.fs.unit.get_energy_density_terms('Sol')

    # Calculate scaling factors now that constraints/expressions are built
    iscale.calculate_scaling_factors(solid_prop)

    # Test scaling
    assert (pytest.approx(1e-3, abs=1e-2) ==
            iscale.get_scaling_factor(
                solid_prop.fs.unit.dens_mass_particle))

    for i, c in solid_prop.fs.unit.material_flow_terms.items():
        assert (pytest.approx(1e-2, abs=1e-2) ==
                iscale.get_scaling_factor(c))
    for i, c in solid_prop.fs.unit.material_density_terms.items():
        assert (pytest.approx(1e-1, abs=1e-2) ==
                iscale.get_scaling_factor(c))
    for i, c in solid_prop.fs.unit.energy_density_terms.items():
        assert (pytest.approx(1e-6, abs=1e-2) ==
                iscale.get_scaling_factor(c))
    for i, c in solid_prop.fs.unit.enthalpy_flow_terms.items():
        assert (pytest.approx(1e-9, abs=1e-2) ==
                iscale.get_scaling_factor(c))

    assert (pytest.approx(1e-2, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                solid_prop.fs.unit.density_particle_constraint))
    assert (pytest.approx(1e-5, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                solid_prop.fs.unit.density_skeletal_constraint))
    for i, c in solid_prop.fs.unit.cp_shomate_eqn.items():
        assert (pytest.approx(1e-6, abs=1e-2) ==
                iscale.get_constraint_transform_applied_scaling_factor(c))
    assert (pytest.approx(1e-6, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                solid_prop.fs.unit.mixture_heat_capacity_eqn))
    for i, c in solid_prop.fs.unit.enthalpy_shomate_eqn.items():
        assert (pytest.approx(1e-6, abs=1e-2) ==
                iscale.get_constraint_transform_applied_scaling_factor(c))
    assert (pytest.approx(1e-6, abs=1e-2) ==
            iscale.get_constraint_transform_applied_scaling_factor(
                solid_prop.fs.unit.mixture_enthalpy_eqn))


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(solid_prop):
    initialization_tester(
            solid_prop)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(solid_prop):

    assert hasattr(solid_prop.fs.unit, "dens_mass_skeletal")
    assert hasattr(solid_prop.fs.unit, "cp_mass")
    assert hasattr(solid_prop.fs.unit, "enth_mass")

    results = solver.solve(solid_prop)

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(solid_prop):
    assert (pytest.approx(3251.75, abs=1e-2) ==
            solid_prop.fs.unit.dens_mass_skeletal.value)
    assert (pytest.approx(1, abs=1e-2) ==
            solid_prop.fs.unit.cp_mass.value)
    assert (pytest.approx(0.0039, abs=1e-2) ==
            solid_prop.fs.unit.enth_mass.value)


@pytest.mark.component
def test_units_consistent(solid_prop):

    # Construct property methods to build the constraints

    assert hasattr(solid_prop.fs.unit, "flow_mass")
    assert hasattr(solid_prop.fs.unit, "particle_porosity")
    assert hasattr(solid_prop.fs.unit, "temperature")
    assert hasattr(solid_prop.fs.unit, "mass_frac_comp")
    assert hasattr(solid_prop.fs.unit, "dens_mass_skeletal")
    assert hasattr(solid_prop.fs.unit, "dens_mass_particle")
    assert hasattr(solid_prop.fs.unit, "cp_mol_comp")
    assert hasattr(solid_prop.fs.unit, "cp_mass")
    assert hasattr(solid_prop.fs.unit, "enth_mass")
    assert hasattr(solid_prop.fs.unit, "enth_mol_comp")

    # Call flow and density methods to construct flow and density expressions
    for i in solid_prop.fs.unit._params.component_list:
        solid_prop.fs.unit.get_material_flow_terms('Sol', i)
        solid_prop.fs.unit.get_material_density_terms('Sol', i)
    solid_prop.fs.unit.get_enthalpy_flow_terms('Sol')
    solid_prop.fs.unit.get_energy_density_terms('Sol')

    assert_units_consistent(solid_prop)
