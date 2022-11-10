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

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Var,
    value,
    units as pyunits,
)
from pyomo.contrib.incidence_analysis import (
    solve_strongly_connected_components,
)
from pyomo.common.collections import ComponentMap
import pyomo.common.unittest as unittest
from pyomo.common.unittest import TestCase
from pyomo.util.subsystems import ParamSweeper
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_large_residuals,
    fixed_variables_set,
    activated_constraints_set,
)

from idaes.core.solvers import get_solver
from idaes.core.util import scaling as iscale

from idaes.models_extra.gas_solid_contactors.properties.oxygen_iron_OC_oxidation.gas_phase_thermo import (
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
    m.fs.unit[0].mole_frac_comp["CO2"].fix(0.0004)
    m.fs.unit[0].mole_frac_comp["H2O"].fix(0.0093)
    m.fs.unit[0].mole_frac_comp["O2"].fix(0.2095)
    m.fs.unit[0].mole_frac_comp["N2"].fix(0.7808)

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
    m.fs.unit[0].mole_frac_comp["CO2"].fix(0.0004)
    m.fs.unit[0].mole_frac_comp["H2O"].fix(0.0093)
    m.fs.unit[0].mole_frac_comp["O2"].fix(0.2095)
    m.fs.unit[0].mole_frac_comp["N2"].fix(0.7808)

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
    assert pytest.approx(1e-3, rel=1e-5) == iscale.get_scaling_factor(
        gas_prop.fs.unit[0].dens_mol
    )

    for i, c in gas_prop.fs.unit[0].material_flow_terms.items():
        assert pytest.approx(1e-2, rel=1e-5) == iscale.get_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].material_density_terms.items():
        assert pytest.approx(1e-2, rel=1e-5) == iscale.get_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].energy_density_terms.items():
        assert pytest.approx(1e-9, rel=1e-5) == iscale.get_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].enthalpy_flow_terms.items():
        assert pytest.approx(1e-9, rel=1e-5) == iscale.get_scaling_factor(c)

    assert pytest.approx(
        1e2, rel=1e-5
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].mw_eqn
    )
    assert pytest.approx(
        1e-5, rel=1e-5
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].ideal_gas
    )
    for i, c in gas_prop.fs.unit[0].comp_conc_eqn.items():
        assert pytest.approx(
            1e-2, rel=1e-5
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].visc_d_constraint.items():
        assert pytest.approx(
            1e5, rel=1e-5
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    for i, c in gas_prop.fs.unit[0].diffusion_comp_constraint.items():
        assert pytest.approx(
            1e5, rel=1e-5
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    assert pytest.approx(
        1, rel=1e-5
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].therm_cond_constraint
    )
    for i, c in gas_prop.fs.unit[0].cp_shomate_eqn.items():
        assert pytest.approx(
            1e-6, rel=1e-5
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    assert pytest.approx(
        1e-6, rel=1e-5
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].mixture_heat_capacity_eqn
    )
    assert pytest.approx(
        1e-6, rel=1e-5
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].cp_mass_basis
    )
    for i, c in gas_prop.fs.unit[0].enthalpy_shomate_eqn.items():
        assert pytest.approx(
            1e-6, rel=1e-5
        ) == iscale.get_constraint_transform_applied_scaling_factor(c)
    assert pytest.approx(
        1e-6, rel=1e-5
    ) == iscale.get_constraint_transform_applied_scaling_factor(
        gas_prop.fs.unit[0].mixture_enthalpy_eqn
    )
    assert pytest.approx(
        1e-5, rel=1e-5
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
    assert hasattr(gas_prop.fs.unit[0], "cp_mol")
    assert hasattr(gas_prop.fs.unit[0], "enth_mol")

    results = solver.solve(gas_prop)

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(gas_prop):
    assert pytest.approx(1, rel=1e-5) == gas_prop.fs.unit[0].mw.value
    assert pytest.approx(1, rel=1e-5) == gas_prop.fs.unit[0].cp_mol.value
    assert pytest.approx(1, rel=1e-5) == gas_prop.fs.unit[0].enth_mol.value


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


# -----------------------------------------------------------------------------


@pytest.mark.unit
def test_state_vars():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = GasPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block()

    assert isinstance(m.fs.state.flow_mol, Var)
    assert isinstance(m.fs.state.temperature, Var)
    assert isinstance(m.fs.state.pressure, Var)
    assert isinstance(m.fs.state.mole_frac_comp, Var)

    assert isinstance(m.fs.state.sum_component_eqn, Constraint)

    assert len(list(m.component_data_objects(Var))) == 7
    assert len(list(m.component_data_objects(Constraint))) == 1

    for name, var in m.fs.state.define_state_vars().items():
        # State vars should be included in the metadata with no method
        assert name in m.fs.properties._metadata._properties
        assert m.fs.properties._metadata._properties[name]["method"] is None


@pytest.mark.unit
def test_indexed_state_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block([1, 2, 3])

    assert len(list(m.component_data_objects(Var))) == 21
    assert len(list(m.component_data_objects(Constraint))) == 3

    for i, state in m.fs.state.items():
        assert isinstance(state.flow_mol, Var)
        assert isinstance(state.temperature, Var)
        assert isinstance(state.pressure, Var)
        assert isinstance(state.mole_frac_comp, Var)

        assert isinstance(state.sum_component_eqn, Constraint)

        for name, var in state.define_state_vars().items():
            # State vars should be included in the metadata with no method
            assert name in m.fs.properties._metadata._properties
            assert m.fs.properties._metadata._properties[name]["method"] is None


@pytest.mark.unit
def test_property_construction_ordered():
    """
    The purpose of this test is to make sure each property can be
    constructed properly. Further, because we know that these
    variables form a DAG, we check that the addition of each new
    property adds one variable and constraint, i.e. that the
    matching defined in the parameter block is in topological order.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block()

    nvar = 7
    ncon = 1
    assert len(list(m.component_data_objects(Var))) == nvar
    assert len(list(m.component_data_objects(Constraint))) == ncon

    matching = [
        ("mw", "mw_eqn"),
        ("dens_mol", "ideal_gas"),
        ("dens_mol_comp", "comp_conc_eqn"),
        ("dens_mass", "dens_mass_basis"),
        ("visc_d", "visc_d_constraint"),
        ("diffusion_comp", "diffusion_comp_constraint"),
        ("therm_cond", "therm_cond_constraint"),
        ("cp_mol_comp", "cp_shomate_eqn"),
        ("cp_mol", "mixture_heat_capacity_eqn"),
        ("cp_mass", "cp_mass_basis"),
        ("enth_mol_comp", "enthalpy_shomate_eqn"),
        ("enth_mol", "mixture_enthalpy_eqn"),
        ("entr_mol", "entropy_correlation"),
    ]

    # Make sure the matching captures all the non-state vars.
    state_vars = m.fs.state.define_state_vars()
    n_state_vars = len(state_vars)
    n_vars = len(m.fs.properties._metadata._properties)
    assert len(matching) == n_vars - n_state_vars

    # Make sure we can add constraints one at a time, and only
    # add one additional variable. This is specific to this
    # property package.
    for varname, conname in matching:
        assert varname not in state_vars
        var = getattr(m.fs.state, varname)
        con = getattr(m.fs.state, conname)
        dim = len(var)
        nvar += dim
        ncon += dim
        assert dim == len(con)
        assert isinstance(var, Var)
        assert isinstance(con, Constraint)
        assert len(list(m.component_data_objects(Var))) == nvar
        assert len(list(m.component_data_objects(Constraint))) == ncon


@pytest.mark.unit
def test_property_construction_unordered():
    """
    The purpose of this test is to make sure "chained construction"
    of properties happens properly. I.e. if we add a property that
    requires several other properties, these are added as well.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block()

    nvar = 7
    ncon = 1
    assert len(list(m.component_data_objects(Var))) == nvar
    assert len(list(m.component_data_objects(Constraint))) == ncon

    # Access this property to force the construction of prerequisite
    # properties.
    m.fs.state.cp_mass

    # Make sure dependent properties have been properly added.
    assert hasattr(m.fs.state, "cp_mol")
    assert hasattr(m.fs.state, "cp_mol_comp")
    assert hasattr(m.fs.state, "mw")
    nvar += len(m.fs.state.cp_mol_comp)
    nvar += len(m.fs.state.cp_mol)
    nvar += len(m.fs.state.cp_mass)
    nvar += len(m.fs.state.mw)
    ncon += len(m.fs.state.cp_mol_comp)
    ncon += len(m.fs.state.cp_mol)
    ncon += len(m.fs.state.cp_mass)
    ncon += len(m.fs.state.mw)
    assert len(list(m.component_data_objects(Var))) == nvar
    assert len(list(m.component_data_objects(Constraint))) == ncon


@pytest.mark.unit
class TestProperties(TestCase):
    """
    The purpose of these tests is to ensure that the property package's
    state block can be used to solve for the specified properties, and
    that the values produced are what we expect. Some of these values are
    computed by hand, while others (viscosity, enthalpy, conductivity, etc.)
    are computed from these constraints we are testing.
    These tests pretty much just make sure the calculated values don't change.
    If the property equations change, and the values change as a result,
    the tolerances may be relaxed. If the new values are significantly
    different, this should be investigated.
    """

    def _make_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = GasPhaseParameterBlock()
        m.fs.state = m.fs.properties.build_state_block(defined_state=True)
        for var in m.fs.state.define_state_vars().values():
            var.fix()
        return m

    def test_mw(self):
        m = self._make_model()
        state = m.fs.state

        # Define a somewhat arbitrary set of values we'd like to use to
        # test molecular weight.
        n_scenario = 7
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scenario,
            "temperature": [300.0 * pyunits.K] * n_scenario,
            "pressure": [1.0 * pyunits.bar] * n_scenario,
            "mole_frac_comp[O2]": [1.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.0],
            "mole_frac_comp[N2]": [0.0, 0.5, 0.25, 1.0, 0.0, 0.0, 0.0],
            "mole_frac_comp[H2O]": [0.0, 0.0, 0.25, 0.0, 1.0, 0.0, 0.5],
            "mole_frac_comp[CO2]": [0.0, 0.0, 0.25, 0.0, 0.0, 1.0, 0.5],
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )
        target_values = {
            "mw": [
                32.0 * pyunits.g / pyunits.mol,
                30.0 * pyunits.g / pyunits.mol,
                30.5 * pyunits.g / pyunits.mol,
                28.0 * pyunits.g / pyunits.mol,
                18.0 * pyunits.g / pyunits.mol,
                44.0 * pyunits.g / pyunits.mol,
                31.0 * pyunits.g / pyunits.mol,
            ]
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct mw and all prerequisites
        state.mw

        param_sweeper = ParamSweeper(
            n_scenario,
            state_values,
            output_values=target_values,
        )
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Check that state block has been been solved correctly
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs have been set properly
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(val, abs=1e-8)

                # Check that the state block computes the property values
                # we expect
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(val, abs=1e-8)

    def test_dens_mol(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 8
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scen,
            "temperature": [
                200.0 * pyunits.K,
                300.0 * pyunits.K,
                600.0 * pyunits.K,
                900.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
            ],
            "pressure": [
                1.0 * pyunits.bar,
                1.0 * pyunits.bar,
                1.0 * pyunits.bar,
                1.0 * pyunits.bar,
                1.0 * pyunits.bar,
                0.5 * pyunits.bar,
                1.5 * pyunits.bar,
                2.0 * pyunits.bar,
            ],
            "mole_frac_comp[O2]": [0.25] * n_scen,
            "mole_frac_comp[N2]": [0.25] * n_scen,
            "mole_frac_comp[H2O]": [0.25] * n_scen,
            "mole_frac_comp[CO2]": [0.25] * n_scen,
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )

        u = pyunits.mol / pyunits.m**3
        target_values = {
            "dens_mol": [
                # Calculated by P*1e5/(T*8.314462618)
                60.136 * u,
                40.091 * u,
                20.045 * u,
                13.364 * u,
                10.023 * u,
                5.011 * u,
                15.034 * u,
                20.045 * u,
            ]
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct dens_mol and all prerequisites
        state.dens_mol
        assert_units_consistent(state.ideal_gas)

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Check that state block equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are what we expect
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(val, rel=1e-3)

                # Check that state block performs the calculations we expect
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(val, rel=1e-3)

    def test_dens_mol_comp(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 5
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scen,
            "temperature": [1200.0 * pyunits.K] * n_scen,
            "pressure": [1.0 * pyunits.bar] * n_scen,
            "mole_frac_comp[O2]": [1.0, 0.0, 0.0, 0.0, 0.25],
            "mole_frac_comp[N2]": [0.0, 1.0, 0.0, 0.0, 0.25],
            "mole_frac_comp[H2O]": [0.0, 0.0, 1.0, 0.0, 0.25],
            "mole_frac_comp[CO2]": [0.0, 0.0, 0.0, 1.0, 0.25],
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )

        u = pyunits.mol / pyunits.m**3
        target_values = {
            "dens_mol": [10.023 * u] * n_scen,
            "dens_mol_comp[O2]": [
                10.023 * u,
                0.0 * u,
                0.0 * u,
                0.0 * u,
                0.25 * 10.023 * u,
            ],
            "dens_mol_comp[N2]": [
                0.0 * u,
                10.023 * u,
                0.0 * u,
                0.0 * u,
                0.25 * 10.023 * u,
            ],
            "dens_mol_comp[H2O]": [
                0.0 * u,
                0.0 * u,
                10.023 * u,
                0.0 * u,
                0.25 * 10.023 * u,
            ],
            "dens_mol_comp[CO2]": [
                0.0 * u,
                0.0 * u,
                0.0 * u,
                10.023 * u,
                0.25 * 10.023 * u,
            ],
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct dens_mol_comp and all prerequisites
        state.dens_mol_comp

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs have been set properly
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    # ^ Problem converting units when value is zero
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure property values are what we expect
                for var, val in target.items():
                    # val = value(pyunits.convert(val, var.get_units()))
                    # ^ Problem converting units when value is zero
                    assert var.value == pytest.approx(value(val), rel=1e-3)

    def test_visc_d(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 8
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scen,
            "temperature": [
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                300.0 * pyunits.K,
                600.0 * pyunits.K,
                900.0 * pyunits.K,
                1200.0 * pyunits.K,
            ],
            "pressure": [1.0 * pyunits.bar] * n_scen,
            "mole_frac_comp[O2]": [1.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[N2]": [0.0, 1.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[H2O]": [0.0, 0.0, 1.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[CO2]": [0.0, 0.0, 0.0, 1.0, 0.25, 0.25, 0.25, 0.25],
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )

        target_values = {
            "visc_d": [
                # These values were copied after a solve.
                # TODO: Cross-reference with another source.
                5.534440949133228e-05 * pyunits.Pa * pyunits.s,
                4.67667824296429e-05 * pyunits.Pa * pyunits.s,
                4.6232771210527155e-05 * pyunits.Pa * pyunits.s,
                4.512867970060493e-05 * pyunits.Pa * pyunits.s,
                1.6181534595116313e-05 * pyunits.Pa * pyunits.s,
                2.866222939903063e-05 * pyunits.Pa * pyunits.s,
                3.909320395131273e-05 * pyunits.Pa * pyunits.s,
                4.838841106600266e-05 * pyunits.Pa * pyunits.s,
            ]
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct visc_d and all prerequisites
        state.visc_d

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

    def test_diffusion_comp(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 11
        bar = pyunits.bar
        K = pyunits.K
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scen,
            "temperature": [
                1200.0 * K,
                1200.0 * K,
                1200.0 * K,
                1200.0 * K,
                300.0 * K,
                600.0 * K,
                900.0 * K,
                1200.0 * K,
                1200.0 * K,
                1200.0 * K,
                1200.0 * K,
            ],
            "pressure": [
                1.0 * bar,
                1.0 * bar,
                1.0 * bar,
                1.0 * bar,
                1.0 * bar,
                1.0 * bar,
                1.0 * bar,
                1.0 * bar,
                0.5 * bar,
                1.5 * bar,
                2.0 * bar,
            ],
            "mole_frac_comp[O2]": [
                # Note that diffusivity is not defined for a pure
                # component in itself (zero gradient, zero net diffusion)
                0.90,
                0.025,
                0.025,
                0.025,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
            "mole_frac_comp[N2]": [
                0.025,
                0.90,
                0.025,
                0.025,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
            "mole_frac_comp[H2O]": [
                0.025,
                0.025,
                0.90,
                0.025,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
            "mole_frac_comp[CO2]": [
                0.025,
                0.025,
                0.025,
                0.90,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
            ],
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )

        cm = pyunits.cm
        s = pyunits.s
        target_values = {
            # These values look reasonable
            # TODO: Verify with external source.
            "diffusion_comp[O2]": [
                3.11792621830951e-5 * cm**2 / s,
                2.456227751888218e-5 * cm**2 / s,
                3.034740091620132e-5 * cm**2 / s,
                1.9479388404541838e-5 * cm**2 / s,
                0.206691259894311e-5 * cm**2 / s,
                0.6952237580376006e-5 * cm**2 / s,
                1.4134625566198689e-5 * cm**2 / s,
                2.3384446637321363e-5 * cm**2 / s,
                4.676889327464272e-5 * cm**2 / s,
                1.5589631091547582e-5 * cm**2 / s,
                1.1692223318660684e-5 * cm**2 / s,
            ],
            "diffusion_comp[N2]": [
                2.457518754629481e-5 * cm**2 / s,
                3.1383309495574956e-5 * cm**2 / s,
                3.0334118458311523e-5 * cm**2 / s,
                1.9790010245966736e-5 * cm**2 / s,
                0.2080439152537244e-5 * cm**2 / s,
                0.6997735302088169e-5 * cm**2 / s,
                1.422712718932098e-5 * cm**2 / s,
                2.353748212168125e-5 * cm**2 / s,
                4.70749642433625e-5 * cm**2 / s,
                1.5691654747787498e-5 * cm**2 / s,
                1.176874106084062e-5 * cm**2 / s,
            ],
            "diffusion_comp[H2O]": [
                3.0845168350215713e-5 * cm**2 / s,
                3.0811129521516416e-5 * cm**2 / s,
                3.719143107603082e-5 * cm**2 / s,
                2.5051274838003996e-5 * cm**2 / s,
                0.24654668546150185e-5 * cm**2 / s,
                0.8292808959890481e-5 * cm**2 / s,
                1.6860147281349098e-5 * cm**2 / s,
                2.7893573307023156e-5 * cm**2 / s,
                5.578714661404631e-5 * cm**2 / s,
                1.8595715538015434e-5 * cm**2 / s,
                1.3946786653511578e-5 * cm**2 / s,
            ],
            "diffusion_comp[CO2]": [
                1.929322600571961e-5 * cm**2 / s,
                1.9589681584041365e-5 * cm**2 / s,
                2.4422863072776697e-5 * cm**2 / s,
                2.7098612589802444e-5 * cm**2 / s,
                0.17964011927809195e-5 * cm**2 / s,
                0.6042349293467888e-5 * cm**2 / s,
                1.2284727588198259e-5 * cm**2 / s,
                2.0323959442351844e-5 * cm**2 / s,
                4.064791888470369e-5 * cm**2 / s,
                1.3549306294901227e-5 * cm**2 / s,
                1.0161979721175922e-5 * cm**2 / s,
            ],
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct diffusion_comp and all prerequisites
        state.diffusion_comp

        assert_units_consistent(state.diffusion_comp_constraint)

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

    def test_therm_cond(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 8
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scen,
            "temperature": [
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                300.0 * pyunits.K,
                600.0 * pyunits.K,
                900.0 * pyunits.K,
                1200.0 * pyunits.K,
            ],
            "pressure": [1.0 * pyunits.bar] * n_scen,
            "mole_frac_comp[O2]": [1.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[N2]": [0.0, 1.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[H2O]": [0.0, 0.0, 1.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[CO2]": [0.0, 0.0, 0.0, 1.0, 0.25, 0.25, 0.25, 0.25],
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )

        u = pyunits.kJ / pyunits.m / pyunits.K / pyunits.s
        target_values = {
            "therm_cond": [
                8.490673044994837e-05 * u,
                7.803113104821915e-05 * u,
                0.0001245121534187936 * u,
                7.844692969560201e-05 * u,
                2.1981943936613706e-05 * u,
                4.567583423706824e-05 * u,
                6.946515568649932e-05 * u,
                9.30078254960681e-05 * u,
            ],
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct therm_cond and all prerequisites
        state.therm_cond

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

    def test_cp(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 8
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scen,
            "temperature": [
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                300.0 * pyunits.K,
                600.0 * pyunits.K,
                900.0 * pyunits.K,
                1200.0 * pyunits.K,
            ],
            "pressure": [1.0 * pyunits.bar] * n_scen,
            "mole_frac_comp[O2]": [1.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[N2]": [0.0, 1.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[H2O]": [0.0, 0.0, 1.0, 0.0, 0.25, 0.25, 0.25, 0.25],
            "mole_frac_comp[CO2]": [0.0, 0.0, 0.0, 1.0, 0.25, 0.25, 0.25, 0.25],
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )

        u = pyunits.kJ / pyunits.mol / pyunits.K
        kJkgK = pyunits.kJ / pyunits.kg / pyunits.K
        target_values = {
            "cp_mol_comp[O2]": [
                0.03566421043844448 * u,
                0.03566421043844444 * u,
                0.03566421043844444 * u,
                0.03566421043844444 * u,
                0.02408660519211111 * u,
                0.031970683705777776 * u,
                0.03435676292601234 * u,
                0.03566421043844444 * u,
            ],
            "cp_mol_comp[N2]": [
                0.03372177593533332 * u,
                0.03372177593533333 * u,
                0.03372177593533333 * u,
                0.03372177593533333 * u,
                0.03059729435133333 * u,
                0.030104019077333333 * u,
                0.03208929344525926 * u,
                0.03372177593533333 * u,
            ],
            "cp_mol_comp[H2O]": [
                0.0437510227322222 * u,
                0.04375102273222223 * u,
                0.04375102273222223 * u,
                0.04375102273222223 * u,
                0.03359738794555556 * u,
                0.036317861208888885 * u,
                0.039997715202839505 * u,
                0.04375102273222223 * u,
            ],
            "cp_mol_comp[CO2]": [
                0.05634605443600005 * u,
                0.056346054436000007 * u,
                0.056346054436000007 * u,
                0.056346054436000007 * u,
                0.037217621149000006 * u,
                0.047317934392 * u,
                0.053001289534111116 * u,
                0.056346054436000007 * u,
            ],
            "cp_mol": [
                0.03566421043844448 * u,
                0.03372177593533333 * u,
                0.04375102273222223 * u,
                0.056346054436000007 * u,
                0.0313747271595 * u,
                0.036427624596 * u,
                0.03986126527705556 * u,
                0.25
                * (
                    # Component values at 1200 K have been computed.
                    0.03566421043844444 * u
                    + 0.03372177593533333 * u
                    + 0.04375102273222223 * u
                    + 0.056346054436000007 * u
                ),
            ],
            "cp_mass": [
                1.1145065762013922 * kJkgK,
                1.2043491405476134 * kJkgK,
                2.4306123740123473 * kJkgK,
                1.280592146272725 * kJkgK,
                1.0286795790000056 * kJkgK,
                1.194348347409837 * kJkgK,
                1.3069267303952685 * kJkgK,
                1.3892054388688537 * kJkgK,
            ],
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct cp_mass and all prerequisites.
        # This constructs cp_mol and cp_mol_comp as well.
        state.cp_mass

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

    def test_enth(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 7
        state_values = {
            "flow_mol": [1.0 * pyunits.mol / pyunits.s] * n_scen,
            "temperature": [
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                1200.0 * pyunits.K,
                # We do not test enthalpy at 300 K. The Shomate equation
                # we use is not valid below 500 K.
                # 300.0*pyunits.K,
                600.0 * pyunits.K,
                900.0 * pyunits.K,
                1200.0 * pyunits.K,
            ],
            "pressure": [1.0 * pyunits.bar] * n_scen,
            "mole_frac_comp[O2]": [1.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25],
            "mole_frac_comp[N2]": [0.0, 1.0, 0.0, 0.0, 0.25, 0.25, 0.25],
            "mole_frac_comp[H2O]": [0.0, 0.0, 1.0, 0.0, 0.25, 0.25, 0.25],
            "mole_frac_comp[CO2]": [0.0, 0.0, 0.0, 1.0, 0.25, 0.25, 0.25],
        }
        state_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in state_values.items()
        )

        u = pyunits.kJ / pyunits.mol
        target_values = {
            "enth_mol_comp[O2]": [
                29.760175857866656 * u,
                29.760175857866656 * u,
                29.760175857866656 * u,
                29.760175857866656 * u,
                # 0.5175085434916653*u,
                9.248259058533336 * u,
                19.241674269713887 * u,
                29.760175857866656 * u,
            ],
            "enth_mol_comp[N2]": [
                28.1081423656 * u,
                28.1081423656 * u,
                28.1081423656 * u,
                28.1081423656 * u,
                # -0.021818752399997976*u,
                8.8939164816 * u,
                18.223311732266666 * u,
                28.1081423656 * u,
            ],
            "enth_mol_comp[H2O]": [
                34.505905041333335 * u,
                34.505905041333335 * u,
                34.505905041333335 * u,
                34.505905041333335 * u,
                # 0.06267505633334736*u,
                10.500564354666665 * u,
                21.939189237444452 * u,
                34.505905041333335 * u,
            ],
            "enth_mol_comp[CO2]": [
                44.474410900800024 * u,
                44.474410900800024 * u,
                44.474410900800024 * u,
                44.474410900800024 * u,
                # 0.06585135367498651*u,
                12.906441898799983 * u,
                28.031785067675003 * u,
                44.474410900800024 * u,
            ],
            "enth_mol": [
                29.760175857866656 * u,
                28.1081423656 * u,
                34.505905041333335 * u,
                44.474410900800024 * u,
                # 0.15605405027500296*u,
                10.387295448399996 * u,
                21.858990076775 * u,
                0.25
                * (
                    29.760175857866656 * u
                    + 28.1081423656 * u
                    + 34.505905041333335 * u
                    + 44.474410900800024 * u
                ),
            ],
        }
        target_values = ComponentMap(
            (state.find_component(name), values)
            for name, values in target_values.items()
        )

        # Construct enth_mol and all prerequisites, including enth_mol_comp
        state.enth_mol

        param_sweeper = ParamSweeper(n_scen, state_values, output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), rel=1e-3)


if __name__ == "__main__":
    unittest.main()
