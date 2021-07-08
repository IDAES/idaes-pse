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

from pyomo.environ import (
        ConcreteModel,
        TerminationCondition,
        SolverStatus,
        Var,
        Constraint,
        )

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.testing import (get_default_solver,
                                     initialization_tester)

from idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation. \
    solid_phase_thermo import SolidPhaseParameterBlock

# Get default solver for testing
solver = get_default_solver()


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
    m.fs.unit.mass_frac_comp["Fe2O3"].fix(0.244)
    m.fs.unit.mass_frac_comp["Fe3O4"].fix(0.202)
    m.fs.unit.mass_frac_comp["Al2O3"].fix(0.554)

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
    assert results.solver.termination_condition == \
        TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok


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


def test_state_vars():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = SolidPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block()

    assert isinstance(m.fs.state.flow_mass, Var)
    assert isinstance(m.fs.state.temperature, Var)
    assert isinstance(m.fs.state.particle_porosity, Var)
    assert isinstance(m.fs.state.mass_frac_comp, Var)

    assert isinstance(m.fs.state.sum_component_eqn, Constraint)

    assert len(list(m.fs.state.component_data_objects(Var))) == 6
    assert len(list(m.component_data_objects(Constraint))) == 1

    for name, var in m.fs.state.define_state_vars().items():
        assert name in m.fs.properties._metadata._properties


def test_indexed_state_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = SolidPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block([1,2,3])

    assert len([v for v in m.component_data_objects(Var) if not v.fixed]) == 18
    assert len(list(m.component_data_objects(Constraint))) == 3

    for i, state in m.fs.state.items():
        assert isinstance(state.flow_mass, Var)
        assert isinstance(state.temperature, Var)
        assert isinstance(state.particle_porosity, Var)
        assert isinstance(state.mass_frac_comp, Var)

        assert isinstance(state.sum_component_eqn, Constraint)


def test_property_construction_ordered():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = SolidPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block()

    # If we construct properties in this order, variables and constraints
    # will be added one at a time.
    matching = [
            ("dens_mass_skeletal", "density_skeletal_constraint"),
            ("dens_mass_particle", "density_particle_constraint"),
            ("cp_mol_comp", "cp_shomate_eqn"),
            ("cp_mass", "mixture_heat_capacity_eqn"),
            ("enth_mol_comp", "enthalpy_shomate_eqn"),
            ("enth_mass", "mixture_enthalpy_eqn"),
            ]

    state_vars = m.fs.state.define_state_vars()
    n_state_vars = len(state_vars)
    n_vars = len(m.fs.properties._metadata._properties)
    assert len(matching) == n_vars - n_state_vars

    nvar = len(list(m.fs.state.component_data_objects(Var)))
    ncon = len(list(m.component_data_objects(Constraint)))

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
        assert len(list(m.fs.state.component_data_objects(Var))) == nvar
        assert len(list(m.component_data_objects(Constraint))) == ncon
