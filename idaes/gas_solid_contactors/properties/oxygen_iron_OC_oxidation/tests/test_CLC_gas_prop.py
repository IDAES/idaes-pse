##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for CLC gas phase thermo state block; tests for construction and solve
Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import (
        ConcreteModel,
        Constraint,
        TerminationCondition,
        SolverStatus,
        Var,
        Reference,
        value,
        units as pyunits,
        )
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.common.unittest import TestCase
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.util.subsystems import ParamSweeper

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import (
        degrees_of_freedom,
        number_large_residuals,
        large_residuals_set,
        )

from idaes.core.util.testing import (get_default_solver,
                                     initialization_tester)

from idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation. \
    gas_phase_thermo import GasPhaseParameterBlock

# Get default solver for testing
solver = get_default_solver()


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
    m.fs.unit.pressure.fix(1.60)
    m.fs.unit.mole_frac_comp["CO2"].fix(0.0004)
    m.fs.unit.mole_frac_comp["H2O"].fix(0.0093)
    m.fs.unit.mole_frac_comp["O2"].fix(0.2095)
    m.fs.unit.mole_frac_comp["N2"].fix(0.7808)

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


@pytest.mark.unit
def test_setInputs_state_block(gas_prop):
    assert degrees_of_freedom(gas_prop.fs.unit) == 0


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
    assert hasattr(gas_prop.fs.unit, "cp_mol")
    assert hasattr(gas_prop.fs.unit, "enth_mol")

    results = solver.solve(gas_prop)

    # Check for optimal solution
    assert results.solver.termination_condition == \
        TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(gas_prop):
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.mw.value)
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.cp_mol.value)
    assert (pytest.approx(1, abs=1e-2) ==
            gas_prop.fs.unit.enth_mol.value)
# -----------------------------------------------------------------------------


def make_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = GasPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block(
            default={
                "parameters": m.fs.properties,
                "defined_state": True,
                })

    m.fs.state.flow_mol.set_value(1.0*pyunits.mol/pyunits.s)
    m.fs.state.temperature.set_value(450.0*pyunits.K)
    m.fs.state.pressure.set_value(1.60*pyunits.bar)
    m.fs.state.mole_frac_comp["CO2"].set_value(0.0004)
    m.fs.state.mole_frac_comp["H2O"].set_value(0.0093)
    m.fs.state.mole_frac_comp["O2"].set_value(0.2095)
    m.fs.state.mole_frac_comp["N2"].set_value(0.7808)

    return m


def test_state_vars():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

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


def test_indexed_state_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block([1,2,3])

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


def test_property_construction_ordered():
    """
    The purpose of this test is to make sure each property can be
    constructed properly. Further, because we know that these
    variables form a DAG, we check that the addition of each new
    property adds one variable and constraint, i.e. that the
    matching defined in the parameter block is in topological order.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = GasPhaseParameterBlock()
    m.fs.state = m.fs.properties.build_state_block()

    nvar = 7
    ncon = 1
    assert len(list(m.component_data_objects(Var))) == nvar
    assert len(list(m.component_data_objects(Constraint))) == ncon

    # Make sure the matching captures all the non-state vars.
    state_vars = m.fs.state.define_state_vars()
    n_state_vars = len(state_vars)
    n_vars = len(m.fs.properties._metadata._properties)
    assert len(m.fs.properties.matching) == n_vars - n_state_vars

    # Make sure we can add constraints one at a time, and only
    # add one additional variable. This is specific to this
    # property package.
    for varname, conname in m.fs.properties.matching:
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


def test_property_construction_unordered():
    """
    The purpose of this test is to make sure "chained construction"
    of properties happens properly. I.e. if we add a property that
    requires several other properties, these are added as well.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
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


def _create_subsystem_block(variables, constraints):
    block = ConcreteModel()
    block.vars = Reference(variables)
    block.cons = Reference(constraints)
    var_set = ComponentSet(variables)
    other_vars = []
    for con in constraints:
        for var in identify_variables(con.body, include_fixed=False):
            # Fixed vars will not be picked up by writer, so no point
            # including them. If we use this block for other purposes
            # than solving, this may change.
            if var not in var_set:
                other_vars.append(var)
                var_set.add(var)
    block.other_vars = Reference(other_vars)
    return block


def _solve_strongly_connected_components(block):
    variables = [var for var in block.component_data_objects(Var)
            if not var.fixed]
    constraints = [con for con in block.component_data_objects(Constraint)
            if con.active]
    assert len(variables) == len(constraints)
    igraph = IncidenceGraphInterface()
    var_block_map, con_block_map = igraph.block_triangularize(
            variables=variables,
            constraints=constraints,
            )
    blocks = set(var_block_map.values())
    n_blocks = len(blocks)
    var_blocks = [[] for b in range(n_blocks)]
    con_blocks = [[] for b in range(n_blocks)]
    for var, b in var_block_map.items():
        var_blocks[b].append(var)
    for con, b in con_block_map.items():
        con_blocks[b].append(con)
    for b_vars, b_cons in zip(var_blocks, con_blocks):
        assert len(b_vars) == len(b_cons)
        if len(b_vars) == 1:
            var = b_vars[0]
            con = b_cons[0]
            calculate_variable_from_constraint(var, con)
        else:
            _temp = _create_subsystem_block(b_vars, b_cons)
            _temp.other_vars[:].fix()
            solver.solve(_temp)
            _temp.other_vars[:].unfix()


class TestProperties(TestCase):
    """
    The purpose of these tests is to ensure that the property package's
    state block can be used to solve for the specified properties, and
    that the values produced are what we expect.
    """
    def _make_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = GasPhaseParameterBlock()
        m.fs.state = m.fs.properties.build_state_block(default={
            # defined_state True because we want to ensure the state
            # block is square.
            "defined_state": True,
            })
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
                "flow_mol": [1.0*pyunits.mol/pyunits.s]*n_scenario,
                "temperature": [300.0*pyunits.K]*n_scenario,
                "pressure": [1.0*pyunits.bar]*n_scenario,
                "mole_frac_comp[O2]":  [1.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.0],
                "mole_frac_comp[N2]":  [0.0, 0.5, 0.25, 1.0, 0.0, 0.0, 0.0],
                "mole_frac_comp[H2O]": [0.0, 0.0, 0.25, 0.0, 1.0, 0.0, 0.5],
                "mole_frac_comp[CO2]": [0.0, 0.0, 0.25, 0.0, 0.0, 1.0, 0.5],
                }
        state_values = ComponentMap((state.find_component(name), values)
                for name, values in state_values.items())
        target_values = {
                "mw": [
                    32.0*pyunits.g/pyunits.mol,
                    30.0*pyunits.g/pyunits.mol,
                    30.5*pyunits.g/pyunits.mol,
                    28.0*pyunits.g/pyunits.mol,
                    18.0*pyunits.g/pyunits.mol,
                    44.0*pyunits.g/pyunits.mol,
                    31.0*pyunits.g/pyunits.mol,
                    ]
                }
        target_values = ComponentMap((state.find_component(name), values)
                for name, values in target_values.items())

        # Construct mw and all prerequisites
        state.mw

        param_sweeper = ParamSweeper(
                n_scenario,
                state_values,
                output_values=target_values,
                )
        with param_sweeper:
            for inputs, target in param_sweeper:
                _solve_strongly_connected_components(state)

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
                "flow_mol": [1.0*pyunits.mol/pyunits.s]*n_scen,
                "temperature": [
                    200.0*pyunits.K,
                    300.0*pyunits.K,
                    600.0*pyunits.K,
                    900.0*pyunits.K,
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    ],
                "pressure": [
                    1.0*pyunits.bar,
                    1.0*pyunits.bar,
                    1.0*pyunits.bar,
                    1.0*pyunits.bar,
                    1.0*pyunits.bar,
                    0.5*pyunits.bar,
                    1.5*pyunits.bar,
                    2.0*pyunits.bar,
                    ],
                "mole_frac_comp[O2]":  [0.25]*n_scen,
                "mole_frac_comp[N2]":  [0.25]*n_scen,
                "mole_frac_comp[H2O]": [0.25]*n_scen,
                "mole_frac_comp[CO2]": [0.25]*n_scen,
                }
        state_values = ComponentMap((state.find_component(name), values)
                for name, values in state_values.items())

        u = pyunits.mol/pyunits.m**3
        target_values = {
                "dens_mol": [
                    # Calculated by P*1e5/(T*8.314462618)
                    60.136*u, 40.091*u, 20.045*u, 13.364*u, 10.023*u,
                    5.011*u, 15.034*u, 20.045*u,
                    ]
                }
        target_values = ComponentMap((state.find_component(name), values)
                for name, values in target_values.items())

        # Construct dens_mol and all prerequisites
        state.dens_mol

        param_sweeper = ParamSweeper(n_scen, state_values,
                output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                _solve_strongly_connected_components(state)

                # Check that state block equations have been convered
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are what we expect
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(val, abs=1e-3)

                # Check that state block performs the calculations we expect
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(val, abs=1e-3)

    def test_dens_mol_comp(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 5
        state_values = {
                "flow_mol": [1.0*pyunits.mol/pyunits.s]*n_scen,
                "temperature": [1200.0*pyunits.K]*n_scen,
                "pressure": [1.0*pyunits.bar]*n_scen,
                "mole_frac_comp[O2]":  [1.0, 0.0, 0.0, 0.0, 0.25],
                "mole_frac_comp[N2]":  [0.0, 1.0, 0.0, 0.0, 0.25],
                "mole_frac_comp[H2O]": [0.0, 0.0, 1.0, 0.0, 0.25],
                "mole_frac_comp[CO2]": [0.0, 0.0, 0.0, 1.0, 0.25],
                }
        state_values = ComponentMap((state.find_component(name), values)
                for name, values in state_values.items())

        u = pyunits.mol/pyunits.m**3
        target_values = {
                "dens_mol": [10.023*u]*n_scen,
                "dens_mol_comp[O2]": [
                    10.023*u, 0.0*u, 0.0*u, 0.0*u, 0.25*10.023*u,
                    ],
                "dens_mol_comp[N2]": [
                    0.0*u, 10.023*u, 0.0*u, 0.0*u, 0.25*10.023*u,
                    ],
                "dens_mol_comp[H2O]": [
                    0.0*u, 0.0*u, 10.023*u, 0.0*u, 0.25*10.023*u,
                    ],
                "dens_mol_comp[CO2]": [
                    0.0*u, 0.0*u, 0.0*u, 10.023*u, 0.25*10.023*u,
                    ],
                }
        target_values = ComponentMap((state.find_component(name), values)
                for name, values in target_values.items())

        # Construct dens_mol_comp and all prerequisites
        state.dens_mol_comp

        param_sweeper = ParamSweeper(n_scen, state_values,
                output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                _solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs have been set properly
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    # ^ Problem converting units when value is zero
                    assert var.value == pytest.approx(value(val), abs=1e-3)

                # Make sure property values are what we expect
                for var, val in target.items():
                    #val = value(pyunits.convert(val, var.get_units()))
                    # ^ Problem converting units when value is zero
                    assert var.value == pytest.approx(value(val), abs=1e-3)

    def test_visc_d(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 8
        state_values = {
                "flow_mol": [1.0*pyunits.mol/pyunits.s]*n_scen,
                "temperature": [
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    300.0*pyunits.K,
                    600.0*pyunits.K,
                    900.0*pyunits.K,
                    1200.0*pyunits.K,
                    ],
                "pressure": [1.0*pyunits.bar]*n_scen,
                "mole_frac_comp[O2]":  [
                    1.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
                "mole_frac_comp[N2]":  [
                    0.0, 1.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
                "mole_frac_comp[H2O]": [
                    0.0, 0.0, 1.0, 0.0, 0.25, 0.25, 0.25, 0.25],
                "mole_frac_comp[CO2]": [
                    0.0, 0.0, 0.0, 1.0, 0.25, 0.25, 0.25, 0.25],
                }
        state_values = ComponentMap((state.find_component(name), values)
                for name, values in state_values.items())

        target_values = {"visc_d": [
            # These values were copied after a solve.
            # TODO: Cross-reference with another source.
            5.534440949133228e-05*pyunits.Pa*pyunits.s,
            4.67667824296429e-05*pyunits.Pa*pyunits.s,
            4.6232771210527155e-05*pyunits.Pa*pyunits.s,
            4.512867970060493e-05*pyunits.Pa*pyunits.s,
            1.6181534595116313e-05*pyunits.Pa*pyunits.s,
            2.866222939903063e-05*pyunits.Pa*pyunits.s,
            3.909320395131273e-05*pyunits.Pa*pyunits.s,
            4.838841106600266e-05*pyunits.Pa*pyunits.s,
            ]}
        target_values = ComponentMap((state.find_component(name), values)
                for name, values in target_values.items())

        # Construct visc_d and all prerequisites
        state.visc_d

        param_sweeper = ParamSweeper(n_scen, state_values,
                output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                _solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), abs=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), abs=1e-3)

    def test_diffusion_comp(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 11
        bar = pyunits.bar
        K = pyunits.K
        state_values = {
                "flow_mol": [1.0*pyunits.mol/pyunits.s]*n_scen,
                "temperature": [
                    1200.0*K, 1200.0*K, 1200.0*K, 1200.0*K,
                    300.0*K, 600.0*K, 900.0*K, 1200.0*K,
                    1200.0*K, 1200.0*K, 1200.0*K,
                    ],
                "pressure": [
                    1.0*bar, 1.0*bar, 1.0*bar, 1.0*bar,
                    1.0*bar, 1.0*bar, 1.0*bar, 1.0*bar,
                    0.5*bar, 1.5*bar, 2.0*bar,
                    ],
                "mole_frac_comp[O2]":  [
                    # Note that diffusivity is not defined for a pure
                    # component in itself (zero gradient, zero net diffusion)
                    0.90, 0.025, 0.025, 0.025, 0.25, 0.25, 0.25, 0.25,
                    0.25, 0.25, 0.25],
                "mole_frac_comp[N2]":  [
                    0.025, 0.90, 0.025, 0.025, 0.25, 0.25, 0.25, 0.25,
                    0.25, 0.25, 0.25],
                "mole_frac_comp[H2O]": [
                    0.025, 0.025, 0.90, 0.025, 0.25, 0.25, 0.25, 0.25,
                    0.25, 0.25, 0.25],
                "mole_frac_comp[CO2]": [
                    0.025, 0.025, 0.025, 0.90, 0.25, 0.25, 0.25, 0.25,
                    0.25, 0.25, 0.25],
                }
        state_values = ComponentMap((state.find_component(name), values)
                for name, values in state_values.items())

        cm = pyunits.cm
        s = pyunits.s
        target_values = {
                # These values look reasonable
                # TODO: Verify with external source.
                "diffusion_comp[O2]": [
                    3.11792621830951*cm**2/s,
                    2.456227751888218*cm**2/s,
                    3.034740091620132*cm**2/s,
                    1.9479388404541838*cm**2/s,
                    0.206691259894311*cm**2/s,
                    0.6952237580376006*cm**2/s,
                    1.4134625566198689*cm**2/s,
                    2.3384446637321363*cm**2/s,
                    4.676889327464272*cm**2/s,
                    1.5589631091547582*cm**2/s,
                    1.1692223318660684*cm**2/s,
                    ],
                "diffusion_comp[N2]": [
                    2.457518754629481*cm**2/s,
                    3.1383309495574956*cm**2/s,
                    3.0334118458311523*cm**2/s,
                    1.9790010245966736*cm**2/s,
                    0.2080439152537244*cm**2/s,
                    0.6997735302088169*cm**2/s,
                    1.422712718932098*cm**2/s,
                    2.353748212168125*cm**2/s,
                    4.70749642433625*cm**2/s,
                    1.5691654747787498*cm**2/s,
                    1.176874106084062*cm**2/s,
                    ],
                "diffusion_comp[H2O]": [
                    3.0845168350215713*cm**2/s,
                    3.0811129521516416*cm**2/s,
                    3.719143107603082*cm**2/s,
                    2.5051274838003996*cm**2/s,
                    0.24654668546150185*cm**2/s,
                    0.8292808959890481*cm**2/s,
                    1.6860147281349098*cm**2/s,
                    2.7893573307023156*cm**2/s,
                    5.578714661404631*cm**2/s,
                    1.8595715538015434*cm**2/s,
                    1.3946786653511578*cm**2/s,
                    ],
                "diffusion_comp[CO2]": [
                    1.929322600571961*cm**2/s,
                    1.9589681584041365*cm**2/s,
                    2.4422863072776697*cm**2/s,
                    2.7098612589802444*cm**2/s,
                    0.17964011927809195*cm**2/s,
                    0.6042349293467888*cm**2/s,
                    1.2284727588198259*cm**2/s,
                    2.0323959442351844*cm**2/s,
                    4.064791888470369*cm**2/s,
                    1.3549306294901227*cm**2/s,
                    1.0161979721175922*cm**2/s,
                    ],
                }
        target_values = ComponentMap((state.find_component(name), values)
                for name, values in target_values.items())

        # Construct diffusion_comp and all prerequisites
        state.diffusion_comp

        param_sweeper = ParamSweeper(n_scen, state_values,
                output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                _solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), abs=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), abs=1e-3)

    def test_therm_cond(self):
        m = self._make_model()
        state = m.fs.state

        n_scen = 8
        state_values = {
                "flow_mol": [1.0*pyunits.mol/pyunits.s]*n_scen,
                "temperature": [
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    1200.0*pyunits.K,
                    300.0*pyunits.K,
                    600.0*pyunits.K,
                    900.0*pyunits.K,
                    1200.0*pyunits.K,
                    ],
                "pressure": [1.0*pyunits.bar]*n_scen,
                "mole_frac_comp[O2]":  [
                    1.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
                "mole_frac_comp[N2]":  [
                    0.0, 1.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25],
                "mole_frac_comp[H2O]": [
                    0.0, 0.0, 1.0, 0.0, 0.25, 0.25, 0.25, 0.25],
                "mole_frac_comp[CO2]": [
                    0.0, 0.0, 0.0, 1.0, 0.25, 0.25, 0.25, 0.25],
                }
        state_values = ComponentMap((state.find_component(name), values)
                for name, values in state_values.items())

        u = pyunits.kJ/pyunits.m/pyunits.K/pyunits.s
        target_values = {
                "therm_cond": [
                    8.490673044994837e-05*u,
                    7.803113104821915e-05*u,
                    0.0001245121534187936*u,
                    7.844692969560201e-05*u,
                    2.1981943936613706e-05*u,
                    4.567583423706824e-05*u,
                    6.946515568649932e-05*u,
                    9.30078254960681e-05*u,
                    ],
                }
        target_values = ComponentMap((state.find_component(name), values)
                for name, values in target_values.items())

        # Construct therm_cond and all prerequisites
        state.therm_cond

        param_sweeper = ParamSweeper(n_scen, state_values,
                output_values=target_values)
        with param_sweeper:
            for inputs, target in param_sweeper:
                _solve_strongly_connected_components(state)

                # Make sure property equations have been converged
                assert number_large_residuals(state, tol=1e-8) == 0

                # Sanity check that inputs are properly set
                for var, val in inputs.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), abs=1e-3)

                # Make sure properties have been calculated as expected
                for var, val in target.items():
                    val = value(pyunits.convert(val, var.get_units()))
                    assert var.value == pytest.approx(value(val), abs=1e-3)


if __name__ == "__main__":
    test_state_vars()
    test_indexed_state_block()
    test_property_construction_ordered()
    test_property_construction_unordered()
    test = TestProperties()
    test.test_mw()
    test.test_dens_mol()
    test.test_dens_mol_comp()
    test.test_visc_d()
    test.test_diffusion_comp()
    test.test_therm_cond()
