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
        units as pyunits,
        )

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom

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


if __name__ == "__main__":
    test_state_vars()
    m = make_model()
    print(degrees_of_freedom(m))
    for var in m.component_data_objects(Var):
        print(var.name, var.fixed)

    for con in m.component_data_objects(Constraint, active=True):
        print(con.name)

    # At this point, the model is very uninteresting. We have only queried
    # variables that have "method": None. This means that no constraints
    # have been added to the state block.
    import pdb; pdb.set_trace()
