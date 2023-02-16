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
Tests for CLC heterogeneous reaction block; tests for construction and solve
Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, Var
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_phase_thermo import (
    GasPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_phase_thermo import (
    SolidPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions import (
    HeteroReactionParameterBlock,
)


# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.fixture(scope="class")
def rxn_prop():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Set up thermo props and reaction props
    m.fs.solid_properties = SolidPhaseParameterBlock()
    m.fs.solid_state_block = m.fs.solid_properties.build_state_block(
        [0], parameters=m.fs.solid_properties, defined_state=True
    )

    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.gas_state_block = m.fs.gas_properties.build_state_block(
        [0], parameters=m.fs.gas_properties, defined_state=True
    )

    m.fs.reactions = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_properties,
        gas_property_package=m.fs.gas_properties,
    )
    m.fs.unit = m.fs.reactions.reaction_block_class(
        [0],
        parameters=m.fs.reactions,
        solid_state_block=m.fs.solid_state_block,
        gas_state_block=m.fs.gas_state_block,
        has_equilibrium=False,
    )

    # Fix required variables to make reaction model square
    # (gas mixture and component densities,
    # solid particle porosity, density and component fractions)
    m.fs.gas_state_block[0].dens_mol.fix(10)
    m.fs.gas_state_block[0].dens_mol_comp.fix(10)
    m.fs.solid_state_block[0].particle_porosity.fix(0.27)
    m.fs.solid_state_block[0].mass_frac_comp["Fe2O3"].fix(0.45)
    m.fs.solid_state_block[0].mass_frac_comp["Fe3O4"].fix(1e-9)
    m.fs.solid_state_block[0].mass_frac_comp["Al2O3"].fix(0.55)
    m.fs.solid_state_block[0].dens_mass_skeletal.fix(1)
    m.fs.solid_state_block[0].temperature.fix(1183.15)  # K

    return m


@pytest.mark.unit
def test_build_reaction_block(rxn_prop):
    assert isinstance(rxn_prop.fs.unit[0].k_rxn, Var)
    assert isinstance(rxn_prop.fs.unit[0].OC_conv, Var)
    assert isinstance(rxn_prop.fs.unit[0].OC_conv_temp, Var)
    assert isinstance(rxn_prop.fs.unit[0].reaction_rate, Var)


@pytest.mark.unit
def test_setInputs_reaction_block(rxn_prop):
    assert degrees_of_freedom(rxn_prop.fs.unit[0]) == 0


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(rxn_prop):
    orig_fixed_vars = fixed_variables_set(rxn_prop)
    orig_act_consts = activated_constraints_set(rxn_prop)

    rxn_prop.fs.unit.initialize()

    assert degrees_of_freedom(rxn_prop) == 0

    fin_fixed_vars = fixed_variables_set(rxn_prop)
    fin_act_consts = activated_constraints_set(rxn_prop)

    assert len(fin_act_consts) == len(orig_act_consts)
    assert len(fin_fixed_vars) == len(orig_fixed_vars)

    for c in fin_act_consts:
        assert c in orig_act_consts
    for v in fin_fixed_vars:
        assert v in orig_fixed_vars


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(rxn_prop):

    assert hasattr(rxn_prop.fs.unit[0], "k_rxn")
    assert hasattr(rxn_prop.fs.unit[0], "OC_conv")
    assert hasattr(rxn_prop.fs.unit[0], "reaction_rate")

    results = solver.solve(rxn_prop.fs.unit[0])

    # Check for optimal solution
    assert check_optimal_termination(results)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(rxn_prop):
    assert pytest.approx(1, abs=1e-2) == rxn_prop.fs.unit[0].k_rxn["R1"].value
    assert pytest.approx(0, abs=1e-2) == rxn_prop.fs.unit[0].OC_conv.value
    assert pytest.approx(0, abs=1e-2) == rxn_prop.fs.unit[0].reaction_rate["R1"].value


@pytest.mark.component
def test_units_consistent(rxn_prop):
    assert_units_consistent(rxn_prop)
