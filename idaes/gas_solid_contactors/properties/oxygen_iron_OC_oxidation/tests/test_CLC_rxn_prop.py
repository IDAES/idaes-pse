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

from pyomo.environ import (
        ConcreteModel,
        TerminationCondition,
        SolverStatus,
        Var,
        Constraint,
        value,
        units as pyunits,
        )

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.testing import (get_default_solver,
                                     initialization_tester)

from idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation. \
    gas_phase_thermo import GasPhaseParameterBlock
from idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation. \
    solid_phase_thermo import SolidPhaseParameterBlock
from idaes.gas_solid_contactors.properties.oxygen_iron_OC_oxidation. \
    hetero_reactions import HeteroReactionParameterBlock


# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
@pytest.fixture(scope="class")
def rxn_prop():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo props and reaction props
    m.fs.solid_properties = SolidPhaseParameterBlock()
    m.fs.solid_state_block = m.fs.solid_properties.build_state_block(
        default={"parameters": m.fs.solid_properties,
                 "defined_state": True})

    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.gas_state_block = m.fs.gas_properties.build_state_block(
        default={"parameters": m.fs.gas_properties,
                 "defined_state": True})

    m.fs.reactions = HeteroReactionParameterBlock(
                default={"solid_property_package": m.fs.solid_properties,
                         "gas_property_package": m.fs.gas_properties})
    m.fs.unit = m.fs.reactions.reaction_block_class(
            default={"parameters": m.fs.reactions,
                     "solid_state_block": m.fs.solid_state_block,
                     "gas_state_block": m.fs.gas_state_block,
                     "has_equilibrium": False})

    # Fix required variables to make reaction model square
    # (gas mixture and component densities,
    # solid particle porosity, density and component fractions)
    m.fs.gas_state_block.dens_mol.fix(10)
    m.fs.gas_state_block.dens_mol_comp.fix(10)
    m.fs.solid_state_block.particle_porosity.fix(0.23)
    m.fs.solid_state_block.mass_frac_comp["Fe2O3"].fix(0.244)
    m.fs.solid_state_block.mass_frac_comp["Fe3O4"].fix(0.202)
    m.fs.solid_state_block.mass_frac_comp["Al2O3"].fix(0.554)
    m.fs.solid_state_block.dens_mass_skeletal.fix(3125)
    m.fs.solid_state_block.temperature.fix(1183.15)  # K

    return m


@pytest.mark.unit
def test_build_reaction_block(rxn_prop):
    assert isinstance(rxn_prop.fs.unit.k_rxn, Var)
    assert isinstance(rxn_prop.fs.unit.OC_conv, Var)
    assert isinstance(rxn_prop.fs.unit.OC_conv_temp, Var)
    assert isinstance(rxn_prop.fs.unit.reaction_rate, Var)


@pytest.mark.unit
def test_setInputs_reaction_block(rxn_prop):
    assert degrees_of_freedom(rxn_prop.fs.unit) == 0


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(rxn_prop):
    initialization_tester(
            rxn_prop)


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(rxn_prop):

    assert hasattr(rxn_prop.fs.unit, "k_rxn")
    assert hasattr(rxn_prop.fs.unit, "OC_conv")
    assert hasattr(rxn_prop.fs.unit, "reaction_rate")

    results = solver.solve(rxn_prop.fs.unit)

    # Check for optimal solution
    assert results.solver.termination_condition == \
        TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solution(rxn_prop):
    assert (pytest.approx(1, abs=1e-2) ==
            rxn_prop.fs.unit.k_rxn['R1'].value)
    assert (pytest.approx(0, abs=1e-2) ==
            rxn_prop.fs.unit.OC_conv.value)
    assert (pytest.approx(0, abs=1e-2) ==
            rxn_prop.fs.unit.reaction_rate['R1'].value)


def test_state_vars():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.solid_properties = SolidPhaseParameterBlock()
    m.fs.solid_state = m.fs.solid_properties.build_state_block(default={
        "parameters": m.fs.solid_properties,
        "defined_state": True,
        })

    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.gas_state = m.fs.gas_properties.build_state_block(default={
        "parameters": m.fs.gas_properties,
        "defined_state": True,
        })

    m.fs.reaction_properties = HeteroReactionParameterBlock(default={
        "solid_property_package": m.fs.solid_properties,
        "gas_property_package": m.fs.gas_properties,
        })
    m.fs.reaction_block = m.fs.reaction_properties.reaction_block_class(
            default={
                "parameters": m.fs.reaction_properties,
                "solid_state_block": m.fs.solid_state,
                "gas_state_block": m.fs.gas_state,
                "has_equilibrium": False,
            })

    for var in m.fs.gas_state.define_state_vars().values():
        var.fix()
    for var in m.fs.solid_state.define_state_vars().values():
        var.fix()

    # Note that these checks are necessary to trigger the construction of
    # the reaction block variables
    assert isinstance(m.fs.reaction_block.k_rxn, Var)
    assert isinstance(m.fs.reaction_block.OC_conv, Var)
    assert isinstance(m.fs.reaction_block.OC_conv_temp, Var)
    assert isinstance(m.fs.reaction_block.reaction_rate, Var)

    assert degrees_of_freedom(m) == 0

    rxn_vars = list(m.fs.reaction_block.component_data_objects(Var))
    rxn_cons = list(m.fs.reaction_block.component_data_objects(Constraint))
    assert len(rxn_vars) == 4
    assert len(rxn_cons) == 4


if __name__ == "__main__":
    test_state_vars()
