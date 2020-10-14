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
Tests for CLC heterogeneous reaction block; tests for construction and solve
Author: Chinedu Okoli
"""

import pytest

from pyomo.environ import (ConcreteModel,
                           Var)

from idaes.core import FlowsheetBlock

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.core.util.testing import (get_default_solver,
                                     initialization_tester)

from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import GasPhaseThermoParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    solid_phase_thermo import SolidPhaseThermoParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    hetero_reactions import HeteroReactionParameterBlock


# Get default solver for testing
solver = get_default_solver()


# -----------------------------------------------------------------------------
@pytest.fixture(scope="class")
def rxn_prop():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo props and reaction props
    m.fs.solid_properties = SolidPhaseThermoParameterBlock()
    m.fs.solid_state_block = m.fs.solid_properties.build_state_block(
        default={"parameters": m.fs.solid_properties,
                 "defined_state": True})

    m.fs.gas_properties = GasPhaseThermoParameterBlock()
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
    m.fs.solid_state_block.particle_porosity.fix(0.27)
    m.fs.solid_state_block.mass_frac_comp["Fe2O3"].fix(0.45)
    m.fs.solid_state_block.mass_frac_comp["Fe3O4"].fix(1e-9)
    m.fs.solid_state_block.mass_frac_comp["Al2O3"].fix(0.55)
    m.fs.solid_state_block.dens_mass_skeletal.fix(1)

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
