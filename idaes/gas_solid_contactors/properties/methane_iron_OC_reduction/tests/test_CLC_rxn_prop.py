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
import sys
import os
import pytest

from pyomo.environ import (ConcreteModel, SolverFactory, TerminationCondition,
                           SolverStatus)

from idaes.core.util.model_statistics import degrees_of_freedom

# Access parent directory (property_packages_subfolder) of the current dir
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

from solid_phase_thermo import Solid_Phase_Thermo_ParameterBlock
from gas_phase_thermo import Gas_Phase_Thermo_ParameterBlock
from hetero_reactions import HeteroReactionParameterBlock

# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
else:
    solver = None

# -----------------------------------------------------------------------------
m = ConcreteModel()

# Set up thermo props and reaction props
m.solid_properties = Solid_Phase_Thermo_ParameterBlock()
m.solid_state_block = m.solid_properties.state_block_class(
    default={"parameters": m.solid_properties,
             "defined_state": True})

m.gas_properties = Gas_Phase_Thermo_ParameterBlock()
m.gas_state_block = m.gas_properties.state_block_class(
    default={"parameters": m.gas_properties,
             "defined_state": True})

m.reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m.solid_properties,
                     "gas_property_package": m.gas_properties})
m.reaction_block = m.reactions.reaction_block_class(
        default={"parameters": m.reactions,
                 "solid_state_block": m.solid_state_block,
                 "gas_state_block": m.gas_state_block,
                 "has_equilibrium": False})


def test_build_reaction_block():
    assert hasattr(m.reaction_block, "k_rxn")
    assert hasattr(m.reaction_block, "OC_conv")
    assert hasattr(m.reaction_block, "OC_conv_temp")
    assert hasattr(m.reaction_block, "reaction_rate")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve():
    # Fix required variables to make reaction model square
    # (gas density, solid density, temperature and component fractions)
    m.reaction_block.gas_state_ref.dens_mole_comp_vap['CH4'].fix(10)
    m.reaction_block.solid_state_ref.dens_mass_sol.fix(
            m.reaction_block.solid_state_ref.dens_mass_sol.value)
    m.reaction_block.solid_state_ref.temperature.fix(1183.15)
    m.reaction_block.solid_state_ref.mass_frac["Fe2O3"].fix(0.45)
    m.reaction_block.solid_state_ref.mass_frac["Fe3O4"].fix(1e-9)
    m.reaction_block.solid_state_ref.mass_frac["Al2O3"].fix(0.55)

    m.reaction_block.initialize()
    results = solver.solve(m.reaction_block, tee=False)
    assert degrees_of_freedom(m.reaction_block) == 0

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok
