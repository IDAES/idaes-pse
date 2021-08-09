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
Tests for CLC gas phase thermo state block; tests for construction and solves
Author: Chinedu Okoli
"""
import pytest
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    SolverStatus

from Models.solid_phase_thermo import Solid_Phase_Thermo_ParameterBlock
from Models.gas_phase_thermo import Gas_Phase_Thermo_ParameterBlock
from Models.hetero_reactions import HeteroReactionParameterBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
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
        default={"property_package": m.solid_properties})
m.reaction_block = m.reactions.reaction_block_class(
        default={"parameters": m.reactions,
                 "state_block": m.solid_state_block,
                 "gas_state_block": m.gas_state_block,
                 "has_equilibrium": False})


def test_build_reaction_block():
    assert hasattr(m.reaction_block, "k_rxn")
    assert hasattr(m.reaction_block, "OC_conv")
    assert hasattr(m.reaction_block, "OC_conv_temp")
    m.reaction_block._reaction_rate()
    assert hasattr(m.reaction_block, "reaction_rate")


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve():
    # Fix some input variables
    m.reaction_block.state_ref.temperature.fix(1183.15)
    m.reaction_block.state_ref.mass_frac["Fe2O3"].fix(0.45)
    m.reaction_block.state_ref.mass_frac["Fe3O4"].fix(1e-9)
    m.reaction_block.state_ref.mass_frac["Al2O3"].fix(0.55)

    m.reaction_block.gas_state_ref.dens_mole_comp_vap['CH4'].fix(10)
    m.reaction_block.gas_state_ref.temperature.fix(298.15)
#    m.reaction_block.gas_state_ref.pressure.fix(1.60)
    m.reaction_block.gas_state_ref.mole_frac["CO2"].fix(0.4772)
    m.reaction_block.gas_state_ref.mole_frac["H2O"].fix(0.0646)
    m.reaction_block.gas_state_ref.mole_frac["CH4"].fix(0.4582)

    m.reaction_block.initialize()
    results = solver.solve(m.reaction_block, tee=False)
    assert degrees_of_freedom(m) == 0

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok
