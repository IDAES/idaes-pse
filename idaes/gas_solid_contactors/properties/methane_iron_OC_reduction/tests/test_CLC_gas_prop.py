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
import sys
import os
import pytest

from pyomo.environ import (ConcreteModel, SolverFactory, TerminationCondition,
                           SolverStatus, value)

from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import GasPhaseThermoParameterBlock

# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
else:
    solver = None

# -----------------------------------------------------------------------------
m = ConcreteModel()

# gas properties and state inlet block
m.properties = GasPhaseThermoParameterBlock()
m.state_block = m.properties.state_block_class(
    default={"parameters": m.properties,
             "defined_state": True})


def test_build_inlet_state_block():
    assert hasattr(m.state_block, "mw_gas")
    assert hasattr(m.state_block, "dens_mole_vap")
    assert hasattr(m.state_block, "dens_mole_comp_vap")
    assert hasattr(m.state_block, "dens_mass_vap")
    assert hasattr(m.state_block, "cp_mol_comp")
    assert hasattr(m.state_block, "cp_mol")
    assert hasattr(m.state_block, "cp_mass")
    assert hasattr(m.state_block, "visc_d")


def test_setInputs_state_block():

    m.state_block.flow_mol.fix(1)
    m.state_block.temperature.fix(450)
    m.state_block.pressure.fix(1.60)
    m.state_block.mole_frac["CO2"].fix(0.4772)
    m.state_block.mole_frac["H2O"].fix(0.0646)
    m.state_block.mole_frac["CH4"].fix(0.4582)

    assert degrees_of_freedom(m.state_block) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve():
    # Initialize and solve block
    m.state_block.initialize()
    results = solver.solve(m.state_block, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    # Check if values are as expected
    assert value(m.state_block.mole_frac['CO2']) == \
        pytest.approx(0.4772, abs=1e-3)
    assert value(m.state_block.mole_frac['CH4']) == \
        pytest.approx(0.4582, abs=1e-3)
    assert degrees_of_freedom(m) == 0
