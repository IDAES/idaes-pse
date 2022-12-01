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
Simplified Flash Unit Model, only for IAPWS with mixed state.
Phase separator: inlet water/steam mixture is separated into liquid and vapor
streams.

Created: August 21 2020

Created on Thu Aug 18 12:59:50 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.models.properties import iapws95

from idaes.models_extra.power_generation.unit_models.helm.phase_separator import (
    HelmPhaseSeparator,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture
def build_phase_separator():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = HelmPhaseSeparator(property_package=m.fs.properties)
    return m


@pytest.mark.unit
def test_basic_build(build_phase_separator):
    """Make a model and make sure it doesn't throw exception"""
    m = build_phase_separator
    assert degrees_of_freedom(m) == 3
    # Check unit config arguments
    assert len(m.fs.unit.config) == 4
    assert m.fs.unit.config.property_package is m.fs.properties


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(build_phase_separator):
    state_args_water_steam = {
        "flow_mol": 1.5e5,  # mol/s
        "pressure": 1.2e7,  # Pa
        "enth_mol": 28365.2608,
    }  # j/mol
    initialization_tester(
        build_phase_separator, dof=3, state_args_water_steam=state_args_water_steam
    )


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_wflash(build_phase_separator):
    m = build_phase_separator
    state_args_water_steam = {
        "flow_mol": 1.5e5,  # mol/s
        "pressure": 1.2e7,  # Pa
        "enth_mol": 28365.2608,
    }  # j/mol
    m.fs.unit.initialize(state_args_water_steam=state_args_water_steam)
    m.fs.unit.inlet.flow_mol.fix()
    m.fs.unit.inlet.enth_mol.fix()
    m.fs.unit.inlet.pressure.fix()
    solver.solve(m, tee=True)
    Fin = pyo.value(m.fs.unit.mixed_state[0].flow_mol)
    FL = pyo.value(m.fs.unit.liq_state[0].flow_mol)
    FV = pyo.value(m.fs.unit.vap_state[0].flow_mol)
    assert degrees_of_freedom(m) == 0
    assert pytest.approx((FL + FV), abs=1e-3) == Fin

    assert pytest.approx(
        pyo.value(
            m.fs.unit.mixed_state[0].flow_mol * m.fs.unit.mixed_state[0].vapor_frac
        ),
        abs=1e-3,
    ) == pyo.value(m.fs.unit.vap_state[0].flow_mol)

    assert pytest.approx(
        pyo.value(m.fs.unit.mixed_state[0].enth_mol_phase["Vap"]), abs=1e-3
    ) == pyo.value(m.fs.unit.vap_state[0].enth_mol)
