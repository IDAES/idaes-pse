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
Carbon purification unit (CPU) model test

Created on May 18 2022 by A. Noring
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import unit model
from idaes.models_extra.power_generation.unit_models.cpu import CarbonProcessingUnit

from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def model():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.cpu = CarbonProcessingUnit()

    return m


@pytest.mark.unit
def test_basic_build(model):
    """Make the cpu model and make sure it doesn't throw an exception"""
    assert degrees_of_freedom(model) == 6


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_cpu(model):
    # fix inputs
    model.fs.cpu.inlet.flow_mol.fix(1260)
    model.fs.cpu.inlet.mole_frac_comp[0, "Ar"].fix(0.005)
    model.fs.cpu.inlet.mole_frac_comp[0, "CO2"].fix(0.900)
    model.fs.cpu.inlet.mole_frac_comp[0, "H2O"].fix(0.05)
    model.fs.cpu.inlet.mole_frac_comp[0, "N2"].fix(0.025)
    model.fs.cpu.inlet.mole_frac_comp[0, "O2"].fix(0.02)

    optarg = {"tol": 1e-7, "linear_solver": "ma27", "max_iter": 40}
    solver.options = optarg

    model.fs.cpu.initialize(optarg=optarg, release_state=False)

    assert degrees_of_freedom(model) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_cpu(model):
    results = solver.solve(model, tee=True)

    # test outlet flowrates, CO2 purity, heat duty, and work
    assert pytest.approx(1103.3, abs=1e-1) == pyo.value(
        model.fs.cpu.pureco2.flow_mol[0]
    )
    assert pytest.approx(37.1, abs=1e-1) == pyo.value(model.fs.cpu.water.flow_mol[0])
    assert pytest.approx(119.6, abs=1e-1) == pyo.value(model.fs.cpu.vent.flow_mol[0])

    assert pyo.value(model.fs.cpu.pureco2.mole_frac_comp[0, "CO2"]) > 0.99

    assert pytest.approx(35262510, abs=1) == pyo.value(model.fs.cpu.heat_duty[0])
    assert pytest.approx(21402861.6, abs=1) == pyo.value(model.fs.cpu.work[0])

    assert degrees_of_freedom(model) == 0
    assert pyo.check_optimal_termination(results)
