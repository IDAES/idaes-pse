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
Tests for the SOFC deep neural net ROM.

Author: A. Noring
"""

import pytest

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints)
from idaes.core.solvers import get_solver
from idaes.core.util.testing import initialization_tester

from idaes.models_extra.power_generation.properties.sofc.sofc_rom_builder \
    import SofcRom

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_rom():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.sofc_rom = SofcRom()

    return m


@pytest.mark.unit
def test_basic_build(build_rom):

    m = build_rom

    assert degrees_of_freedom(m) == 9
    assert number_variables(m) == 1500
    assert number_total_constraints(m) == 1491

    assert hasattr(m.fs.sofc_rom, "input")
    assert hasattr(m.fs.sofc_rom, "norm_input")
    assert hasattr(m.fs.sofc_rom, "norm_output")
    assert hasattr(m.fs.sofc_rom, "output")

    assert hasattr(m.fs.sofc_rom, "current_density")
    assert hasattr(m.fs.sofc_rom, "fuel_temperature")
    assert hasattr(m.fs.sofc_rom, "internal_reforming")
    assert hasattr(m.fs.sofc_rom, "air_temperature")
    assert hasattr(m.fs.sofc_rom, "air_recirculation")
    assert hasattr(m.fs.sofc_rom, "OTC")
    assert hasattr(m.fs.sofc_rom, "fuel_util")
    assert hasattr(m.fs.sofc_rom, "air_util")
    assert hasattr(m.fs.sofc_rom, "pressure")

    assert hasattr(m.fs.sofc_rom, "anode_outlet_temperature")
    assert hasattr(m.fs.sofc_rom, "cathode_outlet_temperature")
    assert hasattr(m.fs.sofc_rom, "stack_voltage")
    assert hasattr(m.fs.sofc_rom, "max_cell_temperature")
    assert hasattr(m.fs.sofc_rom, "deltaT_cell")


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(build_rom):

    m = build_rom

    m.fs.sofc_rom.current_density.fix(4000)
    m.fs.sofc_rom.fuel_temperature.fix(348.3)
    m.fs.sofc_rom.internal_reforming.fix(0.6)
    m.fs.sofc_rom.air_temperature.fix(617.3)
    m.fs.sofc_rom.air_recirculation.fix(0.5)
    m.fs.sofc_rom.OTC.fix(2.1)
    m.fs.sofc_rom.fuel_util.fix(0.8)
    m.fs.sofc_rom.air_util.fix(0.449)
    m.fs.sofc_rom.pressure.fix(1)

    initialization_tester(m, unit=m.fs.sofc_rom, dof=0)


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(build_rom):
    m = build_rom

    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)

    assert pyo.check_optimal_termination(results)

    assert pytest.approx(703.4, abs=1e-1) == pyo.value(
        m.fs.sofc_rom.anode_outlet_temperature)

    assert pytest.approx(724.3, abs=1e-1) == pyo.value(
        m.fs.sofc_rom.cathode_outlet_temperature)

    assert pytest.approx(0.875, abs=1e-3) == pyo.value(
        m.fs.sofc_rom.stack_voltage)

    assert pytest.approx(748.4, abs=1e-1) == pyo.value(
        m.fs.sofc_rom.max_cell_temperature)

    assert pytest.approx(99.2, abs=1e-1) == pyo.value(
        m.fs.sofc_rom.deltaT_cell)
