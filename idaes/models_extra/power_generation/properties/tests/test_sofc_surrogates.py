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
Tests for the SOFC keras surrogates.

Author: A. Noring
"""

import pytest

pytest.importorskip(
    "tensorflow",
    reason="tensorflow must be installed to enable SOFC keras surrogates tests",
)
# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.testing import initialization_tester

from idaes.models_extra.power_generation.properties.sofc.sofc_keras_surrogate import (
    SofcSurrogate,
)

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_rom():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.sofc = SofcSurrogate()

    return m


@pytest.mark.unit
def test_basic_build(build_rom):

    m = build_rom

    assert degrees_of_freedom(m) == 9

    assert hasattr(m.fs.sofc, "current_density")
    assert hasattr(m.fs.sofc, "fuel_temperature")
    assert hasattr(m.fs.sofc, "internal_reforming")
    assert hasattr(m.fs.sofc, "air_temperature")
    assert hasattr(m.fs.sofc, "air_recirculation")
    assert hasattr(m.fs.sofc, "otc_ratio")
    assert hasattr(m.fs.sofc, "fuel_utilization")
    assert hasattr(m.fs.sofc, "air_utilization")
    assert hasattr(m.fs.sofc, "pressure")

    assert hasattr(m.fs.sofc, "anode_outlet_temperature")
    assert hasattr(m.fs.sofc, "cathode_outlet_temperature")
    assert hasattr(m.fs.sofc, "stack_voltage")
    assert hasattr(m.fs.sofc, "max_cell_temperature")
    assert hasattr(m.fs.sofc, "delta_cell_temperature")


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(build_rom):

    m = build_rom

    m.fs.sofc.current_density.fix(4000)
    m.fs.sofc.fuel_temperature.fix(621.45)
    m.fs.sofc.internal_reforming.fix(0.6)
    m.fs.sofc.air_temperature.fix(890.45)
    m.fs.sofc.air_recirculation.fix(0.5)
    m.fs.sofc.otc_ratio.fix(2.1)
    m.fs.sofc.fuel_utilization.fix(0.8)
    m.fs.sofc.air_utilization.fix(0.449)
    m.fs.sofc.pressure.fix(1)

    initialization_tester(m, unit=m.fs.sofc, dof=0)


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve(build_rom):
    m = build_rom

    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)

    assert pyo.check_optimal_termination(results)

    assert pytest.approx(976.9, abs=1e-1) == pyo.value(
        m.fs.sofc.anode_outlet_temperature[0]
    )

    assert pytest.approx(998.5, abs=1e-1) == pyo.value(
        m.fs.sofc.cathode_outlet_temperature[0]
    )

    assert pytest.approx(0.876, abs=1e-3) == pyo.value(m.fs.sofc.stack_voltage[0])

    assert pytest.approx(1022.1, abs=1e-1) == pyo.value(
        m.fs.sofc.max_cell_temperature[0]
    )

    assert pytest.approx(99.4, abs=1e-1) == pyo.value(
        m.fs.sofc.delta_cell_temperature[0]
    )
