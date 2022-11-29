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
Unit operation model for a waterpipe with pressure drop:

Created on Aug 27, 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models_extra.power_generation.unit_models.waterpipe import WaterPipe

# Import Unit Model Modules
from idaes.models.properties import iapws95

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_unit():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterPipe(
        dynamic=False,
        property_package=m.fs.properties,
        has_holdup=True,
        has_heat_transfer=False,
        has_pressure_change=True,
        water_phase="Liq",
        contraction_expansion_at_end="contraction",
    )

    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a unit model and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 9
    # Check unit config arguments
    assert len(m.fs.unit.config) == 11
    assert m.fs.unit.config.has_pressure_change
    assert isinstance(m.fs.unit.deltaP, pyo.Var)
    assert m.fs.unit.config.property_package is m.fs.properties


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unit(build_unit):
    m = build_unit
    m.fs.unit.diameter.fix(0.04)
    m.fs.unit.length.fix(40)
    m.fs.unit.number_of_pipes.fix(100)
    m.fs.unit.elevation_change.fix(25)
    m.fs.unit.area_ratio.fix(0.5)
    m.fs.unit.fcorrection_dp.fix(1.0)

    state_args = {"flow_mol": 10000, "pressure": 1.3e7, "enth_mol": 18000}

    initialization_tester(build_unit, dof=3, state_args=state_args)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_solve_unit(build_unit):
    m = build_unit

    m.fs.unit.inlet.enth_mol.fix()
    m.fs.unit.inlet.flow_mol.fix()
    m.fs.unit.inlet.pressure.fix()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)

    # check material balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].flow_mol
                - m.fs.unit.control_volume.properties_out[0].flow_mol
            ),
            abs=1e-3,
        )
        == 0
    )

    # pressure drop
    assert pytest.approx(-224074.4039, rel=1e-5) == pyo.value(m.fs.unit.deltaP[0])
    # check energy balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].enth_mol
                - m.fs.unit.control_volume.properties_out[0].enth_mol
            ),
            abs=1e-3,
        )
        == 0
    )


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_pipe_expansion():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterPipe(
        dynamic=False,
        property_package=m.fs.properties,
        has_holdup=True,
        has_heat_transfer=False,
        has_pressure_change=True,
        water_phase="Liq",
        contraction_expansion_at_end="expansion",
    )
    m.fs.unit.diameter.fix(0.04)
    m.fs.unit.length.fix(40)
    m.fs.unit.number_of_pipes.fix(100)
    m.fs.unit.elevation_change.fix(25)
    m.fs.unit.area_ratio.fix(0.5)
    m.fs.unit.fcorrection_dp.fix(1.0)

    state_args = {"flow_mol": 10000, "pressure": 1.3e7, "enth_mol": 18000}

    initialization_tester(m, dof=3, state_args=state_args)
    m.fs.unit.inlet.enth_mol.fix()
    m.fs.unit.inlet.flow_mol.fix()
    m.fs.unit.inlet.pressure.fix()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)

    # check material balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].flow_mol
                - m.fs.unit.control_volume.properties_out[0].flow_mol
            ),
            abs=1e-3,
        )
        == 0
    )

    # pressure drop
    assert pytest.approx(-224306.548, rel=1e-5) == pyo.value(m.fs.unit.deltaP[0])
    # check energy balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].enth_mol
                - m.fs.unit.control_volume.properties_out[0].enth_mol
            ),
            abs=1e-3,
        )
        == 0
    )


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_pipe_noexpansion():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterPipe(
        dynamic=False,
        property_package=m.fs.properties,
        has_holdup=True,
        has_heat_transfer=False,
        has_pressure_change=True,
        water_phase="Liq",
        contraction_expansion_at_end="None",
    )
    m.fs.unit.diameter.fix(0.04)
    m.fs.unit.length.fix(40)
    m.fs.unit.number_of_pipes.fix(100)
    m.fs.unit.elevation_change.fix(25)
    m.fs.unit.fcorrection_dp.fix(1.0)

    state_args = {"flow_mol": 10000, "pressure": 1.3e7, "enth_mol": 18000}

    initialization_tester(m, dof=3, state_args=state_args)
    m.fs.unit.inlet.enth_mol.fix()
    m.fs.unit.inlet.flow_mol.fix()
    m.fs.unit.inlet.pressure.fix()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)

    # check material balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].flow_mol
                - m.fs.unit.control_volume.properties_out[0].flow_mol
            ),
            abs=1e-3,
        )
        == 0
    )

    # pressure drop
    assert pytest.approx(-219382.990, rel=1e-3) == pyo.value(m.fs.unit.deltaP[0])
    # check energy balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].enth_mol
                - m.fs.unit.control_volume.properties_out[0].enth_mol
            ),
            abs=1e-3,
        )
        == 0
    )


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_pipe_vaporphase():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterPipe(
        dynamic=False,
        property_package=m.fs.properties,
        has_holdup=True,
        has_heat_transfer=False,
        has_pressure_change=True,
        water_phase="Vap",
        contraction_expansion_at_end="None",
    )
    m.fs.unit.diameter.fix(0.04)
    m.fs.unit.length.fix(40)
    m.fs.unit.number_of_pipes.fix(100)
    m.fs.unit.elevation_change.fix(25)
    m.fs.unit.fcorrection_dp.fix(1.0)

    state_args = {"flow_mol": 10000, "pressure": 3.6e6, "enth_mol": 53942}

    initialization_tester(m, dof=3, state_args=state_args)
    m.fs.unit.inlet.enth_mol.fix()
    m.fs.unit.inlet.flow_mol.fix()
    m.fs.unit.inlet.pressure.fix()
    assert degrees_of_freedom(m) == 0

    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)

    assert pyo.value(m.fs.unit.control_volume.properties_in[0].vapor_frac) == 1

    # check material balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].flow_mol
                - m.fs.unit.control_volume.properties_out[0].flow_mol
            ),
            abs=1e-3,
        )
        == 0
    )

    # pressure drop
    assert pytest.approx(-539105.588, rel=1e-3) == pyo.value(m.fs.unit.deltaP[0])
    # check energy balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].enth_mol
                - m.fs.unit.control_volume.properties_out[0].enth_mol
            ),
            abs=1e-3,
        )
        == 0
    )
