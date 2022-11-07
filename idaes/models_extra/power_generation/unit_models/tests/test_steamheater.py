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
Unit operation model for a steam heater applicable to platen superheater
and roof superheater, model main equations:

* Heat is given by fire-side boiler model
* Calculate pressure change due to friction
* Calculate slag layer wall temperature
* Consider a layer of metal and a layer of slag

Created on Aug 27, 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.models.properties import iapws95

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
from idaes.models_extra.power_generation.unit_models.steamheater import SteamHeater

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_unit():
    m = pyo.ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.unit = SteamHeater(
        dynamic=False,
        property_package=m.fs.properties,
        has_holdup=True,
        has_heat_transfer=True,
        has_pressure_change=True,
        single_side_only=True,
    )
    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a steam heater model and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 12
    # Check unit config arguments
    assert len(m.fs.unit.config) == 10
    assert m.fs.unit.config.has_heat_transfer
    assert isinstance(m.fs.unit.heat_duty, pyo.Var)
    assert isinstance(m.fs.unit.deltaP, pyo.Var)
    assert m.fs.unit.config.property_package is m.fs.properties


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unit(build_unit):
    m = build_unit
    # platen superheater
    m.fs.unit.diameter_in.fix(0.035)
    m.fs.unit.tube_thickness.fix(0.005)
    m.fs.unit.fin_thickness.fix(0.0035)
    m.fs.unit.slag_thickness[:].fix(0.001)
    m.fs.unit.fin_length.fix(0.00966)
    m.fs.unit.tube_length.fix(40.0)
    m.fs.unit.number_tubes.fix(200.0)
    m.fs.unit.therm_cond_slag.fix(1.3)
    m.fs.unit.control_volume.scaling_factor_holdup_energy = 1e-8
    m.fs.unit.heat_fireside[:].fix(5.5e7)  # initial guess
    hpl = iapws95.htpx(798.15 * pyo.units.K, 24790249.01 * pyo.units.Pa)
    m.fs.unit.inlet[:].flow_mol.fix(24194.177)
    m.fs.unit.inlet[:].enth_mol.fix(hpl)
    m.fs.unit.inlet[:].pressure.fix(24790249.01)

    state_args = {"flow_mol": 24194.177, "pressure": 24790249.01, "enth_mol": hpl}

    initialization_tester(build_unit, dof=0, state_args=state_args)


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
    assert pytest.approx(-282158.4020, rel=1e-5) == pyo.value(m.fs.unit.deltaP[0])
    # check energy balance
    assert (
        pytest.approx(
            pyo.value(
                m.fs.unit.control_volume.properties_in[0].enth_mol
                - m.fs.unit.control_volume.properties_out[0].enth_mol
            ),
            abs=1e-3,
        )
        == -2273.2742
    )

    assert pytest.approx(835.7582, abs=1e-3) == pyo.value(
        m.fs.unit.control_volume.properties_out[0].temperature
    )
