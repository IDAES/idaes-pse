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
Drum model
The drum model consists of three main unit operations
1) a flash model (only for water and steam)
2) a mixer model
3) a tank drum level/accumulation model

Inlet Ports:

* water/steam mixture from water wall
* feedwater_inlet: feedwater from economizer/pipe

Outlet Ports:

* liquid_outlet: liquid to downcomer
* steam_outlet: steam leaving the drum

The drum model receives saturated water from the water wall, and this is
separated at the flash unit, the steam leaves the drum, while the water/liquid
mixes with the feedwater from the economizer (or water pipe model).
Finally, the mixed state liquid (stream leaving the mixer), is used as the
control volume of the drum level model, which computes velocity and pressure
drop. The exit of the drum model is the liquid outlet.


Created on Thu Aug 18 12:59:50 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.util.check_units import assert_units_consistent

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.unit_models.drum import Drum

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_drum():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.unit = Drum(
        property_package=m.fs.prop_water,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # fix inputs
    m.fs.unit.drum_diameter.fix(1.2)
    m.fs.unit.drum_length.fix(15.3256)
    m.fs.unit.number_downcomers.fix(6)
    m.fs.unit.downcomer_diameter.fix(0.38)
    m.fs.unit.drum_level[:].fix(0.6)  # drum level, set as unit radius
    m.fs.unit.heat_duty[:].fix(0.0)  # assume no heat loss
    return m


@pytest.mark.unit
def test_basic_build(build_drum):
    """Make a model and make sure it doesn't throw exception"""
    m = build_drum
    assert degrees_of_freedom(m) == 5
    # Check unit config arguments
    assert len(m.fs.unit.config) == 9
    assert m.fs.unit.config.has_heat_transfer
    assert m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.prop_water


@pytest.mark.integration
def test_units(build_drum):
    assert_units_consistent(build_drum)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_drum(build_drum):
    state_args_water_steam = {
        "flow_mol": 103745 / 3600 * 1000,  # mol/s
        "pressure": 12024201.99,  # Pa
        "enth_mol": 28365.2608,
    }  # j/mol

    state_args_feedwater = {
        "flow_mol": 83193 / 3600 * 1000,
        "pressure": 12024201.99,
        "enth_mol": 22723.907,
    }

    initialization_tester(
        build_drum,
        dof=5,
        state_args_water_steam=state_args_water_steam,
        state_args_feedwater=state_args_feedwater,
    )


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_run_drum(build_drum):
    m = build_drum
    # fix inlets
    m.fs.unit.water_steam_inlet.flow_mol[:].fix()
    m.fs.unit.water_steam_inlet.pressure[:].fix()
    m.fs.unit.water_steam_inlet.enth_mol[:].fix()
    m.fs.unit.feedwater_inlet.flow_mol[:].fix()
    # m.fs.unit.feedwater_inlet.pressure[:].fix()
    m.fs.unit.feedwater_inlet.enth_mol[:].fix()

    optarg = {"tol": 1e-7, "linear_solver": "ma27", "max_iter": 40}
    solver.options = optarg
    # solve model
    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)
    assert degrees_of_freedom(m) == 0
    assert pytest.approx(0.6, abs=1e-3) == pyo.value(m.fs.unit.drum_level[0])
    # energy balance
    assert pytest.approx(0, abs=1e-3) == pyo.value(
        m.fs.unit.water_steam_inlet.flow_mol[0]
        * m.fs.unit.water_steam_inlet.enth_mol[0]
        + m.fs.unit.feedwater_inlet.flow_mol[0] * +m.fs.unit.feedwater_inlet.enth_mol[0]
        - m.fs.unit.steam_outlet.flow_mol[0] * m.fs.unit.steam_outlet.enth_mol[0]
        - m.fs.unit.liquid_outlet.flow_mol[0] * m.fs.unit.liquid_outlet.enth_mol[0]
    )
    # pressure drop
    assert pytest.approx(2261.2171, rel=1e-4) == pyo.value(m.fs.unit.deltaP[0])
    # mass balance
    assert pytest.approx(0, abs=1e-3) == pyo.value(
        m.fs.unit.water_steam_inlet.flow_mol[0]
        + m.fs.unit.feedwater_inlet.flow_mol[0]
        - m.fs.unit.steam_outlet.flow_mol[0]
        - m.fs.unit.liquid_outlet.flow_mol[0]
    )

    assert pytest.approx(
        pyo.value(
            m.fs.unit.water_steam_inlet.flow_mol[0]
            * m.fs.unit.flash.mixed_state[0].vapor_frac
        ),
        abs=1e-3,
    ) == pyo.value(m.fs.unit.steam_outlet.flow_mol[0])
