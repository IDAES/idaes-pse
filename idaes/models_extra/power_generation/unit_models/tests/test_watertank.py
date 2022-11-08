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
Watertank model test
The water tank has only one inlet and one outlet

main assumptions:
1.- Heat loss is a variable given by the user (zero heat loss can be
    specified if adiabatic)
2.- Calculate pressure change due to gravity based on water level
    and contraction to downcomer
3.- Water level is either fixed for steady-state model or calculated for
    dynamic model
4.- Assume enthalpy_in == enthalpy_out + heat loss
5.- Subcooled water from economizer and saturated water from waterwall
    are mixed before entering the tank

Created on Nov 04 2020 by J. Ma, M. Zamarripa, D. Caballero
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.unit_models.watertank import WaterTank

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------
@pytest.fixture(scope="module")
def build_watertank_simple():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterTank(
        property_package=m.fs.prop_water,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # fix inputs for simple tank
    m.fs.unit.tank_cross_sect_area.fix(1.14)  # tank cross sectional area
    m.fs.unit.tank_level[:].fix(0.6)  # tank level
    m.fs.unit.heat_duty[:].fix(0.0)  # assume no heat loss
    return m


@pytest.fixture(scope="module")
def build_watertank_rect():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterTank(
        tank_type="rectangular_tank",
        property_package=m.fs.prop_water,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # fix inputs for horizontal cylindrical tank
    m.fs.unit.tank_width.fix(0.5)  # tank width
    m.fs.unit.tank_length.fix(1.2)  # tank length
    m.fs.unit.tank_level[:].fix(0.6)  # tank level
    m.fs.unit.heat_duty[:].fix(0.0)  # assume no heat loss
    return m


@pytest.fixture(scope="module")
def build_watertank_vert_cylin():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterTank(
        tank_type="vertical_cylindrical_tank",
        property_package=m.fs.prop_water,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # fix inputs for horizontal cylindrical tank
    m.fs.unit.tank_diameter.fix(1.2)  # tank diameter
    m.fs.unit.tank_level[:].fix(0.6)  # tank level
    m.fs.unit.heat_duty[:].fix(0.0)  # assume no heat loss
    return m


@pytest.fixture(scope="module")
def build_watertank_hori_cylin():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.unit = WaterTank(
        tank_type="horizontal_cylindrical_tank",
        property_package=m.fs.prop_water,
        has_holdup=False,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # fix inputs for horizontal cylindrical tank
    m.fs.unit.tank_diameter.fix(1.2)  # tank diameter
    m.fs.unit.tank_length.fix(15.3256)  # tank length
    m.fs.unit.tank_level[:].fix(0.6)  # tank level
    m.fs.unit.heat_duty[:].fix(0.0)  # assume no heat loss
    return m


@pytest.fixture(scope="module")
def tank_models(
    build_watertank_simple,
    build_watertank_rect,
    build_watertank_vert_cylin,
    build_watertank_hori_cylin,
):
    return [
        build_watertank_simple,
        build_watertank_rect,
        build_watertank_vert_cylin,
        build_watertank_hori_cylin,
    ]


@pytest.mark.unit
def test_basic_build(tank_models):
    """Make a turbine model and make sure it doesn't throw exception"""
    for i in tank_models:
        # build the simple tank
        m = i
        assert degrees_of_freedom(m) == 3
        # Check unit config arguments
        assert len(m.fs.unit.config) == 10
        assert m.fs.unit.config.has_heat_transfer
        assert m.fs.unit.config.has_pressure_change
        assert m.fs.unit.config.property_package is m.fs.prop_water


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_watertank(tank_models):
    state_args_feedwater = {
        "flow_mol": 83193 / 3600 * 1000,  # mol/s
        "pressure": 12024201.99,  # Pa
        "enth_mol": 22723.907,
    }  # J/mol
    for i in tank_models:
        initialization_tester(i, dof=3, state_args=state_args_feedwater)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_run_watertank(tank_models):
    # volume_values = [simple tank, rectangular, vertical cylindrical,
    #                  horizontal cylindrical]
    volume_values = [0.684, 0.36, 0.679, 8.666]

    optarg = {"tol": 1e-7, "linear_solver": "ma27", "max_iter": 40}

    solver.options = optarg

    for i in tank_models:

        m = i

        # fix inlets
        m.fs.unit.inlet.flow_mol[:].fix()
        m.fs.unit.inlet.pressure[:].fix()
        m.fs.unit.inlet.enth_mol[:].fix()

        # solve model
        results = solver.solve(m, tee=True)

        # Check for optimal solution
        assert pyo.check_optimal_termination(results)
        assert degrees_of_freedom(m) == 0
        assert pytest.approx(0.6, abs=1e-3) == pyo.value(m.fs.unit.tank_level[0])
        # mass balance
        assert pytest.approx(0, abs=1e-3) == pyo.value(
            m.fs.unit.inlet.flow_mol[0] - m.fs.unit.outlet.flow_mol[0]
        )
        # energy balance
        assert pytest.approx(0, abs=1e-3) == pyo.value(
            m.fs.unit.inlet.flow_mol[0] * m.fs.unit.inlet.enth_mol[0]
            - m.fs.unit.outlet.flow_mol[0] * m.fs.unit.outlet.enth_mol[0]
        )
        # pressure drop
        assert pytest.approx(4410.081, rel=1e-3) == pyo.value(m.fs.unit.deltaP[0])

        # volume
        assert pytest.approx(
            volume_values[tank_models.index(i)], abs=1e-3
        ) == pyo.value(m.fs.unit.volume[0])
