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
The boiler 2D heat exchanger model consist of a cross flow shell and tube
heat exchanger. 1-D Cross Flow Heat Exchanger Model with wall temperatures,
discretization based on tube rows


The model includes shell and tube rigorous heat transfer calculations and
pressure drop calculations for shell side. Note that this model assumes no
phase transitions (if user requires phase transitions, they need a general
model)
Created on Nov 25th, 2020 by Boiler Team (J. Ma, M. Zamarripa)
"""
import pytest

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.environ import units as pyunits

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom

# Import Unit Model Modules
from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger_2D import (
    HeatExchangerCrossFlow2D_Header,
)
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_unit():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()
    # Add a 2D cross-flow heat exchanger with header model to the flowsheet
    m.fs.unit = HeatExchangerCrossFlow2D_Header(
        tube_side={"property_package": m.fs.prop_water, "has_pressure_change": True},
        shell_side={"property_package": m.fs.prop_fluegas, "has_pressure_change": True},
        finite_elements=5,
        flow_type="counter_current",
        tube_arrangement="in-line",
        tube_side_water_phase="Vap",
        has_radiation=True,
        radial_elements=5,
        tube_inner_diameter=0.035,
        tube_thickness=0.0035,
        has_header=True,
        header_radial_elements=5,
        header_inner_diameter=0.3,
        header_wall_thickness=0.03,
    )

    # Primary Superheater (NETL baseline report)
    ITM = 0.0254  # inch to meter conversion
    # m.fs.unit.tube_di.fix((2.5-2*0.165)*ITM)
    # m.fs.unit.tube_thickness.fix(0.165*ITM)
    m.fs.unit.pitch_x.fix(3 * ITM)
    # gas path transverse width 54.78 ft / number of columns
    m.fs.unit.pitch_y.fix(54.78 / 108 * 12 * ITM)
    m.fs.unit.tube_length_seg.fix(53.13 * 12 * ITM)
    m.fs.unit.tube_nseg.fix(20 * 2)
    m.fs.unit.tube_ncol.fix(108)
    m.fs.unit.tube_inlet_nrow.fix(4)
    m.fs.unit.delta_elevation.fix(50)
    m.fs.unit.tube_r_fouling = 0.000176  # (0.001 h-ft^2-F/BTU)
    m.fs.unit.tube_r_fouling = 0.003131  # (0.03131 - 0.1779 h-ft^2-F/BTU)

    # material inputs
    m.fs.unit.therm_cond_wall = 43.0  # Carbon steel 1% C in W/m/K
    m.fs.unit.dens_wall = 7850  # kg/m3 or 0.284 lb/in3
    m.fs.unit.cp_wall = 510.8  # J/kg-K (0.05 to 0.25 % C)
    m.fs.unit.Young_modulus = 2.00e5  # 200 GPa (29,000 ksi)
    m.fs.unit.Possion_ratio = 0.29
    m.fs.unit.coefficient_therm_expansion = 1.3e-5

    m.fs.unit.emissivity_wall.fix(0.7)
    m.fs.unit.fcorrection_htc_tube.fix(1.0)
    m.fs.unit.fcorrection_htc_shell.fix(1.0)
    m.fs.unit.fcorrection_dp_tube.fix(1.0)
    m.fs.unit.fcorrection_dp_shell.fix(1.0)
    m.fs.unit.temperature_ambient.fix(350.0)
    m.fs.unit.head_insulation_thickness.fix(0.025)
    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a model and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 11
    # Check unit config arguments
    assert len(m.fs.unit.config) == 19
    assert m.fs.unit.config.shell_side.has_pressure_change
    assert m.fs.unit.config.tube_side.has_pressure_change
    assert m.fs.unit.tube.config.property_package is m.fs.prop_water


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unit(build_unit):
    m = build_unit
    #   Set inputs
    h = pyo.value(iapws95.htpx(773.15 * pyunits.K, 2.5449e7 * pyunits.Pa))
    print(h)
    m.fs.unit.tube_inlet.flow_mol[0].fix(24678.26)  # mol/s
    m.fs.unit.tube_inlet.enth_mol[0].fix(h)  # J/mol
    m.fs.unit.tube_inlet.pressure[0].fix(2.5449e7)  # Pascals

    # FLUE GAS Inlet relative to a primary superheater
    FGrate = 5109.76  # mol/s equivalent of ~1930.08 klb/hr
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.unit.shell_inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.unit.shell_inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.unit.shell_inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.unit.shell_inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.unit.shell_inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.unit.shell_inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.unit.shell_inlet.temperature[0].fix(1102.335)
    m.fs.unit.shell_inlet.pressure[0].fix(100145)

    initialization_tester(build_unit, dof=0)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_run_unit(build_unit):
    m = build_unit
    assert degrees_of_freedom(m) == 0
    optarg = {"tol": 1e-6, "linear_solver": "ma27", "max_iter": 40}
    solver.options = optarg
    # solve model
    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert pyo.check_optimal_termination(results)

    # energy balance
    assert pytest.approx(0, abs=1e-3) == pyo.value(
        m.fs.unit.tube_inlet.flow_mol[0] * m.fs.unit.tube_inlet.enth_mol[0]
        + sum(
            m.fs.unit.shell_inlet.flow_mol_comp[0, i]
            for i in m.fs.unit.shell.config.property_package.component_list
        )
        * m.fs.unit.shell.properties[0, 0].enth_mol
        - m.fs.unit.tube_outlet.flow_mol[0] * m.fs.unit.tube_outlet.enth_mol[0]
        - sum(
            m.fs.unit.shell_outlet.flow_mol_comp[0, i]
            for i in m.fs.unit.shell.config.property_package.component_list
        )
        * m.fs.unit.shell.properties[0, 1].enth_mol
    )
    # pressure drop
    assert pytest.approx(100134.3247, abs=1e-3) == pyo.value(
        m.fs.unit.shell.properties[0, 1].pressure
    )
    # mass balance
    assert pytest.approx(0, abs=1e-3) == pyo.value(
        sum(
            m.fs.unit.shell_inlet.flow_mol_comp[0, i]
            for i in m.fs.unit.shell.config.property_package.component_list
        )
        + m.fs.unit.tube_inlet.flow_mol[0]
        - m.fs.unit.tube_outlet.flow_mol[0]
        - sum(
            m.fs.unit.shell_outlet.flow_mol_comp[0, i]
            for i in m.fs.unit.shell.config.property_package.component_list
        )
    )
