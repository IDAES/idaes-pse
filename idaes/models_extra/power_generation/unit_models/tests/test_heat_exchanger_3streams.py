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
Tests for 0D Boiler heat exchanger model.

Author: Miguel Zamarripa
"""
import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, value, Param
from idaes.core import FlowsheetBlock

# import ideal flue gas prop pack
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
from idaes.models_extra.power_generation.unit_models.heat_exchanger_3streams import (
    HeatExchangerWith3Streams,
)
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.fixture(scope="module")
def build_unit():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()
    m.fs.unit = HeatExchangerWith3Streams(
        side_1_property_package=m.fs.prop_fluegas,
        side_2_property_package=m.fs.prop_fluegas,
        side_3_property_package=m.fs.prop_fluegas,
        has_heat_transfer=True,
        has_pressure_change=True,
        flow_type_side_2="counter-current",
        flow_type_side_3="counter-current",
    )
    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 30
    # Check unit config arguments
    assert len(m.fs.unit.config) == 15
    assert m.fs.unit.config.has_heat_transfer
    assert m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unit(build_unit):
    m = build_unit
    # Set inputs
    # FLUE GAS Inlet NETL baseline report
    FGrate = 5109.76  # mol/s equivalent of ~1930.08 klb/hr
    m.fs.mole_frac_air = Param(
        m.fs.prop_fluegas.component_list,
        mutable=False,
        initialize={
            "O2": 0.20784,
            "N2": 0.783994,
            "NO": 0.000001,
            "CO2": 0.000337339,
            "H2O": 0.0078267,
            "SO2": 0.000001,
        },
        doc="mole fraction of air species",
    )
    m.fs.unit.side_1_inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.unit.side_1_inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.unit.side_1_inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.unit.side_1_inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.unit.side_1_inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.unit.side_1_inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.unit.side_1_inlet.temperature[0].fix(650.335)
    m.fs.unit.side_1_inlet.pressure[0].fix(100145)

    for i in m.fs.prop_fluegas.component_list:
        m.fs.unit.side_2_inlet.flow_mol_comp[0, i].fix(
            FGrate * 0.6 * m.fs.mole_frac_air[i]
        )
    for i in m.fs.prop_fluegas.component_list:
        m.fs.unit.side_3_inlet.flow_mol_comp[0, i].fix(
            FGrate * 0.4 * m.fs.mole_frac_air[i]
        )
    m.fs.unit.side_2_inlet.temperature[0].fix(324.15)
    m.fs.unit.side_2_inlet.pressure[0].fix(100145)
    m.fs.unit.side_3_inlet.temperature[0].fix(373.15)
    m.fs.unit.side_3_inlet.pressure[0].fix(100145)

    iscale.calculate_scaling_factors(m)

    m.fs.unit.ua_side_2[0].fix(1.5e5)
    m.fs.unit.ua_side_3[0].fix(6.2e5)
    m.fs.unit.frac_heatloss.fix(0.05)
    m.fs.unit.deltaP_side_1[0].fix(-1000)
    m.fs.unit.deltaP_side_2[0].fix(-1000)
    m.fs.unit.deltaP_side_3[0].fix(-1000)
    assert degrees_of_freedom(m) == 0

    initialization_tester(m, dof=0)


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_run_unit(build_unit):
    m = build_unit
    optarg = {"tol": 1e-7, "max_iter": 40}
    solver.options = optarg
    # solve model
    results = solver.solve(m, tee=True)
    # Check for optimal solution
    assert check_optimal_termination(results)
    assert degrees_of_freedom(m) == 0
    assert pytest.approx(434.650, abs=1e-3) == value(
        m.fs.unit.side_1_outlet.temperature[0]
    )
    assert pytest.approx(522.135, abs=1e-3) == value(
        m.fs.unit.side_2_outlet.temperature[0]
    )
    assert pytest.approx(642.115, abs=1e-3) == value(
        m.fs.unit.side_3_outlet.temperature[0]
    )
    # energy balance
    assert pytest.approx(0, abs=1e-3) == value(
        m.fs.unit.side_1.properties_in[0].flow_mol
        * m.fs.unit.side_1.properties_in[0].enth_mol
        - m.fs.unit.side_1.properties_out[0].flow_mol
        * m.fs.unit.side_1.properties_out[0].enth_mol
        + m.fs.unit.heat_duty_side_1[0]
    )
    assert pytest.approx(0, abs=1e-3) == value(
        +m.fs.unit.side_2.properties_in[0].flow_mol
        * m.fs.unit.side_2.properties_in[0].enth_mol
        - m.fs.unit.side_2.properties_out[0].flow_mol
        * m.fs.unit.side_2.properties_out[0].enth_mol
        + m.fs.unit.heat_duty_side_2[0]
    )
    assert pytest.approx(0, abs=1e-3) == value(
        +m.fs.unit.side_3.properties_in[0].flow_mol
        * m.fs.unit.side_3.properties_in[0].enth_mol
        - m.fs.unit.side_3.properties_out[0].flow_mol
        * m.fs.unit.side_3.properties_out[0].enth_mol
        + m.fs.unit.heat_duty_side_3[0]
    )
    # pressure drop
    assert pytest.approx(-1000.0, abs=1e-3) == value(m.fs.unit.deltaP_side_1[0])
