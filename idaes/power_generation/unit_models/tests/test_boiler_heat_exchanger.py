##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for 0D Boiler heat exchanger model.

Author: Miguel Zamarripa
"""
import pytest

from pyomo.environ import (ConcreteModel,
                           TerminationCondition,
                           SolverStatus,
                           value,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from idaes.generic_models.properties import iapws95

# import ideal flue gas prop pack
from idaes.power_generation.properties import FlueGasParameterBlock
# Import Power Plant HX Unit Model
from idaes.power_generation.unit_models.boiler_heat_exchanger import (
        BoilerHeatExchanger, TubeArrangement, DeltaTMethod)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------


@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.unit = BoilerHeatExchanger(default={
        "side_1_property_package": m.fs.prop_steam,
        "side_2_property_package": m.fs.prop_fluegas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": True})

    # Check unit config arguments
    # There are 8 to 10 arguments since you can add a side 1 and 2 config by
    # side_1, side_2, or whatever the user named them
    assert len(m.fs.unit.config) >= 12 and len(m.fs.unit.config) <= 16

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.delta_T_method == \
        DeltaTMethod.counterCurrent


@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_boiler_hx():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.unit = BoilerHeatExchanger(default={
        "side_1_property_package": m.fs.prop_steam,
        "side_2_property_package": m.fs.prop_fluegas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": True})

    #   Set inputs
    h = value(iapws95.htpx(773.15*pyunits.K, 2.5449e7*pyunits.Pa))
    print(h)
    m.fs.unit.side_1_inlet.flow_mol[0].fix(24678.26)   # mol/s
    m.fs.unit.side_1_inlet.enth_mol[0].fix(h)           # J/mol
    m.fs.unit.side_1_inlet.pressure[0].fix(2.5449e7)    # Pascals

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 28.3876e3*0.18  # mol/s equivalent of ~1930.08 klb/hr
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*8.69/100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*14.49/100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*74.34/100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*2.47/100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.0006)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.002)
    m.fs.unit.side_2_inlet.temperature[0].fix(1102.335)
    m.fs.unit.side_2_inlet.pressure[0].fix(100145)

    # Primary Superheater
    ITM = 0.0254  # inch to meter conversion
    m.fs.unit.tube_di.fix((2.5-2*0.165)*ITM)
    m.fs.unit.tube_thickness.fix(0.165*ITM)
    m.fs.unit.pitch_x.fix(3*ITM)
    # gas path transverse width 54.78 ft / number of columns
    m.fs.unit.pitch_y.fix(54.78/108*12*ITM)
    m.fs.unit.tube_length.fix(53.13*12*ITM)
    m.fs.unit.tube_nrow.fix(20*2)
    m.fs.unit.tube_ncol.fix(108)
    m.fs.unit.nrow_inlet.fix(4)
    m.fs.unit.delta_elevation.fix(50)
    m.fs.unit.tube_r_fouling = 0.000176  # (0.001 h-ft^2-F/BTU)
    m.fs.unit.tube_r_fouling = 0.003131  # (0.03131 - 0.1779 h-ft^2-F/BTU)
    if m.fs.unit.config.has_radiation is True:
        m.fs.unit.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.unit.fcorrection_htc.fix(1.5)
    # correction factor for pressure drop calc tube side
    m.fs.unit.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.unit.fcorrection_dp_shell.fix(1.0)

    assert degrees_of_freedom(m) == 0

    m.fs.unit.initialize()

    results = solver.solve(m)
    # Check for optimal solution
    assert results.solver.termination_condition == \
        TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok
    assert value(m.fs.unit.side_1.properties_out[0].temperature) == \
        pytest.approx(588.07, 1)
    assert value(m.fs.unit.side_2.properties_out[0].temperature) == \
        pytest.approx(573.07, 1)


@pytest.mark.skipif(not iapws95.iapws95_available(),
                    reason="IAPWS not available")
@pytest.mark.integration
def test_units():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.unit = BoilerHeatExchanger(default={
        "side_1_property_package": m.fs.prop_steam,
        "side_2_property_package": m.fs.prop_fluegas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": True})

    assert_units_consistent(m)
