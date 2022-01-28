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

from pyomo.environ import (
    ConcreteModel,
    check_optimal_termination,
    SolverStatus,
    value,
    units as pyunits,
)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock

from idaes.generic_models.properties import iapws95

# import ideal flue gas prop pack
from idaes.power_generation.properties import (
    FlueGasParameterBlock,
    FlueGasStateBlock,
)

# Import Power Plant HX Unit Model
from idaes.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger,
    TubeArrangement,
    DeltaTMethod,
    delta_temperature_lmtd_callback,
    delta_temperature_amtd_callback,
    delta_temperature_underwood_callback,
    delta_temperature_underwood_tune_callback,
    HeatExchangerFlowPattern,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util import get_solver
import idaes.core.util.scaling as iscale

from idaes.core.util.testing import PhysicalParameterTestBlock
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()

# -----------------------------------------------------------------------------

@pytest.mark.unit
def test_no_deprecated(caplog):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    caplog.clear()
    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd_callback,
            "tube": {"property_package": m.fs.prop_steam},
            "shell": {"property_package": m.fs.prop_fluegas},
            "has_pressure_change": True,
            "has_holdup": True,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )
    n_warn = 0
    n_depreacted = 0
    for record in caplog.records:
        if record.levelno == idaeslog.WARNING:
            n_warn += 1
    assert n_warn == 0 # 1 DeltaTMethod Enum and 1 for delta_T_method option

@pytest.mark.unit
def test_deprecated_delta_T_method(caplog):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    caplog.clear()
    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd_callback,
            "tube": {"property_package": m.fs.prop_steam},
            "shell": {"property_package": m.fs.prop_fluegas},
            "has_pressure_change": True,
            "has_holdup": True,
            "delta_T_method": HeatExchangerFlowPattern.countercurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )
    n_warn = 0
    n_depreacted = 0
    for record in caplog.records:
        if record.levelno == idaeslog.WARNING:
            n_warn += 1
        if "deprecated" in record.msg:
            n_depreacted += 1
    assert n_warn == 1
    assert n_depreacted == 1
    assert m.fs.unit.config.flow_pattern == HeatExchangerFlowPattern.countercurrent

@pytest.mark.unit
def test_deprecated_delta_T_method_enum1(caplog):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    caplog.clear()
    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd_callback,
            "tube": {"property_package": m.fs.prop_steam},
            "shell": {"property_package": m.fs.prop_fluegas},
            "has_pressure_change": True,
            "has_holdup": True,
            "delta_T_method": DeltaTMethod.counterCurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )
    n_warn = 0
    n_depreacted = 0
    for record in caplog.records:
        if record.levelno == idaeslog.WARNING:
            n_warn += 1
        if "deprecated" in record.msg:
            n_depreacted += 1
    assert n_warn == 3
    assert n_depreacted == 3
    assert m.fs.unit.config.flow_pattern == HeatExchangerFlowPattern.countercurrent

@pytest.mark.unit
def test_deprecated_delta_T_method_enum1(caplog):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    caplog.clear()
    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd_callback,
            "tube": {"property_package": m.fs.prop_steam},
            "shell": {"property_package": m.fs.prop_fluegas},
            "has_pressure_change": True,
            "has_holdup": True,
            "delta_T_method": DeltaTMethod.coCurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )
    n_warn = 0
    n_depreacted = 0
    for record in caplog.records:
        if record.levelno == idaeslog.WARNING:
            n_warn += 1
        if "deprecated" in record.msg:
            n_depreacted += 1
    assert n_warn == 3
    assert n_depreacted == 3
    assert m.fs.unit.config.flow_pattern == HeatExchangerFlowPattern.cocurrent

@pytest.mark.unit
def test_deprecated_prop_pack(caplog):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    caplog.clear()
    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd_callback,
            "side_1_property_package": m.fs.prop_steam,
            "side_2_property_package": m.fs.prop_fluegas,
            "has_pressure_change": True,
            "has_holdup": True,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )
    n_warn = 0
    n_depreacted = 0
    for record in caplog.records:
        if record.levelno == idaeslog.WARNING:
            n_warn += 1
        if "deprecated" in record.msg:
            n_depreacted += 1
    assert n_warn == 2
    assert n_depreacted == 2
    assert isinstance(
        m.fs.unit.config.hot_side_config.property_package,
        FlueGasParameterBlock
    )
    assert isinstance(
        m.fs.unit.config.cold_side_config.property_package,
        iapws95.Iapws95ParameterBlock
    )


def tc(delta_temperature_callback=delta_temperature_underwood_tune_callback):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_callback,
            "tube": {"property_package": m.fs.prop_steam},
            "shell": {"property_package": m.fs.prop_fluegas},
            "has_pressure_change": True,
            "has_holdup": True,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )

    # Check unit config arguments
    assert not m.fs.unit.config.dynamic
    assert m.fs.unit.config.has_holdup
    assert m.fs.unit.config.flow_pattern == HeatExchangerFlowPattern.countercurrent


def th(
    delta_temperature_callback=delta_temperature_underwood_tune_callback,
    tout_1=809.55,
    tout_2=788.53,
):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_callback,
            "tube": {"property_package": m.fs.prop_steam},
            "shell": {"property_package": m.fs.prop_fluegas},
            "has_pressure_change": True,
            "has_holdup": False,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )

    #   Set inputs
    h = value(iapws95.htpx(773.15 * pyunits.K, 2.5449e7 * pyunits.Pa))
    m.fs.unit.side_1_inlet.flow_mol[0].fix(24678.26)  # mol/s
    m.fs.unit.side_1_inlet.enth_mol[0].fix(h)  # J/mol
    m.fs.unit.side_1_inlet.pressure[0].fix(2.5449e7)  # Pascals

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 28.3876e3 * 0.18  # mol/s equivalent of ~1930.08 klb/hr
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.unit.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.unit.side_2_inlet.temperature[0].fix(1102.335)
    m.fs.unit.side_2_inlet.pressure[0].fix(100145)

    # Primary Superheater
    ITM = 0.0254  # inch to meter conversion
    m.fs.unit.tube_di.fix((2.5 - 2 * 0.165) * ITM)
    m.fs.unit.tube_thickness.fix(0.165 * ITM)
    m.fs.unit.pitch_x.fix(3 * ITM)
    # gas path transverse width 54.78 ft / number of columns
    m.fs.unit.pitch_y.fix(54.78 / 108 * 12 * ITM)
    m.fs.unit.tube_length.fix(53.13 * 12 * ITM)
    m.fs.unit.tube_nrow.fix(20 * 2)
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
    iscale.calculate_scaling_factors(m)
    m.fs.unit.initialize()

    results = solver.solve(m)
    # Check for optimal solution
    assert check_optimal_termination(results)
    assert value(m.fs.unit.side_1.properties_out[0].temperature) == pytest.approx(
        tout_1, abs=0.5
    )
    assert value(m.fs.unit.side_2.properties_out[0].temperature) == pytest.approx(
        tout_2, abs=0.5
    )



def tu(delta_temperature_callback=delta_temperature_underwood_tune_callback):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.prop_steam = iapws95.Iapws95ParameterBlock()
    m.fs.prop_fluegas = FlueGasParameterBlock()

    m.fs.unit = BoilerHeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_callback,
            "tube": {"property_package": m.fs.prop_steam},
            "shell": {"property_package": m.fs.prop_fluegas},
            "has_pressure_change": True,
            "has_holdup": False,
            "flow_pattern": HeatExchangerFlowPattern.countercurrent,
            "tube_arrangement": TubeArrangement.inLine,
            "side_1_water_phase": "Liq",
            "has_radiation": True,
        }
    )

    assert_units_consistent(m)


@pytest.mark.unit
def test_config_am():
    tc(delta_temperature_amtd_callback)


@pytest.mark.unit
def test_config_lm():
    tc(delta_temperature_lmtd_callback)


@pytest.mark.unit
def test_config_uw():
    tc(delta_temperature_underwood_callback)


@pytest.mark.unit
def test_config_uwt():
    tc(delta_temperature_underwood_tune_callback)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_boiler_hx_am():
    # arithmetic mean is pretty far off
    th(delta_temperature_amtd_callback, tout_1=817.7, tout_2=720)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_boiler_hx_lm():
    th(delta_temperature_lmtd_callback)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_boiler_hx_uw():
    th(delta_temperature_underwood_callback)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_boiler_hx_uwt():
    th(delta_temperature_underwood_tune_callback)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.integration
def test_units_am():
    tu(delta_temperature_amtd_callback)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.integration
def test_units_lm():
    tu(delta_temperature_lmtd_callback)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.integration
def test_units_uw():
    tu(delta_temperature_underwood_callback)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.integration
def test_units_uwt():
    tu(delta_temperature_underwood_tune_callback)
