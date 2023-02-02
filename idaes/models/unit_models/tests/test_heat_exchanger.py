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
Tests for 0D heat exchanger models.

Author: John Eslick
"""
import pytest
import pandas

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    value,
    Var,
    units as pyunits,
)
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from idaes.models.unit_models.heat_exchanger import (
    delta_temperature_lmtd_callback,
    delta_temperature_lmtd2_callback,
    delta_temperature_amtd_callback,
    delta_temperature_underwood_callback,
    delta_temperature_lmtd_smooth_callback,
    HeatExchanger,
    HeatExchangerFlowPattern,
)
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties import iapws95
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.examples.BT_PR import configuration
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError


# Imports to assemble BT-PR with different units
from idaes.core import LiquidPhase, VaporPhase, Component
from idaes.models.properties.modular_properties.state_definitions import FTPx
from idaes.models.properties.modular_properties.eos.ceos import Cubic, CubicType
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    LogBubbleDew,
)
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity
import idaes.models.properties.modular_properties.pure.RPP4 as RPP
from idaes.core.util.exceptions import ConfigurationError

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_bad_option():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(KeyError):
        m.fs.unit = HeatExchanger(**{"I'm a bad option": "hot"})


@pytest.mark.unit
def test_bad_option2():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="cold_side_name cannot be 'hot_side'."
    ):
        m.fs.unit = HeatExchanger(cold_side_name="hot_side")


@pytest.mark.unit
def test_bad_option3():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="hot_side_name cannot be 'cold_side'."
    ):
        m.fs.unit = HeatExchanger(hot_side_name="cold_side")


@pytest.mark.unit
def test_bad_option4():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        ConfigurationError, match="cold_side_name cannot be 'cold_side'."
    ):
        m.fs.unit = HeatExchanger(cold_side_name="cold_side")


@pytest.mark.unit
def test_bad_option5():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(ConfigurationError, match="hot_side_name cannot be 'hot_side'."):
        m.fs.unit = HeatExchanger(hot_side_name="hot_side")


@pytest.mark.unit
def test_hot_and_cold_names_same():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    with pytest.raises(
        NameError,
        match="HeatExchanger hot and cold side cannot have the same name 'shell'.",
    ):
        m.fs.unit = HeatExchanger(hot_side_name="shell", cold_side_name="shell")


@pytest.mark.unit
def test_hot_side_name_clash():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    with pytest.raises(
        ValueError,
        match="fs.unit could not assign hot side alias "
        "build as an attribute of that name already "
        "exists.",
    ):
        m.fs.unit = HeatExchanger(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            hot_side_name="build",
        )


@pytest.mark.unit
def test_cold_side_name_clash():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    with pytest.raises(
        ValueError,
        match="fs.unit could not assign cold side alias "
        "build as an attribute of that name already "
        "exists.",
    ):
        m.fs.unit = HeatExchanger(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            cold_side_name="build",
        )


@pytest.mark.unit
def test_user_names():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = HeatExchanger(
        hot_side_name="shell",
        cold_side_name="tube",
        shell={"property_package": m.fs.properties},
        tube={"property_package": m.fs.properties},
    )

    assert m.fs.unit.config.hot_side.property_package is m.fs.properties
    assert m.fs.unit.config.cold_side.property_package is m.fs.properties

    assert m.fs.unit.shell is m.fs.unit.hot_side
    assert m.fs.unit.tube is m.fs.unit.cold_side

    assert m.fs.unit.shell_inlet is m.fs.unit.hot_side_inlet
    assert m.fs.unit.tube_inlet is m.fs.unit.cold_side_inlet
    assert m.fs.unit.shell_outlet is m.fs.unit.hot_side_outlet
    assert m.fs.unit.tube_outlet is m.fs.unit.cold_side_outlet


@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = HeatExchanger(
        hot_side={"property_package": m.fs.properties},
        cold_side={"property_package": m.fs.properties},
    )

    # Check unit config arguments
    # There are 8 to 10 arguments since you can add a side 1 and 2 config by
    # side_1, side_2, or whatever the user named them
    assert len(m.fs.unit.config) == 8

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert isinstance(m.fs.unit.config.hot_side, ConfigBlock)
    assert isinstance(m.fs.unit.config.cold_side, ConfigBlock)
    assert (
        m.fs.unit.config.delta_temperature_callback is delta_temperature_lmtd_callback
    )
    assert m.fs.unit.config.flow_pattern == HeatExchangerFlowPattern.countercurrent

    # Check hot_side config
    assert len(m.fs.unit.config.hot_side) == 7
    assert (
        m.fs.unit.config.hot_side.material_balance_type
        == MaterialBalanceType.useDefault
    )
    assert m.fs.unit.config.hot_side.energy_balance_type == EnergyBalanceType.useDefault
    assert (
        m.fs.unit.config.hot_side.momentum_balance_type
        == MomentumBalanceType.pressureTotal
    )
    assert not m.fs.unit.config.hot_side.has_phase_equilibrium
    assert not m.fs.unit.config.hot_side.has_pressure_change
    assert m.fs.unit.config.hot_side.property_package is m.fs.properties

    # Check cold_side config
    assert len(m.fs.unit.config.cold_side) == 7
    assert (
        m.fs.unit.config.cold_side.material_balance_type
        == MaterialBalanceType.useDefault
    )
    assert (
        m.fs.unit.config.cold_side.energy_balance_type == EnergyBalanceType.useDefault
    )
    assert (
        m.fs.unit.config.cold_side.momentum_balance_type
        == MomentumBalanceType.pressureTotal
    )
    assert not m.fs.unit.config.cold_side.has_phase_equilibrium
    assert not m.fs.unit.config.cold_side.has_pressure_change
    assert m.fs.unit.config.cold_side.property_package is m.fs.properties


def basic_model(cb=delta_temperature_lmtd_callback):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = iapws95.Iapws95ParameterBlock()

    m.fs.unit = HeatExchanger(
        hot_side={"property_package": m.fs.properties},
        cold_side={"property_package": m.fs.properties},
        delta_temperature_callback=cb,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
    )
    #   Set inputs
    m.fs.unit.hot_side_inlet.flow_mol[0].fix(100)
    m.fs.unit.hot_side_inlet.enth_mol[0].fix(4000)
    m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

    m.fs.unit.cold_side_inlet.flow_mol[0].fix(100)
    m.fs.unit.cold_side_inlet.enth_mol[0].fix(3500)
    m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

    m.fs.unit.area.fix(1000)
    m.fs.unit.overall_heat_transfer_coefficient.fix(100)

    assert degrees_of_freedom(m) == 0
    m.fs.unit.initialize()
    return m


def basic_model2(cb=delta_temperature_lmtd_callback):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = iapws95.Iapws95ParameterBlock()

    m.fs.unit = HeatExchanger(
        hot_side={"property_package": m.fs.properties},
        cold_side={"property_package": m.fs.properties},
        delta_temperature_callback=cb,
        flow_pattern=HeatExchangerFlowPattern.cocurrent,
    )
    #   Set inputs
    m.fs.unit.hot_side_inlet.flow_mol[0].fix(100)
    m.fs.unit.hot_side_inlet.enth_mol[0].fix(4000)
    m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

    m.fs.unit.cold_side_inlet.flow_mol[0].fix(100)
    m.fs.unit.cold_side_inlet.enth_mol[0].fix(3500)
    m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

    m.fs.unit.area.fix(100)
    m.fs.unit.overall_heat_transfer_coefficient.fix(100)

    assert degrees_of_freedom(m) == 0
    m.fs.unit.initialize()
    return m


def basic_model3(cb=delta_temperature_lmtd_callback):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = iapws95.Iapws95ParameterBlock()

    m.fs.unit = HeatExchanger(
        hot_side={"property_package": m.fs.properties},
        cold_side={"property_package": m.fs.properties},
        delta_temperature_callback=cb,
        flow_pattern=HeatExchangerFlowPattern.crossflow,
    )
    #   Set inputs
    m.fs.unit.hot_side_inlet.flow_mol[0].fix(100)
    m.fs.unit.hot_side_inlet.enth_mol[0].fix(4000)
    m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

    m.fs.unit.cold_side_inlet.flow_mol[0].fix(100)
    m.fs.unit.cold_side_inlet.enth_mol[0].fix(3500)
    m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

    m.fs.unit.area.fix(1000)
    m.fs.unit.overall_heat_transfer_coefficient.fix(100)
    m.fs.unit.crossflow_factor.fix(1.0)
    assert degrees_of_freedom(m) == 0
    m.fs.unit.initialize()
    return m


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_lmtd_smooth_cb():
    m = basic_model(delta_temperature_lmtd_smooth_callback)
    results = solver.solve(m)
    assert check_optimal_termination(results)
    # hot in end
    assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(0.464879, rel=1e-3)
    # hot out end
    assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(
        0.465069, rel=1e-3
    )
    assert value(m.fs.unit.heat_duty[0]) == pytest.approx(46497.44)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_lmtd2_cb():
    m = basic_model(delta_temperature_lmtd2_callback)
    results = solver.solve(m)
    assert check_optimal_termination(results)
    # hot in end
    assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(0.464879, rel=1e-3)
    # hot out end
    assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(
        0.465069, rel=1e-3
    )
    assert value(m.fs.unit.heat_duty[0]) == pytest.approx(46497.44)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_lmtd2_co_cb():
    m = basic_model2(delta_temperature_lmtd2_callback)
    results = solver.solve(m)
    assert check_optimal_termination(results)
    # hot in end
    assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(6.63771, rel=1e-3)
    # hot out end
    assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(0.46658, rel=1e-3)
    assert value(m.fs.unit.heat_duty[0]) == pytest.approx(23242.69)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_lmtd2_cross_cb():
    m = basic_model3(delta_temperature_lmtd2_callback)
    results = solver.solve(m)
    assert check_optimal_termination(results)
    # hot in end
    assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(0.464879, rel=1e-3)
    # hot out end
    assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(
        0.465069, rel=1e-3
    )
    assert value(m.fs.unit.heat_duty[0]) == pytest.approx(46497.44)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_lmtd3_cb():
    m = basic_model(delta_temperature_lmtd2_callback)
    results = solver.solve(m)
    assert check_optimal_termination(results)
    # hot in end
    assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(0.464879, rel=1e-3)
    # hot out end
    assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(
        0.465069, rel=1e-3
    )
    assert value(m.fs.unit.heat_duty[0]) == pytest.approx(46497.44)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_amtd_cb():
    # since the delta T at both ends, AMTD ends up about the same as LMTD
    m = basic_model(delta_temperature_amtd_callback)
    results = solver.solve(m)
    assert check_optimal_termination(results)
    # hot in end
    assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(0.464879, rel=1e-3)
    # hot out end
    assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(
        0.465069, rel=1e-3
    )
    assert value(m.fs.unit.heat_duty[0]) == pytest.approx(46497.44)


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_underwood_cb():
    m = basic_model(delta_temperature_underwood_callback)
    results = solver.solve(m)
    assert check_optimal_termination(results)
    # hot in end
    assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(0.464879, rel=1e-3)
    # hot out end
    assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(
        0.465069, rel=1e-3
    )
    assert value(m.fs.unit.heat_duty[0]) == pytest.approx(46497.44)


# -----------------------------------------------------------------------------
class TestBTX_cocurrent(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(valid_phase="Liq")

        m.fs.unit = HeatExchanger(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_pattern=HeatExchangerFlowPattern.cocurrent,
        )

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.cold_side_inlet.temperature[0].fix(300)  # K
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.area.fix(1)
        m.fs.unit.overall_heat_transfer_coefficient.fix(100)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "hot_side_inlet")
        assert len(btx.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_inlet")
        assert len(btx.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "hot_side_outlet")
        assert len(btx.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_outlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_outlet")
        assert len(btx.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_outlet, "pressure")

        assert isinstance(btx.fs.unit.overall_heat_transfer_coefficient, Var)
        assert isinstance(btx.fs.unit.area, Var)
        assert not hasattr(btx.fs.unit, "crossflow_factor")
        assert isinstance(btx.fs.unit.heat_duty, Var)
        assert isinstance(btx.fs.unit.delta_temperature_in, Var)
        assert isinstance(btx.fs.unit.delta_temperature_out, Var)
        assert isinstance(btx.fs.unit.unit_heat_balance, Constraint)
        assert isinstance(btx.fs.unit.delta_temperature, (Var, Expression))
        assert isinstance(btx.fs.unit.heat_transfer_equation, Constraint)

        assert number_variables(btx) == 50
        assert number_total_constraints(btx) == 38
        assert number_unused_variables(btx) == 0

    @pytest.mark.integration
    def test_units(self, btx):
        assert_units_equivalent(
            btx.fs.unit.overall_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(btx.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(btx.fs.unit.delta_temperature_in, pyunits.K)
        assert_units_equivalent(btx.fs.unit.delta_temperature_out, pyunits.K)

        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "HX Area": btx.fs.unit.area,
                "Heat Duty": btx.fs.unit.heat_duty[0],
                "HX Coefficient": btx.fs.unit.overall_heat_transfer_coefficient[0],
            },
            "exprs": {
                "Delta T Driving": btx.fs.unit.delta_temperature[0],
                "Delta T In": btx.fs.unit.delta_temperature_in[0],
                "Delta T Out": btx.fs.unit.delta_temperature_out[0],
            },
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "flow_mol": getattr(pyunits.pint_registry, "mole/second"),
                    "mole_frac_comp benzene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "mole_frac_comp toluene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "temperature": getattr(pyunits.pint_registry, "kelvin"),
                    "pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "Hot Side Inlet": {
                    "flow_mol": 5.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 365,
                    "pressure": 101325.0,
                },
                "Hot Side Outlet": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 298.15,
                    "pressure": 101325.0,
                },
                "Cold Side Inlet": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 300,
                    "pressure": 101325.0,
                },
                "Cold Side Outlet": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 298.15,
                    "pressure": 101325.0,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(5, abs=1e-3) == value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(359.5, abs=1e-1) == value(
            btx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(1, abs=1e-3) == value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(329.9, abs=1e-1) == value(
            btx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.hot_side_inlet.flow_mol[0]
                    - btx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    btx.fs.unit.cold_side_inlet.flow_mol[0]
                    - btx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.hot_side.properties_in[0].enth_mol_phase["Liq"]
                - btx.fs.unit.hot_side.properties_out[0].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.cold_side.properties_in[0].enth_mol_phase["Liq"]
                - btx.fs.unit.cold_side.properties_out[0].enth_mol_phase["Liq"]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# -----------------------------------------------------------------------------
class TestBTX_cocurrent_alt_name(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(valid_phase="Liq")

        m.fs.unit = HeatExchanger(
            hot_side_name="hot",
            cold_side_name="cold",
            hot={"property_package": m.fs.properties},
            cold={"property_package": m.fs.properties},
            flow_pattern=HeatExchangerFlowPattern.cocurrent,
        )

        m.fs.unit.hot_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.hot_inlet.temperature[0].fix(365)  # K
        m.fs.unit.hot_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.hot_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.hot_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.cold_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.cold_inlet.temperature[0].fix(300)  # K
        m.fs.unit.cold_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.cold_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.cold_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.area.fix(1)
        m.fs.unit.overall_heat_transfer_coefficient.fix(100)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "hot_inlet")
        assert len(btx.fs.unit.hot_inlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_inlet, "temperature")
        assert hasattr(btx.fs.unit.hot_inlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_inlet")
        assert len(btx.fs.unit.cold_inlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_inlet, "temperature")
        assert hasattr(btx.fs.unit.cold_inlet, "pressure")

        assert hasattr(btx.fs.unit, "hot_outlet")
        assert len(btx.fs.unit.hot_outlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_outlet, "temperature")
        assert hasattr(btx.fs.unit.hot_outlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_outlet")
        assert len(btx.fs.unit.cold_outlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_outlet, "temperature")
        assert hasattr(btx.fs.unit.cold_outlet, "pressure")

        assert isinstance(btx.fs.unit.overall_heat_transfer_coefficient, Var)
        assert isinstance(btx.fs.unit.area, Var)
        assert not hasattr(btx.fs.unit, "crossflow_factor")
        assert isinstance(btx.fs.unit.heat_duty, Var)
        assert isinstance(btx.fs.unit.delta_temperature_in, Var)
        assert isinstance(btx.fs.unit.delta_temperature_out, Var)
        assert isinstance(btx.fs.unit.unit_heat_balance, Constraint)
        assert isinstance(btx.fs.unit.delta_temperature, (Var, Expression))
        assert isinstance(btx.fs.unit.heat_transfer_equation, Constraint)

        assert number_variables(btx) == 50
        assert number_total_constraints(btx) == 38
        assert number_unused_variables(btx) == 0

    @pytest.mark.integration
    def test_units(self, btx):
        assert_units_equivalent(
            btx.fs.unit.overall_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(btx.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(btx.fs.unit.delta_temperature_in, pyunits.K)
        assert_units_equivalent(btx.fs.unit.delta_temperature_out, pyunits.K)

        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "HX Area": btx.fs.unit.area,
                "Heat Duty": btx.fs.unit.heat_duty[0],
                "HX Coefficient": btx.fs.unit.overall_heat_transfer_coefficient[0],
            },
            "exprs": {
                "Delta T Driving": btx.fs.unit.delta_temperature[0],
                "Delta T In": btx.fs.unit.delta_temperature_in[0],
                "Delta T Out": btx.fs.unit.delta_temperature_out[0],
            },
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "flow_mol": getattr(pyunits.pint_registry, "mole/second"),
                    "mole_frac_comp benzene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "mole_frac_comp toluene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "temperature": getattr(pyunits.pint_registry, "kelvin"),
                    "pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "hot Inlet": {
                    "flow_mol": 5.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 365,
                    "pressure": 101325.0,
                },
                "hot Outlet": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 298.15,
                    "pressure": 101325.0,
                },
                "cold Inlet": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 300,
                    "pressure": 101325.0,
                },
                "cold Outlet": {
                    "flow_mol": 1.0,
                    "mole_frac_comp benzene": 0.5,
                    "mole_frac_comp toluene": 0.5,
                    "temperature": 298.15,
                    "pressure": 101325.0,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(5, abs=1e-3) == value(btx.fs.unit.hot_outlet.flow_mol[0])
        assert pytest.approx(359.5, abs=1e-1) == value(
            btx.fs.unit.hot_outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx.fs.unit.hot_outlet.pressure[0]
        )

        assert pytest.approx(1, abs=1e-3) == value(btx.fs.unit.cold_outlet.flow_mol[0])
        assert pytest.approx(329.9, abs=1e-1) == value(
            btx.fs.unit.cold_outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx.fs.unit.cold_outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.hot_inlet.flow_mol[0]
                    - btx.fs.unit.hot_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    btx.fs.unit.cold_inlet.flow_mol[0]
                    - btx.fs.unit.cold_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            btx.fs.unit.hot_outlet.flow_mol[0]
            * (
                btx.fs.unit.hot.properties_in[0].enth_mol_phase["Liq"]
                - btx.fs.unit.hot.properties_out[0].enth_mol_phase["Liq"]
            )
        )
        cold_side = value(
            btx.fs.unit.cold_outlet.flow_mol[0]
            * (
                btx.fs.unit.cold.properties_in[0].enth_mol_phase["Liq"]
                - btx.fs.unit.cold.properties_out[0].enth_mol_phase["Liq"]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS_countercurrent(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = HeatExchanger(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
        )

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(100)
        m.fs.unit.hot_side_inlet.enth_mol[0].fix(6000)
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(100)
        m.fs.unit.cold_side_inlet.enth_mol[0].fix(5000)
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

        m.fs.unit.area.fix(1000)
        m.fs.unit.overall_heat_transfer_coefficient.fix(100)

        return m

    @pytest.fixture(scope="class")
    def iapws_underwood(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = HeatExchanger(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            delta_temperature_callback=delta_temperature_underwood_callback,
            flow_pattern=HeatExchangerFlowPattern.countercurrent,
        )

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(100)
        m.fs.unit.hot_side_inlet.enth_mol[0].fix(6000)
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(100)
        m.fs.unit.cold_side_inlet.enth_mol[0].fix(5000)
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)

        m.fs.unit.area.fix(1000)
        m.fs.unit.overall_heat_transfer_coefficient.fix(100)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, iapws):
        assert len(iapws.fs.unit.hot_side_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.hot_side_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "hot_side_outlet")
        assert len(iapws.fs.unit.hot_side_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.hot_side_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.hot_side_outlet, "pressure")

        assert len(iapws.fs.unit.cold_side_inlet.vars) == 3
        assert hasattr(iapws.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(iapws.fs.unit.cold_side_inlet, "enth_mol")
        assert hasattr(iapws.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(iapws.fs.unit, "cold_side_outlet")
        assert len(iapws.fs.unit.cold_side_outlet.vars) == 3
        assert hasattr(iapws.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.cold_side_outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.cold_side_outlet, "pressure")

        assert isinstance(iapws.fs.unit.overall_heat_transfer_coefficient, Var)
        assert isinstance(iapws.fs.unit.area, Var)
        assert not hasattr(iapws.fs.unit, "crossflow_factor")
        assert isinstance(iapws.fs.unit.heat_duty, Var)
        assert isinstance(iapws.fs.unit.delta_temperature_in, Var)
        assert isinstance(iapws.fs.unit.delta_temperature_out, Var)
        assert isinstance(iapws.fs.unit.unit_heat_balance, Constraint)
        assert isinstance(iapws.fs.unit.delta_temperature, (Expression, Var))
        assert isinstance(iapws.fs.unit.heat_transfer_equation, Constraint)

        assert number_variables(iapws) == 18
        assert number_total_constraints(iapws) == 10
        assert number_unused_variables(iapws) == 0

    @pytest.mark.integration
    def test_units(self, iapws):
        assert_units_equivalent(
            iapws.fs.unit.overall_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(iapws.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(iapws.fs.unit.delta_temperature_in, pyunits.K)
        assert_units_equivalent(iapws.fs.unit.delta_temperature_out, pyunits.K)

        assert_units_consistent(iapws)

    @pytest.mark.unit
    def test_dof(self, iapws):
        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.unit
    def test_dof_alt_name1(self, iapws):
        iapws.fs.unit.hot_side_inlet.flow_mol[0].fix(100)
        iapws.fs.unit.hot_side_inlet.enth_mol[0].fix(6000)
        iapws.fs.unit.hot_side_inlet.pressure[0].fix(101325)

        iapws.fs.unit.cold_side_inlet.flow_mol[0].fix(100)
        iapws.fs.unit.cold_side_inlet.enth_mol[0].fix(5000)
        iapws.fs.unit.cold_side_inlet.pressure[0].fix(101325)

        iapws.fs.unit.area.fix(1000)
        iapws.fs.unit.overall_heat_transfer_coefficient.fix(100)

        assert degrees_of_freedom(iapws) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "HX Area": iapws.fs.unit.area,
                "Heat Duty": iapws.fs.unit.heat_duty[0],
                "HX Coefficient": iapws.fs.unit.overall_heat_transfer_coefficient[0],
            },
            "exprs": {
                "Delta T Driving": iapws.fs.unit.delta_temperature[0],
                "Delta T In": iapws.fs.unit.delta_temperature_in[0],
                "Delta T Out": iapws.fs.unit.delta_temperature_out[0],
            },
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, iapws):
        stable = iapws.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "Molar Flow": getattr(pyunits.pint_registry, "mole/second"),
                    "Mass Flow": getattr(pyunits.pint_registry, "kg/second"),
                    "T": getattr(pyunits.pint_registry, "K"),
                    "P": getattr(pyunits.pint_registry, "Pa"),
                    "Vapor Fraction": getattr(pyunits.pint_registry, "dimensionless"),
                    "Molar Enthalpy": getattr(pyunits.pint_registry, "J/mole"),
                },
                "Hot Side Inlet": {
                    "Molar Flow": 100,
                    "Mass Flow": 1.8015,
                    "T": 352.67,
                    "P": 101325,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 6000.0,
                },
                "Hot Side Outlet": {
                    "Molar Flow": 1,
                    "Mass Flow": 1.8015e-2,
                    "T": 270.4877112932641,
                    "P": 11032305.8275,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 0.01102138712926277,
                },
                "Cold Side Inlet": {
                    "Molar Flow": 100,
                    "Mass Flow": 1.8015,
                    "T": 339.43,
                    "P": 101325,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 5000.0,
                },
                "Cold Side Outlet": {
                    "Molar Flow": 1,
                    "Mass Flow": 1.8015e-2,
                    "T": 270.4877112932641,
                    "P": 11032305.8275,
                    "Vapor Fraction": 0,
                    "Molar Enthalpy": 0.01102138712926277,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_underwood(self, iapws_underwood):
        initialization_tester(iapws_underwood)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_underwood(self, iapws_underwood):
        results = solver.solve(iapws_underwood)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(100, abs=1e-5) == value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(100, abs=1e-5) == value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
        )

        assert pytest.approx(5070.219, abs=1e0) == value(
            iapws.fs.unit.hot_side_outlet.enth_mol[0]
        )
        assert pytest.approx(5929.781, abs=1e0) == value(
            iapws.fs.unit.cold_side_outlet.enth_mol[0]
        )

        assert pytest.approx(101325, abs=1e2) == value(
            iapws.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(101325, abs=1e2) == value(
            iapws.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, iapws):
        assert (
            abs(
                value(
                    iapws.fs.unit.hot_side_inlet.flow_mol[0]
                    - iapws.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    iapws.fs.unit.cold_side_inlet.flow_mol[0]
                    - iapws.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            iapws.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                iapws.fs.unit.hot_side_inlet.enth_mol[0]
                - iapws.fs.unit.hot_side_outlet.enth_mol[0]
            )
        )
        cold_side = value(
            iapws.fs.unit.cold_side_outlet.flow_mol[0]
            * (
                iapws.fs.unit.cold_side_inlet.enth_mol[0]
                - iapws.fs.unit.cold_side_outlet.enth_mol[0]
            )
        )
        assert abs(hot_side + cold_side) <= 1e-6


# -----------------------------------------------------------------------------
# @pytest.mark.skip(reason="Solutions vary with differnt versions of solver.")
class TestSaponification_crossflow(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = HeatExchanger(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_pattern=HeatExchangerFlowPattern.crossflow,
        )

        m.fs.unit.hot_side_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.hot_side_inlet.temperature[0].fix(320)
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.hot_side_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.cold_side_inlet.flow_vol[0].fix(1e-3)
        m.fs.unit.cold_side_inlet.temperature[0].fix(300)
        m.fs.unit.cold_side_inlet.pressure[0].fix(101325)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.cold_side_inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.area.fix(1000)
        m.fs.unit.overall_heat_transfer_coefficient.fix(100)
        m.fs.unit.crossflow_factor.fix(0.6)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert len(sapon.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.hot_side_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(sapon.fs.unit.hot_side_inlet, "pressure")

        assert len(sapon.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.hot_side_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(sapon.fs.unit.hot_side_outlet, "pressure")

        assert len(sapon.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(sapon.fs.unit.cold_side_inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(sapon.fs.unit.cold_side_inlet, "pressure")

        assert len(sapon.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(sapon.fs.unit.cold_side_outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(sapon.fs.unit.cold_side_outlet, "pressure")

        assert isinstance(sapon.fs.unit.overall_heat_transfer_coefficient, Var)
        assert isinstance(sapon.fs.unit.area, Var)
        assert isinstance(sapon.fs.unit.crossflow_factor, Var)
        assert isinstance(sapon.fs.unit.heat_duty, Var)
        assert isinstance(sapon.fs.unit.delta_temperature_in, Var)
        assert isinstance(sapon.fs.unit.delta_temperature_out, Var)
        assert isinstance(sapon.fs.unit.unit_heat_balance, Constraint)
        assert isinstance(sapon.fs.unit.delta_temperature, (Expression, Var))
        assert isinstance(sapon.fs.unit.heat_transfer_equation, Constraint)

        assert number_variables(sapon) == 39
        assert number_total_constraints(sapon) == 20
        assert number_unused_variables(sapon) == 0

    @pytest.mark.integration
    def test_units(self, sapon):
        assert_units_equivalent(
            sapon.fs.unit.overall_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(sapon.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(sapon.fs.unit.delta_temperature_in, pyunits.K)
        assert_units_equivalent(sapon.fs.unit.delta_temperature_out, pyunits.K)

        assert_units_consistent(sapon)

    @pytest.mark.unit
    def test_dof(self, sapon):
        assert degrees_of_freedom(sapon) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "HX Area": sapon.fs.unit.area,
                "Heat Duty": sapon.fs.unit.heat_duty[0],
                "HX Coefficient": sapon.fs.unit.overall_heat_transfer_coefficient[0],
                "Crossflow Factor": sapon.fs.unit.crossflow_factor[0],
            },
            "exprs": {
                "Delta T Driving": sapon.fs.unit.delta_temperature[0],
                "Delta T In": sapon.fs.unit.delta_temperature_in[0],
                "Delta T Out": sapon.fs.unit.delta_temperature_out[0],
            },
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, sapon):
        stable = sapon.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "Volumetric Flowrate": getattr(
                        pyunits.pint_registry, "m**3/second"
                    ),
                    "Molar Concentration H2O": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration NaOH": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration EthylAcetate": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration SodiumAcetate": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Molar Concentration Ethanol": getattr(
                        pyunits.pint_registry, "mole/m**3"
                    ),
                    "Temperature": getattr(pyunits.pint_registry, "K"),
                    "Pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "Hot Side Inlet": {
                    "Volumetric Flowrate": 1e-3,
                    "Molar Concentration H2O": 55388,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 0,
                    "Molar Concentration Ethanol": 0,
                    "Temperature": 320,
                    "Pressure": 1.0132e05,
                },
                "Hot Side Outlet": {
                    "Volumetric Flowrate": 1.00,
                    "Molar Concentration H2O": 100.00,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 100.00,
                    "Molar Concentration Ethanol": 100.00,
                    "Temperature": 298.15,
                    "Pressure": 1.0132e05,
                },
                "Cold Side Inlet": {
                    "Volumetric Flowrate": 1e-3,
                    "Molar Concentration H2O": 55388,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 0,
                    "Molar Concentration Ethanol": 0,
                    "Temperature": 300,
                    "Pressure": 1.0132e05,
                },
                "Cold Side Outlet": {
                    "Volumetric Flowrate": 1.00,
                    "Molar Concentration H2O": 100.00,
                    "Molar Concentration NaOH": 100.00,
                    "Molar Concentration EthylAcetate": 100.00,
                    "Molar Concentration SodiumAcetate": 100.00,
                    "Molar Concentration Ethanol": 100.00,
                    "Temperature": 298.15,
                    "Pressure": 1.0132e05,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(1e-3, abs=1e-6) == value(
            sapon.fs.unit.hot_side_outlet.flow_vol[0]
        )
        assert pytest.approx(1e-3, abs=1e-6) == value(
            sapon.fs.unit.cold_side_outlet.flow_vol[0]
        )

        assert pytest.approx(55388.0, rel=1e-3) == value(
            sapon.fs.unit.hot_side_outlet.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.hot_side_outlet.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.hot_side_outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.hot_side_outlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.hot_side_outlet.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(55388.0, rel=1e-3) == value(
            sapon.fs.unit.cold_side_outlet.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.cold_side_outlet.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.cold_side_outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.cold_side_outlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.cold_side_outlet.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(301.3, abs=1e-1) == value(
            sapon.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(318.7, abs=1e-1) == value(
            sapon.fs.unit.cold_side_outlet.temperature[0]
        )

        assert pytest.approx(101325, abs=1e2) == value(
            sapon.fs.unit.hot_side_outlet.pressure[0]
        )
        assert pytest.approx(101325, abs=1e2) == value(
            sapon.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        hot_side = value(
            sapon.fs.unit.hot_side_outlet.flow_vol[0]
            * sapon.fs.properties.dens_mol
            * sapon.fs.properties.cp_mol
            * (
                sapon.fs.unit.hot_side_inlet.temperature[0]
                - sapon.fs.unit.hot_side_outlet.temperature[0]
            )
        )
        cold_side = value(
            sapon.fs.unit.cold_side_outlet.flow_vol[0]
            * sapon.fs.properties.dens_mol
            * sapon.fs.properties.cp_mol
            * (
                sapon.fs.unit.cold_side_inlet.temperature[0]
                - sapon.fs.unit.cold_side_outlet.temperature[0]
            )
        )
        assert abs(hot_side + cold_side) <= 1e0


# -----------------------------------------------------------------------------
class TestBT_Generic_cocurrent(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        # As we lack other example prop packs with units, take the generic
        # BT-PR package and change the base units
        configuration2 = {
            # Specifying components
            "components": {
                "benzene": {
                    "type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (78.1136e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (48.9e5, pyunits.Pa),  # [1]
                        "temperature_crit": (562.2, pyunits.K),  # [1]
                        "omega": 0.212,  # [1]
                        "cp_mol_ig_comp_coeff": {
                            "A": (-3.392e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (4.739e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-3.017e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (7.130e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "enth_mol_form_vap_comp_ref": (
                            82.9e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            -269,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-6.98273, None),  # [1]
                            "B": (1.33213, None),
                            "C": (-2.62863, None),
                            "D": (-3.33399, None),
                        },
                    },
                },
                "toluene": {
                    "type": Component,
                    "enth_mol_ig_comp": RPP,
                    "entr_mol_ig_comp": RPP,
                    "pressure_sat_comp": RPP,
                    "phase_equilibrium_form": {("Vap", "Liq"): log_fugacity},
                    "parameter_data": {
                        "mw": (92.1405e-3, pyunits.kg / pyunits.mol),  # [1]
                        "pressure_crit": (41e5, pyunits.Pa),  # [1]
                        "temperature_crit": (591.8, pyunits.K),  # [1]
                        "omega": 0.263,  # [1]
                        "cp_mol_ig_comp_coeff": {
                            "A": (-2.435e1, pyunits.J / pyunits.mol / pyunits.K),  # [1]
                            "B": (5.125e-1, pyunits.J / pyunits.mol / pyunits.K**2),
                            "C": (-2.765e-4, pyunits.J / pyunits.mol / pyunits.K**3),
                            "D": (4.911e-8, pyunits.J / pyunits.mol / pyunits.K**4),
                        },
                        "enth_mol_form_vap_comp_ref": (
                            50.1e3,
                            pyunits.J / pyunits.mol,
                        ),  # [3]
                        "entr_mol_form_vap_comp_ref": (
                            -321,
                            pyunits.J / pyunits.mol / pyunits.K,
                        ),  # [3]
                        "pressure_sat_comp_coeff": {
                            "A": (-7.28607, None),  # [1]
                            "B": (1.38091, None),
                            "C": (-2.83433, None),
                            "D": (-2.79168, None),
                        },
                    },
                },
            },
            # Specifying phases
            "phases": {
                "Liq": {
                    "type": LiquidPhase,
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
                "Vap": {
                    "type": VaporPhase,
                    "equation_of_state": Cubic,
                    "equation_of_state_options": {"type": CubicType.PR},
                },
            },
            # Set base units of measurement
            "base_units": {
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.t,
                "amount": pyunits.mol,
                "temperature": pyunits.degR,
            },
            # Specifying state definition
            "state_definition": FTPx,
            "state_bounds": {
                "flow_mol": (0, 100, 1000, pyunits.mol / pyunits.s),
                "temperature": (273.15, 300, 500, pyunits.K),
                "pressure": (5e4, 1e5, 1e6, pyunits.Pa),
            },
            "pressure_ref": (101325, pyunits.Pa),
            "temperature_ref": (298.15, pyunits.K),
            # Defining phase equilibria
            "phases_in_equilibrium": [("Vap", "Liq")],
            "phase_equilibrium_state": {("Vap", "Liq"): SmoothVLE},
            "bubble_dew_method": LogBubbleDew,
            "parameter_data": {
                "PR_kappa": {
                    ("benzene", "benzene"): 0.000,
                    ("benzene", "toluene"): 0.000,
                    ("toluene", "benzene"): 0.000,
                    ("toluene", "toluene"): 0.000,
                }
            },
        }

        m.fs.properties = GenericParameterBlock(**configuration)
        m.fs.properties2 = GenericParameterBlock(**configuration2)

        m.fs.unit = HeatExchanger(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties2},
            flow_pattern=HeatExchangerFlowPattern.cocurrent,
        )

        m.fs.unit.hot_side_inlet.flow_mol[0].fix(5)  # mol/s
        m.fs.unit.hot_side_inlet.temperature[0].fix(365)  # K
        m.fs.unit.hot_side_inlet.pressure[0].fix(101325)  # Pa
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.hot_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.cold_side_inlet.flow_mol[0].fix(1)  # mol/s
        m.fs.unit.cold_side_inlet.temperature[0].fix(540)  # degR
        m.fs.unit.cold_side_inlet.pressure[0].fix(101.325)  # kPa
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.cold_side_inlet.mole_frac_comp[0, "toluene"].fix(0.5)

        m.fs.unit.area.fix(1)
        m.fs.unit.overall_heat_transfer_coefficient.fix(100)

        m.fs.unit.cold_side.scaling_factor_pressure = 1

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "hot_side_inlet")
        assert len(btx.fs.unit.hot_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_inlet")
        assert len(btx.fs.unit.cold_side_inlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_inlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_inlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_inlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_inlet, "pressure")

        assert hasattr(btx.fs.unit, "hot_side_outlet")
        assert len(btx.fs.unit.hot_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.hot_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.hot_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.hot_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.hot_side_outlet, "pressure")

        assert hasattr(btx.fs.unit, "cold_side_outlet")
        assert len(btx.fs.unit.cold_side_outlet.vars) == 4
        assert hasattr(btx.fs.unit.cold_side_outlet, "flow_mol")
        assert hasattr(btx.fs.unit.cold_side_outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.cold_side_outlet, "temperature")
        assert hasattr(btx.fs.unit.cold_side_outlet, "pressure")

        assert isinstance(btx.fs.unit.overall_heat_transfer_coefficient, Var)
        assert isinstance(btx.fs.unit.area, Var)
        assert not hasattr(btx.fs.unit, "crossflow_factor")
        assert isinstance(btx.fs.unit.heat_duty, Var)
        assert isinstance(btx.fs.unit.delta_temperature_in, Var)
        assert isinstance(btx.fs.unit.delta_temperature_out, Var)
        assert isinstance(btx.fs.unit.unit_heat_balance, Constraint)
        assert isinstance(btx.fs.unit.delta_temperature, (Var, Expression))
        assert isinstance(btx.fs.unit.heat_transfer_equation, Constraint)

        assert number_variables(btx) == 190
        assert number_total_constraints(btx) == 118
        assert number_unused_variables(btx) == 20

    @pytest.mark.integration
    def test_units(self, btx):
        assert_units_equivalent(
            btx.fs.unit.overall_heat_transfer_coefficient,
            pyunits.W / pyunits.m**2 / pyunits.K,
        )
        assert_units_equivalent(btx.fs.unit.area, pyunits.m**2)
        assert_units_equivalent(btx.fs.unit.delta_temperature_in, pyunits.K)
        assert_units_equivalent(btx.fs.unit.delta_temperature_out, pyunits.K)

        assert_units_consistent(btx)

    @pytest.mark.unit
    def test_dof(self, btx):
        assert degrees_of_freedom(btx) == 0

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "HX Area": btx.fs.unit.area,
                "Heat Duty": btx.fs.unit.heat_duty[0],
                "HX Coefficient": btx.fs.unit.overall_heat_transfer_coefficient[0],
            },
            "exprs": {
                "Delta T Driving": btx.fs.unit.delta_temperature[0],
                "Delta T In": btx.fs.unit.delta_temperature_in[0],
                "Delta T Out": btx.fs.unit.delta_temperature_out[0],
            },
        }

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = pandas.DataFrame.from_dict(
            {
                "Units": {
                    "Total Molar Flowrate": getattr(
                        pyunits.pint_registry, "mole/second"
                    ),
                    "Total Mole Fraction benzene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "Total Mole Fraction toluene": getattr(
                        pyunits.pint_registry, "dimensionless"
                    ),
                    "Temperature": getattr(pyunits.pint_registry, "kelvin"),
                    "Pressure": getattr(pyunits.pint_registry, "Pa"),
                },
                "Hot Side Inlet": {
                    "Total Molar Flowrate": 5.0,
                    "Total Mole Fraction benzene": 0.5,
                    "Total Mole Fraction toluene": 0.5,
                    "Temperature": 365,
                    "Pressure": 101325.0,
                },
                "Hot Side Outlet": {
                    "Total Molar Flowrate": 100.0,
                    "Total Mole Fraction benzene": 0.5,
                    "Total Mole Fraction toluene": 0.5,
                    "Temperature": 300,
                    "Pressure": 1e5,
                },
                "Cold Side Inlet": {
                    "Total Molar Flowrate": 1.0,
                    "Total Mole Fraction benzene": 0.5,
                    "Total Mole Fraction toluene": 0.5,
                    "Temperature": 300,
                    "Pressure": 101325.0,
                },
                "Cold Side Outlet": {
                    "Total Molar Flowrate": 100.0,
                    "Total Mole Fraction benzene": 0.5,
                    "Total Mole Fraction toluene": 0.5,
                    "Temperature": 300,
                    "Pressure": 1e5,
                },
            }
        )

        pandas.testing.assert_frame_equal(stable, expected, rtol=1e-4, atol=1e-4)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(5, abs=1e-3) == value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
        )
        assert pytest.approx(359.4, abs=1e-1) == value(
            btx.fs.unit.hot_side_outlet.temperature[0]
        )
        assert pytest.approx(101325, abs=1e-3) == value(
            btx.fs.unit.hot_side_outlet.pressure[0]
        )

        assert pytest.approx(1, abs=1e-3) == value(
            btx.fs.unit.cold_side_outlet.flow_mol[0]
        )
        assert pytest.approx(596.9, abs=1e-1) == value(
            btx.fs.unit.cold_side_outlet.temperature[0]
        )
        assert pytest.approx(101.325, abs=1e-3) == value(
            btx.fs.unit.cold_side_outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, btx):
        assert (
            abs(
                value(
                    btx.fs.unit.hot_side_inlet.flow_mol[0]
                    - btx.fs.unit.hot_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    btx.fs.unit.cold_side_inlet.flow_mol[0]
                    - btx.fs.unit.cold_side_outlet.flow_mol[0]
                )
            )
            <= 1e-6
        )

        hot_side = value(
            btx.fs.unit.hot_side_outlet.flow_mol[0]
            * (
                btx.fs.unit.hot_side.properties_in[0].enth_mol
                - btx.fs.unit.hot_side.properties_out[0].enth_mol
            )
        )
        cold_side = pyunits.convert_value(
            value(
                btx.fs.unit.cold_side_outlet.flow_mol[0]
                * (
                    btx.fs.unit.cold_side.properties_in[0].enth_mol
                    - btx.fs.unit.cold_side.properties_out[0].enth_mol
                )
            ),
            from_units=pyunits.kJ / pyunits.s,
            to_units=pyunits.J / pyunits.s,
        )
        assert abs(hot_side + cold_side) <= 1e-6

    @pytest.mark.component
    def test_initialization_error(self, btx):
        btx.fs.unit.hot_side_outlet.flow_mol[0].fix(20)

        with pytest.raises(InitializationError):
            btx.fs.unit.initialize()
