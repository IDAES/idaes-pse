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
Tests for 0D lumped capacitance heat exchanger model

Author: Rusty Gentile, John Eslick, Andrew Lee
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    value,
    Var,
    Param,
    units as pyunits,
)

from pyomo.dae import DerivativeVar
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from pyomo.core.base.units_container import InconsistentUnitsError

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError, IdaesError

from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester

from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from idaes.models.properties import iapws95
from idaes.models.unit_models import (
    HeatExchangerLumpedCapacitance,
    HeatExchangerFlowPattern,
)

from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback

# Get default solver for testing
solver = get_solver()

# Number of steps for transient simulations
TIME_STEPS = 50


class TestHXRegression(object):
    """
    Regression tests taken from the base 0D HeatExchanger model
    """

    @pytest.mark.unit
    def test_config(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            hot_side_name="shell",
            cold_side_name="tube",
            shell={"property_package": m.fs.properties},
            tube={"property_package": m.fs.properties},
            dynamic_heat_balance=False,
        )

        # Check unit config arguments
        # There are 10 to 12 arguments since you can add a side 1 and 2 config by
        # side_1, side_2, or whatever the user named them, plus
        # "dynamic_heat_balance"
        assert len(m.fs.unit.config) >= 10 and len(m.fs.unit.config) <= 12

        assert not m.fs.unit.config.dynamic
        assert not m.fs.unit.config.has_holdup
        assert isinstance(m.fs.unit.config.shell, ConfigBlock)
        assert isinstance(m.fs.unit.config.tube, ConfigBlock)
        assert (
            m.fs.unit.config.delta_temperature_callback
            is delta_temperature_lmtd_callback
        )
        assert m.fs.unit.config.flow_pattern == HeatExchangerFlowPattern.countercurrent

        # Check shell config
        assert len(m.fs.unit.config.shell) == 7
        assert (
            m.fs.unit.config.shell.material_balance_type
            == MaterialBalanceType.useDefault
        )
        assert (
            m.fs.unit.config.shell.energy_balance_type == EnergyBalanceType.useDefault
        )
        assert (
            m.fs.unit.config.shell.momentum_balance_type
            == MomentumBalanceType.pressureTotal
        )
        assert not m.fs.unit.config.shell.has_phase_equilibrium
        assert not m.fs.unit.config.shell.has_pressure_change
        assert m.fs.unit.config.shell.property_package is m.fs.properties

        # Check tube config
        assert len(m.fs.unit.config.tube) == 7
        assert (
            m.fs.unit.config.tube.material_balance_type
            == MaterialBalanceType.useDefault
        )
        assert m.fs.unit.config.tube.energy_balance_type == EnergyBalanceType.useDefault
        assert (
            m.fs.unit.config.tube.momentum_balance_type
            == MomentumBalanceType.pressureTotal
        )
        assert not m.fs.unit.config.tube.has_phase_equilibrium
        assert not m.fs.unit.config.tube.has_pressure_change
        assert m.fs.unit.config.tube.property_package is m.fs.properties

    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            flow_pattern=HeatExchangerFlowPattern.crossflow,
            dynamic_heat_balance=False,
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

        # Modified from the original:
        # m.fs.unit.overall_heat_transfer_coefficient.fix(100)
        m.fs.unit.ua_hot_side.fix(200 * 1000)
        m.fs.unit.ua_cold_side.fix(200 * 1000)

        m.fs.unit.crossflow_factor.fix(0.6)

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_sapon_solution(self, sapon):

        results = solver.solve(sapon)
        check_optimal_termination(results)

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


class TestHXLCGeneric(object):
    @pytest.fixture()
    def static_flowsheet_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.properties = iapws95.Iapws95ParameterBlock()
        return m

    @pytest.fixture()
    def dynamic_flowsheet_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 1], time_units=pyunits.s)
        m.fs.properties = iapws95.Iapws95ParameterBlock()
        return m

    @pytest.fixture()
    def unconstrained_model(self, dynamic_flowsheet_model):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 1], time_units=pyunits.s)
        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            hot_side_name="shell",
            cold_side_name="tube",
            hot_side={"property_package": m.fs.properties},
            cold_side={"property_package": m.fs.properties},
            dynamic=False,
        )

        return m

    @pytest.fixture()
    def model(self, unconstrained_model):

        m = unconstrained_model
        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=1, wrt=m.fs.time, scheme="BACKWARD")

        m.fs.unit.area.fix(1000)
        m.fs.unit.heat_capacity = 100

        m.fs.unit.shell_inlet.flow_mol[:].fix(100)
        m.fs.unit.shell_inlet.pressure[:].fix(101325)
        m.fs.unit.shell_inlet.enth_mol[:].fix(4000)

        m.fs.unit.tube_inlet.flow_mol[:].fix(100)
        m.fs.unit.tube_inlet.pressure[:].fix(101325)
        m.fs.unit.tube_inlet.enth_mol[:].fix(3000)

        m.fs.unit.ua_cold_side[:].fix(100)
        m.fs.unit.ua_hot_side[:].fix(100)

        return m

    @pytest.mark.unit
    @pytest.mark.xfail(raises=InconsistentUnitsError)
    def test_units(self, model):
        # Note: using the discretizer makes the units of measure on the time
        # derivative term inconsistent...
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_units_unconstrained(self, unconstrained_model):

        # ...but without the discretizer, the units are fine
        assert_units_consistent(unconstrained_model)

        m = unconstrained_model
        assert_units_equivalent(m.fs.unit.ua_hot_side, pyunits.W / pyunits.K)
        assert_units_equivalent(m.fs.unit.ua_cold_side, pyunits.W / pyunits.K)
        assert_units_equivalent(m.fs.unit.ua_hot_side_to_wall, pyunits.W / pyunits.K)
        assert_units_equivalent(m.fs.unit.temperature_wall, pyunits.K)
        assert_units_equivalent(m.fs.unit.heat_capacity_wall, pyunits.J / pyunits.K)
        assert_units_equivalent(
            m.fs.unit.thermal_resistance_wall, pyunits.K / pyunits.W
        )
        assert_units_equivalent(
            m.fs.unit.thermal_fouling_hot_side, pyunits.K / pyunits.W
        )
        assert_units_equivalent(
            m.fs.unit.thermal_fouling_cold_side, pyunits.K / pyunits.W
        )
        assert_units_equivalent(m.fs.unit.dT_wall_dt, pyunits.K / pyunits.s)

    @pytest.mark.unit
    def test_base_vars(self, model):
        assert isinstance(model.fs.unit.area, Var)
        assert isinstance(model.fs.unit.overall_heat_transfer_coefficient, Var)
        assert isinstance(model.fs.unit.heat_duty, Var)
        assert isinstance(model.fs.unit.delta_temperature_in, Var)
        assert isinstance(model.fs.unit.delta_temperature_out, Var)

    @pytest.mark.unit
    def test_new_vars(self, model):
        assert isinstance(model.fs.unit.dT_wall_dt, DerivativeVar)
        assert isinstance(model.fs.unit.temperature_wall, Var)
        assert isinstance(model.fs.unit.ua_cold_side, Var)
        assert isinstance(model.fs.unit.ua_hot_side, Var)
        assert isinstance(model.fs.unit.ua_hot_side_to_wall, Var)
        assert isinstance(model.fs.unit.heat_capacity_wall, Param)
        assert isinstance(model.fs.unit.thermal_resistance_wall, Param)
        assert isinstance(model.fs.unit.thermal_fouling_hot_side, Param)
        assert isinstance(model.fs.unit.thermal_fouling_cold_side, Param)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):

        model.fs.unit.initialize()
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.unit
    def test_static_flowsheet(self, static_flowsheet_model):
        m = static_flowsheet_model
        with pytest.raises(
            ConfigurationError,
            match="dynamic heat balance cannot be used with a "
            "steady-state flowsheet",
        ):
            m.fs.unit = HeatExchangerLumpedCapacitance(
                hot_side_name="shell",
                cold_side_name="tube",
                shell={"property_package": m.fs.properties},
                tube={"property_package": m.fs.properties},
                dynamic=False,
                dynamic_heat_balance=True,
            )

    @pytest.mark.unit
    def test_static_heat_balance(self, static_flowsheet_model):
        m = static_flowsheet_model
        m.fs.unit = HeatExchangerLumpedCapacitance(
            hot_side_name="shell",
            cold_side_name="tube",
            shell={"property_package": m.fs.properties},
            tube={"property_package": m.fs.properties},
            dynamic=False,
            dynamic_heat_balance=False,
        )

        with pytest.raises(
            IdaesError,
            match="heat holdup term cannot be activated when "
            "`dynamic_heat_balance=False`",
        ):
            m.fs.unit.activate_dynamic_heat_eq()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "HX Area": model.fs.unit.area,
                "Heat Duty": model.fs.unit.heat_duty[0],
                "HX Coefficient": model.fs.unit.overall_heat_transfer_coefficient[0],
            },
            "exprs": {
                "Delta T Driving": model.fs.unit.delta_temperature[0],
                "Delta T In": model.fs.unit.delta_temperature_in[0],
                "Delta T Out": model.fs.unit.delta_temperature_out[0],
            },
        }
