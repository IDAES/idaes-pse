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

from idaes.models.properties import swco2, iapws95
from idaes.power_generation.properties import FlueGasParameterBlock
from idaes.models.unit_models import (
    HeatExchangerLumpedCapacitance,
    HeatExchangerFlowPattern,
)

from idaes.models.unit_models.heat_exchanger import delta_temperature_lmtd_callback
import numpy as np

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
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            default={
                "shell": {"property_package": m.fs.properties},
                "tube": {"property_package": m.fs.properties},
                "dynamic_heat_balance": False,
            }
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

    @pytest.fixture()
    def basic_model(self):

        cb = delta_temperature_lmtd_callback
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = iapws95.Iapws95ParameterBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            default={
                "shell": {"property_package": m.fs.properties},
                "tube": {"property_package": m.fs.properties},
                "delta_temperature_callback": cb,
                "flow_pattern": HeatExchangerFlowPattern.countercurrent,
                "dynamic_heat_balance": False,
            }
        )

        #   Set inputs
        m.fs.unit.inlet_1.flow_mol[0].fix(100)
        m.fs.unit.inlet_1.enth_mol[0].fix(4000)
        m.fs.unit.inlet_1.pressure[0].fix(101325)

        m.fs.unit.inlet_2.flow_mol[0].fix(100)
        m.fs.unit.inlet_2.enth_mol[0].fix(3500)
        m.fs.unit.inlet_2.pressure[0].fix(101325)

        m.fs.unit.area.fix(1000)

        # Modified from the original:
        # m.fs.unit.overall_heat_transfer_coefficient.fix(100)
        m.fs.unit.ua_hot_side.fix(200 * 1000)
        m.fs.unit.ua_cold_side.fix(200 * 1000)

        assert degrees_of_freedom(m) == 0
        m.fs.unit.get_costing()
        m.fs.unit.initialize()
        return m

    @pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.unit
    def test_costing(self, basic_model):

        m = basic_model  # (delta_temperature_lmtd_callback)

        assert m.fs.unit.costing.purchase_cost.value == pytest.approx(529738.6793, 1e-5)

        assert_units_consistent(m.fs.unit.costing)

        results = solver.solve(m)

        # Check for optimal solution
        assert check_optimal_termination(results)

        # Check Solution

        # hot in end
        assert value(m.fs.unit.delta_temperature_in[0]) == pytest.approx(
            0.464879, rel=1e-3
        )
        # hot out end
        assert value(m.fs.unit.delta_temperature_out[0]) == pytest.approx(
            0.465069, rel=1e-3
        )
        assert value(m.fs.unit.heat_duty[0]) == pytest.approx(46497.44)

        # Costing
        assert m.fs.unit.costing.purchase_cost.value == pytest.approx(529738.6793, 1e-5)

    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            default={
                "shell": {"property_package": m.fs.properties},
                "tube": {"property_package": m.fs.properties},
                "flow_pattern": HeatExchangerFlowPattern.crossflow,
                "dynamic_heat_balance": False,
            }
        )

        m.fs.unit.inlet_1.flow_vol[0].fix(1e-3)
        m.fs.unit.inlet_1.temperature[0].fix(320)
        m.fs.unit.inlet_1.pressure[0].fix(101325)
        m.fs.unit.inlet_1.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet_1.conc_mol_comp[0, "Ethanol"].fix(0.0)

        m.fs.unit.inlet_2.flow_vol[0].fix(1e-3)
        m.fs.unit.inlet_2.temperature[0].fix(300)
        m.fs.unit.inlet_2.pressure[0].fix(101325)
        m.fs.unit.inlet_2.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
        m.fs.unit.inlet_2.conc_mol_comp[0, "Ethanol"].fix(0.0)

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
            sapon.fs.unit.outlet_1.flow_vol[0]
        )
        assert pytest.approx(1e-3, abs=1e-6) == value(
            sapon.fs.unit.outlet_2.flow_vol[0]
        )

        assert pytest.approx(55388.0, rel=1e-3) == value(
            sapon.fs.unit.outlet_1.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.outlet_1.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.outlet_1.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.outlet_1.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.outlet_1.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(55388.0, rel=1e-3) == value(
            sapon.fs.unit.outlet_2.conc_mol_comp[0, "H2O"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.outlet_2.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100.0, rel=1e-3) == value(
            sapon.fs.unit.outlet_2.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.outlet_2.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(0.0, abs=1e-3) == value(
            sapon.fs.unit.outlet_2.conc_mol_comp[0, "Ethanol"]
        )

        assert pytest.approx(301.3, abs=1e-1) == value(
            sapon.fs.unit.outlet_1.temperature[0]
        )
        assert pytest.approx(318.7, abs=1e-1) == value(
            sapon.fs.unit.outlet_2.temperature[0]
        )

        assert pytest.approx(101325, abs=1e2) == value(
            sapon.fs.unit.outlet_1.pressure[0]
        )
        assert pytest.approx(101325, abs=1e2) == value(
            sapon.fs.unit.outlet_2.pressure[0]
        )


class TestHXLCGeneric(object):
    @pytest.fixture()
    def static_flowsheet_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(default={"dynamic": False})
        m.fs.properties = iapws95.Iapws95ParameterBlock()
        return m

    @pytest.fixture()
    def dynamic_flowsheet_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(
            default={"dynamic": True, "time_set": [0, 1], "time_units": pyunits.s}
        )
        m.fs.properties = iapws95.Iapws95ParameterBlock()
        return m

    @pytest.fixture()
    def unconstrained_model(self, dynamic_flowsheet_model):
        m = dynamic_flowsheet_model
        m.fs.unit = HeatExchangerLumpedCapacitance(
            default={
                "shell": {"property_package": m.fs.properties},
                "tube": {"property_package": m.fs.properties},
                "dynamic": False,
            }
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
                default={
                    "shell": {"property_package": m.fs.properties},
                    "tube": {"property_package": m.fs.properties},
                    "dynamic": False,
                    "dynamic_heat_balance": True,
                }
            )

    @pytest.mark.unit
    def test_static_heat_balance(self, static_flowsheet_model):
        m = static_flowsheet_model
        m.fs.unit = HeatExchangerLumpedCapacitance(
            default={
                "shell": {"property_package": m.fs.properties},
                "tube": {"property_package": m.fs.properties},
                "dynamic": False,
                "dynamic_heat_balance": False,
            }
        )

        with pytest.raises(
            IdaesError,
            match="heat holdup term cannot be activated when "
            "`dynamic_heat_balance=False`",
        ):
            m.fs.unit.activate_dynamic_heat_eq()

    @pytest.mark.unit
    def test_build_valid_configs(self, static_flowsheet_model, dynamic_flowsheet_model):
        def build_unit_model(m, dyn, dyn_hb):
            return HeatExchangerLumpedCapacitance(
                default={
                    "shell": {"property_package": m.fs.properties},
                    "tube": {"property_package": m.fs.properties},
                    "dynamic": dyn,
                    "dynamic_heat_balance": dyn_hb,
                }
            )

        # Static model
        ms = static_flowsheet_model
        ms.fs.unit = build_unit_model(ms, False, False)

        # Dynamic model
        md = dynamic_flowsheet_model
        md.fs.u1 = build_unit_model(ms, False, False)
        md.fs.u2 = build_unit_model(ms, True, False)
        md.fs.u3 = build_unit_model(ms, False, True)
        md.fs.u4 = build_unit_model(ms, True, True)
        md.discretizer = TransformationFactory("dae.finite_difference")
        md.discretizer.apply_to(md, nfe=10, wrt=md.fs.time, scheme="BACKWARD")

        assert True


class TestHXLCTransientSCO2(object):
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(
            default={
                "dynamic": True,
                "time_set": [0, 300, 600, 900, 1200, 1500],
                "time_units": pyunits.s,
            }
        )

        m.fs.sco2 = swco2.SWCO2ParameterBlock()
        m.fs.fluegas = FlueGasParameterBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            default={
                "delta_temperature_callback": delta_temperature_lmtd_callback,
                "hot_side_name": "tube",
                "cold_side_name": "shell",
                "tube": {"property_package": m.fs.sco2, "has_pressure_change": True},
                "shell": {
                    "property_package": m.fs.fluegas,
                    "has_pressure_change": False,
                },
                "flow_pattern": HeatExchangerFlowPattern.crossflow,
                "dynamic": False,
            }
        )

        m.discretizer = TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=TIME_STEPS - 1, wrt=m.fs.time, scheme="BACKWARD")

        # Cold-side boundary conditions
        shell_inlet_temperature = 288.15
        shell_flow = 44004.14222
        shell_area = 690073.9153
        shell_hconv = 22
        m.fs.unit.ua_cold_side[:].fix(shell_area * shell_hconv)
        m.fs.unit.shell_inlet.flow_mol_comp[:, "H2O"].fix(0.01027 * shell_flow)
        m.fs.unit.shell_inlet.flow_mol_comp[:, "CO2"].fix(0.000411592 * shell_flow)
        m.fs.unit.shell_inlet.flow_mol_comp[:, "N2"].fix(0.780066026 * shell_flow)
        m.fs.unit.shell_inlet.flow_mol_comp[:, "O2"].fix(0.209252382 * shell_flow)
        m.fs.unit.shell_inlet.flow_mol_comp[:, "NO"].fix(0)
        m.fs.unit.shell_inlet.flow_mol_comp[:, "SO2"].fix(0)
        m.fs.unit.shell_inlet.temperature[:].fix(shell_inlet_temperature)
        m.fs.unit.shell_inlet.pressure[:].fix(101325)

        # Hot-side boundary conditions
        tube_inlet_temperature = 384.35
        tube_inlet_pressure = 7653000
        tube_outlet_pressure = 7500000
        tube_flow = 13896.84163
        tube_area = 19542.2771
        tube_hconv = 1000
        tube_mass = 1160 * 322
        m.fs.unit.ua_hot_side[:].fix(tube_area * tube_hconv)
        tube_inlet_enthalpy = swco2.htpx(
            T=tube_inlet_temperature * pyunits.K, P=tube_inlet_pressure * pyunits.Pa
        )
        m.fs.unit.tube_inlet.flow_mol[:].fix(tube_flow)
        m.fs.unit.tube_inlet.pressure[:].fix(tube_inlet_pressure)
        m.fs.unit.tube_inlet.enth_mol[:].fix(tube_inlet_enthalpy)
        m.fs.unit.tube_outlet.pressure[:].fix(tube_outlet_pressure)

        # Area has no effect on the results here so long as it's not zero
        m.fs.unit.area.fix(1)
        m.fs.unit.crossflow_factor.fix(0.8)
        m.fs.unit.heat_capacity_wall = tube_mass * 466

        m.fs.unit.tube_outlet[:].enth_mol.setub(tube_inlet_enthalpy)
        m.fs.unit.shell_outlet.temperature[:].setlb(shell_inlet_temperature)

        return m

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_steady_state(self, model):
        """
        Make sure temperatures are all the same for steady state simulation
        """
        model.fs.unit.initialize()
        solver.solve(model)

        # Expected temperatures
        exp_sco2 = np.ones(TIME_STEPS) * 305.2
        exp_air = np.ones(TIME_STEPS) * 370.4
        exp_wall = np.ones(TIME_STEPS) * 339.7

        self.check_temperatures(
            model, np.array(model.fs.time), exp_sco2, exp_air, exp_wall
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.integration
    def test_steps(self, model):
        """
        Add step changes in air temperature similar to the example in the docs
        """
        model.fs.unit.initialize()

        # Add disturbances
        for t in model.fs.time:
            if 300 <= t < 600:
                model.fs.unit.shell_inlet.temperature[t].fix(288.15 - 10)
            elif 600 <= t < 900:
                model.fs.unit.shell_inlet.temperature[t].fix(288.15)
            elif 900 <= t < 1200:
                model.fs.unit.shell_inlet.temperature[t].fix(288.15 + 10)
            elif t >= 1200:
                model.fs.unit.shell_inlet.temperature[t].fix(288.15)

        # Transient solution
        solver.solve(model)

        times = [0, 300, 600, 900, 1200, 1500]
        sco2_exp = [305.2, 304.9, 305.1, 306.5, 305.7, 305.2]
        air_exp = [370.4, 373.1, 370.3, 365.9, 370.7, 370.4]
        wall_exp = [339.4, 338.7, 339.1, 340.7, 339.9, 339.4]

        self.check_temperatures(model, times, sco2_exp, air_exp, wall_exp)

    @pytest.mark.static
    def check_temperatures(
        self,
        model,
        times,
        expected_sco2_temps,
        expected_air_temps,
        expected_wall_temps,
        tol=1e-3,
    ):

        for i, t in enumerate(times):
            assert value(
                model.fs.unit.tube.properties_out[t].temperature
            ) == pytest.approx(expected_sco2_temps[i], rel=tol)
            assert value(
                model.fs.unit.shell.properties_out[t].temperature
            ) == pytest.approx(expected_air_temps[i], rel=tol)
            assert value(model.fs.unit.temperature_wall[t]) == pytest.approx(
                expected_wall_temps[i], rel=tol
            )
