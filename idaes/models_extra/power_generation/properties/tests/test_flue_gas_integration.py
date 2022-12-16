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
    ConcreteModel,
    TransformationFactory,
    value,
    units as pyunits,
)

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

from idaes.core import FlowsheetBlock

from idaes.models.properties import swco2
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
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


class TestHXLCTransientSCO2(object):
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(
            dynamic=True, time_set=[0, 300, 600, 900, 1200, 1500], time_units=pyunits.s
        )

        m.fs.sco2 = swco2.SWCO2ParameterBlock()
        m.fs.fluegas = FlueGasParameterBlock()

        m.fs.unit = HeatExchangerLumpedCapacitance(
            delta_temperature_callback=delta_temperature_lmtd_callback,
            hot_side_name="tube",
            cold_side_name="shell",
            tube={"property_package": m.fs.sco2, "has_pressure_change": True},
            shell={"property_package": m.fs.fluegas, "has_pressure_change": False},
            flow_pattern=HeatExchangerFlowPattern.crossflow,
            dynamic=False,
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
