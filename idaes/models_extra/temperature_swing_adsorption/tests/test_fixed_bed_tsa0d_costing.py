#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for the direct air capture costing model.
"""

__author__ = "Alex Noring, Ryan Hughes"

import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, value

from idaes.core import FlowsheetBlock
from idaes.models_extra.temperature_swing_adsorption import (
    FixedBedTSA0D,
    FixedBedTSA0DInitializer,
    Adsorbent,
    TransformationScheme,
    SteamCalculationType,
)
from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
from idaes.models_extra.temperature_swing_adsorption.costing.dac_costing import (
    get_dac_costing,
    print_dac_costing,
    dac_costing_summary,
    get_dac_costing_data,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@pytest.fixture(scope="class")
def base_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.compressor_props = FlueGasParameterBlock(components=["N2", "CO2"])
    m.fs.steam_props = iapws95.Iapws95ParameterBlock()

    m.fs.unit = FixedBedTSA0D(
        adsorbent=Adsorbent.zeolite_13x,
        number_of_beds=600,
        transformation_method="dae.collocation",
        transformation_scheme=TransformationScheme.lagrangeRadau,
        finite_elements=20,
        collocation_points=6,
        compressor=True,
        compressor_properties=m.fs.compressor_props,
        steam_calculation=SteamCalculationType.rigorous,
        steam_properties=m.fs.steam_props,
    )

    m.fs.unit.inlet.flow_mol_comp[0, "H2O"].fix(0)
    m.fs.unit.inlet.flow_mol_comp[0, "CO2"].fix(40)
    m.fs.unit.inlet.flow_mol_comp[0, "N2"].fix(99960)
    m.fs.unit.inlet.flow_mol_comp[0, "O2"].fix(0)
    m.fs.unit.inlet.temperature.fix(303.15)
    m.fs.unit.inlet.pressure.fix(100000)

    m.fs.unit.temperature_desorption.fix(470)
    m.fs.unit.temperature_adsorption.fix(310)
    m.fs.unit.temperature_heating.fix(500)
    m.fs.unit.temperature_cooling.fix(300)
    m.fs.unit.bed_diameter.fix(4)
    m.fs.unit.bed_height.fix(8)
    m.fs.unit.compressor.unit.efficiency_isentropic.fix(0.8)

    iscale.calculate_scaling_factors(m)

    initializer = FixedBedTSA0DInitializer()
    initializer.initialize(m.fs.unit)

    return m


class TestCostingCaseConfigs:
    @pytest.mark.unit
    def test_get_dac_costing_data_invalid_case(self):
        with pytest.raises(ConfigurationError, match="costing case not defined"):
            get_dac_costing_data("not_a_valid_case")

    @pytest.mark.unit
    def test_get_dac_costing_data_electric_boiler(self):
        data = get_dac_costing_data("electric_boiler")

        assert isinstance(data, dict)
        assert len(data) > 0

    @pytest.mark.unit
    def test_get_dac_costing_data_retrofit_ngcc(self):
        data = get_dac_costing_data("retrofit_ngcc")

        assert isinstance(data, dict)
        assert len(data) > 0

    @pytest.mark.integration
    def test_get_dac_costing_invalid_case(self, base_model):
        with pytest.raises(ConfigurationError, match="costing case not defined"):
            get_dac_costing(base_model.fs.unit, "not_a_valid_case")


@pytest.mark.integration
class TestElectricBoilerCosting:
    @pytest.fixture(scope="class")
    def model(self, base_model):
        get_dac_costing(base_model.fs.unit, "electric_boiler")

        return base_model

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_build(self, model):
        assert hasattr(model.fs, "costing")

        cost_attr_list = ["total_TPC", "total_fixed_OM_cost", "total_variable_OM_cost"]
        for cost_attr in cost_attr_list:
            assert hasattr(model.fs.costing, cost_attr)

        assert number_variables(model) == 3008
        assert number_total_constraints(model) == 2977
        assert number_unused_variables(model) == 12

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_EB_cost_accounts(self, model):
        assert model.fs.raw_water_system.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["3.2", "3.4", "9.5", "14.6"]

        assert model.fs.steam_system.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["3.1", "3.3", "3.5"]

        assert model.fs.cooling_tower.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["9.1"]

        assert model.fs.water_discharge_system.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["3.7"]

        assert model.fs.cooling_water_system.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"]

        assert model.fs.electric_systems.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == [
            "11.1",
            "11.2",
            "11.3",
            "11.4",
            "11.5",
            "11.6",
            "11.7",
            "11.8",
            "11.9",
            "12.4",
            "12.5",
            "12.6",
            "12.7",
            "12.8",
            "12.9",
            "13.1",
            "13.2",
            "13.3",
            "14.4",
            "14.7",
            "14.8",
            "14.9",
            "14.10",
        ]

        assert model.fs.electric_boiler.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.9"]

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_shared_DAC_cost_accounts(self, model):
        assert model.fs.sorbent_makeup.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == [
            "1.5",
            "1.6",
            "1.7",
            "1.8",
            "1.9",
            "2.5",
            "2.6",
            "2.9",
            "10.6",
            "10.7",
            "10.9",
        ]

        assert model.fs.vessels.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.1"]

        assert model.fs.product_compression.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.2"]

        assert model.fs.compressor_aftercooler.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.3"]

        assert model.fs.duct_dampers.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.4"]

        assert model.fs.feed_fans.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.5"]

        assert model.fs.desorption_gas_handling.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.6"]

        assert model.fs.steam_distribution.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.7"]

        assert model.fs.controls_equipment.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.8"]

        assert model.fs.CO2_storage_vessel.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.10"]

        assert model.fs.CO2_dryer.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["15.11"]

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, model):
        solver.options.bound_push = 1e-6
        results = solver.solve(model)
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, model):

        assert pytest.approx(75.1400, abs=1e-4) == value(
            model.fs.costing.annualized_cost
        )
        assert pytest.approx(27.3995, abs=1e-4) == value(
            model.fs.costing.total_fixed_OM_cost
        )
        assert pytest.approx(333.336, abs=1e-3) == value(
            model.fs.costing.total_variable_OM_cost[0]
        )
        assert pytest.approx(0.118581, abs=1e-6) == value(
            model.fs.costing.cost_of_capture
        )

    @pytest.mark.ui
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_report(self, model):
        print_dac_costing(model.fs.unit)
        dac_costing_summary(model.fs.unit)


@pytest.mark.integration
class TestRetrofitNgccCosting:
    @pytest.fixture(scope="class")
    def model(self, base_model):
        get_dac_costing(base_model.fs.unit, "retrofit_ngcc")

        return base_model

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_build(self, model):
        assert hasattr(model.fs, "costing")

        cost_attr_list = ["total_TPC", "total_fixed_OM_cost", "total_variable_OM_cost"]
        for cost_attr in cost_attr_list:
            assert hasattr(model.fs.costing, cost_attr)

        assert number_variables(model) == 2962
        assert number_total_constraints(model) == 2931
        assert number_unused_variables(model) == 12

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_Ngcc_cost_accounts(self, model):
        assert model.fs.gas_flow_to_dac.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["7.3"]

        assert model.fs.steam_flow_system.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == ["8.4"]

        assert model.fs.electric_systems.costing.config.costing_method_arguments[
            "cost_accounts"
        ] == [
            "11.2",
            "11.3",
            "11.4",
            "11.5",
            "11.6",
            "12.1",
            "12.2",
            "12.3",
            "12.4",
            "12.5",
            "12.6",
            "12.7",
            "12.8",
            "12.9",
        ]

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, model):
        solver.options.bound_push = 1e-6
        results = solver.solve(model)
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, model):

        assert pytest.approx(77.3557, abs=1e-4) == value(
            model.fs.costing.annualized_cost
        )
        assert pytest.approx(28.2422, abs=1e-4) == value(
            model.fs.costing.total_fixed_OM_cost
        )
        assert pytest.approx(263.751, abs=1e-3) == value(
            model.fs.costing.total_variable_OM_cost[0]
        )
        assert pytest.approx(0.101344, abs=1e-6) == value(
            model.fs.costing.cost_of_capture
        )
