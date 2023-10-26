#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for the direct air capture costing model.
"""

__author__ = "Alex Noring"

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
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestElectricBoilerCosting:
    @pytest.fixture(scope="class")
    def model(self):
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

        get_dac_costing(m.fs.unit, "electric_boiler")

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_build(self, model):
        assert hasattr(model.fs, "costing")
        assert hasattr(model.fs.costing, "total_TPC")
        assert hasattr(model.fs.costing, "total_fixed_OM_cost")
        assert hasattr(model.fs.costing, "total_variable_OM_cost")

        assert number_variables(model) == 3007
        assert number_total_constraints(model) == 2976
        assert number_unused_variables(model) == 12

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        solver.options.bound_push = 1e-6
        results = solver.solve(model)
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # CO2 rich outlet
        assert pytest.approx(79.1955, abs=1e-4) == value(
            model.fs.costing.annualized_cost
        )
        assert pytest.approx(27.1399, abs=1e-4) == value(
            model.fs.costing.total_fixed_OM_cost
        )
        assert pytest.approx(1090.57, abs=1e-2) == value(
            model.fs.costing.total_variable_OM_cost[0]
        )
        assert pytest.approx(0.317544, abs=1e-6) == value(
            model.fs.costing.cost_of_capture
        )

    @pytest.mark.ui
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_report(self, model):
        print_dac_costing(model.fs.unit)
        dac_costing_summary(model.fs.unit)


class TestRetrofitNgccCosting:
    @pytest.fixture(scope="class")
    def model(self):
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

        get_dac_costing(m.fs.unit, "retrofit_ngcc")

        return m

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_build(self, model):
        assert hasattr(model.fs, "costing")
        assert hasattr(model.fs.costing, "total_TPC")
        assert hasattr(model.fs.costing, "total_fixed_OM_cost")
        assert hasattr(model.fs.costing, "total_variable_OM_cost")

        assert number_variables(model) == 2961
        assert number_total_constraints(model) == 2930
        assert number_unused_variables(model) == 12

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        solver.options.bound_push = 1e-6
        results = solver.solve(model)
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # CO2 rich outlet
        assert pytest.approx(81.1308, abs=1e-4) == value(
            model.fs.costing.annualized_cost
        )
        assert pytest.approx(27.8940, abs=1e-4) == value(
            model.fs.costing.total_fixed_OM_cost
        )
        assert pytest.approx(1020.94, abs=1e-2) == value(
            model.fs.costing.total_variable_OM_cost[0]
        )
        assert pytest.approx(0.300183, abs=1e-6) == value(
            model.fs.costing.cost_of_capture
        )
