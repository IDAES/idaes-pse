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
Tests for fixed bed TSA 0D model.
"""

__author__ = "Alex Noring"

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Expression,
    Var,
    value
)
from pyomo.util.check_units import assert_units_consistent


from idaes.core import FlowsheetBlock
from idaes.models_extra.temperature_swing_adsorption.fixed_bed_tsa0d import \
    FixedBedTSA0D
from idaes.models_extra.temperature_swing_adsorption.util import (
    plot_tsa_profiles, tsa_summary
)

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestTsaZeolite:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.unit = FixedBedTSA0D(
            adsorbent="Zeolite-13X",
            calculate_beds=False,
            number_of_beds=120,
            transformation_method="dae.collocation",
            transformation_scheme="LAGRANGE-RADAU",
            finite_elements=20,
            collocation_points=6,
            compressor=False,
            steam_calculation=None)

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
        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        assert hasattr(model.fs.unit, "c_ref")
        assert not hasattr(model.fs.unit, "nL_inf")
        assert not hasattr(model.fs.unit, "qm0")
        assert not hasattr(model.fs.unit, "compressor")
        assert not hasattr(model.fs.unit, "flow_mass_steam")
        assert not hasattr(model.fs.unit, "steam_heater")
        assert isinstance(model.fs.unit.velocity_in, Expression)

        assert number_variables(model) == 2838
        assert number_total_constraints(model) == 2815
        assert number_unused_variables(model) == 11

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # CO2 rich outlet
        assert pytest.approx(2.34462, abs=1e-3) == value(
            model.fs.unit.co2_rich_stream.flow_mol_comp[0, "CO2"]
        )
        assert pytest.approx(62.07595, abs=1e-3) == value(
            model.fs.unit.co2_rich_stream.flow_mol_comp[0, "N2"]
        )
        assert pytest.approx(470, abs=1e-3) == value(
            model.fs.unit.co2_rich_stream.temperature[0]
        )
        assert pytest.approx(100000, abs=1e-3) == value(
            model.fs.unit.co2_rich_stream.pressure[0]
        )
        # N2 rich outlet
        assert pytest.approx(37.65537, abs=1e-3) == value(
            model.fs.unit.n2_rich_stream.flow_mol_comp[0, "CO2"]
        )
        assert pytest.approx(99897.92404, abs=1e-3) == value(
            model.fs.unit.n2_rich_stream.flow_mol_comp[0, "N2"]
        )
        assert pytest.approx(310, abs=1e-3) == value(
            model.fs.unit.n2_rich_stream.temperature[0]
        )
        assert pytest.approx(100000, abs=1e-3) == value(
            model.fs.unit.n2_rich_stream.pressure[0]
        )
        # air velocity
        assert pytest.approx(2.44178, abs=1e-3) == value(
            model.fs.unit.velocity_in
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.fs.unit.report()
        tsa_summary(model.fs.unit)
        plot_tsa_profiles(model.fs.unit)


class TestTsaMgmof:
    # also testing calculate beds and simple steam calcs
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.unit = FixedBedTSA0D(
            adsorbent="mmen-Mg-MOF-74",
            calculate_beds=True,
            number_of_beds=1,
            transformation_method="dae.collocation",
            transformation_scheme="LAGRANGE-RADAU",
            finite_elements=20,
            collocation_points=6,
            compressor=False,
            steam_calculation="simplified")

        m.fs.unit.inlet.flow_mol_comp[0, "H2O"].fix(0)
        m.fs.unit.inlet.flow_mol_comp[0, "CO2"].fix(0.00960*0.12)
        m.fs.unit.inlet.flow_mol_comp[0, "N2"].fix(0.00960*0.88)
        m.fs.unit.inlet.flow_mol_comp[0, "O2"].fix(0)
        m.fs.unit.inlet.temperature.fix(300)
        m.fs.unit.inlet.pressure.fix(1e5)

        m.fs.unit.temperature_desorption.fix(430)
        m.fs.unit.temperature_adsorption.fix(310)
        m.fs.unit.temperature_heating.fix(440)
        m.fs.unit.temperature_cooling.fix(300)
        m.fs.unit.bed_diameter.fix(3/100)
        m.fs.unit.bed_height.fix(1.2)
        m.fs.unit.pressure_drop.fix(-4000)
        iscale.calculate_scaling_factors(m)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        assert not hasattr(model.fs.unit, "c_ref")
        assert hasattr(model.fs.unit, "nL_inf")
        assert not hasattr(model.fs.unit, "qm0")
        assert not hasattr(model.fs.unit, "compressor")
        assert hasattr(model.fs.unit, "flow_mass_steam")
        assert not hasattr(model.fs.unit, "steam_heater")
        assert isinstance(model.fs.unit.velocity_in, Var)

        assert number_variables(model) == 2840
        assert number_total_constraints(model) == 2816
        assert number_unused_variables(model) == 11

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # CO2 rich outlet
        assert pytest.approx(1.1486e-3, abs=1e-6) == value(
            model.fs.unit.co2_rich_stream.flow_mol_comp[0, "CO2"]
        )
        assert pytest.approx(1.6541e-5, abs=1e-8) == value(
            model.fs.unit.co2_rich_stream.flow_mol_comp[0, "N2"]
        )
        assert pytest.approx(430, abs=1e-1) == value(
            model.fs.unit.co2_rich_stream.temperature[0]
        )
        assert pytest.approx(100000, abs=1e-1) == value(
            model.fs.unit.co2_rich_stream.pressure[0]
        )
        # N2 rich outlet
        assert pytest.approx(3.3229e-6, abs=1e-9) == value(
            model.fs.unit.n2_rich_stream.flow_mol_comp[0, "CO2"]
        )
        assert pytest.approx(8.4314e-3, abs=1e-6) == value(
            model.fs.unit.n2_rich_stream.flow_mol_comp[0, "N2"]
        )
        assert pytest.approx(310, abs=1e-1) == value(
            model.fs.unit.n2_rich_stream.temperature[0]
        )
        assert pytest.approx(100000, abs=1e-1) == value(
            model.fs.unit.n2_rich_stream.pressure[0]
        )
        # air velocity
        assert pytest.approx(0.41747, abs=1e-4) == value(
            model.fs.unit.velocity_in
        )


class TestTsaPolystyrene:
    # also testing compressor and rigorous steam calcs
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.unit = FixedBedTSA0D(
            adsorbent="polystyrene-amine",
            calculate_beds=False,
            number_of_beds=1,
            transformation_method="dae.collocation",
            transformation_scheme="LAGRANGE-RADAU",
            finite_elements=20,
            collocation_points=6,
            compressor=True,
            steam_calculation="rigorous")

        m.fs.unit.inlet.flow_mol_comp[0, "H2O"].fix(0)
        m.fs.unit.inlet.flow_mol_comp[0, "CO2"].fix(0.00960*0.12)
        m.fs.unit.inlet.flow_mol_comp[0, "N2"].fix(0.00960*0.88)
        m.fs.unit.inlet.flow_mol_comp[0, "O2"].fix(0)
        m.fs.unit.inlet.temperature.fix(300)
        m.fs.unit.inlet.pressure.fix(1e5)

        m.fs.unit.temperature_desorption.fix(445)
        m.fs.unit.temperature_adsorption.fix(440)
        m.fs.unit.temperature_heating.fix(450)
        m.fs.unit.temperature_cooling.fix(435)
        m.fs.unit.bed_diameter.fix(10/100)
        m.fs.unit.bed_height.fix(1.2)
        m.fs.unit.compressor.unit.efficiency_isentropic.fix(0.8)
        iscale.calculate_scaling_factors(m)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):
        assert not hasattr(model.fs.unit, "c_ref")
        assert not hasattr(model.fs.unit, "nL_inf")
        assert hasattr(model.fs.unit, "qm0")
        assert hasattr(model.fs.unit, "compressor")
        assert hasattr(model.fs.unit, "flow_mass_steam")
        assert hasattr(model.fs.unit, "steam_heater")
        assert isinstance(model.fs.unit.velocity_in, Expression)

        assert number_variables(model) == 2869
        assert number_total_constraints(model) == 2842
        assert number_unused_variables(model) == 11

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
        results = solver.solve(model)
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # CO2 rich outlet
        assert pytest.approx(1.04099e-5, abs=1e-8) == value(
            model.fs.unit.co2_rich_stream.flow_mol_comp[0, "CO2"]
        )
        assert pytest.approx(6.24804e-5, abs=1e-8) == value(
            model.fs.unit.co2_rich_stream.flow_mol_comp[0, "N2"]
        )
        assert pytest.approx(445, abs=1e-1) == value(
            model.fs.unit.co2_rich_stream.temperature[0]
        )
        assert pytest.approx(100000, abs=1e-1) == value(
            model.fs.unit.co2_rich_stream.pressure[0]
        )
        # N2 rich outlet
        assert pytest.approx(1.14159e-3, abs=1e-6) == value(
            model.fs.unit.n2_rich_stream.flow_mol_comp[0, "CO2"]
        )
        assert pytest.approx(8.3855e-3, abs=1e-6) == value(
            model.fs.unit.n2_rich_stream.flow_mol_comp[0, "N2"]
        )
        assert pytest.approx(440, abs=1e-1) == value(
            model.fs.unit.n2_rich_stream.temperature[0]
        )
        assert pytest.approx(100000, abs=1e-1) == value(
            model.fs.unit.n2_rich_stream.pressure[0]
        )
        # air velocity
        assert pytest.approx(0.06388, abs=1e-4) == value(
            model.fs.unit.velocity_in
        )
