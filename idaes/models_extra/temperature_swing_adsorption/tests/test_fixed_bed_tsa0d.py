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
Tests for conventional fixed bed TSA 0D model.
Author: Alex Noring
"""
import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, Expression, Var, value
from pyomo.util.check_units import assert_units_consistent


from idaes.core import FlowsheetBlock
from idaes.models_extra.temperature_swing_adsorption.fixed_bed_tsa0d import FixedBedTSA0D
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)

from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver
import idaes.core.util.scaling as iscale

import idaes.logger as idaeslog

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

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        initialization_tester(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)
        assert check_optimal_termination(results)

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


class TestTsaMgmof:
    # also testing calculate beds and simple steam calcs
    pass


class TestTsaPolystyrene:
    # also testing compressor and rigorous steam calcs
    pass


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
    compressor=False,
    steam_calculation=None)

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
#m.fs.unit.velocity_in.fix(0.42)
iscale.calculate_scaling_factors(m)


"""
# build model
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.unit = FixedBedTSA0D(
    adsorbent="mmen-Mg-MOF-74",
    calculate_beds=False,
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
#m.fs.unit.velocity_in.fix(0.42)
iscale.calculate_scaling_factors(m)

#m.fs.unit.initialize()

# test_build
# make sure the zeolite vars exist
assert not hasattr(m.fs.unit, "c_ref")
assert hasattr(m.fs.unit, "nL_inf")
assert not hasattr(m.fs.unit, "qm0")
# make sure calculate beds is right
assert isinstance(m.fs.unit.velocity_in, Var)
# make sure there is no compressor or heater
assert not hasattr(m.fs.unit, "compressor")
assert hasattr(m.fs.unit, "flow_mass_steam")
assert not hasattr(m.fs.unit, "steam_heater")

assert number_variables(m) == 2840
assert number_total_constraints(m) == 2816
assert number_unused_variables(m) == 11

# test_dof
assert degrees_of_freedom(m) == 0

# test_units
#assert_units_consistent(m)

# test_initialize
initialization_tester(m)

# test_solve
results = solver.solve(m)
assert check_optimal_termination(results)

# test_solution
# CO2 rich outlet
assert pytest.approx(1.1486e-3, abs=1e-6) == value(
    m.fs.unit.co2_rich_stream.flow_mol_comp[0, "CO2"]
)
assert pytest.approx(1.6541e-5, abs=1e-8) == value(
    m.fs.unit.co2_rich_stream.flow_mol_comp[0, "N2"]
)
assert pytest.approx(430, abs=1e-1) == value(
    m.fs.unit.co2_rich_stream.temperature[0]
)
assert pytest.approx(100000, abs=1e-1) == value(
    m.fs.unit.co2_rich_stream.pressure[0]
)
# N2 rich outlet
assert pytest.approx(3.3229e-6, abs=1e-9) == value(
    m.fs.unit.n2_rich_stream.flow_mol_comp[0, "CO2"]
)
assert pytest.approx(8.4314e-3, abs=1e-6) == value(
    m.fs.unit.n2_rich_stream.flow_mol_comp[0, "N2"]
)
assert pytest.approx(310, abs=1e-1) == value(
    m.fs.unit.n2_rich_stream.temperature[0]
)
assert pytest.approx(100000, abs=1e-1) == value(
    m.fs.unit.n2_rich_stream.pressure[0]
)
# air velocity
assert pytest.approx(2.44178, abs=1e-3) == value(
    m.fs.unit.velocity_in
)"""


# what problems do we have?
# inconsistent units on mg-mof and polystyrine
# had some initial difficulties initializing mg-mof but now its working
# Problem with initialization when calculate_beds=True, velocity_in gets unfixed during initialization and pressure drop gets fixed, so it fails test initialization
