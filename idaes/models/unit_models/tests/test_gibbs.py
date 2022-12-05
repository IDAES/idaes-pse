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
Tests for ControlVolumeBlockData.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    value,
    units,
)
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import FlowsheetBlock, EnergyBalanceType, MomentumBalanceType
from idaes.models.unit_models.gibbs_reactor import GibbsReactor
from idaes.models.properties.activity_coeff_models.methane_combustion_ideal import (
    MethaneParameterBlock as MethaneCombustionParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    fixed_variables_set,
    activated_constraints_set,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import ConfigurationError


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = GibbsReactor(property_package=m.fs.properties)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 9

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.inert_species == []

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 4

    assert not hasattr(m.fs.unit, "inert_species_balance")


@pytest.mark.unit
def test_inerts():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

    assert isinstance(m.fs.unit.inert_species_balance, Constraint)
    assert len(m.fs.unit.inert_species_balance) == 2
    assert m.fs.unit.inert_species_balance[0, "p1", "c1"] != Constraint.Skip
    assert m.fs.unit.inert_species_balance[0, "p2", "c1"] != Constraint.Skip

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 2


@pytest.mark.unit
def test_inerts_dependent_w_multi_phase():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    # Change elemental composition to introduce dependency
    m.fs.properties.element_comp = {
        "c1": {"H": 0, "He": 0, "Li": 3},
        "c2": {"H": 4, "He": 5, "Li": 0},
    }

    m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

    assert isinstance(m.fs.unit.inert_species_balance, Constraint)
    assert len(m.fs.unit.inert_species_balance) == 2
    assert m.fs.unit.inert_species_balance[0, "p1", "c1"] != Constraint.Skip
    assert m.fs.unit.inert_species_balance[0, "p2", "c1"] != Constraint.Skip

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 2


@pytest.mark.unit
def test_inerts_dependent_w_single_phase():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    # Set phase list to only have 1 phase
    m.fs.properties.phase_list = ["p1"]
    # Change elemental composition to introduce dependency
    m.fs.properties.element_comp = {
        "c1": {"H": 0, "He": 0, "Li": 3},
        "c2": {"H": 4, "He": 5, "Li": 0},
    }

    m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

    assert isinstance(m.fs.unit.inert_species_balance, Constraint)
    assert len(m.fs.unit.inert_species_balance) == 0
    assert (0, "p1", "c1") not in m.fs.unit.inert_species_balance

    assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
    assert len(m.fs.unit.gibbs_minimization) == 1


@pytest.mark.unit
def test_invalid_inert():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    with pytest.raises(
        ConfigurationError,
        match="fs.unit invalid component in inert_species "
        "argument. foo is not in the property package "
        "component list.",
    ):
        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties, inert_species=["foo"]
        )


# -----------------------------------------------------------------------------
class TestMethane(object):
    @pytest.fixture(scope="class")
    def methane(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = MethaneCombustionParameterBlock()

        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_mol[0].fix(230.0)
        m.fs.unit.inlet.mole_frac_comp[0, "H2"].fix(0.0435)
        m.fs.unit.inlet.mole_frac_comp[0, "N2"].fix(0.6522)
        m.fs.unit.inlet.mole_frac_comp[0, "O2"].fix(0.1739)
        m.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "CH4"].fix(0.1304)
        m.fs.unit.inlet.mole_frac_comp[0, "CO"].fix(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "NH3"].fix(1e-5)
        m.fs.unit.inlet.temperature[0].fix(1500.0)
        m.fs.unit.inlet.pressure[0].fix(101325.0)

        m.fs.unit.outlet.temperature[0].fix(2844.38)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, methane):
        assert hasattr(methane.fs.unit, "inlet")
        assert len(methane.fs.unit.inlet.vars) == 4
        assert hasattr(methane.fs.unit.inlet, "flow_mol")
        assert hasattr(methane.fs.unit.inlet, "mole_frac_comp")
        assert hasattr(methane.fs.unit.inlet, "temperature")
        assert hasattr(methane.fs.unit.inlet, "pressure")

        assert hasattr(methane.fs.unit, "outlet")
        assert len(methane.fs.unit.outlet.vars) == 4
        assert hasattr(methane.fs.unit.outlet, "flow_mol")
        assert hasattr(methane.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(methane.fs.unit.outlet, "temperature")
        assert hasattr(methane.fs.unit.outlet, "pressure")

        assert hasattr(methane.fs.unit, "gibbs_minimization")
        assert hasattr(methane.fs.unit, "heat_duty")
        assert hasattr(methane.fs.unit, "deltaP")

        assert number_variables(methane) == 80
        assert number_total_constraints(methane) == 67
        assert number_unused_variables(methane) == 0

    @pytest.mark.component
    def test_units(self, methane):
        assert_units_consistent(methane)
        assert_units_equivalent(methane.fs.unit.heat_duty[0], units.W)
        assert_units_equivalent(methane.fs.unit.deltaP[0], units.Pa)

    @pytest.mark.unit
    def test_dof(self, methane):
        assert degrees_of_freedom(methane) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_temperature(self, methane):
        initialization_tester(
            methane,
            optarg={"tol": 1e-6},
            state_args={
                "temperature": 2844.38,
                "pressure": 101325.0,
                "flow_mol": 251.05,
                "mole_frac_comp": {
                    "CH4": 1e-5,
                    "CO": 0.0916,
                    "CO2": 0.0281,
                    "H2": 0.1155,
                    "H2O": 0.1633,
                    "N2": 0.5975,
                    "NH3": 1e-5,
                    "O2": 0.0067,
                },
            },
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_temperature(self, methane):
        results = solver.solve(methane)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution_temperature(self, methane):
        assert pytest.approx(250.06, abs=1e-2) == value(
            methane.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.0974, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.0226, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.1030, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.1769, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.5999, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.0002, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(methane.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            methane.fs.unit.outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation_temperature(self, methane):
        assert (
            abs(
                value(
                    methane.fs.unit.inlet.flow_mol[0]
                    * methane.fs.unit.control_volume.properties_in[0].enth_mol_phase[
                        "Vap"
                    ]
                    - methane.fs.unit.outlet.flow_mol[0]
                    * methane.fs.unit.control_volume.properties_out[0].enth_mol_phase[
                        "Vap"
                    ]
                    + methane.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize_duty(self, methane):
        methane.fs.unit.outlet.temperature[0].unfix()
        methane.fs.unit.heat_duty.fix(-7454077)
        assert degrees_of_freedom(methane) == 0

        orig_fixed_vars = fixed_variables_set(methane)
        orig_act_consts = activated_constraints_set(methane)

        methane.fs.unit.initialize(
            optarg={"tol": 1e-6},
            state_args={
                "temperature": 2844.38,
                "pressure": 101325.0,
                "flow_mol": 251.05,
                "mole_frac_comp": {
                    "CH4": 1e-5,
                    "CO": 0.0916,
                    "CO2": 0.0281,
                    "H2": 0.1155,
                    "H2O": 0.1633,
                    "N2": 0.5975,
                    "NH3": 1e-5,
                    "O2": 0.0067,
                },
            },
        )

        assert degrees_of_freedom(methane) == 0

        fin_fixed_vars = fixed_variables_set(methane)
        fin_act_consts = activated_constraints_set(methane)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve_heat_duty(self, methane):
        solver.options["tol"] = 1e-9
        solver.options["nlp_scaling_method"] = "user-scaling"

        results = solver.solve(methane, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution_duty(self, methane):
        assert pytest.approx(250.06, abs=1e-1) == value(
            methane.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.103863, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.016123, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.096080, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.183897, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.600030, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(2.86950e-06, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(methane.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            methane.fs.unit.outlet.pressure[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation_duty(self, methane):
        assert (
            abs(
                value(
                    methane.fs.unit.inlet.flow_mol[0]
                    * methane.fs.unit.control_volume.properties_in[0].enth_mol_phase[
                        "Vap"
                    ]
                    - methane.fs.unit.outlet.flow_mol[0]
                    * methane.fs.unit.control_volume.properties_out[0].enth_mol_phase[
                        "Vap"
                    ]
                    + methane.fs.unit.heat_duty[0]
                )
            )
            <= 1e-4
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, methane):
        perf_dict = methane.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Heat Duty": methane.fs.unit.heat_duty[0],
                "Pressure Change": methane.fs.unit.deltaP[0],
            }
        }
