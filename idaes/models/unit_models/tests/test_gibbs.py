#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for Gibbs reactor.

Author: Andrew Lee
"""
import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ConcreteModel,
    Constraint,
    Suffix,
    TransformationFactory,
    value,
    units,
)

from idaes.core import FlowsheetBlock, EnergyBalanceType, MomentumBalanceType
from idaes.models.unit_models.gibbs_reactor import GibbsReactor, GibbsReactorScaler
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
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.scaling import AutoScaler, CustomScalerBase, set_scaling_factor
from idaes.core.util.scaling import (
    jacobian_cond,
    extreme_jacobian_rows,
    extreme_jacobian_columns,
)

# Natural gas property package for integration testing
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models_extra.power_generation.properties.natural_gas_PR import get_prop


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver(solver="ipopt_v2", writer_config={"scale_model": True})


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

    assert m.fs.unit.default_scaler is GibbsReactorScaler


class TestGibbsInerts:
    @pytest.mark.unit
    def test_no_inerts(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = GibbsReactor(property_package=m.fs.properties)

        for e in m.fs.unit.lagrange_set:
            assert e in ["H", "He", "Li"]

        assert not hasattr(m.fs.unit, "inert_species_balance")

    @pytest.mark.unit
    def test_inerts(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

        for e in m.fs.unit.lagrange_set:
            assert e in ["H", "He", "Li"]

        assert isinstance(m.fs.unit.inert_species_balance, Constraint)
        assert len(m.fs.unit.inert_species_balance) == 2
        assert m.fs.unit.inert_species_balance[0, "p1", "c1"] != Constraint.Skip
        assert m.fs.unit.inert_species_balance[0, "p2", "c1"] != Constraint.Skip

        assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
        assert len(m.fs.unit.gibbs_minimization) == 2

    @pytest.mark.unit
    def test_inerts_dependent_w_multi_phase(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()
        # Change elemental composition to introduce dependency
        m.fs.properties.element_comp = {
            "c1": {"H": 0, "He": 0, "Li": 3},
            "c2": {"H": 4, "He": 5, "Li": 0},
        }

        m.fs.unit = GibbsReactor(property_package=m.fs.properties, inert_species=["c1"])

        for e in m.fs.unit.lagrange_set:
            assert e in ["H", "He"]

        assert isinstance(m.fs.unit.inert_species_balance, Constraint)
        assert len(m.fs.unit.inert_species_balance) == 2
        assert m.fs.unit.inert_species_balance[0, "p1", "c1"] != Constraint.Skip
        assert m.fs.unit.inert_species_balance[0, "p2", "c1"] != Constraint.Skip

        assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
        assert len(m.fs.unit.gibbs_minimization) == 2

    @pytest.mark.unit
    def test_inerts_dependent_w_single_phase(self):
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

        for e in m.fs.unit.lagrange_set:
            assert e in ["H", "He"]

        assert isinstance(m.fs.unit.inert_species_balance, Constraint)
        assert len(m.fs.unit.inert_species_balance) == 0
        assert (0, "p1", "c1") not in m.fs.unit.inert_species_balance

        assert isinstance(m.fs.unit.gibbs_minimization, Constraint)
        assert len(m.fs.unit.gibbs_minimization) == 1

    @pytest.mark.unit
    def test_invalid_inert(self):
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

    @pytest.mark.integration
    def test_natural_gas_w_inerts(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        ng_config = get_prop(
            components=[
                "H2",
                "CO",
                "H2O",
                "CO2",
                "CH4",
                "C2H6",
                "C3H8",
                "C4H10",
                "N2",
                "O2",
                "Ar",
            ]
        )
        m.fs.ng_props = GenericParameterBlock(**ng_config)
        m.fs.unit = GibbsReactor(
            has_heat_transfer=True,
            has_pressure_change=True,
            inert_species=["N2", "Ar"],
            property_package=m.fs.ng_props,
        )

        for e in m.fs.unit.lagrange_set:
            assert e in ["H", "C", "O"]


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

        # Fix some bounds to avoid potential log(0)
        # TODO: This really should be fixed in the property package, but breaks other tests
        m.fs.unit.control_volume.properties_out[0].pressure.setlb(1000)
        m.fs.unit.control_volume.properties_out[0].mole_frac_phase_comp.setlb(1e-12)

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
    def test_structural_issues(self, methane):
        dt = DiagnosticsToolbox(methane)
        dt.assert_no_structural_warnings()

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
        scaler = AutoScaler()
        scaler.scale_variables_by_magnitude(methane)
        scaler.scale_constraints_by_jacobian_norm(methane)

        results = solver.solve(methane)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_verify_scaling_temperature(self, methane):
        assert jacobian_cond(methane, scaled=False) == pytest.approx(5.703e17, rel=1e-3)
        assert jacobian_cond(methane, scaled=True) == pytest.approx(2511, abs=1)

        assert len(extreme_jacobian_rows(methane, scaled=True)) == 0
        assert len(extreme_jacobian_columns(methane, scaled=True)) == 0

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
                    (
                        methane.fs.unit.inlet.flow_mol[0]
                        * methane.fs.unit.control_volume.properties_in[
                            0
                        ].enth_mol_phase["Vap"]
                        - methane.fs.unit.outlet.flow_mol[0]
                        * methane.fs.unit.control_volume.properties_out[
                            0
                        ].enth_mol_phase["Vap"]
                        + methane.fs.unit.heat_duty[0]
                    )
                    / methane.fs.unit.heat_duty[0]
                )
            )
            <= 1e-6
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, methane):
        scaling = TransformationFactory("core.scale_model")
        sm = scaling.create_using(methane, rename=False)

        dt = DiagnosticsToolbox(sm)
        dt.assert_no_numerical_warnings()

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
        scaler = AutoScaler(overwrite=True)
        scaler.scale_variables_by_magnitude(methane)
        scaler.scale_constraints_by_jacobian_norm(methane)

        results = solver.solve(methane)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_verify_scaling_duty(self, methane):
        assert jacobian_cond(methane, scaled=True) == pytest.approx(9191, abs=1)

        assert len(extreme_jacobian_rows(methane, scaled=True)) == 0
        assert len(extreme_jacobian_columns(methane, scaled=True)) == 0

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution_duty(self, methane):
        methane.fs.unit.outlet.display()
        assert pytest.approx(250.06, abs=1e-1) == value(
            methane.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.097367, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.022591, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.10301, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.17691, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.59989, abs=1e-4) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            methane.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.00024128, abs=1e-4) == value(
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


# TODO: Replace once scaling deployed to property package
class PropertyScaler(CustomScalerBase):
    def variable_scaling_routine(self, model, overwrite):
        pass

    def constraint_scaling_routine(self, model, overwrite):
        for c in model.component_data_objects(ctype=Constraint, descend_into=True):
            self.scale_constraint_by_nominal_value(
                c, scheme="inverse_sum", overwrite=overwrite
            )


class TestInitializers:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = MethaneCombustionParameterBlock()

        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_mol[0].set_value(230.0)
        m.fs.unit.inlet.mole_frac_comp[0, "H2"].set_value(0.0435)
        m.fs.unit.inlet.mole_frac_comp[0, "N2"].set_value(0.6522)
        m.fs.unit.inlet.mole_frac_comp[0, "O2"].set_value(0.1739)
        m.fs.unit.inlet.mole_frac_comp[0, "CO2"].set_value(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "CH4"].set_value(0.1304)
        m.fs.unit.inlet.mole_frac_comp[0, "CO"].set_value(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "H2O"].set_value(1e-5)
        m.fs.unit.inlet.mole_frac_comp[0, "NH3"].set_value(1e-5)
        m.fs.unit.inlet.temperature[0].set_value(1500.0)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.outlet.temperature[0].fix(2844.38)
        m.fs.unit.deltaP.fix(0)

        # Apply scaling - Best guesses for unit model and outlet state conditions
        set_scaling_factor(m.fs.unit.control_volume.heat[0.0], 1e-6)

        set_scaling_factor(m.fs.unit.control_volume.properties_out[0.0].flow_mol, 1e-2)
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].flow_mol_phase, 1e-2
        )  # Only 1 phase, so we "know" this
        # N2 is inert, so will be order 0.1, assume CH4 and H2 are near-totally consumed, assume most O2 consumed
        # Assume moderate amounts of CO2 and H2O, small amounts of CO, trace NH3 NH3
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2"], 1e4
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["N2"], 10
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["O2"], 1e2
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO2"], 10
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CH4"], 1e4
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO"], 1e3
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2O"], 10
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["NH3"], 1e4
        )
        set_scaling_factor(
            m.fs.unit.control_volume.properties_out[0.0].temperature, 1e-3
        )
        set_scaling_factor(m.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5)

        scaler = GibbsReactorScaler()
        scaler.scale_model(
            m.fs.unit,
            submodel_scalers={
                "control_volume.properties_in": PropertyScaler,
                "control_volume.properties_out": PropertyScaler,
            },
        )

        return m

    @pytest.mark.component
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer(
            writer_config={"scale_model": True}
        )
        initializer.initialize(
            model.fs.unit,
        )

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(250.06, abs=1e-2) == value(
            model.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.0974, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.0226, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.1030, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.1769, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.5999, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.0002, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(model.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )

        assert not model.fs.unit.inlet.flow_mol[0].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "N2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "O2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CH4"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "NH3"].fixed
        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed

    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(
            constraint_tolerance=2e-5,
            block_solver_writer_config={"linear_presolve": False, "scale_model": True},
        )
        initializer.initialize(
            model.fs.unit,
            initial_guesses={
                "control_volume.properties_out[0].pressure": 101325.0,
                "control_volume.properties_out[0].flow_mol": 251.05,
                "control_volume.properties_out[0].mole_frac_comp[CH4]": 1e-5,
                "control_volume.properties_out[0].mole_frac_comp[CO]": 0.0916,
                "control_volume.properties_out[0].mole_frac_comp[CO2]": 0.0281,
                "control_volume.properties_out[0].mole_frac_comp[H2]": 0.1155,
                "control_volume.properties_out[0].mole_frac_comp[H2O]": 0.1633,
                "control_volume.properties_out[0].mole_frac_comp[N2]": 0.59478,
                "control_volume.properties_out[0].mole_frac_comp[NH3]": 1e-5,
                "control_volume.properties_out[0].mole_frac_comp[O2]": 0.0067,
            },
        )

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(250.06, abs=1e-2) == value(
            model.fs.unit.outlet.flow_mol[0]
        )
        assert pytest.approx(0.0, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CH4"]
        )
        assert pytest.approx(0.0974, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO"]
        )
        assert pytest.approx(0.0226, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "CO2"]
        )
        assert pytest.approx(0.1030, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2"]
        )
        assert pytest.approx(0.1769, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "H2O"]
        )
        assert pytest.approx(0.5999, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "N2"]
        )
        assert pytest.approx(0.0, abs=1e-5) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "NH3"]
        )
        assert pytest.approx(0.0002, abs=1e-4) == value(
            model.fs.unit.outlet.mole_frac_comp[0, "O2"]
        )
        assert pytest.approx(-7454077, abs=1e2) == value(model.fs.unit.heat_duty[0])
        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )

        assert not model.fs.unit.inlet.flow_mol[0].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "N2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "O2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CH4"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "CO"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.mole_frac_comp[0, "NH3"].fixed
        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed
