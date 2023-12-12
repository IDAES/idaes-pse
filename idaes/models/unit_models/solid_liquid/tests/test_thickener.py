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
Tests for thickener unit model.
Authors: Andrew Lee
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    units,
    value,
    Var,
)
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent, assert_units_equivalent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
    declare_process_block_class,
    MaterialFlowBasis,
    PhysicalParameterBlock,
    StateBlockData,
    StateBlock,
    LiquidPhase,
    Component,
)
from idaes.models.unit_models.solid_liquid import Thickener0D
from idaes.models.unit_models.separator import (
    Separator,
    SeparatorData,
    EnergySplittingType,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.solvers import get_solver
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@declare_process_block_class("TestParameterBlock")
class TestParameterData(PhysicalParameterBlock):
    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = TestStateBlock

        # Add Phase objects
        self.Liq = LiquidPhase()

        # Add Component objects
        self.a = Component()
        self.b = Component()

    @classmethod
    def define_metadata(cls, obj):
        obj.add_default_units(
            {
                "time": units.s,
                "length": units.m,
                "mass": units.kg,
                "amount": units.mol,
                "temperature": units.K,
            }
        )


@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class TestStateBlockData(StateBlockData):
    def build(self):
        super().build()

        # Create state variables
        self.flow_mass = Var(
            initialize=1.0,
            units=units.kg / units.s,
        )
        self.pressure = Var(
            initialize=101325.0,
            bounds=(1e3, 1e6),
            units=units.Pa,
        )
        self.temperature = Var(
            initialize=298.15,
            bounds=(298.15, 323.15),
            units=units.K,
        )
        self.mass_frac_comp = Var(
            self.params.component_list,
            initialize=0.1,
            units=units.dimensionless,
        )

        self.dens_mass = Param(
            initialize=1000,
            units=units.kg / units.m**3,
        )

        if self.config.defined_state is False:
            self.summ_mass_frac_eqn = Constraint(
                expr=sum(self.mass_frac_comp[j] for j in self.component_list) == 1
            )

    def get_material_flow_terms(b, p, j):
        return b.flow_mass * b.mass_frac_comp[j]

    def get_enthalpy_flow_terms(b, p):
        return (
            b.flow_mass * 42e3 * units.J / units.kg * (b.temperature - 273.15 * units.K)
        )

    def default_material_balance_type(self):
        return MaterialBalanceType.componentTotal

    def define_state_vars(b):
        return {
            "flow_mass": b.flow_mass,
            "mass_frac_comp": b.mass_frac_comp,
            "temperature": b.temperature,
            "pressure": b.pressure,
        }

    def get_material_flow_basis(b):
        return MaterialFlowBasis.mass


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_beta_logger(caplog):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = TestParameterBlock()

    m.fs.unit = Thickener0D(
        solid_property_package=m.fs.properties,
        liquid_property_package=m.fs.properties,
    )
    expected = (
        "The Thickener0D model is currently a beta capability and will "
        "likely change in the next release as a more predictive version is "
        "developed."
    )

    assert expected in caplog.text


class TestThickener0DBasic:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = TestParameterBlock()

        m.fs.unit = Thickener0D(
            solid_property_package=m.fs.properties,
            liquid_property_package=m.fs.properties,
        )

        m.fs.unit.solid_inlet.flow_mass.fix(1.33)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "a"].fix(0.2)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "b"].fix(0.8)
        m.fs.unit.solid_inlet.temperature.fix(303.15)
        m.fs.unit.solid_inlet.pressure.fix(101325.0)

        m.fs.unit.liquid_inlet.flow_mass.fix(6.65)
        m.fs.unit.liquid_inlet.mass_frac_comp[0, "a"].fix(0.4)
        m.fs.unit.liquid_inlet.mass_frac_comp[0, "b"].fix(0.6)
        m.fs.unit.liquid_inlet.temperature.fix(320)
        m.fs.unit.liquid_inlet.pressure.fix(2e5)

        # Operating conditions based on Example 5.1 (pg 199)
        # Coulson & Richardson's Chemical Engineering Vol 2 (4th Ed.)
        m.fs.unit.liquid_recovery.fix(0.7)

        # Thickener specific parameters
        m.fs.unit.height_clarification.fix(1.5)
        m.fs.unit.settling_velocity_pinch.fix(0.94e-4)
        m.fs.unit.liquid_solid_pinch.fix(3.7)
        m.fs.unit.settling_time.fix(200)

        return m

    @pytest.mark.unit
    def test_config(self, model):
        # Check unit config arguments
        assert len(model.fs.unit.config) == 9

        assert not model.fs.unit.config.dynamic
        assert not model.fs.unit.config.has_holdup
        assert (
            model.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
        )
        assert (
            model.fs.unit.config.energy_split_basis
            == EnergySplittingType.equal_temperature
        )
        assert (
            model.fs.unit.config.momentum_balance_type
            == MomentumBalanceType.pressureTotal
        )

        assert model.fs.unit.config.solid_property_package is model.fs.properties
        assert model.fs.unit.config.liquid_property_package is model.fs.properties

        assert model.fs.unit.default_initializer is BlockTriangularizationInitializer

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, model):

        assert isinstance(model.fs.unit.solid_inlet, Port)
        assert len(model.fs.unit.solid_inlet.vars) == 4
        assert hasattr(model.fs.unit.solid_inlet, "flow_mass")
        assert hasattr(model.fs.unit.solid_inlet, "mass_frac_comp")
        assert hasattr(model.fs.unit.solid_inlet, "temperature")
        assert hasattr(model.fs.unit.solid_inlet, "pressure")

        assert isinstance(model.fs.unit.solid_outlet, Port)
        assert len(model.fs.unit.solid_outlet.vars) == 4
        assert hasattr(model.fs.unit.solid_outlet, "flow_mass")
        assert hasattr(model.fs.unit.solid_outlet, "mass_frac_comp")
        assert hasattr(model.fs.unit.solid_outlet, "temperature")
        assert hasattr(model.fs.unit.solid_outlet, "pressure")

        assert isinstance(model.fs.unit.liquid_inlet, Port)
        assert len(model.fs.unit.liquid_inlet.vars) == 4
        assert hasattr(model.fs.unit.liquid_inlet, "flow_mass")
        assert hasattr(model.fs.unit.liquid_inlet, "mass_frac_comp")
        assert hasattr(model.fs.unit.liquid_inlet, "temperature")
        assert hasattr(model.fs.unit.liquid_inlet, "pressure")

        assert isinstance(model.fs.unit.retained_liquid_outlet, Port)
        assert len(model.fs.unit.retained_liquid_outlet.vars) == 4
        assert hasattr(model.fs.unit.retained_liquid_outlet, "flow_mass")
        assert hasattr(model.fs.unit.retained_liquid_outlet, "mass_frac_comp")
        assert hasattr(model.fs.unit.retained_liquid_outlet, "temperature")
        assert hasattr(model.fs.unit.retained_liquid_outlet, "pressure")

        assert isinstance(model.fs.unit.recovered_liquid_outlet, Port)
        assert len(model.fs.unit.recovered_liquid_outlet.vars) == 4
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "flow_mass")
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "mass_frac_comp")
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "temperature")
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "pressure")

        assert isinstance(model.fs.unit.split, SeparatorData)
        assert isinstance(model.fs.unit.liquid_recovery, Var)

        assert isinstance(model.fs.unit.area, Var)
        assert isinstance(model.fs.unit.height, Var)
        assert isinstance(model.fs.unit.height_clarification, Var)
        assert isinstance(model.fs.unit.settling_velocity_pinch, Var)
        assert isinstance(model.fs.unit.liquid_solid_pinch, Var)
        assert isinstance(model.fs.unit.liquid_solid_underflow, Var)
        assert isinstance(model.fs.unit.settling_time, Var)

        assert isinstance(model.fs.unit.underflow_sl_constraint, Constraint)
        assert isinstance(model.fs.unit.cross_sectional_area_constraint, Constraint)
        assert isinstance(model.fs.unit.height_constraint, Constraint)

        assert number_variables(model) == 29
        assert number_total_constraints(model) == 14
        # These are the solid properties, as they do not appear in constraints
        assert number_unused_variables(model) == 4

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

        assert_units_equivalent(model.fs.unit.height_clarification, units.m)
        assert_units_equivalent(
            model.fs.unit.settling_velocity_pinch, units.m / units.s
        )
        assert_units_equivalent(model.fs.unit.settling_time, units.s)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.report_structural_issues()
        dt.display_underconstrained_set()
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model, tee=True)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # Solid outlet
        assert pytest.approx(101325.0, rel=1e-8) == value(
            model.fs.unit.solid_outlet.pressure[0]
        )
        assert pytest.approx(303.15, rel=1e-8) == value(
            model.fs.unit.solid_outlet.temperature[0]
        )
        assert pytest.approx(0.2, rel=1e-8) == value(
            model.fs.unit.solid_outlet.mass_frac_comp[0, "a"]
        )
        assert pytest.approx(0.8, rel=1e-8) == value(
            model.fs.unit.solid_outlet.mass_frac_comp[0, "b"]
        )

        assert pytest.approx(1.33, rel=1e-8) == value(
            model.fs.unit.solid_outlet.flow_mass[0]
        )

        # Retained liquid
        assert pytest.approx(2e5, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.pressure[0]
        )
        assert pytest.approx(320, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.temperature[0]
        )
        assert pytest.approx(0.4, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.mass_frac_comp[0, "a"]
        )
        assert pytest.approx(0.6, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.mass_frac_comp[0, "b"]
        )
        assert pytest.approx(6.65 * 0.3, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.flow_mass[0]
        )

        # Recovered liquid
        assert pytest.approx(2e5, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.pressure[0]
        )
        assert pytest.approx(320, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.temperature[0]
        )
        assert pytest.approx(0.4, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.mass_frac_comp[0, "a"]
        )
        assert pytest.approx(0.6, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.mass_frac_comp[0, "b"]
        )
        assert pytest.approx(6.65 * 0.7, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.flow_mass[0]
        )

        # Thickener Vars
        assert pytest.approx(1.5, rel=1e-6) == value(
            model.fs.unit.liquid_solid_underflow[0]
        )
        assert pytest.approx(31.12766, rel=1e-6) == value(model.fs.unit.area)
        model.fs.unit.height.display()
        assert pytest.approx(1.536318, rel=1e-6) == value(model.fs.unit.height)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Area": model.fs.unit.area,
                "Height": model.fs.unit.height,
                "Liquid Recovery": model.fs.unit.liquid_recovery[0],
                "Underflow L/S": model.fs.unit.liquid_solid_underflow[0],
                "Pinch L/S": model.fs.unit.liquid_solid_pinch[0],
                "Critical Settling Velocity": model.fs.unit.settling_velocity_pinch[0],
                "Settling Time": model.fs.unit.settling_time[0],
            }
        }
