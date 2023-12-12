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
from math import isnan
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
    SolidPhase,
    Component,
)
from idaes.models.unit_models.solid_liquid import Thickener0D
from idaes.core.solvers import get_solver
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@declare_process_block_class("SolidParameterBlock")
class SolidParameterData(PhysicalParameterBlock):
    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = SolidStateBlock

        # Add Phase objects
        self.sol = SolidPhase()

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


@declare_process_block_class("SolidStateBlock", block_class=StateBlock)
class SolidStateBlockData(StateBlockData):
    def build(self):
        super().build()

        # Create state variables
        self.flow_mass = Var(
            initialize=1.0,
            units=units.kg / units.s,
        )
        self.flow_vol = Var(
            initialize=1.0,
            units=units.m**3 / units.s,
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
            self.sum_mass_frac_eqn = Constraint(
                expr=sum(self.mass_frac_comp[j] for j in self.component_list) == 1
            )

        self.volumetric_flow = Constraint(
            expr=self.flow_vol * 2400 * units.kg * units.m**-3 == self.flow_mass
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


@declare_process_block_class("LiquidParameterBlock")
class LiquidParameterData(PhysicalParameterBlock):
    def build(self):
        """
        Callable method for Block construction.
        """
        super().build()

        self._state_block_class = LiquidStateBlock

        # Add Phase objects
        self.liq = LiquidPhase()

        # Add Component objects
        self.c = Component()
        self.d = Component()

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


@declare_process_block_class("LiquidStateBlock", block_class=StateBlock)
class LiquidStateBlockData(StateBlockData):
    def build(self):
        super().build()

        # Create state variables
        self.flow_mass = Var(
            initialize=1.0,
            units=units.kg / units.s,
        )
        self.flow_vol = Var(
            initialize=1.0,
            units=units.m**3 / units.s,
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
            self.sum_mass_frac_eqn = Constraint(
                expr=sum(self.mass_frac_comp[j] for j in self.component_list) == 1
            )

        self.volumetric_flow = Constraint(
            expr=self.flow_vol * 1000 * units.kg * units.m**-3 == self.flow_mass
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
class TestThickener0DBasic:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.solid = SolidParameterBlock()
        m.fs.liquid = LiquidParameterBlock()

        m.fs.unit = Thickener0D(
            solid_property_package=m.fs.solid,
            liquid_property_package=m.fs.liquid,
        )

        m.fs.unit.solid_inlet.flow_mass.fix(7e-6 * 2400)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "a"].fix(0.2)
        m.fs.unit.solid_inlet.mass_frac_comp[0, "b"].fix(0.8)
        m.fs.unit.solid_inlet.temperature.fix(303.15)
        m.fs.unit.solid_inlet.pressure.fix(101325.0)

        m.fs.unit.liquid_inlet.flow_mass.fix(0.3 / 0.7 * 7e-6 * 1000)
        m.fs.unit.liquid_inlet.mass_frac_comp[0, "c"].fix(0.3)
        m.fs.unit.liquid_inlet.mass_frac_comp[0, "d"].fix(0.7)
        m.fs.unit.liquid_inlet.temperature.fix(310)
        m.fs.unit.liquid_inlet.pressure.fix(1.5e5)

        m.fs.unit.solid_underflow.flow_mass.fix(5e-6 * 2400)

        # Parameters
        m.fs.unit.area.fix(1.2)
        m.fs.unit.v0.fix(1.18e-4)
        m.fs.unit.v1.fix(1e-5)
        m.fs.unit.C.fix(5)
        m.fs.unit.solid_fraction_max.fix(0.9)

        return m

    @pytest.mark.component
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
        dt.report_structural_issues()
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
        assert value(model.fs.unit.solid_fraction_feed[0]) == pytest.approx(
            0.7, rel=1e-5
        )
        assert value(model.fs.unit.solid_fraction_underflow[0]) == pytest.approx(
            0.801065, rel=1e-5
        )
        assert value(model.fs.unit.solid_fraction_overflow[0]) == pytest.approx(
            0.447338, rel=1e-5
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Area": model.fs.unit.area,
                "Liquid Recovery": model.fs.unit.liquid_split.split_fraction[
                    0, "overflow"
                ],
                "Feed Solid Fraction": model.fs.unit.solid_fraction_feed[0],
                "Underflow Solid Fraction": model.fs.unit.solid_fraction_underflow[0],
            },
            "params": {
                "v0": model.fs.unit.v0,
                "v1": model.fs.unit.v1,
                "C": model.fs.unit.C,
                "solid_fraction_max": model.fs.unit.solid_fraction_max,
            },
        }

    @pytest.mark.ui
    @pytest.mark.component
    def test_get_stream_table_contents(self, model):
        stable = model.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "flow_mass": getattr(units.pint_registry, "kg/second"),
                "mass_frac_comp a": getattr(units.pint_registry, "dimensionless"),
                "mass_frac_comp b": getattr(units.pint_registry, "dimensionless"),
                "temperature": getattr(units.pint_registry, "K"),
                "pressure": getattr(units.pint_registry, "Pa"),
                "mass_frac_comp c": getattr(units.pint_registry, "dimensionless"),
                "mass_frac_comp d": getattr(units.pint_registry, "dimensionless"),
            },
            "Feed Solid": {
                "flow_mass": pytest.approx(0.0168, rel=1e-4),
                "mass_frac_comp a": pytest.approx(0.2, rel=1e-4),
                "mass_frac_comp b": pytest.approx(0.8, rel=1e-4),
                "temperature": pytest.approx(303.15, rel=1e-4),
                "pressure": pytest.approx(101325, rel=1e-4),
                "mass_frac_comp c": float("nan"),
                "mass_frac_comp d": float("nan"),
            },
            "Feed Liquid": {
                "flow_mass": pytest.approx(0.003, rel=1e-4),
                "mass_frac_comp a": float("nan"),
                "mass_frac_comp b": float("nan"),
                "temperature": pytest.approx(310, rel=1e-4),
                "pressure": pytest.approx(1.5e5, rel=1e-4),
                "mass_frac_comp c": pytest.approx(0.3, rel=1e-4),
                "mass_frac_comp d": pytest.approx(0.7, rel=1e-4),
            },
            "Underflow Solid": {
                "flow_mass": pytest.approx(0.012, rel=1e-4),
                "mass_frac_comp a": pytest.approx(0.2, rel=1e-4),
                "mass_frac_comp b": pytest.approx(0.8, rel=1e-4),
                "temperature": pytest.approx(303.15, rel=1e-4),
                "pressure": pytest.approx(101325, rel=1e-4),
                "mass_frac_comp c": float("nan"),
                "mass_frac_comp d": float("nan"),
            },
            "Underflow Liquid": {
                "flow_mass": pytest.approx(0.0012417, rel=1e-4),
                "mass_frac_comp a": float("nan"),
                "mass_frac_comp b": float("nan"),
                "temperature": pytest.approx(310, rel=1e-4),
                "pressure": pytest.approx(1.5e5, rel=1e-4),
                "mass_frac_comp c": pytest.approx(0.3, rel=1e-4),
                "mass_frac_comp d": pytest.approx(0.7, rel=1e-4),
            },
            "Overflow Solid": {
                "flow_mass": pytest.approx(0.0048, rel=1e-4),
                "mass_frac_comp a": pytest.approx(0.2, rel=1e-4),
                "mass_frac_comp b": pytest.approx(0.8, rel=1e-4),
                "temperature": pytest.approx(303.15, rel=1e-4),
                "pressure": pytest.approx(101325, rel=1e-4),
                "mass_frac_comp c": float("nan"),
                "mass_frac_comp d": float("nan"),
            },
            "Overflow Liquid": {
                "flow_mass": pytest.approx(0.0017583, rel=1e-4),
                "mass_frac_comp a": float("nan"),
                "mass_frac_comp b": float("nan"),
                "temperature": pytest.approx(310, rel=1e-4),
                "pressure": pytest.approx(1.5e5, rel=1e-4),
                "mass_frac_comp c": pytest.approx(0.3, rel=1e-4),
                "mass_frac_comp d": pytest.approx(0.7, rel=1e-4),
            },
        }

        stb = stable.to_dict()
        for k, d in stb.items():
            if k == "Units":
                assert d == expected["Units"]
            else:
                for i, dd in d.items():
                    if not isnan(dd):
                        assert dd == expected[k][i]
                    else:
                        if "Liquid" in k:
                            assert i in ["mass_frac_comp a", "mass_frac_comp b"]
                        else:
                            assert i in ["mass_frac_comp c", "mass_frac_comp d"]

    @pytest.mark.unit
    def test_deprecate_initialize(self, model):
        with pytest.raises(
            NotImplementedError,
            match="The Thickener0D unit model does not support the old initialization API. "
            "Please use the new API \(InitializerObjects\) instead.",
        ):
            model.fs.unit.initialize()
