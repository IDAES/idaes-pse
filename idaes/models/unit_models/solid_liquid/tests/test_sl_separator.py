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
Tests for solid-liquid separator unit model.
Authors: Andrew Lee
"""

import pytest

from pyomo.environ import check_optimal_termination, ConcreteModel, units, value, Var
from pyomo.network import Port
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models.solid_liquid import SLSeparator
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
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
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


# -----------------------------------------------------------------------------
class TestSLSeparatorBasic:
    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = SLSeparator(
            solid_property_package=m.fs.properties,
            liquid_property_package=m.fs.properties,
        )

        m.fs.unit.solid_inlet.flow_vol.fix(10)
        m.fs.unit.solid_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.solid_inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.solid_inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.solid_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(1e-6)
        m.fs.unit.solid_inlet.conc_mol_comp[0, "Ethanol"].fix(1e-6)
        m.fs.unit.solid_inlet.temperature.fix(303.15)
        m.fs.unit.solid_inlet.pressure.fix(101325.0)

        m.fs.unit.liquid_inlet.flow_vol.fix(20)
        m.fs.unit.liquid_inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.liquid_inlet.conc_mol_comp[0, "NaOH"].fix(1e-6)
        m.fs.unit.liquid_inlet.conc_mol_comp[0, "EthylAcetate"].fix(1e-6)
        m.fs.unit.liquid_inlet.conc_mol_comp[0, "SodiumAcetate"].fix(50.0)
        m.fs.unit.liquid_inlet.conc_mol_comp[0, "Ethanol"].fix(50.0)
        m.fs.unit.liquid_inlet.temperature.fix(320)
        m.fs.unit.liquid_inlet.pressure.fix(2e5)

        m.fs.unit.liquid_recovery.fix(0.7)

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
        assert hasattr(model.fs.unit.solid_inlet, "flow_vol")
        assert hasattr(model.fs.unit.solid_inlet, "conc_mol_comp")
        assert hasattr(model.fs.unit.solid_inlet, "temperature")
        assert hasattr(model.fs.unit.solid_inlet, "pressure")

        assert isinstance(model.fs.unit.solid_outlet, Port)
        assert len(model.fs.unit.solid_outlet.vars) == 4
        assert hasattr(model.fs.unit.solid_outlet, "flow_vol")
        assert hasattr(model.fs.unit.solid_outlet, "conc_mol_comp")
        assert hasattr(model.fs.unit.solid_outlet, "temperature")
        assert hasattr(model.fs.unit.solid_outlet, "pressure")

        assert isinstance(model.fs.unit.liquid_inlet, Port)
        assert len(model.fs.unit.liquid_inlet.vars) == 4
        assert hasattr(model.fs.unit.liquid_inlet, "flow_vol")
        assert hasattr(model.fs.unit.liquid_inlet, "conc_mol_comp")
        assert hasattr(model.fs.unit.liquid_inlet, "temperature")
        assert hasattr(model.fs.unit.liquid_inlet, "pressure")

        assert isinstance(model.fs.unit.retained_liquid_outlet, Port)
        assert len(model.fs.unit.retained_liquid_outlet.vars) == 4
        assert hasattr(model.fs.unit.retained_liquid_outlet, "flow_vol")
        assert hasattr(model.fs.unit.retained_liquid_outlet, "conc_mol_comp")
        assert hasattr(model.fs.unit.retained_liquid_outlet, "temperature")
        assert hasattr(model.fs.unit.retained_liquid_outlet, "pressure")

        assert isinstance(model.fs.unit.recovered_liquid_outlet, Port)
        assert len(model.fs.unit.recovered_liquid_outlet.vars) == 4
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "flow_vol")
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "conc_mol_comp")
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "temperature")
        assert hasattr(model.fs.unit.recovered_liquid_outlet, "pressure")

        assert isinstance(model.fs.unit.split, SeparatorData)
        assert isinstance(model.fs.unit.liquid_recovery, Var)

        assert number_variables(model) == 34
        assert number_total_constraints(model) == 17
        # These are the solid properties, as they do not appear in constraints
        assert number_unused_variables(model) == 8

    @pytest.mark.component
    def test_units(self, model):
        assert_units_consistent(model)

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model) == 0

    @pytest.mark.component
    def test_structural_issues(self, model):
        dt = DiagnosticsToolbox(model)
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
        assert pytest.approx(100, rel=1e-8) == value(
            model.fs.unit.solid_outlet.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(100, rel=1e-8) == value(
            model.fs.unit.solid_outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(1e-6, rel=1e-8) == value(
            model.fs.unit.solid_outlet.conc_mol_comp[0, "Ethanol"]
        )
        assert pytest.approx(1e-6, rel=1e-8) == value(
            model.fs.unit.solid_outlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(10, rel=1e-8) == value(
            model.fs.unit.solid_outlet.flow_vol[0]
        )

        # Retained liquid
        assert pytest.approx(2e5, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.pressure[0]
        )
        assert pytest.approx(320, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.temperature[0]
        )
        assert pytest.approx(1e-6, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(1e-6, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(50, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.conc_mol_comp[0, "Ethanol"]
        )
        assert pytest.approx(50, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(20 * 0.3, rel=1e-8) == value(
            model.fs.unit.retained_liquid_outlet.flow_vol[0]
        )

        # Recovered liquid
        assert pytest.approx(2e5, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.pressure[0]
        )
        assert pytest.approx(320, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.temperature[0]
        )
        assert pytest.approx(1e-6, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.conc_mol_comp[0, "NaOH"]
        )
        assert pytest.approx(1e-6, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.conc_mol_comp[0, "EthylAcetate"]
        )
        assert pytest.approx(50, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.conc_mol_comp[0, "Ethanol"]
        )
        assert pytest.approx(50, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.conc_mol_comp[0, "SodiumAcetate"]
        )
        assert pytest.approx(20 * 0.7, rel=1e-8) == value(
            model.fs.unit.recovered_liquid_outlet.flow_vol[0]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solid_conservation(self, model):
        assert (
            abs(
                value(
                    model.fs.unit.solid_inlet.flow_vol[0]
                    - model.fs.unit.solid_outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.solid_inlet.flow_vol[0]
                    * sum(
                        model.fs.unit.solid_inlet.conc_mol_comp[0, j]
                        for j in model.fs.properties.component_list
                    )
                    - model.fs.unit.solid_outlet.flow_vol[0]
                    * sum(
                        model.fs.unit.solid_outlet.conc_mol_comp[0, j]
                        for j in model.fs.properties.component_list
                    )
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    (
                        model.fs.unit.solid_inlet.flow_vol[0]
                        * model.fs.properties.dens_mol
                        * model.fs.properties.cp_mol
                        * (
                            model.fs.unit.solid_inlet.temperature[0]
                            - model.fs.properties.temperature_ref
                        )
                    )
                    - (
                        model.fs.unit.solid_outlet.flow_vol[0]
                        * model.fs.properties.dens_mol
                        * model.fs.properties.cp_mol
                        * (
                            model.fs.unit.solid_outlet.temperature[0]
                            - model.fs.properties.temperature_ref
                        )
                    )
                )
            )
            <= 1e-3
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_liquid_conservation(self, model):
        assert (
            abs(
                value(
                    model.fs.unit.liquid_inlet.flow_vol[0]
                    - model.fs.unit.retained_liquid_outlet.flow_vol[0]
                    - model.fs.unit.recovered_liquid_outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    model.fs.unit.liquid_inlet.flow_vol[0]
                    * sum(
                        model.fs.unit.liquid_inlet.conc_mol_comp[0, j]
                        for j in model.fs.properties.component_list
                    )
                    - model.fs.unit.retained_liquid_outlet.flow_vol[0]
                    * sum(
                        model.fs.unit.retained_liquid_outlet.conc_mol_comp[0, j]
                        for j in model.fs.properties.component_list
                    )
                    - model.fs.unit.recovered_liquid_outlet.flow_vol[0]
                    * sum(
                        model.fs.unit.recovered_liquid_outlet.conc_mol_comp[0, j]
                        for j in model.fs.properties.component_list
                    )
                )
            )
            <= 1e-6
        )

        assert (
            abs(
                value(
                    (
                        model.fs.unit.liquid_inlet.flow_vol[0]
                        * model.fs.properties.dens_mol
                        * model.fs.properties.cp_mol
                        * (
                            model.fs.unit.liquid_inlet.temperature[0]
                            - model.fs.properties.temperature_ref
                        )
                    )
                    - (
                        model.fs.unit.retained_liquid_outlet.flow_vol[0]
                        * model.fs.properties.dens_mol
                        * model.fs.properties.cp_mol
                        * (
                            model.fs.unit.retained_liquid_outlet.temperature[0]
                            - model.fs.properties.temperature_ref
                        )
                    )
                    - (
                        model.fs.unit.recovered_liquid_outlet.flow_vol[0]
                        * model.fs.properties.dens_mol
                        * model.fs.properties.cp_mol
                        * (
                            model.fs.unit.recovered_liquid_outlet.temperature[0]
                            - model.fs.properties.temperature_ref
                        )
                    )
                )
            )
            <= 1e-3
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, model):
        perf_dict = model.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Liquid Recovery": model.fs.unit.liquid_recovery[0],
            }
        }

    @pytest.mark.ui
    @pytest.mark.component
    def test_get_stream_table_contents(self, model):
        stable = model.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Volumetric Flowrate": getattr(units.pint_registry, "m**3/second"),
                "Molar Concentration H2O": getattr(units.pint_registry, "mole/m**3"),
                "Molar Concentration NaOH": getattr(units.pint_registry, "mole/m**3"),
                "Molar Concentration EthylAcetate": getattr(
                    units.pint_registry, "mole/m**3"
                ),
                "Molar Concentration SodiumAcetate": getattr(
                    units.pint_registry, "mole/m**3"
                ),
                "Molar Concentration Ethanol": getattr(
                    units.pint_registry, "mole/m**3"
                ),
                "Temperature": getattr(units.pint_registry, "K"),
                "Pressure": getattr(units.pint_registry, "Pa"),
            },
            "Solid Inlet": {
                "Volumetric Flowrate": pytest.approx(10, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(5.5388e4, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(1e-6, abs=1e-4),
                "Molar Concentration Ethanol": pytest.approx(1e-6, rel=1e-4),
                "Temperature": pytest.approx(303.15, rel=1e-4),
                "Pressure": pytest.approx(101325, rel=1e-4),
            },
            "Liquid Inlet": {
                "Volumetric Flowrate": pytest.approx(20, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(5.5388e4, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(1e-6, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(1e-6, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(50, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(50, rel=1e-4),
                "Temperature": pytest.approx(320, rel=1e-4),
                "Pressure": pytest.approx(2e5, rel=1e-4),
            },
            "Solid Outlet": {
                "Volumetric Flowrate": pytest.approx(10, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(5.5388e4, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(100, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(100, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(1e-6, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(1e-6, rel=1e-4),
                "Temperature": pytest.approx(303.15, rel=1e-4),
                "Pressure": pytest.approx(101325, rel=1e-4),
            },
            "Liquid in Solids Outlet": {
                "Volumetric Flowrate": pytest.approx(20 * 0.3, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(5.5388e4, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(1e-6, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(1e-6, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(50, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(50, rel=1e-4),
                "Temperature": pytest.approx(320, rel=1e-4),
                "Pressure": pytest.approx(2e5, rel=1e-4),
            },
            "Recovered Liquid Outlet": {
                "Volumetric Flowrate": pytest.approx(20 * 0.7, rel=1e-4),
                "Molar Concentration H2O": pytest.approx(5.5388e4, rel=1e-4),
                "Molar Concentration NaOH": pytest.approx(1e-6, rel=1e-4),
                "Molar Concentration EthylAcetate": pytest.approx(1e-6, rel=1e-4),
                "Molar Concentration SodiumAcetate": pytest.approx(50, rel=1e-4),
                "Molar Concentration Ethanol": pytest.approx(50, rel=1e-4),
                "Temperature": pytest.approx(320, rel=1e-4),
                "Pressure": pytest.approx(2e5, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected
