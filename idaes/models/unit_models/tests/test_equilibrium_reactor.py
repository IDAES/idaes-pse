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
Tests for equilibrium reactor unit model.
Authors: Andrew Lee
"""

from math import exp
import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ComponentMap,
    ConcreteModel,
    Suffix,
    TransformationFactory,
    units,
    value,
)

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models.equilibrium_reactor import (
    EquilibriumReactor,
    EquilibriumReactorScaler,
    EquilibriumReactorScalerLegacy,
)
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.scaling import (
    get_jacobian,
    jacobian_cond,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
    initialization_tester,
)
from idaes.core.solvers import get_solver
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.scaling import set_scaling_factor

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver("ipopt_v2")


# -----------------------------------------------------------------------------
@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.reactions = ReactionParameterTestBlock(property_package=m.fs.properties)

    m.fs.unit = EquilibriumReactor(
        property_package=m.fs.properties, reaction_package=m.fs.reactions
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 15

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert m.fs.unit.config.has_rate_reactions
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.has_equilibrium_reactions
    assert not m.fs.unit.config.has_phase_equilibrium
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.reaction_package is m.fs.reactions

    assert m.fs.unit.default_initializer is SingleControlVolumeUnitInitializer
    assert m.fs.unit.default_scaler is EquilibriumReactorScaler


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = EquilibriumReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol.fix(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(1e-8)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(1e-8)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, sapon):
        assert hasattr(sapon.fs.unit, "inlet")
        assert len(sapon.fs.unit.inlet.vars) == 4
        assert hasattr(sapon.fs.unit.inlet, "flow_vol")
        assert hasattr(sapon.fs.unit.inlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.inlet, "temperature")
        assert hasattr(sapon.fs.unit.inlet, "pressure")

        assert hasattr(sapon.fs.unit, "outlet")
        assert len(sapon.fs.unit.outlet.vars) == 4
        assert hasattr(sapon.fs.unit.outlet, "flow_vol")
        assert hasattr(sapon.fs.unit.outlet, "conc_mol_comp")
        assert hasattr(sapon.fs.unit.outlet, "temperature")
        assert hasattr(sapon.fs.unit.outlet, "pressure")

        assert hasattr(sapon.fs.unit, "rate_reaction_constraint")
        assert hasattr(sapon.fs.unit, "heat_duty")
        assert hasattr(sapon.fs.unit, "deltaP")

        assert number_variables(sapon) == 26
        assert number_total_constraints(sapon) == 16
        assert number_unused_variables(sapon) == 0

    @pytest.mark.component
    def test_structural_issues(self, sapon):
        dt = DiagnosticsToolbox(sapon)
        dt.assert_no_structural_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, sapon):
        initialization_tester(sapon)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, sapon):
        results = solver.solve(sapon)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(101325.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(304.32, abs=1e-2) == value(
            sapon.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(0.01, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0] - sapon.fs.unit.outlet.flow_vol[0]
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.inlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                    - sapon.fs.unit.outlet.flow_vol[0]
                    * sum(
                        sapon.fs.unit.outlet.conc_mol_comp[0, j]
                        for j in sapon.fs.properties.component_list
                    )
                )
            )
            <= 1e-6
        )
        assert (
            abs(
                value(
                    sapon.fs.unit.inlet.flow_vol[0]
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.inlet.temperature[0]
                        - sapon.fs.properties.temperature_ref
                    )
                    - sapon.fs.unit.outlet.flow_vol[0]
                    * sapon.fs.properties.cp_mol
                    * (
                        sapon.fs.unit.outlet.temperature[0]
                        - sapon.fs.properties.temperature_ref
                    )
                    + sapon.fs.unit.heat_duty[0]
                )
            )
            <= 1e-1
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, sapon):
        dt = DiagnosticsToolbox(sapon)
        dt.assert_no_numerical_warnings()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, sapon):
        perf_dict = sapon.fs.unit._get_performance_contents()

        assert perf_dict == {
            "vars": {
                "Heat Duty": sapon.fs.unit.heat_duty[0],
                "Pressure Change": sapon.fs.unit.deltaP[0],
            }
        }


class TestInitializers:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = EquilibriumReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol[0].set_value(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].set_value(0.0)

        m.fs.unit.inlet.temperature[0].set_value(303.15)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.component
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer()
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert value(model.fs.unit.outlet.flow_vol[0]) == pytest.approx(1e-3, rel=1e-6)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "H2O"]) == pytest.approx(
            55388, rel=1e-5
        )
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "NaOH"]) == pytest.approx(
            0.01611, abs=1e-5
        )
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        ) == pytest.approx(0.01611, abs=1e-5)
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]
        ) == pytest.approx(99.984, rel=1e-5)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]) == pytest.approx(
            99.984, rel=1e-5
        )
        assert value(model.fs.unit.outlet.temperature[0]) == pytest.approx(
            304.32, rel=1e-2
        )
        assert value(model.fs.unit.outlet.pressure[0]) == pytest.approx(
            101325, rel=1e-6
        )

        assert not model.fs.unit.inlet.flow_vol[0].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fixed

        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed

    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(
            constraint_tolerance=2e-5,
            block_solver_writer_config={"linear_presolve": False},
        )
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert value(model.fs.unit.outlet.flow_vol[0]) == pytest.approx(1e-3, rel=1e-6)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "H2O"]) == pytest.approx(
            55388, rel=1e-5
        )
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "NaOH"]) == pytest.approx(
            0.0026102, rel=2e-5
        )
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        ) == pytest.approx(0.0026102, rel=2e-5)
        assert value(
            model.fs.unit.outlet.conc_mol_comp[0, "SodiumAcetate"]
        ) == pytest.approx(99.997, rel=2e-5)
        assert value(model.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]) == pytest.approx(
            99.997, rel=2e-5
        )
        assert value(model.fs.unit.outlet.temperature[0]) == pytest.approx(
            304.32, rel=1e-2
        )
        assert value(model.fs.unit.outlet.pressure[0]) == pytest.approx(
            101325, rel=1e-6
        )

        assert not model.fs.unit.inlet.flow_vol[0].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fixed

        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed


class DummyScaler:
    def __init__(self, **kwargs):
        pass

    def variable_scaling_routine(self, model, **kwargs):
        model._dummy_var_scaler = True

    def constraint_scaling_routine(self, model, **kwargs):
        model._dummy_con_scaler = True


class TestEquilibriumReactorScalerLegacy:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = EquilibriumReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol[0].set_value(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].set_value(0.0)

        m.fs.unit.inlet.temperature[0].set_value(303.15)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.component
    def test_variable_scaling_routine(self, model):
        scaler = EquilibriumReactorScalerLegacy()

        scaler.variable_scaling_routine(model.fs.unit)

        # Inlet state
        sfx_in = model.fs.unit.control_volume.properties_in[0].scaling_factor
        assert isinstance(sfx_in, Suffix)
        assert len(sfx_in) == 8
        assert sfx_in[
            model.fs.unit.control_volume.properties_in[0].flow_vol
        ] == pytest.approx(1e2, rel=1e-8)
        assert sfx_in[
            model.fs.unit.control_volume.properties_in[0].pressure
        ] == pytest.approx(1e-5, rel=1e-8)
        assert sfx_in[
            model.fs.unit.control_volume.properties_in[0].temperature
        ] == pytest.approx(1 / 310.65, rel=1e-8)
        for k, v in model.fs.unit.control_volume.properties_in[0].conc_mol_comp.items():
            if k == "H2O":
                assert sfx_in[v] == pytest.approx(1e-4, rel=1e-8)
            else:
                assert sfx_in[v] == pytest.approx(1e-2, rel=1e-8)

        # Outlet state - should be the same as the inlet
        sfx_out = model.fs.unit.control_volume.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        assert len(sfx_out) == 8
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0].flow_vol
        ] == pytest.approx(1e2, rel=1e-8)
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0].pressure
        ] == pytest.approx(1e-5, rel=1e-8)
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0].temperature
        ] == pytest.approx(1 / 310.65, rel=1e-8)
        for k, v in model.fs.unit.control_volume.properties_out[
            0
        ].conc_mol_comp.items():
            if k == "H2O":
                assert sfx_out[v] == pytest.approx(1e-4, rel=1e-8)
            else:
                assert sfx_out[v] == pytest.approx(1e-2, rel=1e-8)

        # Reaction block
        sfx_rxn = model.fs.unit.control_volume.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        assert len(sfx_rxn) == 2
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0].k_rxn
        ] == pytest.approx(
            1 / (3.132e6 * exp(-43000 / (8.31446262 * 310.65))), rel=1e-8
        )
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0].reaction_rate["R1"]
        ] == pytest.approx(1e2, rel=1e-8)

        # Check that unit model has scaling factors
        sfx_cv = model.fs.unit.control_volume.scaling_factor
        assert isinstance(sfx_cv, Suffix)
        assert len(sfx_cv) == 2
        assert sfx_cv[model.fs.unit.control_volume.heat[0]] == pytest.approx(
            1.917448e-05, rel=1e-3
        )
        assert sfx_cv[model.fs.unit.control_volume.deltaP[0]] == pytest.approx(
            1e-4, rel=1e-3
        )

        # No unit level variables to scale, so no suffix
        assert not hasattr(model.fs.unit, "scaling_factor")

    @pytest.mark.component
    def test_variable_scaling_routine_submodel_scaler(self, model):
        scaler = EquilibriumReactorScalerLegacy()

        scaler_map = ComponentMap()
        scaler_map[model.fs.unit.control_volume.properties_in] = DummyScaler
        scaler_map[model.fs.unit.control_volume.properties_out] = DummyScaler
        scaler_map[model.fs.unit.control_volume.reactions] = DummyScaler

        scaler.variable_scaling_routine(
            model.fs.unit,
            submodel_scalers=scaler_map,
        )

        # Should call DummyScaler submethod for each submodel
        # Should add _dummy_var_scaler = True to all submodels
        assert model.fs.unit.control_volume.properties_in[0]._dummy_var_scaler
        assert model.fs.unit.control_volume.properties_out[0]._dummy_var_scaler
        assert model.fs.unit.control_volume.reactions[0]._dummy_var_scaler

    @pytest.mark.component
    def test_constraint_scaling_routine(self, model):
        scaler = EquilibriumReactorScalerLegacy()

        scaler.constraint_scaling_routine(model.fs.unit)

        # Check that sub-models have suffixes - we will assume they are right at this point
        # No constraints on the inlet properties, so no scaling suffix generated
        # sfx_in = model.fs.unit.control_volume.properties_in[0].scaling_factor

        sfx_out = model.fs.unit.control_volume.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        assert len(sfx_out) == 1
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0.0].conc_water_eqn
        ] == pytest.approx(1e-4, rel=1e-8)

        sfx_rxn = model.fs.unit.control_volume.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        assert len(sfx_rxn) == 2
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0.0].arrhenius_eqn
        ] == pytest.approx(1, rel=1e-8)
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0.0].rate_expression["R1"]
        ] == pytest.approx(1e-4, rel=1e-8)

        # Check that unit model has scaling factors
        sfx_cv = model.fs.unit.control_volume.scaling_factor
        assert isinstance(model.fs.unit.control_volume.scaling_factor, Suffix)
        assert len(sfx_cv) == 12
        assert sfx_cv[
            model.fs.unit.control_volume.enthalpy_balances[0.0]
        ] == pytest.approx(8.03894082e-10, rel=1e-8)
        assert sfx_cv[
            model.fs.unit.control_volume.pressure_balance[0.0]
        ] == pytest.approx(9.86923267e-6, rel=1e-8)
        for c in model.fs.unit.control_volume.material_balances.values():
            assert sfx_cv[c] == pytest.approx(1e-2, rel=1e-8)
        for (
            c
        ) in (
            model.fs.unit.control_volume.rate_reaction_stoichiometry_constraint.values()
        ):
            assert sfx_cv[c] == pytest.approx(1, rel=1e-8)

        sfx_unit = model.fs.unit.scaling_factor
        assert isinstance(sfx_unit, Suffix)
        assert len(sfx_unit) == 1
        assert sfx_unit[
            model.fs.unit.rate_reaction_constraint[0.0, "R1"]
        ] == pytest.approx(1, rel=1e-8)

    @pytest.mark.component
    def test_constraint_scaling_routine_submodel_scaler(self, model):
        scaler = EquilibriumReactorScalerLegacy()

        scaler_map = ComponentMap()
        scaler_map[model.fs.unit.control_volume.properties_in] = DummyScaler
        scaler_map[model.fs.unit.control_volume.properties_out] = DummyScaler
        scaler_map[model.fs.unit.control_volume.reactions] = DummyScaler

        scaler.constraint_scaling_routine(
            model.fs.unit,
            submodel_scalers=scaler_map,
        )

        # Should call DummyScaler submethod for each submodel
        # Should add _dummy_con_scaler = True to all submodels
        assert model.fs.unit.control_volume.properties_in[0]._dummy_con_scaler
        assert model.fs.unit.control_volume.properties_out[0]._dummy_con_scaler
        assert model.fs.unit.control_volume.reactions[0]._dummy_con_scaler

    @pytest.mark.component
    def test_scale_model(self, model):
        scaler = EquilibriumReactorScalerLegacy()

        scaler.scale_model(model.fs.unit)

        # Inlet state
        sfx_in = model.fs.unit.control_volume.properties_in[0].scaling_factor
        assert isinstance(sfx_in, Suffix)
        assert len(sfx_in) == 8
        assert sfx_in[
            model.fs.unit.control_volume.properties_in[0].flow_vol
        ] == pytest.approx(1e2, rel=1e-8)
        assert sfx_in[
            model.fs.unit.control_volume.properties_in[0].pressure
        ] == pytest.approx(1e-5, rel=1e-8)
        assert sfx_in[
            model.fs.unit.control_volume.properties_in[0].temperature
        ] == pytest.approx(1 / 310.65, rel=1e-8)
        for k, v in model.fs.unit.control_volume.properties_in[0].conc_mol_comp.items():
            if k == "H2O":
                assert sfx_in[v] == pytest.approx(1e-4, rel=1e-8)
            else:
                assert sfx_in[v] == pytest.approx(1e-2, rel=1e-8)

        # Outlet state - should be the same as the inlet
        sfx_out = model.fs.unit.control_volume.properties_out[0].scaling_factor
        assert isinstance(sfx_out, Suffix)
        assert len(sfx_out) == 9
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0].flow_vol
        ] == pytest.approx(1e2, rel=1e-8)
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0].pressure
        ] == pytest.approx(1e-5, rel=1e-8)
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0].temperature
        ] == pytest.approx(1 / 310.65, rel=1e-8)
        for k, v in model.fs.unit.control_volume.properties_out[
            0
        ].conc_mol_comp.items():
            if k == "H2O":
                assert sfx_out[v] == pytest.approx(1e-4, rel=1e-8)
            else:
                assert sfx_out[v] == pytest.approx(1e-2, rel=1e-8)
        assert sfx_out[
            model.fs.unit.control_volume.properties_out[0.0].conc_water_eqn
        ] == pytest.approx(1e-4, rel=1e-8)

        # Reaction block
        sfx_rxn = model.fs.unit.control_volume.reactions[0].scaling_factor
        assert isinstance(sfx_rxn, Suffix)
        assert len(sfx_rxn) == 4
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0].k_rxn
        ] == pytest.approx(
            1 / (3.132e6 * exp(-43000 / (8.31446262 * 310.65))), rel=1e-8
        )
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0].reaction_rate["R1"]
        ] == pytest.approx(1e2, rel=1e-8)
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0.0].arrhenius_eqn
        ] == pytest.approx(5.4240896, rel=1e-8)
        assert sfx_rxn[
            model.fs.unit.control_volume.reactions[0.0].rate_expression["R1"]
        ] == pytest.approx(5.4240896e-4, rel=1e-8)

        # Check that unit model has scaling factors
        sfx_cv = model.fs.unit.control_volume.scaling_factor
        assert isinstance(sfx_cv, Suffix)
        assert len(sfx_cv) == 14
        assert sfx_cv[model.fs.unit.control_volume.heat[0]] == pytest.approx(
            1.917448e-05, rel=1e-3
        )
        assert sfx_cv[model.fs.unit.control_volume.deltaP[0]] == pytest.approx(
            1e-4, rel=1e-3
        )
        assert sfx_cv[
            model.fs.unit.control_volume.enthalpy_balances[0.0]
        ] == pytest.approx(7.71546823e-08, rel=1e-8)
        assert sfx_cv[
            model.fs.unit.control_volume.pressure_balance[0.0]
        ] == pytest.approx(1e-5, rel=1e-8)
        for (_, _, j), c in model.fs.unit.control_volume.material_balances.items():
            if j == "H2O":
                assert sfx_cv[c] == pytest.approx(1e-2, rel=1e-8)
            else:
                assert sfx_cv[c] == pytest.approx(1, rel=1e-8)
        for (
            _,
            _,
            j,
        ), c in (
            model.fs.unit.control_volume.rate_reaction_stoichiometry_constraint.items()
        ):
            assert sfx_cv[c] == pytest.approx(1, rel=1e-8)

        sfx_unit = model.fs.unit.scaling_factor
        assert isinstance(sfx_unit, Suffix)
        assert len(sfx_unit) == 1
        assert sfx_unit[
            model.fs.unit.rate_reaction_constraint[0.0, "R1"]
        ] == pytest.approx(1e2, rel=1e-8)

    @pytest.mark.integration
    def test_example_case(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.equil = EquilibriumReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_of_reaction=True,
        )

        m.fs.equil.inlet.flow_vol.fix(1.0e-03)
        m.fs.equil.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.equil.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.equil.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.equil.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(1e-8)
        m.fs.equil.inlet.conc_mol_comp[0, "Ethanol"].fix(1e-8)

        m.fs.equil.inlet.temperature.fix(303.15)
        m.fs.equil.inlet.pressure.fix(101325.0)

        initializer = BlockTriangularizationInitializer()
        initializer.initialize(m.fs.equil)

        set_scaling_factor(m.fs.equil.control_volume.properties_in[0].flow_vol, 1)

        scaler = EquilibriumReactorScalerLegacy()
        scaler.scale_model(m.fs.equil)

        m.fs.equil.inlet.flow_vol.fix(1)
        m.fs.equil.inlet.conc_mol_comp[0, "NaOH"].fix(200.0)
        m.fs.equil.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.equil.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(50)
        m.fs.equil.inlet.conc_mol_comp[0, "Ethanol"].fix(1e-8)

        m.fs.equil.inlet.temperature.fix(320)

        solver = get_solver(
            "ipopt_v2", writer_config={"linear_presolve": True, "scale_model": True}
        )
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            1.030e4, rel=1e-3
        )


class TestEquilibriumReactorScaler:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = EquilibriumReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol[0].set_value(1.0e-03)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].set_value(0.0)

        m.fs.unit.inlet.temperature[0].set_value(303.15)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.integration
    def test_example_case(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.equil = EquilibriumReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_equilibrium_reactions=False,
            has_heat_of_reaction=True,
        )

        m.fs.equil.inlet.flow_vol.fix(1.0e-03)
        m.fs.equil.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.equil.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.equil.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.equil.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(1e-8)
        m.fs.equil.inlet.conc_mol_comp[0, "Ethanol"].fix(1e-8)

        m.fs.equil.inlet.temperature.fix(303.15)
        m.fs.equil.inlet.pressure.fix(101325.0)

        initializer = BlockTriangularizationInitializer()
        initializer.initialize(m.fs.equil)

        set_scaling_factor(m.fs.equil.control_volume.properties_in[0].flow_vol, 1)

        scaler = EquilibriumReactorScaler()
        scaler.scale_model(m.fs.equil)

        m.fs.equil.inlet.flow_vol.fix(1)
        m.fs.equil.inlet.conc_mol_comp[0, "NaOH"].fix(200.0)
        m.fs.equil.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.equil.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(50)
        m.fs.equil.inlet.conc_mol_comp[0, "Ethanol"].fix(1e-8)

        m.fs.equil.inlet.temperature.fix(320)

        solver = get_solver(
            "ipopt_v2", writer_config={"linear_presolve": True, "scale_model": True}
        )
        results = solver.solve(m, tee=True)
        assert_optimal_termination(results)

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(m, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(218.88, rel=1e-3)
