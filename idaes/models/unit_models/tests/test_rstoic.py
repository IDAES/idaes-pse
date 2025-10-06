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
Tests for IDAES Stoichiometric reactor.

Author: Chinedu Okoli, Andrew Lee, Ryan Hughes
"""
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    value,
    units,
    Constraint,
    Var,
    TransformationFactory,
)

from pyomo.common.collections import ComponentMap

from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)

from idaes.models.unit_models.stoichiometric_reactor import (
    StoichiometricReactor,
    StoichiometricReactorScaler,
)

from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
)

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import (
    PhysicalParameterTestBlock,
    ReactionParameterTestBlock,
    initialization_tester,
)
from idaes.core.initialization import (
    BlockTriangularizationInitializer,
    SingleControlVolumeUnitInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
from idaes.core.scaling import set_scaling_factor

from idaes.core.util.scaling import (
    get_jacobian,
    jacobian_cond,
)

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

    m.fs.unit = StoichiometricReactor(
        property_package=m.fs.properties, reaction_package=m.fs.reactions
    )

    # Check unit config arguments
    assert len(m.fs.unit.config) == 12

    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.has_heat_of_reaction
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.reaction_package is m.fs.reactions


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

        m.fs.unit = StoichiometricReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol.fix(1)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(1e-8)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(1e-8)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.rate_reaction_extent[0, "R1"].fix(90)
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

        assert hasattr(sapon.fs.unit, "rate_reaction_extent")

        assert number_variables(sapon) == 24
        assert number_total_constraints(sapon) == 13
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
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, sapon):
        assert pytest.approx(101325.0, abs=1e-2) == value(
            sapon.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(304.21, abs=1e-2) == value(
            sapon.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(90, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]
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

        assert pytest.approx(4410000, abs=1e3) == value(
            sapon.fs.unit.control_volume.heat_of_reaction[0]
        )
        assert (
            abs(
                value(
                    (
                        sapon.fs.unit.inlet.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.inlet.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    - (
                        sapon.fs.unit.outlet.flow_vol[0]
                        * sapon.fs.properties.dens_mol
                        * sapon.fs.properties.cp_mol
                        * (
                            sapon.fs.unit.outlet.temperature[0]
                            - sapon.fs.properties.temperature_ref
                        )
                    )
                    + sapon.fs.unit.control_volume.heat_of_reaction[0]
                )
            )
            <= 1e-3
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
                "Reaction Extent [R1]": sapon.fs.unit.rate_reaction_extent[0, "R1"],
                "Heat Duty": sapon.fs.unit.heat_duty[0],
                "Pressure Change": sapon.fs.unit.deltaP[0],
            }
        }


class TestInitializersSapon:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.unit = StoichiometricReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol[0].set_value(1)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].set_value(0.0)

        m.fs.unit.inlet.temperature[0].set_value(303.15)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.rate_reaction_extent[0, "R1"].fix(90)
        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.component
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer()
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(304.21, abs=1e-2) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(90, abs=1e-2) == value(
            model.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]
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
        # Need to limit tolerance on 1x1 solver otherwise it exceeds the iteration limit
        initializer = BlockTriangularizationInitializer(
            constraint_tolerance=2e-5, calculate_variable_options={"eps": 1e-6}
        )
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(101325.0, abs=1e-2) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(304.21, abs=1e-2) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(90, abs=1e-2) == value(
            model.fs.unit.outlet.conc_mol_comp[0, "Ethanol"]
        )

        assert not model.fs.unit.inlet.flow_vol[0].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fixed

        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed


class TestStoichiometricReactorScaler:
    @pytest.mark.unit
    def test_default_scaler(self):
        assert StoichiometricReactor().default_scaler == StoichiometricReactorScaler

    @pytest.fixture
    def sapon_model(self):
        # saponification example
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()
        m.fs.reactions = SaponificationReactionParameterBlock(
            property_package=m.fs.properties
        )

        m.fs.R101 = StoichiometricReactor(
            property_package=m.fs.properties,
            reaction_package=m.fs.reactions,
            has_heat_transfer=True,
            has_heat_of_reaction=True,
            has_pressure_change=True,
        )

        m.fs.R101.inlet.flow_vol[0].fix(1.0e-03 * units.m**3 / units.s)
        m.fs.R101.inlet.conc_mol_comp[0, "H2O"].fix(55388.0 * units.mol / units.m**3)
        m.fs.R101.inlet.conc_mol_comp[0, "NaOH"].fix(100.0 * units.mol / units.m**3)
        m.fs.R101.inlet.conc_mol_comp[0, "EthylAcetate"].fix(
            100.0 * units.mol / units.m**3
        )
        m.fs.R101.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(
            0.0 * units.mol / units.m**3
        )
        m.fs.R101.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0 * units.mol / units.m**3)

        m.fs.R101.inlet.temperature[0].fix(303.15 * units.K)
        m.fs.R101.inlet.pressure[0].fix(101325.0 * units.Pa)

        m.fs.R101.heat_duty.fix(0 * units.W)
        m.fs.R101.deltaP.fix(0 * units.Pa)

        m.fs.R101.conversion = Var(
            initialize=0.90, bounds=(0, 1), units=units.dimensionless
        )

        m.fs.R101.conv_constraint = Constraint(
            expr=m.fs.R101.conversion * m.fs.R101.inlet.conc_mol_comp[0, "NaOH"]
            == (
                m.fs.R101.inlet.conc_mol_comp[0, "NaOH"]
                - m.fs.R101.outlet.conc_mol_comp[0, "NaOH"]
            )
        )

        m.fs.R101.conversion.fix(0.9)

        set_scaling_factor(m.fs.R101.conversion, 10)
        set_scaling_factor(m.fs.R101.conv_constraint, 10)

        submodel_scaler_obj = m.fs.properties.state_block_class.default_scaler()
        submodel_scaler_obj.default_scaling_factors["flow_vol"] = 1e3
        submodel_scaler_obj.default_scaling_factors["conc_mol_comp"] = 0.1
        submodel_scalers = ComponentMap()
        submodel_scalers[m.fs.R101.control_volume.properties_in] = submodel_scaler_obj
        submodel_scalers[m.fs.R101.control_volume.properties_out] = submodel_scaler_obj

        scaler = StoichiometricReactorScaler()
        scaler.scale_model(model=m.fs.R101, submodel_scalers=submodel_scalers)

        return m

    @pytest.mark.component
    def test_sapon_sfs(self, sapon_model):
        assert sapon_model.fs.R101.control_volume.scaling_factor[
            sapon_model.fs.R101.control_volume.heat[0.0]
        ] == pytest.approx(1.917e-5, rel=1e-3)
        assert (
            sapon_model.fs.R101.control_volume.scaling_factor[
                sapon_model.fs.R101.control_volume.rate_reaction_stoichiometry_constraint[
                    0.0, "Liq", "H2O"
                ]
            ]
            == 0.1
        )

        assert (
            sapon_model.fs.R101.control_volume.properties_in[0.0].scaling_factor[
                sapon_model.fs.R101.control_volume.properties_in[0.0].flow_vol
            ]
            == 1e3
        )

        assert (
            sapon_model.fs.R101.control_volume.properties_out[0.0].scaling_factor[
                sapon_model.fs.R101.control_volume.properties_out[0.0].pressure
            ]
            == 1e-5
        )

    @pytest.mark.integration
    def test_sapon_example_case(self, sapon_model):
        # test condition number
        initializer = BlockTriangularizationInitializer()
        initializer.initialize(sapon_model.fs.R101)

        solver = get_solver(
            "ipopt_v2", writer_config={"linear_presolve": True, "scale_model": True}
        )
        results = solver.solve(sapon_model, tee=False)
        assert check_optimal_termination(results)

        # Check condition number to confirm scaling
        sm = TransformationFactory("core.scale_model").create_using(
            sapon_model, rename=False
        )
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(714.68, rel=1e-3)
