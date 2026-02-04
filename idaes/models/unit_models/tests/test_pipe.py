#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for Pipe.

Author: Douglas Allan
"""
import re
import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    TransformationFactory,
    value,
    Var,
)
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
)
from idaes.models.unit_models.pipe import Pipe, PipeScaler

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
    number_derivative_variables,
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
from idaes.core.scaling.util import jacobian_cond


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver("ipopt_v2")


# -----------------------------------------------------------------------------
@pytest.mark.component
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = Pipe(property_package=m.fs.properties)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 15

    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.energy_balance_type == EnergyBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_heat_transfer
    assert not m.fs.unit.config.has_pressure_change
    assert not m.fs.unit.config.has_phase_equilibrium
    assert m.fs.unit.config.property_package is m.fs.properties

    # Make sure that we got rid of all the PFR options
    assert not hasattr(m.fs.unit.config, "reaction_package")
    assert not hasattr(m.fs.unit.config, "reaction_package_args")
    assert not hasattr(m.fs.unit.config, "has_equilibrium_reactions")
    assert not hasattr(m.fs.unit.config, "has_heat_of_reaction")

    assert m.fs.unit.config.length_domain_set == [0.0, 1.0]
    assert m.fs.unit.config.transformation_method == "dae.finite_difference"
    assert m.fs.unit.config.transformation_scheme == "BACKWARD"
    assert m.fs.unit.config.finite_elements == 20
    assert m.fs.unit.config.collocation_points == 3

    assert m.fs.unit.default_scaler is PipeScaler


@pytest.mark.component
def test_no_rate_reactions():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.reactions = ReactionParameterTestBlock(property_package=m.fs.properties)

    with pytest.raises(
        ValueError,
        match=re.escape(
            "key 'reaction_package' not defined for ConfigDict '' and implicit (undefined) keys are not allowed"
        ),
    ):
        m.fs.unit = Pipe(
            property_package=m.fs.properties, reaction_package=m.fs.reactions
        )


# -----------------------------------------------------------------------------
class TestSaponification(object):
    @pytest.fixture(scope="class")
    def sapon(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = Pipe(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol.fix(1.0)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(1e-8)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fix(1e-8)

        m.fs.unit.inlet.temperature.fix(303.15)
        m.fs.unit.inlet.pressure.fix(101325.0)

        m.fs.unit.length.fix(0.5)
        m.fs.unit.area.fix(0.1)

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

        assert isinstance(sapon.fs.unit.area, Var)
        assert isinstance(sapon.fs.unit.length, Var)
        assert isinstance(sapon.fs.unit.volume, Var)
        assert not hasattr(sapon.fs.unit, "performance_eqn")
        assert hasattr(sapon.fs.unit.control_volume, "heat")
        assert hasattr(sapon.fs.unit, "heat_duty")
        assert hasattr(sapon.fs.unit, "deltaP")

        assert number_variables(sapon) == 486
        assert number_total_constraints(sapon) == 427
        assert number_unused_variables(sapon) == 9
        assert number_derivative_variables(sapon) == 0

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
        assert (
            pytest.approx(101325.0, abs=1e-2) == sapon.fs.unit.outlet.pressure[0].value
        )
        assert (
            pytest.approx(303.15, abs=1e-2) == sapon.fs.unit.outlet.temperature[0].value
        )
        assert pytest.approx(100, abs=1e-2) == value(
            sapon.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_conservation(self, sapon):
        assert value(sapon.fs.unit.inlet.flow_vol[0]) == pytest.approx(
            value(sapon.fs.unit.outlet.flow_vol[0]), rel=1e-6
        )
        for j in sapon.fs.properties.component_list:
            assert value(
                sapon.fs.unit.inlet.flow_vol[0]
                * sapon.fs.unit.inlet.conc_mol_comp[0, j]
            ) == pytest.approx(
                value(
                    sapon.fs.unit.outlet.flow_vol[0]
                    * sapon.fs.unit.outlet.conc_mol_comp[0, j]
                ),
                rel=1e-6,
                abs=1e-6,
            )
        h_in = value(
            (
                sapon.fs.unit.inlet.flow_vol[0]
                * sapon.fs.properties.dens_mol
                * sapon.fs.properties.cp_mol
                * (
                    sapon.fs.unit.inlet.temperature[0]
                    - sapon.fs.properties.temperature_ref
                )
            )
        )
        h_out = value(
            sapon.fs.unit.outlet.flow_vol[0]
            * sapon.fs.properties.dens_mol
            * sapon.fs.properties.cp_mol
            * (
                sapon.fs.unit.outlet.temperature[0]
                - sapon.fs.properties.temperature_ref
            )
        )
        assert h_in == pytest.approx(h_out, abs=1e-6)

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

        assert perf_dict == {"vars": {"Area": sapon.fs.unit.area}}

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_scaling(self, sapon):
        unit = sapon.fs.unit

        assert jacobian_cond(unit, scaled=False) == pytest.approx(2.0373e9, rel=1e-3)

        scaler_obj = unit.default_scaler()
        scaler_obj.default_scaling_factors["length"] = 2
        scaler_obj.default_scaling_factors["area"] = 10
        scaler_obj.set_variable_scaling_factor(unit.inlet.flow_vol[0], 1)

        scaler_obj.scale_model(unit)

        assert scaler_obj.get_scaling_factor(unit.length) == 2
        assert scaler_obj.get_scaling_factor(unit.area) == 10

        assert jacobian_cond(unit, scaled=True) == pytest.approx(7.62474e3, rel=1e-3)


class TestInitializers:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = SaponificationParameterBlock()

        m.fs.unit = Pipe(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        m.fs.unit.inlet.flow_vol[0].set_value(1.0)
        m.fs.unit.inlet.conc_mol_comp[0, "H2O"].set_value(55388.0)
        m.fs.unit.inlet.conc_mol_comp[0, "NaOH"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].set_value(100.0)
        m.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].set_value(0.0)
        m.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].set_value(0.0)

        m.fs.unit.inlet.temperature[0].set_value(303.15)
        m.fs.unit.inlet.pressure[0].set_value(101325.0)

        m.fs.unit.length.fix(0.5)
        m.fs.unit.area.fix(0.1)

        m.fs.unit.heat_duty.fix(0)
        m.fs.unit.deltaP.fix(0)

        return m

    @pytest.mark.component
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer()
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert (
            pytest.approx(101325.0, rel=1e-5) == model.fs.unit.outlet.pressure[0].value
        )
        assert (
            pytest.approx(303.15, rel=1e-5) == model.fs.unit.outlet.temperature[0].value
        )
        assert pytest.approx(100, rel=1e-5) == value(
            model.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
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
        initializer = BlockTriangularizationInitializer(constraint_tolerance=1e-5)

        # Need to exclude unused variables for x=0 in post-check
        initializer.initialize(model.fs.unit, exclude_unused_vars=True)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert (
            pytest.approx(101325.0, rel=1e-5) == model.fs.unit.outlet.pressure[0].value
        )
        assert (
            pytest.approx(303.15, rel=1e-5) == model.fs.unit.outlet.temperature[0].value
        )
        assert pytest.approx(100, rel=1e-5) == value(
            model.fs.unit.outlet.conc_mol_comp[0, "EthylAcetate"]
        )

        assert not model.fs.unit.inlet.flow_vol[0].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "H2O"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "NaOH"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "EthylAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "SodiumAcetate"].fixed
        assert not model.fs.unit.inlet.conc_mol_comp[0, "Ethanol"].fixed

        assert not model.fs.unit.inlet.temperature[0].fixed
        assert not model.fs.unit.inlet.pressure[0].fixed
