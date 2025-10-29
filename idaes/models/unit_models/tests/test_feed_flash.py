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
Tests for feed with flash.
Authors: Andrew Lee, Daison Caballero
"""

import pytest

from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    ComponentMap,
    TransformationFactory,
    value,
    units as pyunits,
)

from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models.feed_flash import FeedFlash, FlashType, FeedFlashScaler
from idaes.models.properties import iapws95
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.properties.modular_properties import GenericParameterBlock
from idaes.models.properties.modular_properties.examples.BT_ideal import (
    configuration as BTIdeal_config,
)
from idaes.core.util.model_statistics import (
    number_variables,
    number_total_constraints,
    number_unused_variables,
)
from idaes.core.util.testing import PhysicalParameterTestBlock, initialization_tester
from idaes.core.solvers import get_solver

from idaes.core.initialization import (
    SingleControlVolumeUnitInitializer,
    BlockTriangularizationInitializer,
    InitializationStatus,
)
from idaes.core.util import DiagnosticsToolbox
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

    m.fs.unit = FeedFlash(property_package=m.fs.properties)

    # Check unit config arguments
    assert len(m.fs.unit.config) == 6

    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.flash_type == FlashType.isothermal
    assert m.fs.unit.config.property_package is m.fs.properties

    assert m.fs.unit.default_initializer is SingleControlVolumeUnitInitializer


@pytest.mark.unit
def test_scaler_object():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = PhysicalParameterTestBlock()
    m.fs.unit = FeedFlash(property_package=m.fs.properties)

    assert m.fs.unit.default_scaler is FeedFlashScaler

    scaler_obj = m.fs.unit.default_scaler()
    scaler_obj.scale_model(m.fs.unit)

    assert m.fs.unit.control_volume.properties_in[0].variables_scaled
    assert m.fs.unit.control_volume.properties_in[0].constraints_scaled

    assert m.fs.unit.control_volume.properties_out[0].variables_scaled
    assert m.fs.unit.control_volume.properties_out[0].constraints_scaled


# -----------------------------------------------------------------------------
class TestBTXIdeal(object):
    @pytest.fixture(scope="class")
    def btx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
        )

        m.fs.unit = FeedFlash(property_package=m.fs.properties)

        m.fs.unit.flow_mol.fix(1)
        m.fs.unit.temperature.fix(368)
        m.fs.unit.pressure.fix(101325)
        m.fs.unit.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.mole_frac_comp[0, "toluene"].fix(0.5)

        # Legacy property package, does not bound many variables which triggers
        # warnings for potential evaluation errors.
        # Fixing property package is out of scope for now.
        m.fs.unit.control_volume.properties_in[0.0].temperature_bubble.setlb(300)
        m.fs.unit.control_volume.properties_in[0.0].temperature_bubble.setub(550)
        m.fs.unit.control_volume.properties_in[0.0].temperature_dew.setlb(300)
        m.fs.unit.control_volume.properties_in[0.0].temperature_dew.setub(550)
        m.fs.unit.control_volume.properties_in[0.0]._temperature_equilibrium.setlb(300)
        m.fs.unit.control_volume.properties_in[0.0]._temperature_equilibrium.setub(550)
        m.fs.unit.control_volume.properties_in[0.0].pressure_sat_comp.setlb(1e4)
        m.fs.unit.control_volume.properties_in[0.0].pressure_sat_comp.setub(5e6)

        m.fs.unit.control_volume.properties_out[0.0].temperature_bubble.setlb(300)
        m.fs.unit.control_volume.properties_out[0.0].temperature_bubble.setub(550)
        m.fs.unit.control_volume.properties_out[0.0].temperature_dew.setlb(300)
        m.fs.unit.control_volume.properties_out[0.0].temperature_dew.setub(550)
        m.fs.unit.control_volume.properties_out[0.0]._temperature_equilibrium.setlb(300)
        m.fs.unit.control_volume.properties_out[0.0]._temperature_equilibrium.setub(550)
        m.fs.unit.control_volume.properties_out[0.0].pressure_sat_comp.setlb(1e4)
        m.fs.unit.control_volume.properties_out[0.0].pressure_sat_comp.setub(5e6)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, btx):
        assert hasattr(btx.fs.unit, "flow_mol")
        assert hasattr(btx.fs.unit, "mole_frac_comp")
        assert hasattr(btx.fs.unit, "temperature")
        assert hasattr(btx.fs.unit, "pressure")

        assert hasattr(btx.fs.unit, "outlet")
        assert len(btx.fs.unit.outlet.vars) == 4
        assert hasattr(btx.fs.unit.outlet, "flow_mol")
        assert hasattr(btx.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(btx.fs.unit.outlet, "temperature")
        assert hasattr(btx.fs.unit.outlet, "pressure")

        assert hasattr(btx.fs.unit, "isothermal")

        assert number_variables(btx) == 34
        assert number_total_constraints(btx) == 29
        assert number_unused_variables(btx) == 0

    @pytest.mark.component
    def test_structural_issues(self, btx):
        dt = DiagnosticsToolbox(btx)
        dt.assert_no_structural_warnings()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, btx):
        perf_dict = btx.fs.unit._get_performance_contents()

        assert perf_dict is None

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, btx):
        stable = btx.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "flow_mol": getattr(pyunits.pint_registry, "mole/second"),
                "mole_frac_comp benzene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "mole_frac_comp toluene": getattr(
                    pyunits.pint_registry, "dimensionless"
                ),
                "temperature": getattr(pyunits.pint_registry, "kelvin"),
                "pressure": getattr(pyunits.pint_registry, "Pa"),
            },
            "Outlet": {
                "flow_mol": pytest.approx(1.0, rel=1e-4),
                "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
                "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
                "temperature": pytest.approx(298.15, rel=1e-4),
                "pressure": pytest.approx(101325.0, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, btx):
        initialization_tester(btx)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, btx):
        results = solver.solve(btx)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, btx):
        assert pytest.approx(101325.0, abs=1e3) == value(btx.fs.unit.outlet.pressure[0])
        assert pytest.approx(368.00, abs=1e-0) == value(
            btx.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(1.0, abs=1e-2) == value(btx.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(0.396, abs=1e-3) == value(
            btx.fs.unit.control_volume.properties_out[0].flow_mol_phase["Vap"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, btx):
        dt = DiagnosticsToolbox(btx)
        dt.assert_no_numerical_warnings()


# -----------------------------------------------------------------------------
class TestBTIdealModular(object):
    @pytest.fixture(scope="class")
    def bt_modular(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = GenericParameterBlock(**BTIdeal_config)

        m.fs.unit = FeedFlash(property_package=m.fs.properties)

        m.fs.unit.flow_mol.fix(1.0)
        m.fs.unit.temperature.fix(368.0)
        m.fs.unit.pressure.fix(101325.0)
        m.fs.unit.mole_frac_comp[0, "benzene"].fix(0.5)
        m.fs.unit.mole_frac_comp[0, "toluene"].fix(0.5)

        # Legacy property package, does not bound many variables which triggers
        # warnings for potential evaluation errors.
        # Fixing property package is out of scope for now.
        m.fs.unit.control_volume.properties_in[0.0]._teq.setlb(300)
        m.fs.unit.control_volume.properties_in[0.0]._teq.setub(550)

        m.fs.unit.control_volume.properties_out[0.0]._teq.setlb(300)
        m.fs.unit.control_volume.properties_out[0.0]._teq.setub(550)

        return m

    @pytest.mark.build
    @pytest.mark.unit
    def test_build(self, bt_modular):
        assert hasattr(bt_modular.fs.unit, "flow_mol")
        assert hasattr(bt_modular.fs.unit, "mole_frac_comp")
        assert hasattr(bt_modular.fs.unit, "temperature")
        assert hasattr(bt_modular.fs.unit, "pressure")

        assert hasattr(bt_modular.fs.unit, "outlet")
        assert len(bt_modular.fs.unit.outlet.vars) == 4
        assert hasattr(bt_modular.fs.unit.outlet, "flow_mol")
        assert hasattr(bt_modular.fs.unit.outlet, "mole_frac_comp")
        assert hasattr(bt_modular.fs.unit.outlet, "temperature")
        assert hasattr(bt_modular.fs.unit.outlet, "pressure")

        assert hasattr(bt_modular.fs.unit, "isothermal")

        assert number_variables(bt_modular.fs.unit) == 42
        assert number_total_constraints(bt_modular) == 37
        # assert number_unused_variables(bt_modular) == 0

    @pytest.mark.component
    def test_structural_issues(self, bt_modular):
        dt = DiagnosticsToolbox(bt_modular)
        dt.assert_no_structural_warnings()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, bt_modular):
        perf_dict = bt_modular.fs.unit._get_performance_contents()

        assert perf_dict is None

    # TODO the formatting for modular properties is broken (see #1684).
    # This test can be fixed once that issue is fixed
    # @pytest.mark.ui
    # @pytest.mark.unit
    # def test_get_stream_table_contents(self, bt_modular):
    #     stable = bt_modular.fs.unit._get_stream_table_contents()

    #     expected = {
    #         "Units": {
    #             "flow_mol": getattr(pyunits.pint_registry, "mole/second"),
    #             "mole_frac_comp benzene": getattr(
    #                 pyunits.pint_registry, "dimensionless"
    #             ),
    #             "mole_frac_comp toluene": getattr(
    #                 pyunits.pint_registry, "dimensionless"
    #             ),
    #             "temperature": getattr(pyunits.pint_registry, "K"),
    #             "pressure": getattr(pyunits.pint_registry, "Pa"),
    #         },
    #         "Outlet": {
    #             "flow_mol": pytest.approx(1.0, rel=1e-4),
    #             "mole_frac_comp benzene": pytest.approx(0.5, rel=1e-4),
    #             "mole_frac_comp toluene": pytest.approx(0.5, rel=1e-4),
    #             "temperature": pytest.approx(368.0, rel=1e-4),
    #             "pressure": pytest.approx(101325.0, rel=1e-4),
    #         },
    #     }

    #     assert stable.to_dict() == expected

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, bt_modular):
        initialization_tester(bt_modular)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, bt_modular):
        results = solver.solve(bt_modular)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, bt_modular):
        assert pytest.approx(101325.0, abs=1e-3) == value(
            bt_modular.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(368.0, abs=1e-3) == value(
            bt_modular.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(1.0, abs=1e-3) == value(
            bt_modular.fs.unit.outlet.flow_mol[0]
        )

        assert pytest.approx(0.603, abs=1e-3) == value(
            bt_modular.fs.unit.control_volume.properties_out[0].flow_mol_phase["Liq"]
        )
        assert pytest.approx(0.396, abs=1e-3) == value(
            bt_modular.fs.unit.control_volume.properties_out[0].flow_mol_phase["Vap"]
        )

        assert pytest.approx(0.412, abs=1e-3) == value(
            bt_modular.fs.unit.control_volume.properties_out[0].mole_frac_phase_comp[
                "Liq", "benzene"
            ]
        )
        assert pytest.approx(0.588, abs=1e-3) == value(
            bt_modular.fs.unit.control_volume.properties_out[0].mole_frac_phase_comp[
                "Liq", "toluene"
            ]
        )
        assert pytest.approx(0.634, abs=1e-3) == value(
            bt_modular.fs.unit.control_volume.properties_out[0].mole_frac_phase_comp[
                "Vap", "benzene"
            ]
        )
        assert pytest.approx(0.366, abs=1e-3) == value(
            bt_modular.fs.unit.control_volume.properties_out[0].mole_frac_phase_comp[
                "Vap", "toluene"
            ]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, bt_modular):
        dt = DiagnosticsToolbox(bt_modular)
        dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_scaler_object(self, bt_modular):
        jac, _ = get_jacobian(bt_modular, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            4.5495232e7, rel=1e-3
        )

        property_scaler = bt_modular.fs.unit.control_volume.properties_in[
            0
        ].default_scaler()
        property_scaler.default_scaling_factors["flow_mol_phase"] = 1

        submodel_scalers = ComponentMap()
        submodel_scalers[bt_modular.fs.unit.control_volume.properties_in] = (
            property_scaler
        )
        submodel_scalers[bt_modular.fs.unit.control_volume.properties_out] = (
            property_scaler
        )

        scaler_object = bt_modular.fs.unit.default_scaler()
        scaler_object.scale_model(bt_modular.fs.unit, submodel_scalers=submodel_scalers)

        sm = TransformationFactory("core.scale_model").create_using(
            bt_modular, rename=False
        )
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            7.54763e4, rel=1e-3
        )


# -----------------------------------------------------------------------------
@pytest.mark.iapws
@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
class TestIAPWS(object):
    @pytest.fixture(scope="class")
    def iapws(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock(
            phase_presentation=iapws95.PhaseType.LG
        )

        m.fs.unit = FeedFlash(
            property_package=m.fs.properties, flash_type=FlashType.isenthalpic
        )

        m.fs.unit.flow_mol.fix(100)
        m.fs.unit.enth_mol.fix(24000)
        m.fs.unit.pressure.fix(101325)

        return m

    @pytest.mark.unit
    def test_build(self, iapws):
        assert hasattr(iapws.fs.unit, "flow_mol")
        assert hasattr(iapws.fs.unit, "enth_mol")
        assert hasattr(iapws.fs.unit, "pressure")

        assert hasattr(iapws.fs.unit, "outlet")
        assert len(iapws.fs.unit.outlet.vars) == 3
        assert hasattr(iapws.fs.unit.outlet, "flow_mol")
        assert hasattr(iapws.fs.unit.outlet, "enth_mol")
        assert hasattr(iapws.fs.unit.outlet, "pressure")

        assert hasattr(iapws.fs.unit, "isenthalpic")

        assert number_variables(iapws) == 6
        assert number_total_constraints(iapws) == 3
        assert number_unused_variables(iapws) == 0

    @pytest.mark.component
    def test_structural_issues(self, iapws):
        dt = DiagnosticsToolbox(iapws)
        dt.assert_no_structural_warnings()

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_performance_contents(self, iapws):
        perf_dict = iapws.fs.unit._get_performance_contents()

        assert perf_dict is None

    @pytest.mark.ui
    @pytest.mark.unit
    def test_get_stream_table_contents(self, iapws):
        stable = iapws.fs.unit._get_stream_table_contents()

        expected = {
            "Units": {
                "Molar Flow": getattr(pyunits.pint_registry, "mole/second"),
                "Mass Flow": getattr(pyunits.pint_registry, "kg/second"),
                "T": getattr(pyunits.pint_registry, "K"),
                "P": getattr(pyunits.pint_registry, "Pa"),
                "Vapor Fraction": getattr(pyunits.pint_registry, "dimensionless"),
                "Molar Enthalpy": getattr(pyunits.pint_registry, "J/mole"),
            },
            "Outlet": {
                "Molar Flow": pytest.approx(1, rel=1e-4),
                "Mass Flow": pytest.approx(1.8015e-2, rel=1e-4),
                "T": pytest.approx(270.4877112932641, rel=1e-4),
                "P": pytest.approx(11032305.8275, rel=1e-4),
                "Vapor Fraction": pytest.approx(0, abs=1e-4),
                "Molar Enthalpy": pytest.approx(0.01102138712926277, rel=1e-4),
            },
        }

        assert stable.to_dict() == expected

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, iapws):
        initialization_tester(iapws)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, iapws):
        results = solver.solve(iapws)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, iapws):
        assert pytest.approx(101325.0, abs=1e3) == value(
            iapws.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(24000, abs=1e3) == value(iapws.fs.unit.outlet.enth_mol[0])
        assert pytest.approx(100.0, abs=1e-2) == value(iapws.fs.unit.outlet.flow_mol[0])

        assert pytest.approx(373.12, abs=1e-2) == value(
            iapws.fs.unit.control_volume.properties_out[0].temperature
        )
        assert pytest.approx(0.5953, abs=1e-4) == value(
            iapws.fs.unit.control_volume.properties_out[0].phase_frac["Liq"]
        )

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_numerical_issues(self, iapws):
        dt = DiagnosticsToolbox(iapws)
        dt.assert_no_numerical_warnings()

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_scaler_object(self, iapws):
        jac, _ = get_jacobian(iapws, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            5.76015e6, rel=1e-3
        )

        props_scaler = iapws.fs.unit.control_volume.properties_in.default_scaler()
        props_scaler.default_scaling_factors["flow_mol"] = 1e-2
        submodel_scalers = ComponentMap()
        submodel_scalers[iapws.fs.unit.control_volume.properties_in] = props_scaler
        submodel_scalers[iapws.fs.unit.control_volume.properties_out] = props_scaler

        scaler_object = iapws.fs.unit.default_scaler()
        scaler_object.scale_model(iapws.fs.unit, submodel_scalers=submodel_scalers)

        sm = TransformationFactory("core.scale_model").create_using(iapws, rename=False)
        jac, _ = get_jacobian(sm, scaled=False)
        assert (jacobian_cond(jac=jac, scaled=False)) == pytest.approx(
            72.58140, rel=1e-3
        )


class TestInitializersIAWPS:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = iapws95.Iapws95ParameterBlock(
            phase_presentation=iapws95.PhaseType.LG
        )

        m.fs.unit = FeedFlash(
            property_package=m.fs.properties, flash_type=FlashType.isenthalpic
        )

        m.fs.unit.flow_mol.fix(100)
        m.fs.unit.enth_mol.fix(24000)
        m.fs.unit.pressure.fix(101325)

        return m

    @pytest.mark.integration
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer()
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(101325.0, abs=1e3) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(24000, abs=1e3) == value(model.fs.unit.outlet.enth_mol[0])
        assert pytest.approx(100.0, abs=1e-2) == value(model.fs.unit.outlet.flow_mol[0])

        assert pytest.approx(373.12, abs=1e-2) == value(
            model.fs.unit.control_volume.properties_out[0].temperature
        )
        assert pytest.approx(0.5953, abs=1e-4) == value(
            model.fs.unit.control_volume.properties_out[0].phase_frac["Liq"]
        )

    @pytest.mark.integration
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(101325.0, abs=1e3) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(24000, abs=1e3) == value(model.fs.unit.outlet.enth_mol[0])
        assert pytest.approx(100.0, abs=1e-2) == value(model.fs.unit.outlet.flow_mol[0])

        assert pytest.approx(373.12, abs=1e-2) == value(
            model.fs.unit.control_volume.properties_out[0].temperature
        )
        assert pytest.approx(0.5953, abs=1e-4) == value(
            model.fs.unit.control_volume.properties_out[0].phase_frac["Liq"]
        )


class TestInitializersBT:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal"
        )

        m.fs.unit = FeedFlash(property_package=m.fs.properties)

        m.fs.unit.flow_mol[0].set_value(1)
        m.fs.unit.temperature[0].set_value(368)
        m.fs.unit.pressure[0].set_value(101325)
        m.fs.unit.mole_frac_comp[0, "benzene"].set_value(0.5)
        m.fs.unit.mole_frac_comp[0, "toluene"].set_value(0.5)

        return m

    @pytest.mark.component
    def test_general_hierarchical(self, model):
        initializer = SingleControlVolumeUnitInitializer()
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(101325.0, abs=1e3) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(368.00, abs=1e-0) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(1.0, abs=1e-2) == value(model.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(0.396, abs=1e-3) == value(
            model.fs.unit.control_volume.properties_out[0].flow_mol_phase["Vap"]
        )

        assert not model.fs.unit.flow_mol[0].fixed
        assert not model.fs.unit.temperature[0].fixed
        assert not model.fs.unit.pressure[0].fixed
        assert not model.fs.unit.mole_frac_comp[0, "benzene"].fixed
        assert not model.fs.unit.mole_frac_comp[0, "toluene"].fixed

    @pytest.mark.component
    def test_block_triangularization(self, model):
        initializer = BlockTriangularizationInitializer(constraint_tolerance=2e-5)
        initializer.initialize(model.fs.unit)

        assert initializer.summary[model.fs.unit]["status"] == InitializationStatus.Ok

        assert pytest.approx(101325.0, abs=1e3) == value(
            model.fs.unit.outlet.pressure[0]
        )
        assert pytest.approx(368.00, abs=1e-0) == value(
            model.fs.unit.outlet.temperature[0]
        )
        assert pytest.approx(1.0, abs=1e-2) == value(model.fs.unit.outlet.flow_mol[0])
        assert pytest.approx(0.396, abs=1e-3) == value(
            model.fs.unit.control_volume.properties_out[0].flow_mol_phase["Vap"]
        )

        assert not model.fs.unit.flow_mol[0].fixed
        assert not model.fs.unit.temperature[0].fixed
        assert not model.fs.unit.pressure[0].fixed
        assert not model.fs.unit.mole_frac_comp[0, "benzene"].fixed
        assert not model.fs.unit.mole_frac_comp[0, "toluene"].fixed
