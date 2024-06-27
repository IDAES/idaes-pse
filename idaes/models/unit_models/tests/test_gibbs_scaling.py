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
Tests for Gibbs reactor Scaler.

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

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.gibbs_reactor import GibbsReactor, GibbsReactorScaler
from idaes.models.properties.activity_coeff_models.methane_combustion_ideal import (
    MethaneParameterBlock as MethaneCombustionParameterBlock,
)
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.solvers import get_solver
from idaes.core.util.scaling import (
    jacobian_cond,
)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver(solver="ipopt_v2", writer_config={"scale_model": True})


# -----------------------------------------------------------------------------
@pytest.fixture
def test_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.properties = PhysicalParameterTestBlock()

    m.fs.unit = GibbsReactor(
        property_package=m.fs.properties,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    return m


class DummyScaler:
    def variable_scaling_routine(self, model, overwrite):
        model._dummy_scaler_test = overwrite

    def constraint_scaling_routine(self, model, overwrite):
        model._dummy_scaler_test = overwrite


class TestVariableScaling:
    @pytest.mark.unit
    def test_variable_scaling_no_input(self, test_model):
        scaler = GibbsReactorScaler()

        scaler.variable_scaling_routine(test_model.fs.unit)

        for v in test_model.fs.unit.lagrange_mult.values():
            print(test_model.fs.unit.scaling_factor[v])
            assert test_model.fs.unit.scaling_factor[v] == pytest.approx(
                1 / (8.314 * 500), rel=1e-4
            )

        for v in test_model.fs.unit.control_volume.heat.values():
            assert test_model.fs.unit.control_volume.scaling_factor[v] == pytest.approx(
                1e-6, rel=1e-4
            )

        for v in test_model.fs.unit.control_volume.deltaP.values():
            assert test_model.fs.unit.control_volume.scaling_factor[v] == pytest.approx(
                1e-3, rel=1e-4
            )

    @pytest.mark.unit
    def test_variable_scaling_no_heat_deltaP(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        scaler = GibbsReactorScaler()

        scaler.variable_scaling_routine(m.fs.unit)

        for v in m.fs.unit.lagrange_mult.values():
            assert m.fs.unit.scaling_factor[v] == pytest.approx(
                1 / (8.314 * 500), rel=1e-4
            )

    @pytest.mark.unit
    def test_variable_scaling_inlet_state(self, test_model):
        prop_in = test_model.fs.unit.control_volume.properties_in[0]
        sfx = prop_in.scaling_factor = Suffix(direction=Suffix.EXPORT)
        sfx[prop_in.temperature] = 1e-2
        sfx[prop_in.pressure] = 1e-5
        for j in prop_in.flow_mol_phase_comp.values():
            sfx[j] = 1e-2

        scaler = GibbsReactorScaler()

        scaler.variable_scaling_routine(test_model.fs.unit)

        # Outlet properties should now have scaling factors
        prop_out = test_model.fs.unit.control_volume.properties_out[0]
        assert prop_out.scaling_factor[prop_out.temperature] == 1e-2
        assert prop_out.scaling_factor[prop_out.pressure] == 1e-5
        for j in prop_out.flow_mol_phase_comp.values():
            prop_out.scaling_factor[j] == 1e-2

        for v in test_model.fs.unit.lagrange_mult.values():
            assert test_model.fs.unit.scaling_factor[v] == pytest.approx(
                1 / (8.314 * 100), rel=1e-4
            )

        for v in test_model.fs.unit.control_volume.heat.values():
            assert test_model.fs.unit.control_volume.scaling_factor[v] == pytest.approx(
                1e-6, rel=1e-4
            )

        for v in test_model.fs.unit.control_volume.deltaP.values():
            assert test_model.fs.unit.control_volume.scaling_factor[v] == pytest.approx(
                1e-3, rel=1e-4
            )

    @pytest.mark.unit
    def test_variable_scaling_submodel_scalers(self, test_model):
        scaler = GibbsReactorScaler()

        scaler.variable_scaling_routine(
            test_model.fs.unit,
            submodel_scalers={
                "control_volume.properties_in": DummyScaler(),
                "control_volume.properties_out": DummyScaler(),
            },
        )

        # Check to see if testing attribute was created correctly
        assert not test_model.fs.unit.control_volume.properties_in[0]._dummy_scaler_test
        assert not test_model.fs.unit.control_volume.properties_out[
            0
        ]._dummy_scaler_test


class TestConstraintScaling:
    @pytest.mark.unit
    def test_constraint_scaling_no_inputs(self, test_model):
        scaler = GibbsReactorScaler()

        scaler.constraint_scaling_routine(test_model.fs.unit)

        sfx = test_model.fs.unit.control_volume.scaling_factor

        assert sfx[
            test_model.fs.unit.control_volume.element_balances[0.0, "H"]
        ] == pytest.approx(0.05, rel=1e-5)
        assert sfx[
            test_model.fs.unit.control_volume.element_balances[0.0, "He"]
        ] == pytest.approx(0.0357143, rel=1e-5)
        assert sfx[
            test_model.fs.unit.control_volume.element_balances[0.0, "Li"]
        ] == pytest.approx(0.0277778, rel=1e-5)
        assert sfx[
            test_model.fs.unit.control_volume.enthalpy_balances[0.0]
        ] == pytest.approx(0.2, rel=1e-5)
        assert sfx[
            test_model.fs.unit.control_volume.pressure_balance[0.0]
        ] == pytest.approx(0.333333, rel=1e-5)

        for k, v in test_model.fs.unit.gibbs_minimization.items():
            if k[2] == "c1":
                assert test_model.fs.unit.scaling_factor[v] == pytest.approx(
                    0.142857, rel=1e-5
                )
            else:
                assert test_model.fs.unit.scaling_factor[v] == pytest.approx(
                    0.0625, rel=1e-5
                )

    @pytest.mark.unit
    def test_constraint_scaling_inerts(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
            inert_species=["c1"],
        )

        scaler = GibbsReactorScaler()

        scaler.constraint_scaling_routine(m.fs.unit)

        sfx = m.fs.unit.control_volume.scaling_factor

        assert sfx[
            m.fs.unit.control_volume.element_balances[0.0, "H"]
        ] == pytest.approx(0.05, rel=1e-5)
        assert sfx[
            m.fs.unit.control_volume.element_balances[0.0, "He"]
        ] == pytest.approx(0.0357143, rel=1e-5)
        assert sfx[
            m.fs.unit.control_volume.element_balances[0.0, "Li"]
        ] == pytest.approx(0.0277778, rel=1e-5)
        assert sfx[m.fs.unit.control_volume.enthalpy_balances[0.0]] == pytest.approx(
            0.2, rel=1e-5
        )
        assert sfx[m.fs.unit.control_volume.pressure_balance[0.0]] == pytest.approx(
            0.333333, rel=1e-5
        )

        for k, v in m.fs.unit.gibbs_minimization.items():
            if k[2] == "c1":
                assert m.fs.unit.scaling_factor[v] == pytest.approx(0.142857, rel=1e-5)
            else:
                assert m.fs.unit.scaling_factor[v] == pytest.approx(0.0625, rel=1e-5)

        for k, v in m.fs.unit.inert_species_balance.items():
            assert m.fs.unit.scaling_factor[v] == pytest.approx(0.5, rel=1e-5)

    @pytest.mark.unit
    def test_constraint_scaling_submodel_scalers(self, test_model):
        scaler = GibbsReactorScaler()

        scaler.constraint_scaling_routine(
            test_model.fs.unit,
            submodel_scalers={
                "control_volume.properties_in": DummyScaler(),
                "control_volume.properties_out": DummyScaler(),
            },
        )

        # Check to see if testing attribute was created correctly
        assert not test_model.fs.unit.control_volume.properties_in[0]._dummy_scaler_test
        assert not test_model.fs.unit.control_volume.properties_out[
            0
        ]._dummy_scaler_test


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
