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
import os
import pytest

from pyomo.environ import (
    assert_optimal_termination,
    ComponentMap,
    ConcreteModel,
    Constraint,
    Suffix,
    TransformationFactory,
    units,
    value,
    Var,
)

from idaes.core import FlowsheetBlock
from idaes.models.unit_models.gibbs_reactor import GibbsReactor, GibbsReactorScaler
from idaes.models.properties.activity_coeff_models.methane_combustion_ideal import (
    MethaneParameterBlock as MethaneCombustionParameterBlock,
)
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.scaling import (
    jacobian_cond,
)
from idaes.core.util import from_json, StoreSpec
from idaes.core.scaling import CustomScalerBase, get_scaling_factor, set_scaling_factor


# Get solution json from scaling tests
FILENAME = "gibbs_solution.json"
local_path = os.path.dirname(os.path.realpath(__file__))
fname = os.path.join(local_path, "..", "..", "..", "core", "scaling", "tests", FILENAME)


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


@pytest.mark.unit
class TestVariableScaling:

    def test_variable_scaling_no_input(self, test_model):
        scaler = GibbsReactorScaler()

        scaler.variable_scaling_routine(test_model.fs.unit)

        for v in test_model.fs.unit.lagrange_mult.values():
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

    def test_variable_scaling_no_heat_deltaP(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = PhysicalParameterTestBlock()

        m.fs.unit = GibbsReactor(
            property_package=m.fs.properties,
            has_heat_transfer=False,
            has_pressure_change=False,
        )

        scaler = GibbsReactorScaler()

        scaler.variable_scaling_routine(m.fs.unit)

        for v in m.fs.unit.lagrange_mult.values():
            assert m.fs.unit.scaling_factor[v] == pytest.approx(
                1 / (8.314 * 500), rel=1e-4
            )

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

    def test_variable_scaling_submodel_scalers(self, test_model):
        scaler = GibbsReactorScaler()

        scaler_map = ComponentMap()
        scaler_map[test_model.fs.unit.control_volume.properties_in] = DummyScaler()
        scaler_map[test_model.fs.unit.control_volume.properties_out] = DummyScaler()

        scaler.variable_scaling_routine(
            test_model.fs.unit,
            submodel_scalers=scaler_map,
        )

        # Check to see if testing attribute was created correctly
        assert not test_model.fs.unit.control_volume.properties_in[0]._dummy_scaler_test
        assert not test_model.fs.unit.control_volume.properties_out[
            0
        ]._dummy_scaler_test


@pytest.mark.unit
class TestConstraintScaling:

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
        ] == pytest.approx(0.25, rel=1e-5)
        assert sfx[
            test_model.fs.unit.control_volume.pressure_balance[0.0]
        ] == pytest.approx(5e-6, rel=1e-5)

        for k, v in test_model.fs.unit.gibbs_minimization.items():
            if k[2] == "c1":
                assert test_model.fs.unit.scaling_factor[v] == pytest.approx(
                    1.53846e-3, rel=1e-5
                )
            else:
                assert test_model.fs.unit.scaling_factor[v] == pytest.approx(
                    6.45161e-4, rel=1e-5
                )

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
            0.25, rel=1e-5
        )
        assert sfx[m.fs.unit.control_volume.pressure_balance[0.0]] == pytest.approx(
            5e-6, rel=1e-5
        )

        for k, v in m.fs.unit.gibbs_minimization.items():
            assert m.fs.unit.scaling_factor[v] == pytest.approx(6.45161e-4, rel=1e-5)

        for k, v in m.fs.unit.inert_species_balance.items():
            assert m.fs.unit.scaling_factor[v] == pytest.approx(0.5, rel=1e-5)

    def test_constraint_scaling_submodel_scalers(self, test_model):
        scaler = GibbsReactorScaler()

        scaler_map = ComponentMap()
        scaler_map[test_model.fs.unit.control_volume.properties_in] = DummyScaler()
        scaler_map[test_model.fs.unit.control_volume.properties_out] = DummyScaler()

        scaler.constraint_scaling_routine(
            test_model.fs.unit,
            submodel_scalers=scaler_map,
        )

        # Check to see if testing attribute was created correctly
        assert not test_model.fs.unit.control_volume.properties_in[0]._dummy_scaler_test
        assert not test_model.fs.unit.control_volume.properties_out[
            0
        ]._dummy_scaler_test


# -----------------------------------------------------------------------------
class SMScaler(CustomScalerBase):
    def variable_scaling_routine(self, model, overwrite):
        pass

    def constraint_scaling_routine(self, model, overwrite):
        for c in model.component_data_objects(ctype=Constraint, descend_into=True):
            self.scale_constraint_by_nominal_value(
                c, scheme="inverse_sum", overwrite=overwrite
            )


# TODO: Turn this into a testing harness?
@pytest.mark.integration
class TestMethaneScaling(object):
    @pytest.fixture
    def methane(self):
        model = ConcreteModel()
        model.fs = FlowsheetBlock(dynamic=False)

        model.fs.properties = MethaneCombustionParameterBlock()

        model.fs.unit = GibbsReactor(
            property_package=model.fs.properties,
            has_heat_transfer=True,
            has_pressure_change=True,
        )

        model.fs.unit.inlet.flow_mol[0].fix(230.0)
        model.fs.unit.inlet.mole_frac_comp[0, "H2"].fix(0.0435)
        model.fs.unit.inlet.mole_frac_comp[0, "N2"].fix(0.6522)
        model.fs.unit.inlet.mole_frac_comp[0, "O2"].fix(0.1739)
        model.fs.unit.inlet.mole_frac_comp[0, "CO2"].fix(1e-5)
        model.fs.unit.inlet.mole_frac_comp[0, "CH4"].fix(0.1304)
        model.fs.unit.inlet.mole_frac_comp[0, "CO"].fix(1e-5)
        model.fs.unit.inlet.mole_frac_comp[0, "H2O"].fix(1e-5)
        model.fs.unit.inlet.mole_frac_comp[0, "NH3"].fix(1e-5)
        model.fs.unit.inlet.temperature[0].fix(1500.0)
        model.fs.unit.inlet.pressure[0].fix(101325.0)

        model.fs.unit.outlet.temperature[0].fix(2844.38)
        model.fs.unit.deltaP.fix(0)

        # Set imperfect scaling factors for all variables, representing an initial "best-guess"
        # Feed states are known exactly - set scaling based on these
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].flow_mol, 1 / 230
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].flow_mol_phase, 1 / 230
        )  # Only 1 phase, so we "know" this
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["H2"],
            1 / 0.0435,
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["N2"],
            1 / 0.6522,
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["O2"],
            1 / 0.1739,
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CO2"], 1e5
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CH4"],
            1 / 0.1304,
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["CO"], 1e5
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["H2O"], 1e5
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].mole_frac_comp["NH3"], 1e5
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].temperature, 1 / 1500
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_in[0.0].pressure, 1e-5
        )
        # Assume user does not know anything about enthalpy

        # Best guesses for unit model and outlet state conditions
        set_scaling_factor(model.fs.unit.control_volume.heat[0.0], 1e-6)

        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].flow_mol, 1e-2
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].flow_mol_phase, 1e-2
        )  # Only 1 phase, so we "know" this
        # N2 is inert, so will be order 0.1, assume CH4 and H2 are near-totally consumed, assume most O2 consumed
        # Assume moderate amounts of CO2 and H2O, small amounts of CO, trace NH3 NH3
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2"], 1e4
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["N2"], 10
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["O2"], 1e2
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO2"], 10
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CH4"], 1e4
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["CO"], 1e3
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["H2O"], 10
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].mole_frac_comp["NH3"], 1e4
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].temperature, 1e-3
        )
        set_scaling_factor(
            model.fs.unit.control_volume.properties_out[0.0].pressure, 1e-5
        )

        from_json(model, fname=fname, wts=StoreSpec.value())

        return model

    def test_variable_scaling_only(self, methane):
        unscaled = jacobian_cond(methane, scaled=False)

        scaler_map = ComponentMap()
        scaler_map[methane.fs.unit.control_volume.properties_in] = SMScaler()
        scaler_map[methane.fs.unit.control_volume.properties_out] = SMScaler()

        scaler = GibbsReactorScaler()
        scaler.variable_scaling_routine(
            methane.fs.unit,
            submodel_scalers=scaler_map,
        )

        scaled = jacobian_cond(methane, scaled=True)

        count = 0
        for c in methane.component_data_objects(ctype=[Var], descend_into=True):
            sf = get_scaling_factor(c)
            if sf is None:
                count += 1
        assert count == 52

        assert scaled < unscaled
        assert scaled == pytest.approx(8.908989e16, rel=1e-5)

    def test_constraint_scaling_only(self, methane):
        unscaled = jacobian_cond(methane, scaled=False)

        scaler_map = ComponentMap()
        scaler_map[methane.fs.unit.control_volume.properties_in] = SMScaler()
        scaler_map[methane.fs.unit.control_volume.properties_out] = SMScaler()

        scaler = GibbsReactorScaler()
        scaler.constraint_scaling_routine(
            methane.fs.unit,
            submodel_scalers=scaler_map,
        )

        scaled = jacobian_cond(methane, scaled=True)

        count = 0
        for c in methane.component_data_objects(ctype=[Constraint], descend_into=True):
            sf = get_scaling_factor(c)
            if sf is None:
                count += 1
        assert count == 0

        assert scaled < unscaled
        assert scaled == pytest.approx(9.316e15, rel=1e-2)

    def test_full_scaling(self, methane):
        unscaled = jacobian_cond(methane, scaled=False)

        scaler_map = ComponentMap()
        scaler_map[methane.fs.unit.control_volume.properties_in] = SMScaler()
        scaler_map[methane.fs.unit.control_volume.properties_out] = SMScaler()

        scaler = GibbsReactorScaler()
        scaler.scale_model(
            methane.fs.unit,
            submodel_scalers=scaler_map,
        )

        scaled = jacobian_cond(methane, scaled=True)

        assert scaled < unscaled
        assert scaled == pytest.approx(7.653e15, rel=1e-2)
