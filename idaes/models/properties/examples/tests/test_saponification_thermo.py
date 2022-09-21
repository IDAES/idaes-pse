#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for saponification property package example.
Authors: Andrew Lee
"""

import pytest
from pyomo.environ import ConcreteModel, Constraint, Param, value, Var
from pyomo.util.check_units import assert_units_consistent
from idaes.core import MaterialBalanceType, EnergyBalanceType, MaterialFlowBasis

from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
    SaponificationStateBlock,
)

from idaes.core.solvers import get_solver


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = SaponificationParameterBlock()

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert len(model.params.config) == 1

    @pytest.mark.unit
    def test_build(self, model):
        assert model.params.state_block_class is SaponificationStateBlock

        assert len(model.params.phase_list) == 1
        for i in model.params.phase_list:
            assert i == "Liq"

        assert len(model.params.component_list) == 5
        for i in model.params.component_list:
            assert i in ["H2O", "NaOH", "EthylAcetate", "SodiumAcetate", "Ethanol"]

        assert isinstance(model.params.cp_mol, Param)
        assert value(model.params.cp_mol) == 75.327

        assert isinstance(model.params.dens_mol, Param)
        assert value(model.params.dens_mol) == 55388

        assert isinstance(model.params.pressure_ref, Param)
        assert value(model.params.pressure_ref) == 101325

        assert isinstance(model.params.temperature_ref, Param)
        assert value(model.params.temperature_ref) == 298.15


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = SaponificationParameterBlock()

        model.props = model.params.build_state_block([1])

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert isinstance(model.props[1].flow_vol, Var)
        assert value(model.props[1].flow_vol) == 1

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 298.15

        assert isinstance(model.props[1].conc_mol_comp, Var)
        assert len(model.props[1].conc_mol_comp) == 5
        for i in model.props[1].conc_mol_comp:
            assert value(model.props[1].conc_mol_comp[i]) == 100

        assert isinstance(model.props[1].conc_water_eqn, Constraint)
        assert len(model.props[1].flow_vol) == 1

    @pytest.mark.unit
    def test_build_defined_state(self):
        model = ConcreteModel()
        model.params = SaponificationParameterBlock()

        model.props = model.params.build_state_block([1], defined_state=True)

        assert isinstance(model.props[1].flow_vol, Var)
        assert value(model.props[1].flow_vol) == 1

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 101325

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 298.15

        assert isinstance(model.props[1].conc_mol_comp, Var)
        assert len(model.props[1].conc_mol_comp) == 5
        for i in model.props[1].conc_mol_comp:
            assert value(model.props[1].conc_mol_comp[i]) == 100

        assert not hasattr(model.props[1], "conc_water_eqn")

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert str(model.props[1].get_material_flow_terms(p, j)) == str(
                    model.props[1].flow_vol * model.props[1].conc_mol_comp[j]
                )

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert str(model.props[1].get_enthalpy_flow_terms(p)) == str(
                model.props[1].flow_vol
                * model.props[1].params.dens_mol
                * model.props[1].params.cp_mol
                * (model.props[1].temperature - model.props[1].params.temperature_ref)
            )

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert str(model.props[1].get_material_density_terms(p, j)) == str(
                    model.props[1].conc_mol_comp[j]
                )

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert str(model.props[1].get_energy_density_terms(p)) == str(
                model.props[1].params.dens_mol
                * model.props[1].params.cp_mol
                * (model.props[1].temperature - model.props[1].params.temperature_ref)
            )

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert (
            model.props[1].default_material_balance_type()
            == MaterialBalanceType.componentPhase
        )

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert (
            model.props[1].default_energy_balance_type()
            == EnergyBalanceType.enthalpyTotal
        )

    @pytest.mark.unit
    def test_get_material_flow_basis(self, model):
        assert model.props[1].get_material_flow_basis() == MaterialFlowBasis.molar

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_vol", "conc_mol_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_vol", "conc_mol_comp", "temperature", "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in [
                "Volumetric Flowrate",
                "Molar Concentration",
                "Temperature",
                "Pressure",
            ]

    @pytest.mark.unit
    def test_model_check_none(self, model, caplog):
        assert model.props[1].model_check() is None
        assert "Temperature set below lower bound" not in caplog.text
        assert "Temperature set above upper bound" not in caplog.text
        assert "Pressure set below lower bound" not in caplog.text
        assert "Pressure set above upper bound" not in caplog.text

    @pytest.mark.unit
    def test_model_check_low_T(self, model, caplog):
        model.props[1].temperature.value = 200
        assert model.props[1].model_check() is None
        assert "Temperature set below lower bound" in caplog.text
        assert "Temperature set above upper bound" not in caplog.text
        assert "Pressure set below lower bound" not in caplog.text
        assert "Pressure set above upper bound" not in caplog.text

    @pytest.mark.unit
    def test_model_check_high_T(self, model, caplog):
        model.props[1].temperature.value = 350
        assert model.props[1].model_check() is None
        assert "Temperature set below lower bound" not in caplog.text
        assert "Temperature set above upper bound" in caplog.text
        assert "Pressure set below lower bound" not in caplog.text
        assert "Pressure set above upper bound" not in caplog.text
        # Reset temeprature
        model.props[1].temperature.value = 298.15

    @pytest.mark.unit
    def test_model_check_low_P(self, model, caplog):
        model.props[1].pressure.value = 1e2
        assert model.props[1].model_check() is None
        assert "Temperature set below lower bound" not in caplog.text
        assert "Temperature set above upper bound" not in caplog.text
        assert "Pressure set below lower bound" in caplog.text
        assert "Pressure set above upper bound" not in caplog.text

    @pytest.mark.unit
    def test_model_check_high_P(self, model, caplog):
        model.props[1].pressure.value = 1e7
        assert model.props[1].model_check() is None
        assert "Temperature set below lower bound" not in caplog.text
        assert "Temperature set above upper bound" not in caplog.text
        assert "Pressure set below lower bound" not in caplog.text
        assert "Pressure set above upper bound" in caplog.text
        # Reset pressure
        model.props[1].pressure.value = 101325

    @pytest.mark.unit
    def test_initialize(self, model):
        assert not model.props[1].flow_vol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].conc_mol_comp:
            assert not model.props[1].conc_mol_comp[i].fixed

        model.props.initialize(hold_state=False, outlvl=1)

        assert not model.props[1].flow_vol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].conc_mol_comp:
            assert not model.props[1].conc_mol_comp[i].fixed

    @pytest.mark.unit
    def test_initialize_hold(self, model):
        assert not model.props[1].flow_vol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].conc_mol_comp:
            assert not model.props[1].conc_mol_comp[i].fixed

        flags = model.props.initialize(hold_state=True)

        assert model.props[1].flow_vol.fixed
        assert model.props[1].temperature.fixed
        assert model.props[1].pressure.fixed
        for i in model.props[1].conc_mol_comp:
            assert model.props[1].conc_mol_comp[i].fixed

        model.props.release_state(flags, outlvl=1)

        assert not model.props[1].flow_vol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].conc_mol_comp:
            assert not model.props[1].conc_mol_comp[i].fixed

    @pytest.mark.component
    def check_units(self, model):
        assert_units_consistent(model)
