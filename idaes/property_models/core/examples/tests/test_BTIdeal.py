##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Tests for Ideal + Ideal Liquid (i.e. no activity coefficient) state block;
only tests for construction as parameters need to be provided or estimated
from VLE data to compute the activity coefficients.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Param,
                           Set,
                           SolverFactory,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Var)

from idaes.core import (FlowsheetBlock,
                        MaterialBalanceType,
                        EnergyBalanceType,
                        MaterialFlowBasis)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import get_default_solver

from idaes.property_models.core.state_definitions import FPTx
import idaes.property_models.core.eos.ideal as ideal
from idaes.property_models.core.phase_equil import smooth_VLE
from idaes.property_models.core.generic.bubble_dew import (bubble_temp_ideal,
                                                           dew_temp_ideal,
                                                           bubble_press_ideal,
                                                           dew_press_ideal)

import idaes.property_models.core.pure.Perrys as Perrys
import idaes.property_models.core.pure.RPP as RPP

from idaes.property_models.core.examples.BT_ideal \
    import BTIdealParameterBlock


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = BTIdealParameterBlock()

        return model

    def test_config(self, model):
        assert len(model.params.config) == 15

        assert model.params.config.state_definition == FPTx

        assert model.params.config.state_bounds == {
                "flow_mol": (0, 1000),
                "temperature": (273.15, 450),
                "pressure": (5e4, 1e6)}

        assert model.params.config.equation_of_state == {
                "Vap": ideal,
                "Liq": ideal}

        assert model.params.config.phase_equilibrium_formulation == smooth_VLE

        assert model.params.config.bubble_temperature == bubble_temp_ideal
        assert model.params.config.dew_temperature == dew_temp_ideal
        assert model.params.config.bubble_pressure == bubble_press_ideal
        assert model.params.config.dew_pressure == dew_press_ideal

        assert model.params.config.dens_mol_comp_liq == Perrys
        assert model.params.config.enth_mol_comp_liq == Perrys
        assert model.params.config.enth_mol_comp_ig == RPP
        assert model.params.config.entr_mol_comp_liq == Perrys
        assert model.params.config.entr_mol_comp_ig == RPP
        assert model.params.config.pressure_sat_comp == RPP

    def test_build(self, model):
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]

        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ['benzene',
                         'toluene']

        for i in model.params.phase_comp.values():
            assert i == model.params.component_list

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 2

        assert model.params.phase_equilibrium_list == {
                1: ["benzene", ("Vap", "Liq")],
                2: ["toluene", ("Vap", "Liq")]}

        assert isinstance(model.params.pressure_ref, Param)
        assert isinstance(model.params.temperature_ref, Param)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = BTIdealParameterBlock()

        model.props = model.params.state_block_class(
                [1],
                default={"parameters": model.params,
                         "defined_state": True})

        return model

    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 500
        assert model.props[1].flow_mol.ub == 1000
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 5.25e5
        assert model.props[1].pressure.ub == 1e6
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 361.575
        assert model.props[1].temperature.ub == 450
        assert model.props[1].temperature.lb == 273.15

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 2
        for i in model.props[1].mole_frac_comp:
            assert value(model.props[1].mole_frac_comp[i]) == 0.5

        # Check supporting variables
        assert isinstance(model.props[1].flow_mol_phase, Var)
        assert len(model.props[1].flow_mol_phase) == 2

        assert isinstance(model.props[1].mole_frac_phase_comp, Var)
        assert len(model.props[1].mole_frac_phase_comp) == 4

        assert isinstance(model.props[1].phase_frac, Var)
        assert len(model.props[1].phase_frac) == 2

        assert isinstance(model.props[1].total_flow_balance, Constraint)
        assert len(model.props[1].total_flow_balance) == 1

        assert isinstance(model.props[1].component_flow_balances, Constraint)
        assert len(model.props[1].component_flow_balances) == 2

        assert isinstance(model.props[1].sum_mole_frac, Constraint)
        assert len(model.props[1].sum_mole_frac) == 1

        assert not hasattr(model.props[1], "sum_mole_frac_out")

        assert isinstance(model.props[1].phase_fraction_constraint, Constraint)
        assert len(model.props[1].phase_fraction_constraint) == 2

    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_flow_terms(p, j) == (
                    model.props[1].flow_mol_phase[p] *
                    model.props[1].mole_frac_phase_comp[p, j])

    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_enthalpy_flow_terms(p) == (
                model.props[1].flow_mol_phase[p] *
                model.props[1].enth_mol_phase[p])

    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_density_terms(p, j) == (
                    model.props[1].dens_mol_phase[p] *
                    model.props[1].mole_frac_phase_comp[p, j])

    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_energy_density_terms(p) == (
                model.props[1].dens_mol_phase[p] *
                model.props[1].enth_mol_phase[p])

    def test_default_material_balance_type(self, model):
        assert model.props[1].default_material_balance_type() == \
            MaterialBalanceType.componentTotal

    def test_default_energy_balance_type(self, model):
        assert model.props[1].default_energy_balance_type() == \
            EnergyBalanceType.enthalpyTotal

    def test_get_material_flow_basis(self, model):
        assert model.props[1].get_material_flow_basis() == \
            MaterialFlowBasis.molar

    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    def test_initialize(self, model):
        assert not model.props[1].flow_mol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].mole_frac_comp:
            assert not model.props[1].mole_frac_comp[i].fixed

        model.props.initialize(hold_state=False, outlvl=1)

        assert not model.props[1].flow_mol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].mole_frac_comp:
            assert not model.props[1].mole_frac_comp[i].fixed

    def test_initialize_hold(self, model):
        assert not model.props[1].flow_mol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].mole_frac_comp:
            assert not model.props[1].mole_frac_comp[i].fixed

        flags = model.props.initialize(hold_state=True)

        assert model.props[1].flow_mol.fixed
        assert model.props[1].temperature.fixed
        assert model.props[1].pressure.fixed
        for i in model.props[1].mole_frac_comp:
            assert model.props[1].mole_frac_comp[i].fixed

        model.props.release_state(flags, outlvl=1)

        assert not model.props[1].flow_mol.fixed
        assert not model.props[1].temperature.fixed
        assert not model.props[1].pressure.fixed
        for i in model.props[1].mole_frac_comp:
            assert not model.props[1].mole_frac_comp[i].fixed
