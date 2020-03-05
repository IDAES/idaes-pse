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
Author: Andrew Lee
"""
import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Param,
                           Set,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Var)

from idaes.core import (MaterialBalanceType,
                        EnergyBalanceType,
                        MaterialFlowBasis)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)
from idaes.core.util.testing import get_default_solver

from idaes.generic_models.properties.core.state_definitions import FTPx
import idaes.generic_models.properties.core.eos.ideal as ideal
from idaes.generic_models.properties.core.phase_equil import smooth_VLE
from idaes.generic_models.properties.core.phase_equil.bubble_dew import (
        bubble_temp_ideal,
        dew_temp_ideal,
        bubble_press_ideal,
        dew_press_ideal)

import idaes.generic_models.properties.core.pure.Perrys as Perrys
import idaes.generic_models.properties.core.pure.RPP as RPP

from idaes.generic_models.properties.core.examples.BT_ideal \
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
        assert len(model.params.config) == 19

        assert len(model.params.config.phase_list) == 2
        for i in model.params.config.phase_list:
            assert i in ["Liq", "Vap"]

        assert len(model.params.config.component_list) == 2
        for i in model.params.config.component_list:
            assert i in ['benzene',
                         'toluene']

        assert model.params.config.phase_component_list is None

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
                "flow_mol": (0, 1000),
                "temperature": (273.15, 450),
                "pressure": (5e4, 1e6)}

        assert model.params.config.equation_of_state == {
                "Vap": ideal,
                "Liq": ideal}

        assert model.params.config.phase_equilibrium_formulation == smooth_VLE
        assert len(model.params.config.phase_equilibrium_dict) == 2
        assert model.params.config.phase_equilibrium_dict == {
                1: ["benzene", ("Vap", "Liq")],
                2: ["toluene", ("Vap", "Liq")]}

        assert model.params.config.temperature_bubble == bubble_temp_ideal
        assert model.params.config.temperature_dew == dew_temp_ideal
        assert model.params.config.pressure_bubble == bubble_press_ideal
        assert model.params.config.pressure_dew == dew_press_ideal

        assert model.params.config.dens_mol_liq_comp == Perrys
        assert model.params.config.enth_mol_liq_comp == Perrys
        assert model.params.config.enth_mol_ig_comp == RPP
        assert model.params.config.entr_mol_liq_comp == Perrys
        assert model.params.config.entr_mol_ig_comp == RPP
        assert model.params.config.pressure_sat_comp == RPP

    def test_build(self, model):
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]

        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ['benzene',
                         'toluene']

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

    def test_dof(self, model):
        # Fix state
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(368)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["benzene"].fix(0.5)
        model.props[1].mole_frac_comp["toluene"].fix(0.5)

        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={'tol': 1e-6})

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp["Liq", "benzene"].value == \
            pytest.approx(0.4121, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp["Vap", "benzene"].value == \
            pytest.approx(0.6339, abs=1e-4)
        assert model.props[1].phase_frac["Vap"].value == \
            pytest.approx(0.3961, abs=1e-4)

    @pytest.mark.ui
    def test_report(self, model):
        model.props[1].report()
