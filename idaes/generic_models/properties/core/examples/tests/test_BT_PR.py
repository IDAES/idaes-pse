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
                           Expression,
                           ExternalFunction,
                           Param,
                           Set,
                           SolverStatus,
                           sqrt,
                           TerminationCondition,
                           value,
                           Var)

from idaes.core import (MaterialBalanceType,
                        EnergyBalanceType,
                        MaterialFlowBasis,
                        Component)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)
from idaes.core.util.testing import get_default_solver
from idaes.core.util.constants import Constants as const

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.phase_equil import smooth_VLE

from idaes.generic_models.properties.core.examples.BT_PR \
    import configuration


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()


class TestParamBlock(object):
    def test_build(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

        assert isinstance(model.params.phase_list, Set)
        assert len(model.params.phase_list) == 2
        for i in model.params.phase_list:
            assert i in ["Liq", "Vap"]
        assert model.params.Liq.is_liquid_phase()
        assert model.params.Vap.is_vapor_phase()

        assert isinstance(model.params.component_list, Set)
        assert len(model.params.component_list) == 2
        for i in model.params.component_list:
            assert i in ['benzene',
                         'toluene']
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 4
        for i in model.params._phase_component_set:
            assert i in [("Liq", "benzene"), ("Liq", "toluene"),
                         ("Vap", "benzene"), ("Vap", "toluene")]

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
                "flow_mol": (0, 1000),
                "temperature": (273.15, 450),
                "pressure": (5e4, 1e6)}

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): smooth_VLE}

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 2
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"benzene": ("Vap", "Liq")},
            "PE2": {"toluene": ("Vap", "Liq")}}

        assert model.params.pressure_ref.value == 1e5
        assert model.params.temperature_ref.value == 300

        assert isinstance(model.params.PR_kappa, Var)
        for i in model.params.PR_kappa:
            assert model.params.PR_kappa[i].fixed


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})

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

        # Test cubic components
        assert isinstance(model.props[1].PR_fw, Expression)
        assert len(model.props[1].PR_fw) == 2
        for i in model.params.component_list:
            omega = model.params.get_component(i).omega
            assert value(model.props[1].PR_fw[i]) == value(
                0.37464 + 1.54226*omega - 0.26992*omega**2)

        assert isinstance(model.props[1].PR_a, Expression)
        assert len(model.props[1].PR_a) == 2
        for i in model.params.component_list:
            Tc = model.params.get_component(i).temperature_crit
            Pc = model.params.get_component(i).pressure_crit
            assert pytest.approx(value(model.props[1].PR_a[i]),
                                 rel=1e-5) == value(
                0.45724*((const.gas_constant*Tc)**2/Pc) *
                (((1+model.props[1].PR_fw[i] *
                   (1-sqrt(model.props[1].temperature/Tc)))**2)))

        assert isinstance(model.props[1].PR_b, Expression)
        assert len(model.props[1].PR_b) == 2
        for i in model.params.component_list:
            Tc = model.params.get_component(i).temperature_crit
            Pc = model.params.get_component(i).pressure_crit
            assert value(model.props[1].PR_b[i]) == value(
                0.07780*const.gas_constant*Tc/Pc)

        assert isinstance(model.props[1].PR_am, Expression)
        assert len(model.props[1].PR_am) == 2
        for p in model.params.phase_list:
            assert pytest.approx(value(model.props[1].PR_am[p]),
                                 rel=1e-5) == value(
                sum(sum(model.props[1].mole_frac_phase_comp[p, i] *
                        model.props[1].mole_frac_phase_comp[p, j] *
                        sqrt(model.props[1].PR_a[i]*model.props[1].PR_a[j]) *
                        (1-model.params.PR_kappa[i, j])
                        for j in model.params.component_list)
                    for i in model.params.component_list))

        assert isinstance(model.props[1].PR_bm, Expression)
        assert len(model.props[1].PR_bm) == 2
        for p in model.params.phase_list:
            assert pytest.approx(value(model.props[1].PR_bm[p]),
                                 rel=1e-5) == value(
                sum(model.props[1].mole_frac_phase_comp[p, i] *
                    model.props[1].PR_b[i]
                    for i in model.props[1].params.component_list))

        assert isinstance(model.props[1].PR_A, Expression)
        assert len(model.props[1].PR_A) == 2
        for p in model.params.phase_list:
            assert pytest.approx(value(model.props[1].PR_A[p]),
                                 rel=1e-5) == value(
                    model.props[1].PR_am[p]*model.props[1].pressure /
                    (const.gas_constant*model.props[1].temperature)**2)

        assert isinstance(model.props[1].PR_B, Expression)
        assert len(model.props[1].PR_B) == 2
        for p in model.params.phase_list:
            assert pytest.approx(value(model.props[1].PR_B[p]),
                                 rel=1e-5) == value(
                    model.props[1].PR_bm[p]*model.props[1].pressure /
                    (const.gas_constant*model.props[1].temperature))

        assert isinstance(model.props[1].PR_delta, Expression)
        assert len(model.props[1].PR_delta) == 4
        for p in model.params.phase_list:
            for i in model.params.component_list:
                assert pytest.approx(value(model.props[1].PR_delta[p, i]),
                                     rel=1e-5) == value(
                    2*sqrt(model.props[1].PR_a[i])/model.props[1].PR_am[p] *
                    sum(model.props[1].mole_frac_phase_comp[p, j] *
                        sqrt(model.props[1].PR_a[j]) *
                        (1-model.params.PR_kappa[i, j])
                        for j in model.params.component_list))

        assert isinstance(model.props[1].PR_dadT, Expression)
        assert len(model.props[1].PR_dadT) == 2
        for p in model.params.phase_list:
            assert pytest.approx(value(model.props[1].PR_dadT[p]),
                                 rel=1e-5) == value(
                -((const.gas_constant/2)*sqrt(0.45724) *
                  sum(sum(model.props[1].mole_frac_phase_comp[p, i] *
                          model.props[1].mole_frac_phase_comp[p, j] *
                          (1-model.params.PR_kappa[i, j]) *
                          (model.props[1].PR_fw[j] *
                           sqrt(model.props[1].PR_a[i] *
                                model.params.get_component(j).temperature_crit /
                                model.params.get_component(j).pressure_crit) +
                           model.props[1].PR_fw[i] *
                           sqrt(model.props[1].PR_a[j] *
                                model.params.get_component(i).temperature_crit /
                                model.params.get_component(i).pressure_crit))
                          for j in model.params.component_list)
                      for i in model.params.component_list) /
                  sqrt(model.props[1].temperature)))

        # Test equilibrium state Expressions
        assert isinstance(model.props[1]._PR_a_eq, Expression)
        assert len(model.props[1]._PR_a_eq) == 2
        for i in model.props[1]._PR_a_eq:
            Tc = model.params.get_component(i[2]).temperature_crit
            Pc = model.params.get_component(i[2]).pressure_crit
            assert pytest.approx(value(model.props[1]._PR_a_eq[i]),
                                 rel=1e-5) == value(
                0.45724*((const.gas_constant*Tc)**2/Pc) *
                (((1+model.props[1].PR_fw[i[2]] *
                   (1-sqrt(model.props[1]._teq[i[0], i[1]]/Tc)))**2)))

        assert isinstance(model.props[1]._PR_am_eq, Expression)
        assert len(model.props[1]._PR_am_eq) == 2
        for idx in model.props[1]._PR_am_eq:
            assert pytest.approx(value(model.props[1]._PR_am_eq[idx]),
                                 rel=1e-5) == value(
                sum(sum(model.props[1].mole_frac_phase_comp[idx[2], i] *
                        model.props[1].mole_frac_phase_comp[idx[2], j] *
                        sqrt(model.props[1]._PR_a_eq[idx[0], idx[1], i] *
                             model.props[1]._PR_a_eq[idx[0], idx[1], j]) *
                        (1-model.params.PR_kappa[i, j])
                        for j in model.params.component_list)
                    for i in model.params.component_list))

        assert isinstance(model.props[1]._PR_A_eq, Expression)
        assert len(model.props[1]._PR_A_eq) == 2
        for i in model.props[1]._PR_A_eq:
            assert pytest.approx(value(model.props[1]._PR_A_eq[i]),
                                 rel=1e-5) == value(
                        model.props[1]._PR_am_eq[i] *
                        model.props[1].pressure /
                        (const.gas_constant *
                         model.props[1]._teq[i[0], i[1]])**2)

        assert isinstance(model.props[1]._PR_B_eq, Expression)
        assert len(model.props[1]._PR_B_eq) == 2
        for i in model.props[1]._PR_B_eq:
            assert pytest.approx(value(model.props[1]._PR_B_eq[i]),
                                 rel=1e-5) == value(
                        model.props[1].PR_bm[i[2]] *
                        model.props[1].pressure /
                        (const.gas_constant *
                         model.props[1]._teq[i[0], i[1]]))

        assert isinstance(model.props[1]._PR_delta_eq, Expression)
        assert len(model.props[1]._PR_delta_eq) == 4
        for idx in model.props[1]._PR_delta_eq:
            assert pytest.approx(value(model.props[1]._PR_delta_eq[idx]),
                                 rel=1e-5) == value(
                2*sqrt(model.props[1]._PR_a_eq[idx[0], idx[1], idx[3]]) /
                model.props[1]._PR_am_eq[idx[0], idx[1], idx[2]] *
                sum(model.props[1].mole_frac_phase_comp[idx[2], j] *
                    sqrt(model.props[1]._PR_a_eq[idx[0], idx[1], j]) *
                    (1-model.params.PR_kappa[idx[3], j])
                    for j in model.params.component_list))

        # Check for external function components
        assert isinstance(model.props[1]._ext_func_param, Param)
        assert model.props[1]._ext_func_param.value == 0  # 0 == PR
        assert isinstance(model.props[1].proc_Z_liq, ExternalFunction)
        assert isinstance(model.props[1].proc_Z_vap, ExternalFunction)

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
