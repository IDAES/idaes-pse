##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Author: Andrew Lee, Alejandro Garciadiego
"""
import pytest
from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Set,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           Var,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import (MaterialBalanceType,
                        EnergyBalanceType,
                        MaterialFlowBasis,
                        Component)
from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)
from idaes.core.util.testing import get_default_solver

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)

from idaes.generic_models.properties.core.state_definitions import FTPx
from idaes.generic_models.properties.core.phase_equil import smooth_VLE

from idaes.generic_models.properties.core.examples.ASU_PR \
    import configuration


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_default_solver()

# Test for configuration dictionaries with parameters from Properties of Gases
# and liquids 4th edition
class TestParamBlock(object):
    @pytest.mark.unit
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
        assert len(model.params.component_list) == 3
        for i in model.params.component_list:
            assert i in ['nitrogen',
                         'argon',
                         'oxygen']
            assert isinstance(model.params.get_component(i), Component)

        assert isinstance(model.params._phase_component_set, Set)
        assert len(model.params._phase_component_set) == 6
        for i in model.params._phase_component_set:
            assert i in [("Liq", "nitrogen"), ("Liq", "argon"), ("Liq", "oxygen"),
                         ("Vap", "nitrogen"), ("Vap", "argon"), ("Vap", "oxygen")]

        assert model.params.config.state_definition == FTPx

        assert model.params.config.state_bounds == {
            "flow_mol": (0, 100, 1000, pyunits.mol/pyunits.s),
            "temperature": (10, 300, 350, pyunits.K),
            "pressure": (5e4, 1e5, 1e7, pyunits.Pa)}

        assert model.params.config.phase_equilibrium_state == {
            ("Vap", "Liq"): smooth_VLE}

        assert isinstance(model.params.phase_equilibrium_idx, Set)
        assert len(model.params.phase_equilibrium_idx) == 3
        for i in model.params.phase_equilibrium_idx:
            assert i in ["PE1", "PE2", "PE3"]

        assert model.params.phase_equilibrium_list == {
            "PE1": {"nitrogen": ("Vap", "Liq")},
            "PE2": {"argon": ("Vap", "Liq")},
            "PE3": {"oxygen": ("Vap", "Liq")}}

        assert model.params.pressure_ref.value == 101325
        assert model.params.temperature_ref.value == 298.15

        assert model.params.nitrogen.mw.value == 28.0135E-3
        assert model.params.nitrogen.pressure_crit.value == 34e5
        assert model.params.nitrogen.temperature_crit.value == 126.2

        assert model.params.argon.mw.value == 39.948E-3
        assert model.params.argon.pressure_crit.value == 48.98e5
        assert model.params.argon.temperature_crit.value == 150.86

        assert model.params.oxygen.mw.value == 31.999E-3
        assert model.params.oxygen.pressure_crit.value == 50.43e5
        assert model.params.oxygen.temperature_crit.value == 154.58

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})

        return model

    @pytest.mark.unit
    def test_build(self, model):
        # Check state variable values and bounds
        assert isinstance(model.props[1].flow_mol, Var)
        assert value(model.props[1].flow_mol) == 100
        assert model.props[1].flow_mol.ub == 1000
        assert model.props[1].flow_mol.lb == 0

        assert isinstance(model.props[1].pressure, Var)
        assert value(model.props[1].pressure) == 1e5
        assert model.props[1].pressure.ub == 1e7
        assert model.props[1].pressure.lb == 5e4

        assert isinstance(model.props[1].temperature, Var)
        assert value(model.props[1].temperature) == 300
        assert model.props[1].temperature.ub == 350
        assert model.props[1].temperature.lb == 10

        assert isinstance(model.props[1].mole_frac_comp, Var)
        assert len(model.props[1].mole_frac_comp) == 3
        for i in model.props[1].mole_frac_comp:
            assert value(model.props[1].mole_frac_comp[i]) == 1/3

        # Check supporting variables
        assert isinstance(model.props[1].flow_mol_phase, Var)
        assert len(model.props[1].flow_mol_phase) == 2

        assert isinstance(model.props[1].mole_frac_phase_comp, Var)
        assert len(model.props[1].mole_frac_phase_comp) == 6

        assert isinstance(model.props[1].phase_frac, Var)
        assert len(model.props[1].phase_frac) == 2

        assert isinstance(model.props[1].total_flow_balance, Constraint)
        assert len(model.props[1].total_flow_balance) == 1

        assert isinstance(model.props[1].component_flow_balances, Constraint)
        assert len(model.props[1].component_flow_balances) == 3

        assert isinstance(model.props[1].sum_mole_frac, Constraint)
        assert len(model.props[1].sum_mole_frac) == 1

        assert not hasattr(model.props[1], "sum_mole_frac_out")

        assert isinstance(model.props[1].phase_fraction_constraint, Constraint)
        assert len(model.props[1].phase_fraction_constraint) == 2

        assert_units_consistent(model)

    @pytest.mark.unit
    def test_get_material_flow_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_flow_terms(p, j) == (
                    model.props[1].flow_mol_phase[p] *
                    model.props[1].mole_frac_phase_comp[p, j])

    @pytest.mark.unit
    def test_get_enthalpy_flow_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_enthalpy_flow_terms(p) == (
                model.props[1].flow_mol_phase[p] *
                model.props[1].enth_mol_phase[p])

    @pytest.mark.unit
    def test_get_material_density_terms(self, model):
        for p in model.params.phase_list:
            for j in model.params.component_list:
                assert model.props[1].get_material_density_terms(p, j) == (
                    model.props[1].dens_mol_phase[p] *
                    model.props[1].mole_frac_phase_comp[p, j])

    @pytest.mark.unit
    def test_get_energy_density_terms(self, model):
        for p in model.params.phase_list:
            assert model.props[1].get_energy_density_terms(p) == (
                model.props[1].dens_mol_phase[p] *
                model.props[1].enth_mol_phase[p])

    @pytest.mark.unit
    def test_default_material_balance_type(self, model):
        assert model.props[1].default_material_balance_type() == \
            MaterialBalanceType.componentTotal

    @pytest.mark.unit
    def test_default_energy_balance_type(self, model):
        assert model.props[1].default_energy_balance_type() == \
            EnergyBalanceType.enthalpyTotal

    @pytest.mark.unit
    def test_get_material_flow_basis(self, model):
        assert model.props[1].get_material_flow_basis() == \
            MaterialFlowBasis.molar

    @pytest.mark.unit
    def test_define_state_vars(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    @pytest.mark.unit
    def test_define_port_members(self, model):
        sv = model.props[1].define_state_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["flow_mol",
                         "mole_frac_comp",
                         "temperature",
                         "pressure"]

    @pytest.mark.unit
    def test_define_display_vars(self, model):
        sv = model.props[1].define_display_vars()

        assert len(sv) == 4
        for i in sv:
            assert i in ["Total Molar Flowrate",
                         "Total Mole Fraction",
                         "Temperature",
                         "Pressure"]

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        # Fix state
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(85.00)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["nitrogen"].fix(1/3)
        model.props[1].mole_frac_comp["argon"].fix(1/3)
        model.props[1].mole_frac_comp["oxygen"].fix(1/3)

        assert degrees_of_freedom(model.props[1]) == 0

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
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        # Check phase equilibrium results
        assert model.props[1].mole_frac_phase_comp["Liq", "nitrogen"].value == \
            pytest.approx(0.1739, abs=1e-4)
        assert model.props[1].mole_frac_phase_comp["Vap", "nitrogen"].value == \
            pytest.approx(0.4221, abs=1e-4)
        assert model.props[1].phase_frac["Vap"].value == \
            pytest.approx(0.6422, abs=1e-4)

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.props[1].report()
