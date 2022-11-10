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
Tests for Ideal + Ideal Liquid (i.e. no activity coefficient) state block;
only tests for construction as parameters need to be provided or estimated
from VLE data to compute the activity coefficients.

Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import check_optimal_termination, ConcreteModel, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

solver = get_solver()


# -----------------------------------------------------------------------------
class TestFTPz_LV_inlet:
    @pytest.fixture(scope="class")
    def model(self):
        # Create a flowsheet for test
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties_ideal_vl = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
        )
        m.fs.state_block_ideal_vl = m.fs.properties_ideal_vl.build_state_block(
            [0], defined_state=True
        )

        m.fs.state_block_ideal_vl[0].flow_mol.fix(1)
        m.fs.state_block_ideal_vl[0].temperature.fix(368)
        m.fs.state_block_ideal_vl[0].pressure.fix(101325)
        m.fs.state_block_ideal_vl[0].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state_block_ideal_vl[0].mole_frac_comp["toluene"].fix(0.5)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert len(model.fs.properties_ideal_vl.config) == 4

        assert model.fs.properties_ideal_vl.config.valid_phase == ("Liq", "Vap")
        assert len(model.fs.properties_ideal_vl.phase_list) == 2
        assert model.fs.properties_ideal_vl.phase_list == ["Liq", "Vap"]

        assert model.fs.state_block_ideal_vl[0].config.defined_state
        assert hasattr(model.fs.state_block_ideal_vl[0], "eq_phase_equilibrium")
        assert not hasattr(model.fs.state_block_ideal_vl[0], "eq_activity_coeff")
        assert not hasattr(model.fs.state_block_ideal_vl[0], "eq_mol_frac_out")

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.fs.state_block_ideal_vl[0]) == 0

    @pytest.mark.component
    def test_units_consistent(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.state_block_ideal_vl.initialize()

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert value(
            model.fs.state_block_ideal_vl[0].mole_frac_phase_comp["Liq", "benzene"]
        ) == pytest.approx(0.4121, abs=1e-3)
        assert value(
            model.fs.state_block_ideal_vl[0].mole_frac_phase_comp["Vap", "benzene"]
        ) == pytest.approx(0.6339, abs=1e-3)


class TestFTPz_L_inlet:
    @pytest.fixture(scope="class")
    def model(self):
        # Create a flowsheet for test
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties_ideal_l = BTXParameterBlock(
            valid_phase="Liq", activity_coeff_model="Ideal", state_vars="FTPz"
        )
        m.fs.state_block_ideal_l = m.fs.properties_ideal_l.build_state_block(
            [0], has_phase_equilibrium=False, defined_state=True
        )

        m.fs.state_block_ideal_l[0].flow_mol.fix(1)
        m.fs.state_block_ideal_l[0].temperature.fix(368)
        m.fs.state_block_ideal_l[0].pressure.fix(101325)
        m.fs.state_block_ideal_l[0].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state_block_ideal_l[0].mole_frac_comp["toluene"].fix(0.5)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert len(model.fs.properties_ideal_l.config) == 4

        assert model.fs.properties_ideal_l.config.valid_phase == "Liq"
        assert len(model.fs.properties_ideal_l.phase_list) == 1
        assert model.fs.properties_ideal_l.phase_list == ["Liq"]

        assert model.fs.state_block_ideal_l[0].config.defined_state
        assert not hasattr(model.fs.state_block_ideal_l[0], "eq_phase_equilibrium")
        assert not hasattr(model.fs.state_block_ideal_l[0], "eq_activity_coeff")
        assert not hasattr(model.fs.state_block_ideal_l[0], "eq_mol_frac_out")

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.fs.state_block_ideal_l[0]) == 0

    @pytest.mark.component
    def test_units_consistent(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.state_block_ideal_l.initialize()

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert value(
            model.fs.state_block_ideal_l[0].mole_frac_phase_comp["Liq", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)
        assert value(
            model.fs.state_block_ideal_l[0].mole_frac_phase_comp["Liq", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)


class TestFTPz_V_inlet:
    @pytest.fixture(scope="class")
    def model(self):
        # Create a flowsheet for test
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties_ideal_v = BTXParameterBlock(
            valid_phase="Vap", activity_coeff_model="Ideal", state_vars="FTPz"
        )
        m.fs.state_block_ideal_v = m.fs.properties_ideal_v.build_state_block(
            [0], has_phase_equilibrium=False, defined_state=True
        )

        m.fs.state_block_ideal_v[0].flow_mol.fix(1)
        m.fs.state_block_ideal_v[0].temperature.fix(368)
        m.fs.state_block_ideal_v[0].pressure.fix(101325)
        m.fs.state_block_ideal_v[0].mole_frac_comp["benzene"].fix(0.5)
        m.fs.state_block_ideal_v[0].mole_frac_comp["toluene"].fix(0.5)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert len(model.fs.properties_ideal_v.config) == 4

        assert model.fs.properties_ideal_v.config.valid_phase == "Vap"
        assert len(model.fs.properties_ideal_v.phase_list) == 1
        assert model.fs.properties_ideal_v.phase_list == ["Vap"]

        assert model.fs.state_block_ideal_v[0].config.defined_state
        assert not hasattr(model.fs.state_block_ideal_v[0], "eq_phase_equilibrium")
        assert not hasattr(model.fs.state_block_ideal_v[0], "eq_activity_coeff")
        assert not hasattr(model.fs.state_block_ideal_v[0], "eq_mol_frac_out")

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.fs.state_block_ideal_v[0]) == 0

    @pytest.mark.component
    def test_units_consistent(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.state_block_ideal_v.initialize()

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert value(
            model.fs.state_block_ideal_v[0].mole_frac_phase_comp["Vap", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)
        assert value(
            model.fs.state_block_ideal_v[0].mole_frac_phase_comp["Vap", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)


class TestFTPz_LV_outlet:
    @pytest.fixture(scope="class")
    def model(self):
        # Create a flowsheet for test
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties_ideal_vl = BTXParameterBlock(
            valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
        )
        m.fs.state_block_ideal_vl = m.fs.properties_ideal_vl.build_state_block(
            [0], defined_state=False
        )

        m.fs.state_block_ideal_vl[0].flow_mol.fix(1)
        m.fs.state_block_ideal_vl[0].temperature.fix(368)
        m.fs.state_block_ideal_vl[0].pressure.fix(101325)
        m.fs.state_block_ideal_vl[0].mole_frac_comp["benzene"].fix(0.5)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert len(model.fs.properties_ideal_vl.config) == 4

        assert model.fs.properties_ideal_vl.config.valid_phase == ("Liq", "Vap")
        assert len(model.fs.properties_ideal_vl.phase_list) == 2
        assert model.fs.properties_ideal_vl.phase_list == ["Liq", "Vap"]

        assert not model.fs.state_block_ideal_vl[0].config.defined_state
        assert hasattr(model.fs.state_block_ideal_vl[0], "eq_phase_equilibrium")
        assert not hasattr(model.fs.state_block_ideal_vl[0], "eq_activity_coeff")
        assert hasattr(model.fs.state_block_ideal_vl[0], "eq_mol_frac_out")

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.fs.state_block_ideal_vl[0]) == 0

    @pytest.mark.component
    def test_units_consistent(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.state_block_ideal_vl.initialize()

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert value(
            model.fs.state_block_ideal_vl[0].mole_frac_phase_comp["Liq", "benzene"]
        ) == pytest.approx(0.4121, abs=1e-3)
        assert value(
            model.fs.state_block_ideal_vl[0].mole_frac_phase_comp["Vap", "benzene"]
        ) == pytest.approx(0.6339, abs=1e-3)


class TestFTPz_L_outlet:
    @pytest.fixture(scope="class")
    def model(self):
        # Create a flowsheet for test
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties_ideal_l = BTXParameterBlock(
            valid_phase="Liq", activity_coeff_model="Ideal", state_vars="FTPz"
        )
        m.fs.state_block_ideal_l = m.fs.properties_ideal_l.build_state_block(
            [0], has_phase_equilibrium=False, defined_state=False
        )

        m.fs.state_block_ideal_l[0].flow_mol.fix(1)
        m.fs.state_block_ideal_l[0].temperature.fix(368)
        m.fs.state_block_ideal_l[0].pressure.fix(101325)
        m.fs.state_block_ideal_l[0].mole_frac_comp["benzene"].fix(0.5)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert len(model.fs.properties_ideal_l.config) == 4

        assert model.fs.properties_ideal_l.config.valid_phase == "Liq"
        assert len(model.fs.properties_ideal_l.phase_list) == 1
        assert model.fs.properties_ideal_l.phase_list == ["Liq"]

        assert not model.fs.state_block_ideal_l[0].config.defined_state
        assert not hasattr(model.fs.state_block_ideal_l[0], "eq_phase_equilibrium")
        assert not hasattr(model.fs.state_block_ideal_l[0], "eq_activity_coeff")
        assert hasattr(model.fs.state_block_ideal_l[0], "eq_mol_frac_out")

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.fs.state_block_ideal_l[0]) == 0

    @pytest.mark.component
    def test_units_consistent(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.state_block_ideal_l.initialize()

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert value(
            model.fs.state_block_ideal_l[0].mole_frac_phase_comp["Liq", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)
        assert value(
            model.fs.state_block_ideal_l[0].mole_frac_phase_comp["Liq", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)


class TestFTPz_V_outlet:
    @pytest.fixture(scope="class")
    def model(self):
        # Create a flowsheet for test
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties_ideal_v = BTXParameterBlock(
            valid_phase="Vap", activity_coeff_model="Ideal", state_vars="FTPz"
        )
        m.fs.state_block_ideal_v = m.fs.properties_ideal_v.build_state_block(
            [0], has_phase_equilibrium=False, defined_state=False
        )

        m.fs.state_block_ideal_v[0].flow_mol.fix(1)
        m.fs.state_block_ideal_v[0].temperature.fix(368)
        m.fs.state_block_ideal_v[0].pressure.fix(101325)
        m.fs.state_block_ideal_v[0].mole_frac_comp["benzene"].fix(0.5)

        return m

    @pytest.mark.unit
    def test_build(self, model):
        assert len(model.fs.properties_ideal_v.config) == 4

        assert model.fs.properties_ideal_v.config.valid_phase == "Vap"
        assert len(model.fs.properties_ideal_v.phase_list) == 1
        assert model.fs.properties_ideal_v.phase_list == ["Vap"]

        assert not model.fs.state_block_ideal_v[0].config.defined_state
        assert not hasattr(model.fs.state_block_ideal_v[0], "eq_phase_equilibrium")
        assert not hasattr(model.fs.state_block_ideal_v[0], "eq_activity_coeff")
        assert hasattr(model.fs.state_block_ideal_v[0], "eq_mol_frac_out")

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.fs.state_block_ideal_v[0]) == 0

    @pytest.mark.component
    def test_units_consistent(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.fs.state_block_ideal_v.initialize()

        assert degrees_of_freedom(model) == 0

        fin_fixed_vars = fixed_variables_set(model)
        fin_act_consts = activated_constraints_set(model)

        assert len(fin_act_consts) == len(orig_act_consts)
        assert len(fin_fixed_vars) == len(orig_fixed_vars)

        for c in fin_act_consts:
            assert c in orig_act_consts
        for v in fin_fixed_vars:
            assert v in orig_fixed_vars

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solution(self, model):
        assert value(
            model.fs.state_block_ideal_v[0].mole_frac_phase_comp["Vap", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)
        assert value(
            model.fs.state_block_ideal_v[0].mole_frac_phase_comp["Vap", "benzene"]
        ) == pytest.approx(0.5, abs=1e-3)
