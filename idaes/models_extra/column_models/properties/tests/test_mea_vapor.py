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
Author: Andrew Lee
"""
import pytest
from pyomo.environ import check_optimal_termination, ConcreteModel, value
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    fixed_variables_set,
    activated_constraints_set,
)
from idaes.core.solvers import get_solver

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)


from idaes.models_extra.column_models.properties.MEA_vapor import flue_gas, wet_co2


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestFlueGasBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**flue_gas)

        model.props = model.params.build_state_block([1], defined_state=True)

        model.props[1].calculate_scaling_factors()

        # Fix state
        model.props[1].flow_mol.fix(21.48)
        model.props[1].temperature.fix(317.88)
        model.props[1].pressure.fix(107650)
        model.props[1].mole_frac_comp["CO2"].fix(0.11453)
        model.props[1].mole_frac_comp["H2O"].fix(0.08526)
        model.props[1].mole_frac_comp["N2"].fix(0.73821)
        model.props[1].mole_frac_comp["O2"].fix(0.06200)

        return model

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})

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
        # Values verified with older model
        assert pytest.approx(30.6689467, rel=1e-8) == value(
            model.props[1].cp_mol_phase["Vap"]
        )

        assert pytest.approx(1.77021548e-05, rel=1e-8) == value(
            model.props[1].visc_d_phase["Vap"]
        )

        assert pytest.approx(2.50126426e-2, rel=1e-8) == value(
            model.props[1].therm_cond_phase["Vap"]
        )

        assert pytest.approx(1.75238546e-05, rel=1e-8) == value(
            model.props[1].diffus_phase_comp["Vap", "CO2"]
        )
        assert pytest.approx(2.64401213e-05, rel=1e-8) == value(
            model.props[1].diffus_phase_comp["Vap", "H2O"]
        )
        assert pytest.approx(2.06616153e-05, rel=1e-8) == value(
            model.props[1].diffus_phase_comp["Vap", "N2"]
        )
        assert pytest.approx(2.14369756e-05, rel=1e-8) == value(
            model.props[1].diffus_phase_comp["Vap", "O2"]
        )


class TestWetCO2Block(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**wet_co2)

        model.props = model.params.build_state_block([1], defined_state=True)

        model.props[1].calculate_scaling_factors()

        # Fix state
        model.props[1].flow_mol.fix(17.496)
        model.props[1].temperature.fix(396.6)
        model.props[1].pressure.fix(183430)
        model.props[1].mole_frac_comp["CO2"].fix(0.0145)
        model.props[1].mole_frac_comp["H2O"].fix(0.9855)

        return model

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})

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
        # Values verified with older model
        assert pytest.approx(34.3944239, rel=1e-8) == value(
            model.props[1].cp_mol_phase["Vap"]
        )

        assert pytest.approx(1.36007669e-05, rel=1e-8) == value(
            model.props[1].visc_d_phase["Vap"]
        )

        assert pytest.approx(2.65105253e-2, rel=1e-8) == value(
            model.props[1].therm_cond_phase["Vap"]
        )

        assert pytest.approx(1.90464661e-05, rel=1e-8) == value(
            model.props[1].diffus_phase_comp["Vap", "CO2"]
        )
        assert pytest.approx(1.90464661e-05, rel=1e-8) == value(
            model.props[1].diffus_phase_comp["Vap", "H2O"]
        )
