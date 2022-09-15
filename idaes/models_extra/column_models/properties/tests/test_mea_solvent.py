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

from idaes.models_extra.column_models.properties.MEA_solvent import configuration


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(**configuration)

        model.props = model.params.build_state_block([1], defined_state=True)

        model.props[1].calculate_scaling_factors()

        # Fix state
        model.props[1].flow_mol.fix(83.89)
        model.props[1].temperature.fix(392.5)
        model.props[1].pressure.fix(183700)
        model.props[1].mole_frac_comp["CO2"].fix(0.0326)
        model.props[1].mole_frac_comp["H2O"].fix(0.8589)
        model.props[1].mole_frac_comp["MEA"].fix(0.1085)

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
        assert value(model.props[1].log_k_eq["carbamate"]) == pytest.approx(
            -1.98499794, rel=1e-8
        )
        assert value(model.props[1].log_k_eq["bicarbonate"]) == pytest.approx(
            -7.5776101, rel=1e-8
        )

        assert pytest.approx(1346.056734, rel=1e-8) == value(
            model.props[1].conc_mol_phase_comp_apparent["Liq", "CO2"]
        )
        assert pytest.approx(35464.0530, rel=1e-8) == value(
            model.props[1].conc_mol_phase_comp_apparent["Liq", "H2O"]
        )
        assert pytest.approx(4479.97410, rel=1e-8) == value(
            model.props[1].conc_mol_phase_comp_apparent["Liq", "MEA"]
        )

        assert pytest.approx(92.0060276, rel=1e-8) == value(
            model.props[1].cp_mol_phase["Liq"]
        )

        assert pytest.approx(8842.75903, rel=1e-8) == value(
            model.props[1].henry["Liq", "CO2"]
        )

        assert pytest.approx(194455.113, rel=1e-8) == value(
            model.props[1].pressure_sat_comp["H2O"]
        )
        assert pytest.approx(15510.7449, rel=1e-8) == value(
            model.props[1].pressure_sat_comp["MEA"]
        )

        assert pytest.approx(2.42188900e-05, rel=1e-8) == value(
            model.props[1].vol_mol_phase["Liq"]
        )

        assert pytest.approx(5.58038484e-4, rel=1e-8) == value(
            model.props[1].visc_d_phase["Liq"]
        )

        assert pytest.approx(0.657512692, rel=1e-8) == value(
            model.props[1].mass_frac_phase_comp_apparent["Liq", "H2O"]
        )
        assert pytest.approx(0.281537026, rel=1e-8) == value(
            model.props[1].mass_frac_phase_comp_apparent["Liq", "MEA"]
        )
        assert pytest.approx(6.095028209e-2, rel=1e-8) == value(
            model.props[1].mass_frac_phase_comp_apparent["Liq", "CO2"]
        )

        assert pytest.approx(0.404942934, rel=1e-8) == value(
            model.props[1].therm_cond_phase["Liq"]
        )

        assert pytest.approx(0.0533310121, rel=1e-8) == value(
            model.props[1].surf_tens_phase["Liq"]
        )

        assert pytest.approx(8.2258789e-9, rel=1e-8) == value(
            model.props[1].diffus_phase_comp_true["Liq", "CO2"]
        )
        assert pytest.approx(4.47017415e-09, rel=1e-8) == value(
            model.props[1].diffus_phase_comp_true["Liq", "MEA"]
        )
        assert pytest.approx(2.17984326e-09, rel=1e-8) == value(
            model.props[1].diffus_phase_comp_true["Liq", "MEACOO_-"]
        )
        assert pytest.approx(2.17984326e-09, rel=1e-8) == value(
            model.props[1].diffus_phase_comp_true["Liq", "MEA_+"]
        )
        assert pytest.approx(-41250.348, rel=1e-8) == value(
            model.props[1].enth_mol_phase_comp["Liq", "MEA"]
        )
        assert pytest.approx(-36850.713, rel=1e-8) == value(
            model.props[1].enth_mol_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(-83998.005, rel=1e-8) == value(
            model.props[1].enth_mol_phase_comp["Liq", "CO2"]
        )
