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
from pyomo.environ import (ConcreteModel,
                           SolverStatus,
                           TerminationCondition,
                           value,
                           units as pyunits)
from pyomo.util.check_units import assert_units_consistent

from idaes.core.util.model_statistics import (degrees_of_freedom,
                                              fixed_variables_set,
                                              activated_constraints_set)
from idaes.core.util import get_solver

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterBlock)


from idaes.power_generation.carbon_capture.mea_solvent_system.properties.MEA_solvent \
    import configuration


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.params = GenericParameterBlock(default=configuration)

        model.props = model.params.build_state_block(
                [1],
                default={"defined_state": True})

        model.props[1].calculate_scaling_factors()

        # Fix state
        # Initial test conditions: 0.2 kg MEA/kg H2O, 0.2 mol CO2/mol MEA
        # 1 kg H20 = 55.5mol
        # 0.2 kg MEA/kg H2O = 0.2 kg MEA = 3.27 mol
        # 0.2 mol CO2/mol MEA = 0.655 mol
        model.props[1].flow_mol.fix(1)
        model.props[1].temperature.fix(298.15)
        model.props[1].pressure.fix(101325)
        model.props[1].mole_frac_comp["H2O"].fix(0.93388)
        model.props[1].mole_frac_comp["MEA"].fix(0.05510)
        model.props[1].mole_frac_comp["CO2"].fix(0.01102)

        # Set values of mole_frac_phase_comp to get correct property values
        model.props[1].mole_frac_phase_comp["Liq", "H2O"].set_value(0.93388)
        model.props[1].mole_frac_phase_comp["Liq", "MEA"].set_value(0.05510)
        model.props[1].mole_frac_phase_comp["Liq", "CO2"].set_value(0.01102)

        return model

    @pytest.mark.unit
    def test_dof(self, model):
        assert degrees_of_freedom(model.props[1]) == 0

    @pytest.mark.unit
    def test_properties(self, model):
        model.props[1].visc_d_phase.display()

        assert pytest.approx(79.5139, rel=1e-5) == value(
            model.props[1].cp_mol_phase["Liq"])
        assert pytest.approx(-84000*0.01102, rel=1e-5) == value(
            model.props[1].enth_mol_phase["Liq"])
        assert pytest.approx(2988.034, rel=1e-5) == value(
            model.props[1].henry["Liq", "CO2"])
        assert pytest.approx(3185.403, rel=1e-5) == value(
            model.props[1].pressure_sat_comp["H2O"])
        assert pytest.approx(49.2646, rel=1e-5) == value(
            model.props[1].pressure_sat_comp["MEA"])
        assert pytest.approx(2.032875e-05, rel=1e-5) == value(
            model.props[1].vol_mol_phase["Liq"])
        assert pytest.approx(9.495423e-4, rel=1e-5) == value(
            model.props[1].visc_d_phase["Liq"])

    @pytest.mark.component
    def test_unit_consistency(self, model):
        assert_units_consistent(model)

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
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

    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.component
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert results.solver.status == SolverStatus.ok

    # @pytest.mark.skipif(solver is None, reason="Solver not available")
    # @pytest.mark.component
    # def test_solution(self, model):
    #     # Check phase equilibrium results
    #     assert model.props[1].mole_frac_phase_comp["Liq", "benzene"].value == \
    #         pytest.approx(0.4121, abs=1e-4)
    #     assert model.props[1].mole_frac_phase_comp["Vap", "benzene"].value == \
    #         pytest.approx(0.6339, abs=1e-4)
    #     assert model.props[1].phase_frac["Vap"].value == \
    #         pytest.approx(0.3961, abs=1e-4)

    #     assert value(model.props[1].conc_mol_phase_comp["Vap", "benzene"]) == \
    #         pytest.approx(20.9946, abs=1e-4)

