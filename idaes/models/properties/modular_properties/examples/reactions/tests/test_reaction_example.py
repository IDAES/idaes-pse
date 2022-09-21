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
from pyomo.environ import check_optimal_termination, Block, ConcreteModel, Set, value
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
from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
)

from idaes.models.properties.modular_properties.examples.reactions.reaction_example import (
    thermo_configuration,
    rxn_configuration,
)


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.thermo_params = GenericParameterBlock(**thermo_configuration)

        model.rxn_params = GenericReactionParameterBlock(
            property_package=model.thermo_params, **rxn_configuration
        )

        rate_config = model.rxn_params.config.rate_reactions
        equil_config = model.rxn_params.config.equilibrium_reactions

        assert isinstance(model.rxn_params.rate_reaction_idx, Set)
        assert len(model.rxn_params.rate_reaction_idx) == 1
        assert "R1" in model.rxn_params.rate_reaction_idx

        assert isinstance(model.rxn_params.equilibrium_reaction_idx, Set)
        assert len(model.rxn_params.equilibrium_reaction_idx) == 1
        assert "R2" in model.rxn_params.equilibrium_reaction_idx

        assert isinstance(model.rxn_params.rate_reaction_stoichiometry, dict)
        assert len(model.rxn_params.rate_reaction_stoichiometry) == 4
        for k, v in model.rxn_params.rate_reaction_stoichiometry.items():
            if (k[1], k[2]) in rate_config[k[0]].stoichiometry.keys():
                assert v == rate_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert isinstance(model.rxn_params.equilibrium_reaction_stoichiometry, dict)
        assert len(model.rxn_params.equilibrium_reaction_stoichiometry) == 4
        for k, v in model.rxn_params.equilibrium_reaction_stoichiometry.items():
            if (k[1], k[2]) in equil_config[k[0]].stoichiometry.keys():
                assert v == equil_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert isinstance(model.rxn_params.reaction_R1, Block)
        assert isinstance(model.rxn_params.reaction_R2, Block)

        assert_units_consistent(model)


class TestStateBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.thermo_params = GenericParameterBlock(**thermo_configuration)

        model.rxn_params = GenericReactionParameterBlock(
            property_package=model.thermo_params, **rxn_configuration
        )

        model.props = model.thermo_params.build_state_block([1], defined_state=True)

        model.rxns = model.rxn_params.build_reaction_block(
            [1], state_block=model.props, has_equilibrium=True
        )

        assert_units_consistent(model)

        return model

    @pytest.mark.unit
    def test_dof(self, model):
        # Fix state
        # Note that flow of D cannot be specified due to equlibrium
        model.props[1].flow_mol_comp["A"].fix(1)
        model.props[1].flow_mol_comp["B"].fix(1)
        model.props[1].flow_mol_comp["C"].fix(1)
        model.props[1].temperature.fix(368)
        model.props[1].pressure.fix(101325)

        assert degrees_of_freedom(model) == 0

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.unit
    def test_initialize(self, model):
        orig_fixed_vars = fixed_variables_set(model)
        orig_act_consts = activated_constraints_set(model)

        model.props.initialize(optarg={"tol": 1e-6})
        model.rxns.initialize(optarg={"tol": 1e-6})

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
    @pytest.mark.unit
    def test_solve(self, model):
        results = solver.solve(model)

        # Check for optimal solution
        assert check_optimal_termination(results)

    @pytest.mark.initialize
    @pytest.mark.solver
    @pytest.mark.skipif(solver is None, reason="Solver not available")
    @pytest.mark.unit
    def test_solution(self, model):
        assert model.props[1].flow_mol_comp["D"].value == pytest.approx(
            7.0849, rel=1e-4
        )

        assert model.props[1].mole_frac_phase_comp["Liq", "A"].value == pytest.approx(
            0.09916, rel=1e-4
        )
        assert model.props[1].mole_frac_phase_comp["Liq", "B"].value == pytest.approx(
            0.09916, rel=1e-4
        )
        assert model.props[1].mole_frac_phase_comp["Liq", "C"].value == pytest.approx(
            0.09916, rel=1e-4
        )
        assert model.props[1].mole_frac_phase_comp["Liq", "D"].value == pytest.approx(
            0.70253, rel=1e-4
        )

        assert value(model.rxns[1].k_rxn["R1"]) == pytest.approx(0.72121, rel=1e-4)

        assert value(model.rxns[1].reaction_rate["R1"]) == pytest.approx(
            value(
                model.rxns[1].k_rxn["R1"]
                * model.props[1].mole_frac_phase_comp["Liq", "A"] ** 1
                * model.props[1].mole_frac_phase_comp["Liq", "B"] ** 1
            ),
            rel=1e-6,
        )

        assert value(model.rxns[1].k_eq["R2"]) == pytest.approx(71.4505, rel=1e-4)

        assert value(model.rxns[1].k_eq["R2"]) == pytest.approx(
            value(
                model.props[1].mole_frac_phase_comp["Liq", "D"] ** 1
                * model.props[1].mole_frac_phase_comp["Liq", "B"] ** -1
                * model.props[1].mole_frac_phase_comp["Liq", "C"] ** -1
            ),
            rel=1e-6,
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.rxns[1].report()
