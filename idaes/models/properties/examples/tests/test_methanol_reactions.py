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
from pyomo.environ import (check_optimal_termination,
                           Block,
                           ConcreteModel,
                           Set,
                           value)
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

from idaes.models.properties.examples.methanol_ideal_VLE import \
    config_dict as vap_config_dict
from idaes.models.properties.examples.methanol_reactions import \
    config_dict as rxn_config_dict


# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


class TestParamBlock(object):
    @pytest.mark.unit
    def test_build(self):
        model = ConcreteModel()
        model.thermo_params = GenericParameterBlock(default=vap_config_dict)

        model.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.thermo_params,
                     **rxn_config_dict}
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
        assert len(model.rxn_params.rate_reaction_stoichiometry) == 5
        for k, v in model.rxn_params.rate_reaction_stoichiometry.items():
            if (k[1], k[2]) in rate_config[k[0]].stoichiometry.keys():
                assert v == rate_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert isinstance(model.rxn_params.equilibrium_reaction_stoichiometry, dict)
        assert len(model.rxn_params.equilibrium_reaction_stoichiometry) == 5
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
        model.thermo_params = GenericParameterBlock(default=vap_config_dict)

        model.rxn_params = GenericReactionParameterBlock(
            default={"property_package": model.thermo_params,
                     **rxn_config_dict}
        )

        model.props = model.thermo_params.build_state_block(
            [1], default={"defined_state": True}
        )

        model.rxns = model.rxn_params.build_reaction_block(
            [1], default={"state_block": model.props, "has_equilibrium": False}
        )

        assert_units_consistent(model)

        return model

    @pytest.mark.unit
    def test_dof(self, model):
        # Fix state
        model.props[1].flow_mol.fix(415.44)
        model.props[1].enth_mol.fix(-1.3643e5)
        model.props[1].pressure.fix(5.1e6)
        model.props[1].mole_frac_comp["CH4"].fix(2.2964e-6)
        model.props[1].mole_frac_comp["CO"].fix(0.11438)
        model.props[1].mole_frac_comp["H2"].fix(0.23743)
        model.props[1].mole_frac_comp["CH3OH"].fix(0.64818)

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
        assert model.props[1].flow_mol.value == pytest.approx(
            415.44, rel=1e-4
        )

        assert model.props[1].mole_frac_comp["CH4"].value == pytest.approx(
            2.2964e-6, rel=1e-4
        )
        assert model.props[1].mole_frac_comp["CO"].value == pytest.approx(
            0.11438, rel=1e-4
        )
        assert model.props[1].mole_frac_comp["H2"].value == pytest.approx(
            0.23743, rel=1e-4
        )
        assert model.props[1].mole_frac_comp["CH3OH"].value == pytest.approx(
            0.64818, rel=1e-4
        )

        assert value(model.rxns[1].k_rxn["R1"]) == pytest.approx(1.0, rel=1e-4)

        assert value(model.rxns[1].reaction_rate["R1"]) == pytest.approx(
            value(
                model.rxns[1].k_rxn["R1"]
                * model.props[1].mole_frac_comp["CO"] ** 1
                * model.props[1].mole_frac_comp["H2"] ** 2
            ),
            rel=1e-6,
        )

        assert value(model.rxns[1].k_eq["R2"]) == pytest.approx(71.4505, rel=1e-4)

        assert value(model.rxns[1].k_eq["R2"]) == pytest.approx(
            value(
                model.props[1].mole_frac_comp["CH3OH"] ** 1
                * model.props[1].mole_frac_comp["CO"] ** -1
                * model.props[1].mole_frac_comp["H2"] ** -2
            ),
            rel=1e-6,
        )

    @pytest.mark.ui
    @pytest.mark.unit
    def test_report(self, model):
        model.rxns[1].report()
