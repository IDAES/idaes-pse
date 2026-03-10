#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for saponification property package example.
Authors: Andrew Lee
"""

from math import exp
import pytest

from pyomo.environ import ConcreteModel, Constraint, Param, Suffix, units, value, Var

from idaes.core import MaterialFlowBasis
from idaes.models.properties.examples.saponification_reactions import (
    SaponificationReactionParameterBlock,
    ReactionBlock,
    SaponificationReactionScaler,
)
from idaes.models.properties.examples.saponification_thermo import (
    SaponificationParameterBlock,
)
from idaes.core.solvers import get_solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver("ipopt_v2")


class TestParamBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = SaponificationParameterBlock()
        model.rparams = SaponificationReactionParameterBlock(
            property_package=model.pparams
        )

        return model

    @pytest.mark.unit
    def test_config(self, model):
        assert len(model.rparams.config) == 2

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rparams.reaction_block_class is ReactionBlock

        assert len(model.rparams.rate_reaction_idx) == 1
        for i in model.rparams.rate_reaction_idx:
            assert i == "R1"

        assert len(model.rparams.rate_reaction_stoichiometry) == 5
        for i in model.rparams.rate_reaction_stoichiometry:
            assert i in [
                ("R1", "Liq", "NaOH"),
                ("R1", "Liq", "EthylAcetate"),
                ("R1", "Liq", "SodiumAcetate"),
                ("R1", "Liq", "Ethanol"),
                ("R1", "Liq", "H2O"),
            ]

        assert isinstance(model.rparams.arrhenius, Param)
        assert value(model.rparams.arrhenius) == 3.132e6

        assert isinstance(model.rparams.energy_activation, Param)
        assert value(model.rparams.energy_activation) == 43000

        assert isinstance(model.rparams.dh_rxn, Param)
        assert len(model.rparams.dh_rxn) == 1
        for i in model.rparams.dh_rxn:
            assert value(model.rparams.dh_rxn[i]) == -49000


class TestReactionBlock(object):
    @pytest.fixture(scope="class")
    def model(self):
        model = ConcreteModel()
        model.pparams = SaponificationParameterBlock()
        model.rparams = SaponificationReactionParameterBlock(
            property_package=model.pparams
        )

        model.props = model.pparams.build_state_block([1])

        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        return model

    @pytest.mark.unit
    def test_build(self, model):
        assert model.rxns[1].conc_mol_comp_ref is model.props[1].conc_mol_comp
        assert model.rxns[1].temperature_ref is model.props[1].temperature
        assert model.rxns[1].dh_rxn is model.rparams.dh_rxn

        assert model.rxns[1].default_scaler is SaponificationReactionScaler

    @pytest.mark.unit
    def test_rate_constant(self, model):
        assert isinstance(model.rxns[1].k_rxn, Var)
        assert isinstance(model.rxns[1].arrhenius_eqn, Constraint)

    @pytest.mark.unit
    def test_rxn_rate(self, model):
        assert isinstance(model.rxns[1].reaction_rate, Var)
        assert isinstance(model.rxns[1].rate_expression, Constraint)

    @pytest.mark.unit
    def test_get_reaction_rate_basis(self, model):
        assert model.rxns[1].get_reaction_rate_basis() == MaterialFlowBasis.molar

    @pytest.mark.unit
    def test_model_check(self, model):
        assert model.rxns[1].model_check() is None

    @pytest.mark.unit
    def test_initialize(self, model):
        assert model.rxns.initialize(outlvl=1) is None

    def check_units(self, model):
        units.assert_units_consistent(model)


class TestSaponificationReactionScaler(object):
    @pytest.mark.unit
    def test_variable_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = SaponificationParameterBlock()
        model.rparams = SaponificationReactionParameterBlock(
            property_package=model.pparams
        )

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, SaponificationReactionScaler)

        scaler.variable_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 2
        assert sfx[model.rxns[1].k_rxn] == pytest.approx(
            1 / (3.132e6 * exp(-43000 / (8.31446262 * 298.15))), rel=1e-8
        )
        assert sfx[model.rxns[1].reaction_rate["R1"]] == pytest.approx(1e2, rel=1e-8)

    @pytest.mark.unit
    def test_constraint_scaling_routine(self):
        model = ConcreteModel()
        model.pparams = SaponificationParameterBlock()
        model.rparams = SaponificationReactionParameterBlock(
            property_package=model.pparams
        )

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, SaponificationReactionScaler)

        scaler.constraint_scaling_routine(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        assert len(sfx) == 2
        assert sfx[model.rxns[1].arrhenius_eqn] == pytest.approx(1, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R1"]] == pytest.approx(1e-4, rel=1e-8)

    @pytest.mark.unit
    def test_scale_model(self):
        model = ConcreteModel()
        model.pparams = SaponificationParameterBlock()
        model.rparams = SaponificationReactionParameterBlock(
            property_package=model.pparams
        )

        model.props = model.pparams.build_state_block([1])
        model.rxns = model.rparams.build_reaction_block([1], state_block=model.props)

        # Trigger build of reaction properties
        model.rxns[1].reaction_rate

        scaler = model.rxns[1].default_scaler()
        assert isinstance(scaler, SaponificationReactionScaler)

        scaler.scale_model(model.rxns[1])

        assert isinstance(model.rxns[1].scaling_factor, Suffix)

        sfx = model.rxns[1].scaling_factor
        k_rxn_sf = 1 / (3.132e6 * exp(-43000 / (8.31446262 * 298.15)))
        assert len(sfx) == 4
        assert sfx[model.rxns[1].k_rxn] == pytest.approx(k_rxn_sf, rel=1e-8)
        assert sfx[model.rxns[1].reaction_rate["R1"]] == pytest.approx(1e2, rel=1e-8)
        assert sfx[model.rxns[1].arrhenius_eqn] == pytest.approx(k_rxn_sf, rel=1e-8)
        assert sfx[model.rxns[1].rate_expression["R1"]] == pytest.approx(
            1e-4 * k_rxn_sf, rel=1e-8
        )
