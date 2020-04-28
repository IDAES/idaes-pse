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
Tests for generic reaction package core code

Author: Andrew Lee
"""
import pytest

from pyomo.environ import Block, ConcreteModel, Expression, Set, Var, value

from idaes.generic_models.properties.core.generic.generic_reaction import (
        GenericReactionParameterBlock)
from idaes.generic_models.properties.core.reactions.dh_rxn import \
    constant_dh_rxn

from idaes.core.util.testing import PhysicalParameterTestBlock

from idaes.core.util.exceptions import ConfigurationError


class TestGenericReactionParameterBlock(object):
    def test_rate_build(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        m.rxn_params = GenericReactionParameterBlock(default={
            "property_package": m.params,
            "rate_reactions": {
                "r1": {"stoichiometry": {("p1", "c1"): -1,
                                         ("p1", "c2"): 2},
                       "heat_of_reaction": "foo",
                       "rate_form": "foo"}}})

        rxn_config = m.rxn_params.config.rate_reactions

        assert isinstance(m.rxn_params.rate_reaction_idx, Set)
        assert len(m.rxn_params.rate_reaction_idx) == 1
        assert "r1" in m.rxn_params.rate_reaction_idx

        assert not hasattr(self, "equilibrium_reaction_idx")

        assert isinstance(m.rxn_params.rate_reaction_stoichiometry, dict)
        assert len(m.rxn_params.rate_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.rate_reaction_stoichiometry.items():
            if (k[1], k[2]) in rxn_config[k[0]].stoichiometry.keys():
                assert v == rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "equilibrium_reaction_stoichiometry")

        assert isinstance(m.rxn_params.reaction_idx, Set)
        assert m.rxn_params.reaction_idx == m.rxn_params.rate_reaction_idx

        assert isinstance(m.rxn_params.reaction_r1, Block)

    def test_rate_build_no_stoichiometry(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params rate reaction r1 was not "
                           "provided with a stoichiometry configuration "
                           "argument."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "rate_reactions": {
                    "r1": {"heat_of_reaction": "foo",
                           "rate_form": "foo"}}})

    def test_rate_build_invalid_phase_stoichiometry(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params stoichiometry for rate reaction "
                           "r1 included unrecognised phase p7."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "rate_reactions": {
                    "r1": {"stoichiometry": {("p7", "c1"): -1,
                                             ("p1", "c2"): 2},
                           "heat_of_reaction": "foo",
                           "rate_form": "foo"}}})

    def test_rate_build_invalid_component_stoichiometry(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params stoichiometry for rate reaction "
                           "r1 included unrecognised component c7."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "rate_reactions": {
                    "r1": {"stoichiometry": {("p1", "c7"): -1,
                                             ("p1", "c2"): 2},
                           "heat_of_reaction": "foo",
                           "rate_form": "foo"}}})

    def test_rate_build_no_form(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params rate reaction r1 was not "
                           "provided with a rate_form configuration "
                           "argument."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "rate_reactions": {
                    "r1": {"stoichiometry": {("p1", "c1"): -1,
                                             ("p1", "c2"): 2},
                           "heat_of_reaction": "foo"}}})

    def test_equil_build(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        m.rxn_params = GenericReactionParameterBlock(default={
            "property_package": m.params,
            "equilibrium_reactions": {
                "e1": {"stoichiometry": {("p2", "c1"): -3,
                                         ("p2", "c2"): 4},
                       "heat_of_reaction": "foo",
                       "equilibrium_form": "foo"}}})

        rxn_config = m.rxn_params.config.equilibrium_reactions

        assert isinstance(m.rxn_params.equilibrium_reaction_idx, Set)
        assert len(m.rxn_params.equilibrium_reaction_idx) == 1
        assert "e1" in m.rxn_params.equilibrium_reaction_idx

        assert not hasattr(self, "rate_reaction_idx")

        assert isinstance(m.rxn_params.equilibrium_reaction_stoichiometry,
                          dict)
        assert len(m.rxn_params.equilibrium_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.equilibrium_reaction_stoichiometry.items():
            if (k[1], k[2]) in rxn_config[k[0]].stoichiometry.keys():
                assert v == rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "rate_reaction_stoichiometry")

        assert isinstance(m.rxn_params.reaction_idx, Set)
        assert m.rxn_params.reaction_idx == \
            m.rxn_params.equilibrium_reaction_idx

        assert isinstance(m.rxn_params.reaction_e1, Block)

    def test_equil_build_no_stoichiometry(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params equilibrium reaction e1 was not "
                           "provided with a stoichiometry configuration "
                           "argument."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "equilibrium_reactions": {
                    "e1": {"heat_of_reaction": "foo",
                           "equilibrium_form": "foo"}}})

    def test_equil_build_invalid_phase_stoichiometry(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params stoichiometry for equilibrium "
                           "reaction e1 included unrecognised phase p7."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "equilibrium_reactions": {
                    "e1": {"stoichiometry": {("p7", "c1"): -3,
                                             ("p2", "c2"): 4},
                           "heat_of_reaction": "foo",
                           "equilibrium_form": "foo"}}})

    def test_equil_build_invalid_component_stoichiometry(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params stoichiometry for equilibrium "
                           "reaction e1 included unrecognised component c7."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "equilibrium_reactions": {
                    "e1": {"stoichiometry": {("p2", "c7"): -3,
                                             ("p2", "c2"): 4},
                           "heat_of_reaction": "foo",
                           "equilibrium_form": "foo"}}})

    def test_equil_build_no_form(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        with pytest.raises(ConfigurationError,
                           match="rxn_params equilibrium reaction e1 was not "
                           "provided with a equilibrium_form configuration "
                           "argument."):
            m.rxn_params = GenericReactionParameterBlock(default={
                "property_package": m.params,
                "equilibrium_reactions": {
                    "e1": {"stoichiometry": {("p2", "c1"): -3,
                                             ("p2", "c2"): 4},
                           "heat_of_reaction": "foo"}}})

    def test_rate_and_equil_build(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        m.rxn_params = GenericReactionParameterBlock(default={
            "property_package": m.params,
            "rate_reactions": {
                "r1": {"stoichiometry": {("p1", "c1"): -1,
                                         ("p1", "c2"): 2},
                       "heat_of_reaction": "foo",
                       "rate_form": "foo"}},
            "equilibrium_reactions": {
                "e1": {"stoichiometry": {("p2", "c1"): -3,
                                         ("p2", "c2"): 4},
                       "heat_of_reaction": "foo",
                       "equilibrium_form": "foo"}}})

        r_rxn_config = m.rxn_params.config.rate_reactions

        assert isinstance(m.rxn_params.rate_reaction_idx, Set)
        assert len(m.rxn_params.rate_reaction_idx) == 1
        assert "r1" in m.rxn_params.rate_reaction_idx

        assert not hasattr(self, "equilibrium_reaction_idx")

        assert isinstance(m.rxn_params.rate_reaction_stoichiometry, dict)
        assert len(m.rxn_params.rate_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.rate_reaction_stoichiometry.items():
            if (k[1], k[2]) in r_rxn_config[k[0]].stoichiometry.keys():
                assert v == r_rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "equilibrium_reaction_stoichiometry")

        e_rxn_config = m.rxn_params.config.equilibrium_reactions

        assert isinstance(m.rxn_params.equilibrium_reaction_idx, Set)
        assert len(m.rxn_params.equilibrium_reaction_idx) == 1
        assert "e1" in m.rxn_params.equilibrium_reaction_idx

        assert not hasattr(self, "rate_reaction_idx")

        assert isinstance(m.rxn_params.equilibrium_reaction_stoichiometry,
                          dict)
        assert len(m.rxn_params.equilibrium_reaction_stoichiometry) == 4
        for k, v in m.rxn_params.equilibrium_reaction_stoichiometry.items():
            if (k[1], k[2]) in e_rxn_config[k[0]].stoichiometry.keys():
                assert v == e_rxn_config[k[0]].stoichiometry[k[1], k[2]]
            else:
                assert v == 0

        assert not hasattr(self, "rate_reaction_stoichiometry")

        assert isinstance(m.rxn_params.reaction_idx, Set)
        assert m.rxn_params.reaction_idx == (
            m.rxn_params.rate_reaction_idx |
            m.rxn_params.equilibrium_reaction_idx)
        assert len(m.rxn_params.reaction_idx) == 2

        assert isinstance(m.rxn_params.reaction_r1, Block)
        assert isinstance(m.rxn_params.reaction_e1, Block)

    def test_build_parameters(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()

        m.rxn_params = GenericReactionParameterBlock(default={
            "property_package": m.params,
            "rate_reactions": {
                "r1": {"stoichiometry": {("p1", "c1"): -1,
                                         ("p1", "c2"): 2},
                       "heat_of_reaction": constant_dh_rxn,
                       "rate_form": "foo",
                       "parameter_data": {
                           "dh_rxn_ref": -10000}}},
            "equilibrium_reactions": {
                "e1": {"stoichiometry": {("p2", "c1"): -3,
                                         ("p2", "c2"): 4},
                       "heat_of_reaction": constant_dh_rxn,
                       "equilibrium_form": "foo",
                       "parameter_data": {
                           "dh_rxn_ref": -20000}}}})

        assert isinstance(m.rxn_params.reaction_r1.dh_rxn_ref, Var)
        assert m.rxn_params.reaction_r1.dh_rxn_ref.fixed
        assert m.rxn_params.reaction_r1.dh_rxn_ref.value == -10000
        assert isinstance(m.rxn_params.reaction_e1.dh_rxn_ref, Var)
        assert m.rxn_params.reaction_e1.dh_rxn_ref.fixed
        assert m.rxn_params.reaction_e1.dh_rxn_ref.value == -20000


# -----------------------------------------------------------------------------
class TestGenericReactionBlock(object):
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        # Add a dummy thermo package for validation
        m.params = PhysicalParameterTestBlock()
        m.sblock = m.params.build_state_block([1])

        m.rxn_params = GenericReactionParameterBlock(default={
            "property_package": m.params,
            "rate_reactions": {
                "r1": {"stoichiometry": {("p1", "c1"): -1,
                                         ("p1", "c2"): 2},
                       "heat_of_reaction": constant_dh_rxn,
                       "rate_form": "foo",
                       "parameter_data": {
                           "dh_rxn_ref": -10000}}},
            "equilibrium_reactions": {
                "e1": {"stoichiometry": {("p2", "c1"): -3,
                                         ("p2", "c2"): 4},
                       "heat_of_reaction": constant_dh_rxn,
                       "equilibrium_form": "foo",
                       "parameter_data": {
                           "dh_rxn_ref": -20000}}}})

        m.rblock = m.rxn_params.build_reaction_block(
            [1], default={"state_block": m.sblock})

        return m

    def test_dh_rxn(self, model):
        assert isinstance(model.rblock[1].dh_rxn, Expression)
        assert len(model.rblock[1].dh_rxn) == 2
        model.rblock[1].dh_rxn.display()
        assert value(model.rblock[1].dh_rxn["r1"]) == -10000
        assert value(model.rblock[1].dh_rxn["e1"]) == -20000
