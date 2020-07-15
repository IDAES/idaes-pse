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
Tests for rate forms
"""

import pytest

from pyomo.environ import Block, ConcreteModel, Var
from pyomo.common.config import ConfigBlock

from idaes.generic_models.properties.core.reactions.dh_rxn import *

from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.misc import add_object_reference


@pytest.fixture
def model():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()
    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = Block()

    m.rparams.config = ConfigBlock(implicit=True)

    m.rparams.config.property_package = m.pparams
    m.rparams.config.rate_reactions = ConfigBlock(implicit=True)
    m.rparams.config.rate_reactions.r1 = ConfigBlock(implicit=True)
    m.rparams.config.rate_reactions.r1 = {
        "stoichiometry": {("p1", "c1"): -1,
                          ("p1", "c2"): 2},
        "parameter_data": {}}
    m.rparams.config.equilibrium_reactions = ConfigBlock(implicit=True)
    m.rparams.config.equilibrium_reactions.e1 = ConfigBlock(implicit=True)
    m.rparams.config.equilibrium_reactions.e1 = {
        "stoichiometry": {("p1", "c1"): -1,
                          ("p1", "c2"): 2},
        "parameter_data": {}}

    m.rparams.reaction_r1 = Block()
    m.rparams.reaction_e1 = Block()

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    return m


@pytest.mark.unit
def test_constant_dh_rxn(model):
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "dh_rxn_ref": 1}
    model.rparams.config.equilibrium_reactions.e1.parameter_data = {
        "dh_rxn_ref": 10}

    constant_dh_rxn.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"])
    constant_dh_rxn.build_parameters(
        model.rparams.reaction_e1,
        model.rparams.config.equilibrium_reactions["e1"])

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.dh_rxn_ref, Var)
    assert model.rparams.reaction_r1.dh_rxn_ref.value == 1

    assert isinstance(model.rparams.reaction_e1.dh_rxn_ref, Var)
    assert model.rparams.reaction_e1.dh_rxn_ref.value == 10

    # Check expressions
    rform = constant_dh_rxn.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(model.rparams.reaction_r1.dh_rxn_ref)

    rform = constant_dh_rxn.return_expression(
        model.rxn[1], model.rparams.reaction_e1, "e1", 300)

    assert str(rform) == str(model.rparams.reaction_e1.dh_rxn_ref)
