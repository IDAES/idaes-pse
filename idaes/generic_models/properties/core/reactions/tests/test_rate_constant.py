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

from pyomo.environ import Block, ConcreteModel, exp, Var
from pyomo.common.config import ConfigBlock

from idaes.generic_models.properties.core.reactions.rate_constant import *

from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.misc import add_object_reference
from idaes.core.util.constants import Constants as c


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

    m.rparams.reaction_r1 = Block()

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    return m


@pytest.mark.unit
def test_arrhenius(model):
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": 1,
        "energy_activation": 500}

    arrhenius.build_parameters(
        model.rparams.reaction_r1,
        model.rparams.config.rate_reactions["r1"])

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        model.rparams.reaction_r1.arrhenius_const *
        exp(-model.rparams.reaction_r1.energy_activation /
            (c.gas_constant*300)))
