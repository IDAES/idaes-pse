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

from idaes.generic_models.properties.core.reactions.equilibrium_constant import *

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
    m.rparams.config.equilibrium_reactions = ConfigBlock(implicit=True)
    m.rparams.config.equilibrium_reactions.r1 = ConfigBlock(implicit=True)
    m.rparams.config.equilibrium_reactions.r1 = {
        "stoichiometry": {("p1", "c1"): -1,
                          ("p1", "c2"): 2},
        "parameter_data": {}}

    m.rparams.reaction_r1 = Block()

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].dh_rxn = Var(["r1"], initialize=1)

    return m


@pytest.mark.unit
def test_van_t_hoff(model):
    model.rparams.config.equilibrium_reactions.r1.parameter_data = {
        "k_eq_ref": 1,
        "T_eq_ref": 500}

    van_t_hoff.build_parameters(
        model.rparams.reaction_r1,
        model.rparams.config.equilibrium_reactions["r1"])

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.k_eq_ref, Var)
    assert model.rparams.reaction_r1.k_eq_ref.value == 1

    assert isinstance(model.rparams.reaction_r1.T_eq_ref, Var)
    assert model.rparams.reaction_r1.T_eq_ref.value == 500

    # Check expression
    rform = van_t_hoff.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        model.rparams.reaction_r1.k_eq_ref *
        exp(-(model.rxn[1].dh_rxn["r1"]/c.gas_constant) *
            (1/300 - 1/model.rparams.reaction_r1.T_eq_ref)))
