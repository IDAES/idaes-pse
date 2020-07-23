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
import types

from pyomo.environ import Block, ConcreteModel, Var, units as pyunits, value
from pyomo.common.config import ConfigBlock
from pyomo.util.check_units import assert_units_equivalent

from idaes.generic_models.properties.core.reactions.equilibrium_constant import *

from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.misc import add_object_reference
from idaes.core.util.constants import Constants as c
from idaes.core.property_meta import PropertyClassMetadata


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

    m.meta_object = PropertyClassMetadata()
    m.meta_object.default_units["temperature"] = pyunits.K
    m.meta_object.default_units["mass"] = pyunits.kg
    m.meta_object.default_units["length"] = pyunits.m
    m.meta_object.default_units["time"] = pyunits.s
    m.meta_object.default_units["amount"] = pyunits.mol

    def get_metadata(self):
        return m.meta_object
    m.rparams.get_metadata = types.MethodType(get_metadata, m.rparams)

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].dh_rxn = Var(["r1"],
                          initialize=1,
                          units=pyunits.J/pyunits.mol)

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
        model.rxn[1], model.rparams.reaction_r1, "r1", 300*pyunits.K)

    assert value(rform) == pytest.approx(0.99984, rel=1e-3)
    assert_units_equivalent(rform, None)
