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

from idaes.generic_models.properties.core.reactions.equilibrium_forms import *

from idaes.core import SolidPhase
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.misc import add_object_reference


@pytest.fixture
def model():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = Block()

    m.rparams.config = ConfigBlock(implicit=True)

    m.rparams.config.property_package = m.pparams
    m.rparams.config.equilibrium_reactions = ConfigBlock(implicit=True)
    m.rparams.config.equilibrium_reactions.r1 = ConfigBlock(implicit=True)
    m.rparams.config.equilibrium_reactions.r1 = {
        "stoichiometry": {("p1", "c1"): -1,
                          ("p1", "c2"): 2,
                          ("sol", "c1"): -3,
                          ("sol", "c2"): 4},
        "parameter_data": {}}

    m.rparams.reaction_r1 = Block()

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    return m


@pytest.mark.unit
def test_mole_frac_power_law_equil_no_order(model):
    mole_frac_power_law_equil.build_parameters(
        model.rparams.reaction_r1,
        model.rparams.config.equilibrium_reactions["r1"])

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.reaction_order, Var)
    assert len(model.rparams.reaction_r1.reaction_order) == 6

    assert model.rparams.reaction_r1.reaction_order["p1", "c1"].value == -1
    assert model.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert model.rparams.reaction_r1.reaction_order["p2", "c1"].value == 0
    assert model.rparams.reaction_r1.reaction_order["p2", "c2"].value == 0
    # Solids should have zero order, as they are excluded
    assert model.rparams.reaction_r1.reaction_order["sol", "c1"].value == 0
    assert model.rparams.reaction_r1.reaction_order["sol", "c2"].value == 0

    # Check reaction form
    rform = mole_frac_power_law_equil.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        model.rxn[1].k_eq["r1"] == (
            model.thermo[1].mole_frac_phase_comp["p1", "c1"] **
            model.rparams.reaction_r1.reaction_order["p1", "c1"] *
            model.thermo[1].mole_frac_phase_comp["p1", "c2"] **
            model.rparams.reaction_r1.reaction_order["p1", "c2"]))


@pytest.mark.unit
def test_mole_frac_power_law_equil_with_order(model):
    model.rparams.config.equilibrium_reactions.r1.parameter_data = {
        "reaction_order": {("p1", "c1"): 1,
                           ("p1", "c2"): 2,
                           ("p2", "c1"): 3,
                           ("p2", "c2"): 4,
                           ("sol", "c1"): 5,
                           ("sol", "c2"): 6}}

    mole_frac_power_law_equil.build_parameters(
        model.rparams.reaction_r1,
        model.rparams.config.equilibrium_reactions["r1"])

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.reaction_order, Var)
    assert len(model.rparams.reaction_r1.reaction_order) == 6
    assert model.rparams.reaction_r1.reaction_order["p1", "c1"].value == 1
    assert model.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert model.rparams.reaction_r1.reaction_order["p2", "c1"].value == 3
    assert model.rparams.reaction_r1.reaction_order["p2", "c2"].value == 4
    assert model.rparams.reaction_r1.reaction_order["sol", "c1"].value == 5
    assert model.rparams.reaction_r1.reaction_order["sol", "c2"].value == 6

    # Check reaction form
    rform = mole_frac_power_law_equil.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        model.rxn[1].k_eq["r1"] == (
            model.thermo[1].mole_frac_phase_comp["p1", "c1"] **
            model.rparams.reaction_r1.reaction_order["p1", "c1"] *
            model.thermo[1].mole_frac_phase_comp["p1", "c2"] **
            model.rparams.reaction_r1.reaction_order["p1", "c2"] *
            model.thermo[1].mole_frac_phase_comp["p2", "c1"] **
            model.rparams.reaction_r1.reaction_order["p2", "c1"] *
            model.thermo[1].mole_frac_phase_comp["p2", "c2"] **
            model.rparams.reaction_r1.reaction_order["p2", "c2"] *
            model.thermo[1].mole_frac_phase_comp["sol", "c1"] **
            model.rparams.reaction_r1.reaction_order["sol", "c1"] *
            model.thermo[1].mole_frac_phase_comp["sol", "c2"] **
            model.rparams.reaction_r1.reaction_order["sol", "c2"]))
