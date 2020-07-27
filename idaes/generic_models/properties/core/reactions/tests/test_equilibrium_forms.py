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

from pyomo.environ import Block, ConcreteModel, Var, units as pyunits

from idaes.generic_models.properties.core.generic.generic_reaction import \
    GenericReactionParameterBlock, ConcentrationForm
from idaes.generic_models.properties.core.reactions.equilibrium_forms import *

from idaes.core import SolidPhase
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.misc import add_object_reference


@pytest.mark.unit
def test_power_law_equil_no_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(default={
        "property_package": m.pparams,
        "base_units": {"time": pyunits.s,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "length": pyunits.m,
                       "temperature": pyunits.K},
        "equilibrium_reactions": {
            "r1": {"stoichiometry": {("p1", "c1"): -1,
                                     ("p1", "c2"): 2,
                                     ("sol", "c1"): -3,
                                     ("sol", "c2"): 4},
                   "equilibrium_form": power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction}}})

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    power_law_equil.build_parameters(
        m.rparams.reaction_r1,
        m.rparams.config.equilibrium_reactions["r1"])

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6

    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == -1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 0
    # Solids should have zero order, as they are excluded
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 0
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 0

    # Check reaction form
    rform = power_law_equil.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        m.rxn[1].k_eq["r1"] == (
            m.thermo[1].mole_frac_phase_comp["p1", "c1"] **
            m.rparams.reaction_r1.reaction_order["p1", "c1"] *
            m.thermo[1].mole_frac_phase_comp["p1", "c2"] **
            m.rparams.reaction_r1.reaction_order["p1", "c2"]))


@pytest.mark.unit
def test_power_law_equil_with_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()

    # Add a solid phase for testing
    m.pparams.sol = SolidPhase()

    m.thermo = m.pparams.build_state_block([1])

    # Create a dummy reaction parameter block
    m.rparams = GenericReactionParameterBlock(default={
        "property_package": m.pparams,
        "base_units": {"time": pyunits.s,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "length": pyunits.m,
                       "temperature": pyunits.K},
        "equilibrium_reactions": {
            "r1": {"stoichiometry": {("p1", "c1"): -1,
                                     ("p1", "c2"): 2,
                                     ("sol", "c1"): -3,
                                     ("sol", "c2"): 4},
                   "equilibrium_form": power_law_equil,
                   "concentration_form": ConcentrationForm.moleFraction,
                   "parameter_data": {
                       "reaction_order": {("p1", "c1"): 1,
                                          ("p1", "c2"): 2,
                                          ("p2", "c1"): 3,
                                          ("p2", "c2"): 4,
                                          ("sol", "c1"): 5,
                                          ("sol", "c2"): 6}}}}})

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_eq = Var(["r1"], initialize=1)

    power_law_equil.build_parameters(
        m.rparams.reaction_r1,
        m.rparams.config.equilibrium_reactions["r1"])

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 6
    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == 1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 3
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 4
    assert m.rparams.reaction_r1.reaction_order["sol", "c1"].value == 5
    assert m.rparams.reaction_r1.reaction_order["sol", "c2"].value == 6

    # Check reaction form
    rform = power_law_equil.return_expression(
        m.rxn[1], m.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        m.rxn[1].k_eq["r1"] == (
            m.thermo[1].mole_frac_phase_comp["p1", "c1"] **
            m.rparams.reaction_r1.reaction_order["p1", "c1"] *
            m.thermo[1].mole_frac_phase_comp["p1", "c2"] **
            m.rparams.reaction_r1.reaction_order["p1", "c2"] *
            m.thermo[1].mole_frac_phase_comp["p2", "c1"] **
            m.rparams.reaction_r1.reaction_order["p2", "c1"] *
            m.thermo[1].mole_frac_phase_comp["p2", "c2"] **
            m.rparams.reaction_r1.reaction_order["p2", "c2"] *
            m.thermo[1].mole_frac_phase_comp["sol", "c1"] **
            m.rparams.reaction_r1.reaction_order["sol", "c1"] *
            m.thermo[1].mole_frac_phase_comp["sol", "c2"] **
            m.rparams.reaction_r1.reaction_order["sol", "c2"]))
