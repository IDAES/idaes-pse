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
Tests for rate forms
"""

import pytest

from pyomo.environ import Block, ConcreteModel, Var, units as pyunits

from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.reactions.rate_forms import *

from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.core.util.misc import add_object_reference


@pytest.mark.unit
def test_power_law_rate_no_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()
    m.thermo = m.pparams.build_state_block([1])

    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        rate_reactions={
            "r1": {
                "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_rxn = Var(["r1"], initialize=1)

    power_law_rate.build_parameters(
        m.rparams.reaction_r1, m.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 4
    for i, v in m.rparams.reaction_r1.reaction_order.items():
        try:
            stoic = m.rparams.config.rate_reactions.r1.stoichiometry[i]
        except KeyError:
            stoic = 0

        if stoic < 1:
            assert v.value == -stoic
        else:
            assert v.value == 0

    # Check reaction form
    rform = power_law_rate.return_expression(m.rxn[1], m.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        m.rxn[1].k_rxn["r1"]
        * m.thermo[1].mole_frac_phase_comp["p1", "c1"]
        ** m.rparams.reaction_r1.reaction_order["p1", "c1"]
    )


@pytest.mark.unit
def test_power_law_rate_with_order():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()
    m.thermo = m.pparams.build_state_block([1])

    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        base_units={
            "time": pyunits.s,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
        rate_reactions={
            "r1": {
                "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                "rate_form": power_law_rate,
                "concentration_form": ConcentrationForm.moleFraction,
                "parameter_data": {
                    "reaction_order": {
                        ("p1", "c1"): 1,
                        ("p1", "c2"): 2,
                        ("p2", "c1"): 3,
                        ("p2", "c2"): 4,
                    }
                },
            }
        },
    )

    # Create a dummy state block
    m.rxn = Block([1])
    add_object_reference(
        m.rxn[1], "phase_component_set", m.pparams._phase_component_set
    )
    add_object_reference(m.rxn[1], "params", m.rparams)
    add_object_reference(m.rxn[1], "state_ref", m.thermo[1])

    m.rxn[1].k_rxn = Var(["r1"], initialize=1)

    power_law_rate.build_parameters(
        m.rparams.reaction_r1, m.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(m.rparams.reaction_r1.reaction_order, Var)
    assert len(m.rparams.reaction_r1.reaction_order) == 4
    assert m.rparams.reaction_r1.reaction_order["p1", "c1"].value == 1
    assert m.rparams.reaction_r1.reaction_order["p1", "c2"].value == 2
    assert m.rparams.reaction_r1.reaction_order["p2", "c1"].value == 3
    assert m.rparams.reaction_r1.reaction_order["p2", "c2"].value == 4

    # Check reaction form
    rform = power_law_rate.return_expression(m.rxn[1], m.rparams.reaction_r1, "r1", 300)

    assert str(rform) == str(
        m.rxn[1].k_rxn["r1"]
        * (
            m.thermo[1].mole_frac_phase_comp["p1", "c1"]
            ** m.rparams.reaction_r1.reaction_order["p1", "c1"]
            * m.thermo[1].mole_frac_phase_comp["p1", "c2"]
            ** m.rparams.reaction_r1.reaction_order["p1", "c2"]
            * m.thermo[1].mole_frac_phase_comp["p2", "c1"]
            ** m.rparams.reaction_r1.reaction_order["p2", "c1"]
            * m.thermo[1].mole_frac_phase_comp["p2", "c2"]
            ** m.rparams.reaction_r1.reaction_order["p2", "c2"]
        )
    )
