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

from pyomo.environ import ConcreteModel, Var, units as pyunits, value
from pyomo.util.check_units import assert_units_equivalent

from idaes.models.properties.modular_properties.base.generic_reaction import (
    GenericReactionParameterBlock,
    ConcentrationForm,
)
from idaes.models.properties.modular_properties.reactions.rate_constant import *

from idaes.core import MaterialFlowBasis
from idaes.core.util.testing import PhysicalParameterTestBlock


@pytest.fixture
def model():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()
    m.thermo = m.pparams.build_state_block([1])

    m.rparams = GenericReactionParameterBlock(
        property_package=m.pparams,
        reaction_basis=MaterialFlowBasis.molar,
        rate_reactions={
            "r1": {
                "stoichiometry": {("p1", "c1"): -1, ("p1", "c2"): 2},
                "rate_form": "foo",
                "concentration_form": ConcentrationForm.moleFraction,
            }
        },
        base_units={
            "amount": pyunits.mol,
            "mass": pyunits.kg,
            "time": pyunits.s,
            "length": pyunits.m,
            "temperature": pyunits.K,
        },
    )

    m.rxn = m.rparams.build_reaction_block(
        [1], state_block=m.thermo, has_equilibrium=False
    )

    return m


@pytest.mark.unit
def test_arrhenius_mole_frac(model):
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": 1,
        "energy_activation": 500,
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(rform, pyunits.mol / pyunits.m**3 / pyunits.s)


@pytest.mark.unit
def test_arrhenius_mole_frac_convert(model):
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": (1e-3, pyunits.kmol / pyunits.m**3 / pyunits.s),
        "energy_activation": (0.5, pyunits.kJ / pyunits.mol),
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(rform, pyunits.mol / pyunits.m**3 / pyunits.s)


@pytest.mark.unit
def test_arrhenius_molarity(model):
    model.rparams.config.rate_reactions.r1.concentration_form = (
        ConcentrationForm.molarity
    )
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": 1,
        "energy_activation": 500,
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(rform, 1 / pyunits.s)


@pytest.mark.unit
def test_arrhenius_molarity_convert(model):
    model.rparams.config.rate_reactions.r1.concentration_form = (
        ConcentrationForm.molarity
    )
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": (3600, 1 / pyunits.hr),
        "energy_activation": (0.5, pyunits.kJ / pyunits.mol),
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(rform, 1 / pyunits.s)


@pytest.mark.unit
def test_arrhenius_molality(model):
    model.rparams.config.rate_reactions.r1.concentration_form = (
        ConcentrationForm.molality
    )
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": 1,
        "energy_activation": 500,
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(rform, pyunits.kg / pyunits.m**3 / pyunits.s)


@pytest.mark.unit
def test_arrhenius_molality_convert(model):
    model.rparams.config.rate_reactions.r1.concentration_form = (
        ConcentrationForm.molality
    )
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": (1e3, pyunits.g / pyunits.m**3 / pyunits.s),
        "energy_activation": (0.5, pyunits.kJ / pyunits.mol),
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(rform, pyunits.kg / pyunits.m**3 / pyunits.s)


@pytest.mark.unit
def test_arrhenius_partial_pressure(model):
    model.rparams.config.rate_reactions.r1.concentration_form = (
        ConcentrationForm.partialPressure
    )
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": 1,
        "energy_activation": 500,
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(
        rform, pyunits.mol / pyunits.m**3 / pyunits.s / pyunits.Pa
    )


@pytest.mark.unit
def test_arrhenius_partial_pressure_convert(model):
    model.rparams.config.rate_reactions.r1.concentration_form = (
        ConcentrationForm.partialPressure
    )
    model.rparams.config.rate_reactions.r1.parameter_data = {
        "arrhenius_const": (
            1e-3,
            pyunits.kmol / pyunits.m**3 / pyunits.s / pyunits.Pa,
        ),
        "energy_activation": (0.5, pyunits.kJ / pyunits.mol),
    }

    arrhenius.build_parameters(
        model.rparams.reaction_r1, model.rparams.config.rate_reactions["r1"]
    )

    # Check parameter construction
    assert isinstance(model.rparams.reaction_r1.arrhenius_const, Var)
    assert model.rparams.reaction_r1.arrhenius_const.value == 1

    assert isinstance(model.rparams.reaction_r1.energy_activation, Var)
    assert model.rparams.reaction_r1.energy_activation.value == 500

    # Check expressions
    rform = arrhenius.return_expression(
        model.rxn[1], model.rparams.reaction_r1, "r1", 300 * pyunits.K
    )

    assert value(rform) == pytest.approx(0.81836, rel=1e-3)
    assert_units_equivalent(
        rform, pyunits.mol / pyunits.m**3 / pyunits.s / pyunits.Pa
    )
