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

from pyomo.environ import ConcreteModel, Var, units as pyunits, value
from pyomo.util.check_units import assert_units_equivalent

from idaes.generic_models.properties.core.generic.generic_reaction import \
    GenericReactionParameterBlock, ConcentrationForm
from idaes.generic_models.properties.core.reactions.equilibrium_constant import *
from idaes.core import MaterialFlowBasis
from idaes.core.util.testing import PhysicalParameterTestBlock


@pytest.fixture
def model():
    m = ConcreteModel()

    # # Add a test thermo package for validation
    m.pparams = PhysicalParameterTestBlock()
    m.thermo = m.pparams.build_state_block([1])

    m.rparams = GenericReactionParameterBlock(default={
         "property_package": m.pparams,
         "reaction_basis": MaterialFlowBasis.molar,
         "equilibrium_reactions": {
             "r1": {"stoichiometry": {("p1", "c1"): -1,
                                      ("p1", "c2"): 2},
                    "equilibrium_form": "foo",
                    "concentration_form": ConcentrationForm.moleFraction}},
         "base_units": {"amount": pyunits.mol,
                        "mass": pyunits.kg,
                        "time": pyunits.s,
                        "length": pyunits.m,
                        "temperature": pyunits.K}})

    m.rxn = m.rparams.build_reaction_block([1], default={
        "state_block": m.thermo, "has_equilibrium": False})

    m.rxn[1].dh_rxn = Var(["r1"],
                          initialize=1,
                          units=pyunits.J/pyunits.mol)

    return m


@pytest.mark.unit
def test_van_t_hoff_mole_frac(model):
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


@pytest.mark.unit
def test_van_t_hoff_mole_frac_convert(model):
    model.rparams.config.equilibrium_reactions.r1.parameter_data = {
        "k_eq_ref": (1, None),
        "T_eq_ref": (900, pyunits.degR)}

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


@pytest.mark.unit
def test_van_t_hoff_molarity(model):
    model.rparams.config.equilibrium_reactions.r1.concentration_form = \
        ConcentrationForm.molarity
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
    assert_units_equivalent(rform, pyunits.mol/pyunits.m**3)


@pytest.mark.unit
def test_van_t_hoff_molarity_convert(model):
    model.rparams.config.equilibrium_reactions.r1.concentration_form = \
        ConcentrationForm.molarity
    model.rparams.config.equilibrium_reactions.r1.parameter_data = {
        "k_eq_ref": (1e-3, pyunits.kmol/pyunits.m**3),
        "T_eq_ref": (900, pyunits.degR)}

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
    assert_units_equivalent(rform, pyunits.mol/pyunits.m**3)


@pytest.mark.unit
def test_van_t_hoff_molality(model):
    model.rparams.config.equilibrium_reactions.r1.concentration_form = \
        ConcentrationForm.molality
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
    assert_units_equivalent(rform, pyunits.mol/pyunits.kg)


@pytest.mark.unit
def test_van_t_hoff_molality_convert(model):
    model.rparams.config.equilibrium_reactions.r1.concentration_form = \
        ConcentrationForm.molality
    model.rparams.config.equilibrium_reactions.r1.parameter_data = {
        "k_eq_ref": (1e-3, pyunits.kmol/pyunits.kg),
        "T_eq_ref": (900, pyunits.degR)}

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
    assert_units_equivalent(rform, pyunits.mol/pyunits.kg)


@pytest.mark.unit
def test_van_t_hoff_partial_pressure(model):
    model.rparams.config.equilibrium_reactions.r1.concentration_form = \
        ConcentrationForm.partialPressure
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
    assert_units_equivalent(rform, pyunits.Pa)


@pytest.mark.unit
def test_van_t_hoff_partial_pressure_convert(model):
    model.rparams.config.equilibrium_reactions.r1.concentration_form = \
        ConcentrationForm.partialPressure
    model.rparams.config.equilibrium_reactions.r1.parameter_data = {
        "k_eq_ref": (1e-3, pyunits.kPa),
        "T_eq_ref": (900, pyunits.degR)}

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
    assert_units_equivalent(rform, pyunits.Pa)
