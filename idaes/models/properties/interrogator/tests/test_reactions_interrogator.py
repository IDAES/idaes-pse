#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
Tests for Reaction Interrogator Tool

@author: alee
"""
import pytest

from pyomo.environ import ConcreteModel, units as pyunits

from idaes.core import FlowsheetBlock, LiquidPhase, Solute
from idaes.models.unit_models import CSTR, PFR
from idaes.models.properties.interrogator import (
    PropertyInterrogatorBlock,
    ReactionInterrogatorBlock,
)


@pytest.mark.unit
def test_interrogator_parameter_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    # Check that parameter block has expected attributes
    assert isinstance(m.fs.rxn_params.required_properties, dict)
    assert len(m.fs.rxn_params.required_properties) == 0


@pytest.mark.unit
def test_interrogator_rxn_block_unindexed_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    # Check get_term methods return an unindexed dummy var
    assert m.fs.rxns[0].prop_unindexed is m.fs.rxns[0]._dummy_var

    # Call again to make sure duplicates are skipped in required_properties
    assert m.fs.rxns[0].prop_unindexed is m.fs.rxns[0]._dummy_var

    # Check that get_term calls were logged correctly
    assert m.fs.rxn_params.required_properties == {"prop_unindexed": ["fs.rxns"]}


@pytest.mark.unit
def test_interrogator_rxn_block_phase_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    # Check get_term methods return an unindexed dummy var
    assert m.fs.rxns[0].prop_phase["Liq"] is m.fs.rxns[0]._dummy_var_phase["Liq"]
    assert m.fs.rxns[0].prop_phase["Vap"] is m.fs.rxns[0]._dummy_var_phase["Vap"]

    # Check that get_term calls were logged correctly
    assert m.fs.rxn_params.required_properties == {"prop_phase": ["fs.rxns"]}


@pytest.mark.unit
def test_interrogator_rxn_block_comp_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    # Check get_term methods return an unindexed dummy var
    assert m.fs.rxns[0].prop_comp["A"] is m.fs.rxns[0]._dummy_var_comp["A"]
    assert m.fs.rxns[0].prop_comp["B"] is m.fs.rxns[0]._dummy_var_comp["B"]

    # Check that get_term calls were logged correctly
    assert m.fs.rxn_params.required_properties == {"prop_comp": ["fs.rxns"]}


@pytest.mark.unit
def test_interrogator_rxn_block_phase_comp_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    # Check get_term methods return an unindexed dummy var
    assert (
        m.fs.rxns[0].prop_phase_comp["Liq", "A"]
        is m.fs.rxns[0]._dummy_var_phase_comp["Liq", "A"]
    )
    assert (
        m.fs.rxns[0].prop_phase_comp["Vap", "B"]
        is m.fs.rxns[0]._dummy_var_phase_comp["Vap", "B"]
    )

    # Check that get_term calls were logged correctly
    assert m.fs.rxn_params.required_properties == {"prop_phase_comp": ["fs.rxns"]}


@pytest.mark.unit
def test_interrogator_rxn_block_reaction_rate_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    # Check get_term methods return an unindexed dummy var
    assert m.fs.rxns[0].reaction_rate["R1"] is m.fs.rxns[0]._dummy_reaction_idx["R1"]
    assert m.fs.rxns[0].reaction_rate["R1"] is m.fs.rxns[0]._dummy_reaction_idx["R1"]

    # Check that get_term calls were logged correctly
    assert m.fs.rxn_params.required_properties == {"reaction_rate": ["fs.rxns"]}


@pytest.mark.unit
def test_interrogator_rxn_block_dh_rxn_call():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    # Check get_term methods return an unindexed dummy var
    assert m.fs.rxns[0].dh_rxn["R1"] is m.fs.rxns[0]._dummy_reaction_idx["R1"]
    assert m.fs.rxns[0].dh_rxn["R1"] is m.fs.rxns[0]._dummy_reaction_idx["R1"]

    # Check that get_term calls were logged correctly
    assert m.fs.rxn_params.required_properties == {"dh_rxn": ["fs.rxns"]}


@pytest.mark.unit
def test_interrogator_initialize_method():
    # Initialize method should return an TypeError
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    with pytest.raises(
        TypeError,
        match="Models constructed using the Reaction "
        "Interrogator package cannot be used to solve a "
        "flowsheet. Please rebuild your flowsheet using a "
        "valid reaction package.",
    ):
        m.fs.rxns.initialize()


@pytest.fixture(scope="module")
def model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)

    m.fs.params = PropertyInterrogatorBlock()
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.R01 = CSTR(
        property_package=m.fs.params,
        reaction_package=m.fs.rxn_params,
        has_heat_of_reaction=True,
    )

    m.fs.R02 = PFR(property_package=m.fs.params, reaction_package=m.fs.rxn_params)

    return m


# Test for physical parameters too, as these can not be checked without having
# a reaction package too.
@pytest.mark.component
def test_interrogate_flowsheet(model):
    assert model.fs.params.required_properties == {
        "material flow terms": ["fs.R01", "fs.R02"],
        "enthalpy flow terms": ["fs.R01", "fs.R02"],
        "material density terms": ["fs.R01", "fs.R02"],
        "energy density terms": ["fs.R01", "fs.R02"],
        "pressure": ["fs.R01", "fs.R02"],
    }

    assert model.fs.rxn_params.required_properties == {
        "reaction_rate": ["fs.R01", "fs.R02"],
        "dh_rxn": ["fs.R01"],
    }


@pytest.mark.component
def test_list_required_properties(model):
    prop_list = model.fs.params.list_required_properties()
    assert prop_list == [
        "material density terms",
        "material flow terms",
        "enthalpy flow terms",
        "energy density terms",
        "pressure",
    ]

    rxn_list = model.fs.rxn_params.list_required_properties()
    assert rxn_list == ["dh_rxn", "reaction_rate"]


@pytest.mark.unit
def test_list_models_requiring_property(model):
    for k in model.fs.params.required_properties.keys():
        model_list = model.fs.params.list_models_requiring_property(k)
        assert model_list == ["fs.R01", "fs.R02"]


@pytest.mark.unit
def test_list_properties_required_by_model_by_name(model):
    prop_list = model.fs.rxn_params.list_properties_required_by_model("fs.R01")

    assert prop_list == ["dh_rxn", "reaction_rate"]


@pytest.mark.unit
def test_list_properties_required_by_model_by_object(model):
    prop_list = model.fs.rxn_params.list_properties_required_by_model(model.fs.R01)

    assert prop_list == ["dh_rxn", "reaction_rate"]


@pytest.mark.unit
def test_list_properties_required_by_model_invalid_model(model):
    with pytest.raises(ValueError):
        model.fs.rxn_params.list_properties_required_by_model("foo")


@pytest.mark.unit
def test_print_required_properties(model, capsys):
    model.fs.params.print_required_properties()
    model.fs.rxn_params.print_required_properties()

    captured = capsys.readouterr()
    assert (
        captured.out
        == """
==========================================================================
Property Interrogator Summary

The Flowsheet requires the following properties (times required):

    material density terms                                               2
    material flow terms                                                  2
    enthalpy flow terms                                                  2
    energy density terms                                                 2
    pressure                                                             2

Note: User constraints may require additional properties which are not
reported here.

==========================================================================
Reaction Property Interrogator Summary

The Flowsheet requires the following reaction properties (times required):

    dh_rxn                                                               1
    reaction_rate                                                        2

Note: User constraints may require additional properties which are not
reported here.
"""
    )


@pytest.mark.unit
def test_print_models_requiring_property(model, capsys):
    model.fs.rxn_params.print_models_requiring_property("reaction_rate")

    captured = capsys.readouterr()
    assert (
        captured.out
        == """
The following models in the Flowsheet require reaction_rate:
    fs.R01
    fs.R02
"""
    )


@pytest.mark.unit
def test_print_properties_reqruied_by_model(model, capsys):
    model.fs.rxn_params.print_properties_required_by_model("fs.R01")

    captured = capsys.readouterr()
    assert (
        captured.out
        == """
The following reaction properties are required by model fs.R01:
    dh_rxn
    reaction_rate
"""
    )


@pytest.mark.unit
def test_interrogator_rxn_block_unindexed_call_custom_phase_comp():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.params = PropertyInterrogatorBlock(
        phase_list={"P1": LiquidPhase, "P2": None},
        component_list={"c1": Solute, "c2": None},
    )
    m.fs.rxn_params = ReactionInterrogatorBlock(property_package=m.fs.params)

    m.fs.props = m.fs.params.build_state_block([0])
    m.fs.rxns = m.fs.rxn_params.build_reaction_block([0], state_block=m.fs.props)

    # Check phase and component lists
    assert m.fs.rxn_params.phase_list == ["P1", "P2"]
    assert m.fs.rxn_params.component_list == ["c1", "c2"]

    # Check get_term methods return an unindexed dummy var
    assert m.fs.rxns[0].prop_unindexed is m.fs.rxns[0]._dummy_var

    # Call again to make sure duplicates are skipped in required_properties
    assert m.fs.rxns[0].prop_unindexed is m.fs.rxns[0]._dummy_var

    # Check that get_term calls were logged correctly
    assert m.fs.rxn_params.required_properties == {"prop_unindexed": ["fs.rxns"]}
