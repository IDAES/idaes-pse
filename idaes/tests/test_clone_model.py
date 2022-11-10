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
Created on Thu Mar 18 15:05:41 2021

@author: andrew
"""
import pytest

import pyomo.environ as pe
import idaes.core
import idaes.models.unit_models
import idaes.models.properties.swco2
from pyomo.network import Arc


@pytest.fixture(scope="function")
def model():
    m = pe.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock()
    m.fs.properties = idaes.models.properties.swco2.SWCO2ParameterBlock()
    m.fs.heater = idaes.models.unit_models.Heater(
        dynamic=False, property_package=m.fs.properties, has_pressure_change=True
    )
    m.fs.heater2 = idaes.models.unit_models.Heater(
        dynamic=False, property_package=m.fs.properties, has_pressure_change=True
    )
    m.fs.stream = Arc(source=m.fs.heater.outlet, destination=m.fs.heater2.inlet)

    return m


@pytest.mark.unit
def test_expand_arcs_and_clone(model):
    # Check that port references were attached to model
    assert hasattr(model.fs.heater, "_enth_mol_inlet_ref")
    assert hasattr(model.fs.heater, "_flow_mol_inlet_ref")
    assert hasattr(model.fs.heater, "_pressure_inlet_ref")
    assert hasattr(model.fs.heater, "_enth_mol_outlet_ref")
    assert hasattr(model.fs.heater, "_flow_mol_outlet_ref")
    assert hasattr(model.fs.heater, "_pressure_outlet_ref")

    assert hasattr(model.fs.heater2, "_enth_mol_inlet_ref")
    assert hasattr(model.fs.heater2, "_flow_mol_inlet_ref")
    assert hasattr(model.fs.heater2, "_pressure_inlet_ref")
    assert hasattr(model.fs.heater2, "_enth_mol_outlet_ref")
    assert hasattr(model.fs.heater2, "_flow_mol_outlet_ref")
    assert hasattr(model.fs.heater2, "_pressure_outlet_ref")

    assert model.fs.stream._expanded_block is None

    m2 = pe.TransformationFactory("network.expand_arcs").create_using(model)

    # Check that Arcs were expanded
    assert isinstance(m2.fs.stream._expanded_block, pe.Block)
    assert model.fs.stream._expanded_block is None


@pytest.mark.unit
def test_clone(model):
    m2 = model.clone()
    assert (
        m2.fs.heater.inlet.flow_mol[0]
        is m2.fs.heater.control_volume.properties_in[0].flow_mol
    )

    assert not (
        m2.fs.heater.inlet.flow_mol[0]
        is model.fs.heater.control_volume.properties_in[0].flow_mol
    )

    pe.TransformationFactory("network.expand_arcs").apply_to(m2)

    # Check that Arcs were expanded
    assert isinstance(m2.fs.stream._expanded_block, pe.Block)
    assert model.fs.stream._expanded_block is None
