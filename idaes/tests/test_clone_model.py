#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:05:41 2021

@author: andrew
"""
import pytest

import pyomo.environ as pe
import idaes.core
import idaes.generic_models.unit_models
import idaes.generic_models.properties.swco2
from pyomo.network import Arc


@pytest.fixture(scope="function")
def model():
    m = pe.ConcreteModel()
    m.fs = idaes.core.FlowsheetBlock()
    m.fs.properties = \
        idaes.generic_models.properties.swco2.SWCO2ParameterBlock()
    m.fs.heater = idaes.generic_models.unit_models.Heater(default={
        'dynamic': False,
        'property_package': m.fs.properties,
        'has_pressure_change': True})
    m.fs.heater2 = idaes.generic_models.unit_models.Heater(default={
        'dynamic': False,
        'property_package': m.fs.properties,
        'has_pressure_change': True})
    m.fs.stream = Arc(source=m.fs.heater.outlet,
                      destination=m.fs.heater2.inlet)

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

    m2 = pe.TransformationFactory('network.expand_arcs').create_using(model)

    # Check that Arcs were expanded
    assert isinstance(m2.fs.stream._expanded_block, pe.Block)
    assert model.fs.stream._expanded_block is None


@pytest.mark.unit
def test_clone(model):
    m2 = model.clone()
    assert (m2.fs.heater.inlet.flow_mol[0] is
            m2.fs.heater.control_volume.properties_in[0].flow_mol)

    assert not (m2.fs.heater.inlet.flow_mol[0] is
                model.fs.heater.control_volume.properties_in[0].flow_mol)

    pe.TransformationFactory('network.expand_arcs').apply_to(m2)

    # Check that Arcs were expanded
    assert isinstance(m2.fs.stream._expanded_block, pe.Block)
    assert model.fs.stream._expanded_block is None
