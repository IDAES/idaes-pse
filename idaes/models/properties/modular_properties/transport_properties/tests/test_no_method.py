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
Test no_method method to specify that no method has been assigned to calculate
a physical property for a phase.

Authors: Douglas Allan
"""

import pytest
from sys import modules

from pyomo.environ import ConcreteModel, units as pyunits

from idaes.models.properties.modular_properties.transport_properties import NoMethod

from idaes.core import declare_process_block_class, LiquidPhase, VaporPhase

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.models.properties.modular_properties.eos.ideal import Ideal
from idaes.models.properties.modular_properties.pure import ConstantProperties


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True

    def parameters(self):
        self.parameters_set = True


def define_state(b):
    b.state_defined = True


# Dummy method to avoid errors when setting metadata dict
def set_metadata(b):
    pass


def construct_dummy_model(component_dict):
    m = ConcreteModel()
    components = [comp for comp in component_dict.keys()]
    comp1 = components[0]
    comp2 = components[1]

    m.params = DummyParameterBlock(
        components={
            comp1: {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (component_dict[comp1]["mw"], pyunits.g / pyunits.mol),
                },
            },
            comp2: {
                "visc_d_phase_comp": {"Vap": ConstantProperties},
                "parameter_data": {
                    "mw": (component_dict[comp2]["mw"], pyunits.g / pyunits.mol),
                },
            },
        },
        phases={
            "Vap": {
                "type": VaporPhase,
                "equation_of_state": Ideal,
                "visc_d_phase": NoMethod,
            },
        },
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
        state_definition=modules[__name__],
        pressure_ref=100000.0,
        temperature_ref=300,
    )

    m.props = m.params.state_block_class([1], defined_state=False, parameters=m.params)
    return m


@pytest.mark.unit
def test_NoMethod_visc_d_phase(caplog):
    component_dict = {
        "N2": {"mw": 28.014},
        "CO2": {"mw": 44.009},
    }
    m = construct_dummy_model(component_dict)
    caplog.clear()
    NoMethod.visc_d_phase.build_parameters(m.params.Vap)
    assert len(caplog.records) == 1
    assert (
        caplog.records[0].message
        == "Skipping construction of dynamic viscosity for phase Vap"
    )
    NoMethod.visc_d_phase.return_expression(m.props[1], "Vap")
    assert len(caplog.records) == 1


@pytest.mark.unit
def test_NoMethod_therm_cond_phase(caplog):
    component_dict = {
        "N2": {"mw": 28.014},
        "CO2": {"mw": 44.009},
    }
    m = construct_dummy_model(component_dict)
    caplog.clear()

    NoMethod.therm_cond_phase.build_parameters(m.params.Vap)
    assert len(caplog.records) == 1
    assert (
        caplog.records[0].message
        == "Skipping construction of thermal conductivity for phase Vap"
    )
    NoMethod.therm_cond_phase.return_expression(m.props[1], "Vap")
    assert len(caplog.records) == 1
