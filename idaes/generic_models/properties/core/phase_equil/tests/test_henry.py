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
Library of common forms for phase equilibrium constraints
"""
from pyomo.environ import ConcreteModel, Expression, Var, units as pyunits

from idaes.generic_models.properties.core.generic.generic_property import \
    GenericParameterBlock
from idaes.generic_models.properties.core.state_definitions import FTPx

from idaes.generic_models.properties.core.phase_equil.henry import \
    ConstantH
from idaes.generic_models.properties.core.phase_equil.forms import \
    fugacity
import pytest
from idaes.core.util.exceptions import ConfigurationError


# Dummy EoS
class DummyEoS(object):
    def common(self, *args, **kwargs):
        pass

    def build_parameters(b):
        pass

    def fug_phase_comp_eq(b, p, j, pp):
        return b.x[p, j]

    def log_fug_phase_comp_eq(b, p, j, pp):
        return 42*b.x[p, j]


@pytest.mark.unit
def test_henry_invalid_phase_type():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    with pytest.raises(ConfigurationError,
                       match="params component H2O was marked as a Henry's "
                       "Law component in phase Vap, but this is not a Liquid "
                       "phase."):
        m.params = GenericParameterBlock(default={
            "components": {"H2O": {
                "parameter_data": {"temperature_crit": 647.3},
                "henry_component": {"Vap": ConstantH},
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity}}},
            "phases": {"Liq": {"equation_of_state": DummyEoS},
                       "Vap": {"equation_of_state": DummyEoS}},
            "state_definition": FTPx,
            "pressure_ref": 1e5,
            "temperature_ref": 300,
            "base_units": {"time": pyunits.s,
                           "length": pyunits.m,
                           "mass": pyunits.kg,
                           "amount": pyunits.mol,
                           "temperature": pyunits.K}})


@pytest.mark.unit
def test_henry_invalid_phase_name():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    with pytest.raises(ConfigurationError,
                       match="params component H2O was marked as a Henry's "
                       "Law component in phase foo, but this is not a valid "
                       "phase name."):
        m.params = GenericParameterBlock(default={
            "components": {"H2O": {
                "parameter_data": {"temperature_crit": 647.3},
                "henry_component": {"foo": ConstantH},
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity}}},
            "phases": {"Liq": {"equation_of_state": DummyEoS},
                       "Vap": {"equation_of_state": DummyEoS}},
            "state_definition": FTPx,
            "pressure_ref": 1e5,
            "temperature_ref": 300,
            "base_units": {"time": pyunits.s,
                           "length": pyunits.m,
                           "mass": pyunits.kg,
                           "amount": pyunits.mol,
                           "temperature": pyunits.K}})


@pytest.mark.unit
def test_constant_H():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    m.params = GenericParameterBlock(default={
        "components": {"H2O": {"parameter_data": {"temperature_crit": 647.3,
                                                  "henry_ref": {
                                                      "Liq": 86}},
                               "henry_component": {
                                   "Liq": ConstantH},
                               "phase_equilibrium_form": {
                                   ("Vap", "Liq"): fugacity}}},
        "phases": {"Liq": {"equation_of_state": DummyEoS},
                   "Vap": {"equation_of_state": DummyEoS}},
        "state_definition": FTPx,
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {"time": pyunits.s,
                       "length": pyunits.m,
                       "mass": pyunits.kg,
                       "amount": pyunits.mol,
                       "temperature": pyunits.K}})

    assert isinstance(m.params.H2O.henry_ref_Liq, Var)
    assert m.params.H2O.henry_ref_Liq.value == 86

    m.state = m.params.build_state_block([0])

    H = ConstantH.return_expression(m.state[0], "Liq", "H2O")
    assert H is m.params.H2O.henry_ref_Liq

    assert isinstance(m.state[0].henry, Expression)
    assert len(m.state[0].henry) == 1
    assert m.state[0].henry["Liq", "H2O"]._expr is m.params.H2O.henry_ref_Liq
