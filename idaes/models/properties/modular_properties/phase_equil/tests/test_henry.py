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
Library of common forms for phase equilibrium constraints
"""
from pyomo.environ import ConcreteModel, Expression, value, Var, units as pyunits

from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterBlock,
)
from idaes.models.properties.modular_properties.state_definitions import FTPx

from idaes.models.properties.modular_properties.phase_equil.henry import (
    ConstantH,
    HenryType,
    henry_equilibrium_ratio,
    _raise_henry_type_error,
    get_henry_concentration_term,
)
from idaes.models.properties.modular_properties.phase_equil.forms import fugacity
from pytest import approx
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
        return 42 * b.x[p, j]


@pytest.mark.unit
def test_henry_invalid_phase_type():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    with pytest.raises(
        ConfigurationError,
        match="params component H2O was marked as a Henry's "
        "Law component in phase Vap, but this is not a Liquid "
        "phase.",
    ):
        m.params = GenericParameterBlock(
            components={
                "H2O": {
                    "parameter_data": {"temperature_crit": 647.3},
                    "henry_component": {"Vap": ConstantH},
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                }
            },
            phases={
                "Liq": {"equation_of_state": DummyEoS},
                "Vap": {"equation_of_state": DummyEoS},
            },
            state_definition=FTPx,
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )


@pytest.mark.unit
def test_henry_invalid_phase_name():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    with pytest.raises(
        ConfigurationError,
        match="params component H2O was marked as a Henry's "
        "Law component in phase foo, but this is not a valid "
        "phase name.",
    ):
        m.params = GenericParameterBlock(
            components={
                "H2O": {
                    "parameter_data": {"temperature_crit": 647.3},
                    "henry_component": {"foo": ConstantH},
                    "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
                }
            },
            phases={
                "Liq": {"equation_of_state": DummyEoS},
                "Vap": {"equation_of_state": DummyEoS},
            },
            state_definition=FTPx,
            pressure_ref=100000.0,
            temperature_ref=300,
            base_units={
                "time": pyunits.s,
                "length": pyunits.m,
                "mass": pyunits.kg,
                "amount": pyunits.mol,
                "temperature": pyunits.K,
            },
        )


@pytest.mark.unit
def test_constant_H():
    m = ConcreteModel()

    # Add a dummy var for use in constructing expressions
    m.x = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    m.mole_frac_phase_comp = Var(["Vap", "Liq"], ["H2O"], initialize=1)

    # Create a dummy parameter block
    m.params = GenericParameterBlock(
        components={
            "H2O": {
                "parameter_data": {"temperature_crit": 647.3, "henry_ref": {"Liq": 86}},
                "henry_component": {
                    "Liq": {"method": ConstantH, "type": HenryType.Kpx}
                },
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            }
        },
        phases={
            "Liq": {"equation_of_state": DummyEoS},
            "Vap": {"equation_of_state": DummyEoS},
        },
        state_definition=FTPx,
        pressure_ref=100000.0,
        temperature_ref=300,
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
    )

    assert isinstance(m.params.H2O.henry_ref_Liq, Var)
    assert m.params.H2O.henry_ref_Liq.value == 86

    m.state = m.params.build_state_block([0])

    H = ConstantH.return_expression(m.state[0], "Liq", "H2O")
    assert H is m.params.H2O.henry_ref_Liq

    assert isinstance(m.state[0].henry, Expression)
    assert len(m.state[0].henry) == 1
    assert m.state[0].henry["Liq", "H2O"]._expr is m.params.H2O.henry_ref_Liq


@pytest.mark.unit
def test_invalid_henry_type():
    with pytest.raises(
        ConfigurationError, match="Unrecognized value for HenryType HenryType.Dummy"
    ):
        _raise_henry_type_error(HenryType.Dummy)

    m = ConcreteModel()
    m.henry = {("Liq", "H2O"): 1}
    henry_dict = {"type": HenryType.Dummy}

    def debugMethod():
        return 42

    m.params = GenericParameterBlock(
        components={
            "H2O": {
                "parameter_data": {"temperature_crit": 647.3, "henry_ref": {"Liq": 86}},
                "henry_component": {
                    "Liq": {"method": debugMethod, "type": HenryType.Dummy}
                },
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            }
        },
        phases={
            "Liq": {"equation_of_state": DummyEoS},
            "Vap": {"equation_of_state": DummyEoS},
        },
        state_definition=FTPx,
        pressure_ref=100000.0,
        temperature_ref=300,
        base_units={
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
    )

    with pytest.raises(
        ConfigurationError, match="Unrecognized value for HenryType HenryType.Dummy"
    ):
        get_henry_concentration_term(m, henry_dict)

    with pytest.raises(
        ConfigurationError, match="Unrecognized value for HenryType HenryType.Dummy"
    ):
        henry_equilibrium_ratio(m, "Liq", "H2O")


@pytest.mark.component
def test_equilibrium_ratio():
    config_dict = {
        "components": {
            "H2O": {
                "parameter_data": {"temperature_crit": 647.3, "henry_ref": {"Liq": 7}},
                "henry_component": {"Liq": {"method": ConstantH, "type": None}},
                "phase_equilibrium_form": {("Vap", "Liq"): fugacity},
            }
        },
        "phases": {
            "Liq": {"equation_of_state": DummyEoS},
            "Vap": {"equation_of_state": DummyEoS},
        },
        "state_definition": FTPx,
        "pressure_ref": 1e5,
        "temperature_ref": 300,
        "base_units": {
            "time": pyunits.s,
            "length": pyunits.m,
            "mass": pyunits.kg,
            "amount": pyunits.mol,
            "temperature": pyunits.K,
        },
    }

    henry_types = [HenryType.Hcp, HenryType.Kpc, HenryType.Hxp, HenryType.Kpx]
    expected_split = [3 / (5 * 7), (3 * 7) / 5, 1 / (5 * 7), 7 / 5]

    for htype, split in zip(henry_types, expected_split):
        config_dict["components"]["H2O"]["henry_component"]["Liq"]["type"] = htype
        m = ConcreteModel()
        m.params = GenericParameterBlock(**config_dict)

        m.state = m.params.build_state_block([0])
        m.state[0].mole_frac_phase_comp["Liq", "H2O"].value = 0.5
        m.state[0].dens_mol_phase = Var(
            ["Vap", "Liq"], initialize=0, units=pyunits.mol / (pyunits.m**3)
        )
        m.state[0].dens_mol_phase["Liq"].value = 3
        m.state[0].pressure.value = 5

        assert approx(split, 1e-6) == value(
            henry_equilibrium_ratio(m.state[0], "Liq", "H2O")
        )
        assert (
            pyunits.get_units(henry_equilibrium_ratio(m.state[0], "Liq", "H2O")) is None
        )
