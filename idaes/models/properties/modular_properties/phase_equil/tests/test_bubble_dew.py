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
Tests for smooth VLE formulation

Authors: Andrew Lee
"""

import pytest
from sys import modules

from pyomo.environ import ConcreteModel, Constraint, Set, value, Var, units as pyunits

from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
)
from idaes.core import declare_process_block_class
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS


# Dummy class to use for Psat calls
Psat = {"H2O": 1e5, "EtOH": 5e4}


# Dummy method to avoid errors when setting metadata dict
def set_metadata(b):
    pass


def pressure_sat_comp(b, j, T):
    return Psat[j.local_name]


@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True

    def parameters(self):
        self.parameters_set = True


def define_state(b):
    b.state_defined = True


@pytest.fixture(scope="class")
def frame():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(
        components={
            "H2O": {"pressure_sat_comp": pressure_sat_comp},
            "EtOH": {"pressure_sat_comp": pressure_sat_comp},
        },
        phases={
            "Liq": {"equation_of_state": DummyEoS},
            "Vap": {"equation_of_state": DummyEoS},
        },
        state_definition=modules[__name__],
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
    m.params._pe_pairs = Set(initialize=[("Vap", "Liq")])

    m.props = m.params.build_state_block([1], defined_state=False)

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300)
    m.props[1]._teq = Var(initialize=300)
    m.props[1].mole_frac_comp = Var(m.params.component_list, initialize=0.5)

    return m


class TestBubbleTempIdeal(object):
    @pytest.mark.unit
    def test_build(self, frame):
        frame.props[1].temperature_bubble = Var(frame.params._pe_pairs)
        frame.props[1]._mole_frac_tbub = Var(
            frame.params._pe_pairs, frame.params.component_list, initialize=0.5
        )

        IdealBubbleDew.temperature_bubble(frame.props[1])

        assert isinstance(frame.props[1].eq_temperature_bubble, Constraint)
        assert len(frame.props[1].eq_temperature_bubble) == 1

        assert isinstance(frame.props[1].eq_mole_frac_tbub, Constraint)
        assert len(frame.props[1].eq_mole_frac_tbub) == 2
        for k in frame.props[1].eq_mole_frac_tbub:
            assert k in [("Vap", "Liq", "H2O"), ("Vap", "Liq", "EtOH")]

    @pytest.mark.unit
    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1 / 10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1 / 10

            frame.props[1].pressure = value(
                sum(
                    frame.props[1].mole_frac_comp[j] * Psat[j]
                    for j in frame.params.component_list
                )
            )

            for pp in frame.params._pe_pairs:
                for j in frame.params.component_list:
                    frame.props[1]._mole_frac_tbub[pp[0], pp[1], j].value = value(
                        frame.props[1].mole_frac_comp[j]
                        * Psat[j]
                        / frame.props[1].pressure
                    )

                assert value(
                    frame.props[1].eq_temperature_bubble[pp[0], pp[1]].body
                ) == pytest.approx(0, abs=1e-8)
                for k in frame.params.component_list:
                    assert value(
                        frame.props[1].eq_mole_frac_tbub[pp[0], pp[1], k].body
                    ) == pytest.approx(0, abs=1e-8)


class TestDewTempIdeal(object):
    @pytest.mark.unit
    def test_build(self, frame):
        frame.props[1].temperature_dew = Var(frame.params._pe_pairs)
        frame.props[1]._mole_frac_tdew = Var(
            frame.params._pe_pairs, frame.params.component_list, initialize=0.5
        )

        IdealBubbleDew.temperature_dew(frame.props[1])

        assert isinstance(frame.props[1].eq_temperature_dew, Constraint)
        assert len(frame.props[1].eq_temperature_dew) == 1

        assert isinstance(frame.props[1].eq_mole_frac_tdew, Constraint)
        assert len(frame.props[1].eq_mole_frac_tdew) == 2
        for k in frame.props[1].eq_mole_frac_tdew:
            assert k in [("Vap", "Liq", "H2O"), ("Vap", "Liq", "EtOH")]

    @pytest.mark.unit
    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1 / 10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1 / 10

            frame.props[1].pressure = value(
                1
                / sum(
                    frame.props[1].mole_frac_comp[j] / Psat[j]
                    for j in frame.params.component_list
                )
            )

            for pp in frame.params._pe_pairs:
                for j in frame.params.component_list:
                    frame.props[1]._mole_frac_tdew[pp[0], pp[1], j].value = value(
                        frame.props[1].mole_frac_comp[j]
                        * frame.props[1].pressure
                        / Psat[j]
                    )

                assert value(
                    frame.props[1].eq_temperature_dew[pp[0], pp[1]].body
                ) == pytest.approx(0, abs=1e-8)
                for k in frame.params.component_list:
                    assert value(
                        frame.props[1].eq_mole_frac_tdew[pp[0], pp[1], k].body
                    ) == pytest.approx(0, abs=1e-8)


class TestBubblePresIdeal(object):
    @pytest.mark.unit
    def test_build(self, frame):
        frame.props[1].pressure_bubble = Var(frame.params._pe_pairs)
        frame.props[1]._mole_frac_pbub = Var(
            frame.params._pe_pairs, frame.params.component_list, initialize=0.5
        )

        IdealBubbleDew.pressure_bubble(frame.props[1])

        assert isinstance(frame.props[1].eq_pressure_bubble, Constraint)
        assert len(frame.props[1].eq_pressure_bubble) == 1

        assert isinstance(frame.props[1].eq_mole_frac_pbub, Constraint)
        assert len(frame.props[1].eq_mole_frac_pbub) == 2
        for k in frame.props[1].eq_mole_frac_pbub:
            assert k in [("Vap", "Liq", "H2O"), ("Vap", "Liq", "EtOH")]

    @pytest.mark.unit
    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1 / 10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1 / 10

            for pp in frame.params._pe_pairs:
                frame.props[1].pressure_bubble[pp[0], pp[1]] = value(
                    sum(
                        frame.props[1].mole_frac_comp[j] * Psat[j]
                        for j in frame.params.component_list
                    )
                )

                for j in frame.params.component_list:
                    frame.props[1]._mole_frac_pbub[pp[0], pp[1], j].value = value(
                        frame.props[1].mole_frac_comp[j]
                        * Psat[j]
                        / frame.props[1].pressure_bubble[pp[0], pp[1]]
                    )

                assert value(
                    frame.props[1].eq_pressure_bubble[pp[0], pp[1]].body
                ) == pytest.approx(0, abs=1e-8)
                for k in frame.params.component_list:
                    assert value(
                        frame.props[1].eq_mole_frac_pbub[pp[0], pp[1], k].body
                    ) == pytest.approx(0, abs=1e-8)


class TestDewPressureIdeal(object):
    @pytest.mark.unit
    def test_build(self, frame):
        frame.props[1].pressure_dew = Var(frame.params._pe_pairs)
        frame.props[1]._mole_frac_pdew = Var(
            frame.params._pe_pairs, frame.params.component_list, initialize=0.5
        )

        IdealBubbleDew.pressure_dew(frame.props[1])

        assert isinstance(frame.props[1].eq_pressure_dew, Constraint)
        assert len(frame.props[1].eq_pressure_dew) == 1

        assert isinstance(frame.props[1].eq_mole_frac_pdew, Constraint)
        assert len(frame.props[1].eq_mole_frac_pdew) == 2
        for k in frame.props[1].eq_mole_frac_pdew:
            assert k in [("Vap", "Liq", "H2O"), ("Vap", "Liq", "EtOH")]

    @pytest.mark.unit
    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1 / 10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1 / 10

            for pp in frame.params._pe_pairs:
                frame.props[1].pressure_dew[pp[0], pp[1]] = value(
                    1
                    / sum(
                        frame.props[1].mole_frac_comp[j] / Psat[j]
                        for j in frame.params.component_list
                    )
                )

                for j in frame.params.component_list:
                    frame.props[1]._mole_frac_pdew[pp[0], pp[1], j].value = value(
                        frame.props[1].mole_frac_comp[j]
                        * frame.props[1].pressure_dew[pp[0], pp[1]]
                        / Psat[j]
                    )

                assert value(
                    frame.props[1].eq_pressure_dew[pp[0], pp[1]].body
                ) == pytest.approx(0, abs=1e-8)
                for k in frame.params.component_list:
                    assert value(
                        frame.props[1].eq_mole_frac_pdew[pp[0], pp[1], k].body
                    ) == pytest.approx(0, abs=1e-8)
