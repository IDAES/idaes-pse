##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Tests for smooth VLE formulation

Authors: Andrew Lee
"""

import pytest
import types
from sys import modules

from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Block,
                           Set,
                           value,
                           Var)

from idaes.generic_models.properties.core.phase_equil import bubble_dew
from idaes.core import declare_process_block_class
from idaes.generic_models.properties.core.generic.generic_property import (
        GenericParameterData)
from idaes.generic_models.properties.core.generic.tests import dummy_eos


# Dummy class to use for Psat calls
Psat = {"H2O": 1e5, "EtOH": 5e4}


def pressure_sat_comp(b, j, T):
    return Psat[j]


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
    m.params = DummyParameterBlock(default={
                "component_list": ["H2O", "EtOH"],
                "phase_list": ["Liq", "Vap"],
                "state_definition": modules[__name__],
                "equation_of_state": {"Vap": dummy_eos,
                                      "Liq": dummy_eos},
                "pressure_sat_comp": pressure_sat_comp})

    m.props = m.params.state_block_class([1],
                                         default={"defined_state": False,
                                         "parameters": m.params})

    # Add common variables
    m.props[1].pressure = Var(initialize=101325)
    m.props[1].temperature = Var(initialize=300)
    m.props[1]._teq = Var(initialize=300)
    m.props[1].mole_frac_comp = Var(m.params.component_list,
                                    initialize=0.5)

    return m


class TestBubbleTempIdeal(object):
    def test_build(self, frame):
        frame.props[1].temperature_bubble = Var()
        frame.props[1]._mole_frac_tbub = Var(frame.params.component_list,
                                             initialize=0.5)

        bubble_dew.bubble_temp_ideal(frame.props[1])

        assert isinstance(frame.props[1].eq_temperature_bubble, Constraint)
        assert len(frame.props[1].eq_temperature_bubble) == 1

        assert isinstance(frame.props[1].eq_mole_frac_tbub, Constraint)
        assert len(frame.props[1].eq_mole_frac_tbub) == 2
        for k in frame.props[1].eq_mole_frac_tbub:
            assert k in frame.params.component_list

    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1/10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1/10

            frame.props[1].pressure = value(
                sum(frame.props[1].mole_frac_comp[j] * Psat[j]
                    for j in frame.params.component_list))
            for j in frame.params.component_list:
                frame.props[1]._mole_frac_tbub[j].value = value(
                    frame.props[1].mole_frac_comp[j] * Psat[j] /
                    frame.props[1].pressure)

            assert value(frame.props[1].eq_temperature_bubble.body) == \
                pytest.approx(0, abs=1e-8)
            for k in frame.params.component_list:
                assert value(frame.props[1].eq_mole_frac_tbub[k].body) == \
                    pytest.approx(0, abs=1e-8)


class TestDewTempIdeal(object):
    def test_build(self, frame):
        frame.props[1].temperature_dew = Var()
        frame.props[1]._mole_frac_tdew = Var(frame.params.component_list,
                                             initialize=0.5)

        bubble_dew.dew_temp_ideal(frame.props[1])

        assert isinstance(frame.props[1].eq_temperature_dew, Constraint)
        assert len(frame.props[1].eq_temperature_dew) == 1

        assert isinstance(frame.props[1].eq_mole_frac_tdew, Constraint)
        assert len(frame.props[1].eq_mole_frac_tdew) == 2
        for k in frame.props[1].eq_mole_frac_tdew:
            assert k in frame.params.component_list

    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1/10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1/10

            frame.props[1].pressure = value(
                1 / sum(frame.props[1].mole_frac_comp[j] / Psat[j]
                        for j in frame.params.component_list))
            for j in frame.params.component_list:
                frame.props[1]._mole_frac_tdew[j].value = value(
                    frame.props[1].mole_frac_comp[j] *
                    frame.props[1].pressure / Psat[j])

            assert value(frame.props[1].eq_temperature_dew.body) == \
                pytest.approx(0, abs=1e-8)
            for k in frame.params.component_list:
                assert value(frame.props[1].eq_mole_frac_tdew[k].body) == \
                    pytest.approx(0, abs=1e-8)


class TestBubblePresIdeal(object):
    def test_build(self, frame):
        frame.props[1].pressure_bubble = Var()
        frame.props[1]._mole_frac_pbub = Var(frame.params.component_list,
                                             initialize=0.5)

        bubble_dew.bubble_press_ideal(frame.props[1])

        assert isinstance(frame.props[1].eq_pressure_bubble, Constraint)
        assert len(frame.props[1].eq_pressure_bubble) == 1

        assert isinstance(frame.props[1].eq_mole_frac_pbub, Constraint)
        assert len(frame.props[1].eq_mole_frac_pbub) == 2
        for k in frame.props[1].eq_mole_frac_pbub:
            assert k in frame.params.component_list

    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1/10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1/10

            frame.props[1].pressure_bubble = value(
                sum(frame.props[1].mole_frac_comp[j] * Psat[j]
                    for j in frame.params.component_list))
            for j in frame.params.component_list:
                frame.props[1]._mole_frac_pbub[j].value = value(
                    frame.props[1].mole_frac_comp[j] * Psat[j] /
                    frame.props[1].pressure_bubble)

            assert value(frame.props[1].eq_pressure_bubble.body) == \
                pytest.approx(0, abs=1e-8)
            for k in frame.params.component_list:
                assert value(frame.props[1].eq_mole_frac_pbub[k].body) == \
                    pytest.approx(0, abs=1e-8)


class TestDewPressureIdeal(object):
    def test_build(self, frame):
        frame.props[1].pressure_dew = Var()
        frame.props[1]._mole_frac_pdew = Var(frame.params.component_list,
                                             initialize=0.5)

        bubble_dew.dew_press_ideal(frame.props[1])

        assert isinstance(frame.props[1].eq_pressure_dew, Constraint)
        assert len(frame.props[1].eq_pressure_dew) == 1

        assert isinstance(frame.props[1].eq_mole_frac_pdew, Constraint)
        assert len(frame.props[1].eq_mole_frac_pdew) == 2
        for k in frame.props[1].eq_mole_frac_pdew:
            assert k in frame.params.component_list

    def test_expressions(self, frame):
        for x1 in range(0, 11, 1):
            frame.props[1].mole_frac_comp["H2O"].value = x1/10
            frame.props[1].mole_frac_comp["EtOH"].value = 1 - x1/10

            frame.props[1].pressure_dew = value(
                1 / sum(frame.props[1].mole_frac_comp[j] / Psat[j]
                        for j in frame.params.component_list))
            for j in frame.params.component_list:
                frame.props[1]._mole_frac_pdew[j].value = value(
                    frame.props[1].mole_frac_comp[j] *
                    frame.props[1].pressure_dew / Psat[j])

            assert value(frame.props[1].eq_pressure_dew.body) == \
                pytest.approx(0, abs=1e-8)
            for k in frame.params.component_list:
                assert value(frame.props[1].eq_mole_frac_pdew[k].body) == \
                    pytest.approx(0, abs=1e-8)
