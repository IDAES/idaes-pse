#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for methods for calculating bubble and dew points

Authors: Andrew Lee, Douglas Allan
"""

import pytest
from sys import modules

from pyomo.environ import ConcreteModel, Constraint, Set, value, Var, units as pyunits

from idaes.core.base.phases import PhaseType

from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    IdealBubbleDew,
    IdealBubbleDewScaler,
)
from idaes.core import declare_process_block_class
from idaes.models.properties.modular_properties.base.generic_property import (
    GenericParameterData,
)
from idaes.models.properties.modular_properties.base.tests.dummy_eos import DummyEoS
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.scaling.util import set_scaling_factor


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

    set_scaling_factor(m.props[1].pressure, 1e-5)
    set_scaling_factor(m.props[1].mole_frac_comp["H2O"], 7)
    set_scaling_factor(m.props[1].mole_frac_comp["EtOH"], 11)

    return m


@pytest.fixture(scope="class")
def frame_inert():
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
            "Sold": {"equation_of_state": DummyEoS},
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

    return m


@pytest.fixture(scope="class")
def frame_noncondensable():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(
        components={
            "H2O": {"pressure_sat_comp": pressure_sat_comp},
            "N2": {"valid_phase_types": [PhaseType.vaporPhase]},
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

    set_scaling_factor(m.props[1].pressure, 1e-5)
    set_scaling_factor(m.props[1].mole_frac_comp["H2O"], 7)
    set_scaling_factor(m.props[1].mole_frac_comp["N2"], 11)

    return m


@pytest.fixture(scope="class")
def frame_nonvolatile():
    m = ConcreteModel()

    # Dummy params block
    m.params = DummyParameterBlock(
        components={
            "H2O": {"pressure_sat_comp": pressure_sat_comp},
            "NaCl": {"valid_phase_types": [PhaseType.liquidPhase]},
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

    set_scaling_factor(m.props[1].pressure, 1e-5)
    set_scaling_factor(m.props[1].mole_frac_comp["H2O"], 7)
    set_scaling_factor(m.props[1].mole_frac_comp["NaCl"], 11)

    return m


@pytest.mark.unit
def test_default_scaler():
    assert IdealBubbleDew.default_scaler == IdealBubbleDewScaler


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

    @pytest.mark.unit
    def test_scaling(self, frame):
        blk = frame.props[1]
        assert len(blk.scaling_factor) == 3
        scaler_obj = IdealBubbleDewScaler()

        scaler_obj.variable_scaling_routine(blk)
        # No variables to scale
        assert len(blk.scaling_factor) == 3

        scaler_obj.constraint_scaling_routine(blk)

        assert len(blk.scaling_factor) == 6
        assert blk.scaling_factor[
            blk.eq_mole_frac_tbub["Vap", "Liq", "H2O"]
        ] == pytest.approx(7e-5)
        assert blk.scaling_factor[
            blk.eq_mole_frac_tbub["Vap", "Liq", "EtOH"]
        ] == pytest.approx(11e-5)
        assert blk.scaling_factor[blk.eq_temperature_bubble["Vap", "Liq"]] == 1e-5

    @pytest.mark.unit
    def test_inert_phases(self, frame_inert):
        with pytest.raises(
            ConfigurationError,
            match="Ideal assumption for calculating bubble and/or dew points is only valid "
            "for systems with two phases. Please use LogBubbleDew approach instead.",
        ):
            IdealBubbleDew.temperature_bubble(frame_inert.props[1])

    @pytest.mark.unit
    def test_build_noncondensable(self, frame_noncondensable):
        params = frame_noncondensable.params
        blk = frame_noncondensable.props[1]
        blk.temperature_dew = Var(params._pe_pairs)
        blk.mole_frac_tdew = Var(
            params._pe_pairs, params.component_list, initialize=0.5
        )

        IdealBubbleDew.temperature_bubble(blk)

        # Constraints skipped
        assert isinstance(blk.eq_temperature_bubble, Constraint)
        assert len(blk.eq_temperature_bubble) == 0

        assert isinstance(blk.eq_mole_frac_tbub, Constraint)
        assert len(blk.eq_mole_frac_tbub) == 0

    @pytest.mark.unit
    def test_scale_noncondensable(self, frame_noncondensable):
        blk = frame_noncondensable.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        # No variables to scale
        scaler_obj.variable_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3

        # Skipped constraints don't get scaling factors
        scaler_obj.constraint_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3

    @pytest.mark.unit
    def test_build_nonvolatile(self, frame_nonvolatile):
        blk = frame_nonvolatile.props[1]
        params = frame_nonvolatile.params

        blk.temperature_bubble = Var(params._pe_pairs)
        blk._mole_frac_tbub = Var(
            params._pe_pairs, params.component_list, initialize=0.5
        )

        IdealBubbleDew.temperature_bubble(blk)

        assert isinstance(blk.eq_temperature_bubble, Constraint)
        assert len(blk.eq_temperature_bubble) == 1

        assert isinstance(blk.eq_mole_frac_tbub, Constraint)
        assert len(blk.eq_mole_frac_tbub) == 2
        for k in blk.eq_mole_frac_tbub:
            assert k in [("Vap", "Liq", "H2O"), ("Vap", "Liq", "NaCl")]

    @pytest.mark.unit
    def test_expressions_nonvolatile(self, frame_nonvolatile):
        blk = frame_nonvolatile.props[1]
        params = frame_nonvolatile.params
        for x1 in range(0, 11, 1):
            blk.mole_frac_comp["H2O"].value = x1 / 10
            blk.mole_frac_comp["NaCl"].value = 1 - x1 / 10

            blk.pressure = value(blk.mole_frac_comp["H2O"] * Psat["H2O"])

            for pp in params._pe_pairs:
                blk._mole_frac_tbub[pp[0], pp[1], "H2O"].value = 1
                blk._mole_frac_tbub[pp[0], pp[1], "NaCl"].value = 0

                assert value(
                    blk.eq_temperature_bubble[pp[0], pp[1]].body
                ) == pytest.approx(0, abs=1e-8)
                for k in params.component_list:
                    assert value(
                        blk.eq_mole_frac_tbub[pp[0], pp[1], k].body
                    ) == pytest.approx(0, abs=1e-8)

    @pytest.mark.unit
    def test_scaling_nonvolatile(self, frame_nonvolatile):
        blk = frame_nonvolatile.props[1]
        assert len(blk.scaling_factor) == 3
        scaler_obj = IdealBubbleDewScaler()

        scaler_obj.variable_scaling_routine(blk)
        # No variables to scale
        assert len(blk.scaling_factor) == 3

        scaler_obj.constraint_scaling_routine(blk)

        assert len(blk.scaling_factor) == 6
        assert blk.scaling_factor[
            blk.eq_mole_frac_tbub["Vap", "Liq", "H2O"]
        ] == pytest.approx(7e-5)
        assert blk.scaling_factor[
            blk.eq_mole_frac_tbub["Vap", "Liq", "NaCl"]
        ] == pytest.approx(11)
        assert blk.scaling_factor[blk.eq_temperature_bubble["Vap", "Liq"]] == 1e-5


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

    @pytest.mark.unit
    def test_scaling(self, frame):
        blk = frame.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        scaler_obj.variable_scaling_routine(blk)
        # No variables to scale
        assert len(blk.scaling_factor) == 3

        scaler_obj.constraint_scaling_routine(blk)

        assert len(blk.scaling_factor) == 6
        assert blk.scaling_factor[
            blk.eq_mole_frac_tdew["Vap", "Liq", "H2O"]
        ] == pytest.approx(7e-5)
        assert blk.scaling_factor[
            blk.eq_mole_frac_tdew["Vap", "Liq", "EtOH"]
        ] == pytest.approx(11e-5)
        assert blk.scaling_factor[blk.eq_temperature_dew["Vap", "Liq"]] == 1

    @pytest.mark.unit
    def test_inert_phases(self, frame_inert):
        with pytest.raises(
            ConfigurationError,
            match="Ideal assumption for calculating bubble and/or dew points is only valid "
            "for systems with two phases. Please use LogBubbleDew approach instead.",
        ):
            IdealBubbleDew.temperature_dew(frame_inert.props[1])

    @pytest.mark.unit
    def test_build_nonvolatile(self, frame_nonvolatile):
        params = frame_nonvolatile.params
        blk = frame_nonvolatile.props[1]
        blk.temperature_dew = Var(params._pe_pairs)
        blk.mole_frac_tdew = Var(
            params._pe_pairs, params.component_list, initialize=0.5
        )

        IdealBubbleDew.temperature_dew(blk)

        # Constraints skipped
        assert isinstance(blk.eq_temperature_dew, Constraint)
        assert len(blk.eq_temperature_dew) == 0

        assert isinstance(blk.eq_mole_frac_tdew, Constraint)
        assert len(blk.eq_mole_frac_tdew) == 0

    @pytest.mark.unit
    def test_scale_nonvolatile(self, frame_nonvolatile):
        blk = frame_nonvolatile.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        # No variables to scale
        scaler_obj.variable_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3

        # Skipped constraints don't get scaling factors
        scaler_obj.constraint_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3

    @pytest.mark.unit
    def test_build_noncondensable(self, frame_noncondensable):
        params = frame_noncondensable.params
        blk = frame_noncondensable.props[1]
        blk.temperature_dew = Var(params._pe_pairs)
        blk._mole_frac_tdew = Var(
            params._pe_pairs, params.component_list, initialize=0.5
        )

        IdealBubbleDew.temperature_dew(blk)

        assert isinstance(blk.eq_temperature_dew, Constraint)
        assert len(blk.eq_temperature_dew) == 1

        assert isinstance(blk.eq_mole_frac_tdew, Constraint)
        assert len(blk.eq_mole_frac_tdew) == 2
        for k in blk.eq_mole_frac_tdew:
            assert k in [("Vap", "Liq", "H2O"), ("Vap", "Liq", "N2")]

    @pytest.mark.unit
    def test_expressions_noncondensable(self, frame_noncondensable):
        blk = frame_noncondensable.props[1]
        params = frame_noncondensable.params
        for x1 in range(1, 11, 1):
            blk.mole_frac_comp["H2O"].value = x1 / 10
            blk.mole_frac_comp["N2"].value = 1 - x1 / 10

            blk.pressure = value(Psat["H2O"] / blk.mole_frac_comp["H2O"])

            for pp in params._pe_pairs:
                blk._mole_frac_tdew[pp[0], pp[1], "H2O"].value = value(
                    blk.mole_frac_comp["H2O"] * blk.pressure / Psat["H2O"]
                )
                blk._mole_frac_tdew[pp[0], pp[1], "N2"].value = 0

                assert value(
                    blk.eq_temperature_dew[pp[0], pp[1]].body
                ) == pytest.approx(0, abs=1e-8)
                for k in params.component_list:
                    assert value(
                        blk.eq_mole_frac_tdew[pp[0], pp[1], k].body
                    ) == pytest.approx(0, abs=1e-8)

    @pytest.mark.unit
    def test_scaling_noncondensable(self, frame_noncondensable):
        blk = frame_noncondensable.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        scaler_obj.variable_scaling_routine(blk)
        # No variables to scale
        assert len(blk.scaling_factor) == 3

        scaler_obj.constraint_scaling_routine(blk)

        assert len(blk.scaling_factor) == 6
        assert blk.scaling_factor[
            blk.eq_mole_frac_tdew["Vap", "Liq", "H2O"]
        ] == pytest.approx(7e-5)
        assert blk.scaling_factor[blk.eq_mole_frac_tdew["Vap", "Liq", "N2"]] == 11
        assert blk.scaling_factor[blk.eq_temperature_dew["Vap", "Liq"]] == 1


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

    @pytest.mark.unit
    def test_scaling(self, frame):
        blk = frame.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        scaler_obj.variable_scaling_routine(blk)
        # No variables to scale
        assert len(blk.scaling_factor) == 3

        scaler_obj.constraint_scaling_routine(blk)

        assert len(blk.scaling_factor) == 6
        assert blk.scaling_factor[
            blk.eq_mole_frac_pbub["Vap", "Liq", "H2O"]
        ] == pytest.approx(7e-5)
        assert blk.scaling_factor[
            blk.eq_mole_frac_pbub["Vap", "Liq", "EtOH"]
        ] == pytest.approx(11e-5)
        assert blk.scaling_factor[blk.eq_pressure_bubble["Vap", "Liq"]] == 1e-5

    @pytest.mark.unit
    def test_inert_phases(self, frame_inert):
        with pytest.raises(
            ConfigurationError,
            match="Ideal assumption for calculating bubble and/or dew points is only valid "
            "for systems with two phases. Please use LogBubbleDew approach instead.",
        ):
            IdealBubbleDew.pressure_bubble(frame_inert.props[1])

    # TODO test nonvolatiles---currently broken due to #1665

    @pytest.mark.unit
    def test_build_noncondensable(self, frame_noncondensable):
        params = frame_noncondensable.params
        blk = frame_noncondensable.props[1]
        blk.pressure_bubble = Var(params._pe_pairs)
        blk.mole_frac_pbub = Var(
            params._pe_pairs, params.component_list, initialize=0.5
        )

        IdealBubbleDew.pressure_bubble(blk)

        # Constraints skipped
        assert isinstance(blk.eq_pressure_bubble, Constraint)
        assert len(blk.eq_pressure_bubble) == 0

        assert isinstance(blk.eq_mole_frac_pbub, Constraint)
        assert len(blk.eq_mole_frac_pbub) == 0

    @pytest.mark.unit
    def test_scale_noncondensable(self, frame_noncondensable):
        blk = frame_noncondensable.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        # No variables to scale
        scaler_obj.variable_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3

        # Skipped constraints don't get scaling factors
        scaler_obj.constraint_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3


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

    @pytest.mark.unit
    def test_scaling(self, frame):
        blk = frame.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        scaler_obj.variable_scaling_routine(blk)
        # No variables to scale
        assert len(blk.scaling_factor) == 3

        scaler_obj.constraint_scaling_routine(blk)

        assert len(blk.scaling_factor) == 6
        assert blk.scaling_factor[
            blk.eq_mole_frac_pdew["Vap", "Liq", "H2O"]
        ] == pytest.approx(7e-5)
        assert blk.scaling_factor[
            blk.eq_mole_frac_pdew["Vap", "Liq", "EtOH"]
        ] == pytest.approx(11e-5)
        assert blk.scaling_factor[blk.eq_pressure_dew["Vap", "Liq"]] == 1

    @pytest.mark.unit
    def test_inert_phases(self, frame_inert):
        with pytest.raises(
            ConfigurationError,
            match="Ideal assumption for calculating bubble and/or dew points is only valid "
            "for systems with two phases. Please use LogBubbleDew approach instead.",
        ):
            IdealBubbleDew.pressure_dew(frame_inert.props[1])

    @pytest.mark.unit
    def test_build_nonvolatile(self, frame_nonvolatile):
        params = frame_nonvolatile.params
        blk = frame_nonvolatile.props[1]
        blk.pressure_dew = Var(params._pe_pairs)
        blk.mole_frac_pdew = Var(
            params._pe_pairs, params.component_list, initialize=0.5
        )

        IdealBubbleDew.pressure_dew(blk)

        # Constraints skipped
        assert isinstance(blk.eq_pressure_dew, Constraint)
        assert len(blk.eq_pressure_dew) == 0

        assert isinstance(blk.eq_mole_frac_pdew, Constraint)
        assert len(blk.eq_mole_frac_pdew) == 0

    @pytest.mark.unit
    def test_scale_nonvolatile(self, frame_nonvolatile):
        blk = frame_nonvolatile.props[1]
        assert len(blk.scaling_factor) == 3

        scaler_obj = IdealBubbleDewScaler()

        # No variables to scale
        scaler_obj.variable_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3

        # Skipped constraints don't get scaling factors
        scaler_obj.constraint_scaling_routine(blk)
        assert len(blk.scaling_factor) == 3

    # TODO Test noncondensables. Currently broken due to #1665
    # @pytest.mark.unit
    # def test_build_noncondensable(self, frame_noncondensable):
    #     params = frame_noncondensable.params
    #     blk = frame_noncondensable.props[1]
    #     blk.pressure_dew = Var(params._pe_pairs)
    #     blk._mole_frac_pdew = Var(
    #         params._pe_pairs, params.component_list, initialize=0.5
    #     )

    #     IdealBubbleDew.pressure_dew(blk)

    #     assert isinstance(blk.eq_pressure_dew, Constraint)
    #     assert len(blk.eq_pressure_dew) == 1

    #     assert isinstance(blk.eq_mole_frac_pdew, Constraint)
    #     assert len(blk.eq_mole_frac_pdew) == 2
    #     for k in blk.eq_mole_frac_pdew:
    #         assert k in [("Vap", "Liq", "H2O"), ("Vap", "Liq", "N2")]

    # @pytest.mark.unit
    # def test_expressions_noncondensable(self, frame_noncondensable):
    #     blk = frame_noncondensable.props[1]
    #     params = frame_noncondensable.params
    #     for x1 in range(1, 11, 1):
    #         blk.mole_frac_comp["H2O"].value = x1 / 10
    #         blk.mole_frac_comp["N2"].value = 1 - x1 / 10

    #         for pp in params._pe_pairs:
    #             blk.pressure_dew[pp[0], pp[1]] = value(
    #                 Psat["H2O"] / blk.mole_frac_comp["H2O"]
    #             )

    #             blk._mole_frac_pdew[pp[0], pp[1], "H2O"].value = 1
    #             blk._mole_frac_pdew[pp[0], pp[1], "N2"].value = 0

    #             assert value(
    #                 blk.eq_pressure_dew[pp[0], pp[1]].body
    #             ) == pytest.approx(0, abs=1e-8)
    #             for k in params.component_list:
    #                 assert value(
    #                     blk.eq_mole_frac_pdew[pp[0], pp[1], k].body
    #                 ) == pytest.approx(0, abs=1e-8)

    # @pytest.mark.unit
    # def test_scaling_noncondensable(self, frame_noncondensable):
    #     blk = frame_noncondensable.props[1]
    #     assert len(blk.scaling_factor) == 3

    #     scaler_obj = IdealBubbleDewScaler()

    #     scaler_obj.variable_scaling_routine(blk)
    #     # No variables to scale
    #     assert len(blk.scaling_factor) == 3

    #     scaler_obj.constraint_scaling_routine(blk)

    #     assert len(blk.scaling_factor) == 6
    #     assert blk.scaling_factor[blk.eq_mole_frac_pdew["Vap","Liq","H2O"]] == pytest.approx(7e-5)
    #     assert blk.scaling_factor[blk.eq_mole_frac_pdew["Vap","Liq","N2"]] == 11
    #     assert blk.scaling_factor[blk.eq_pressure_dew["Vap","Liq"]] == 1
