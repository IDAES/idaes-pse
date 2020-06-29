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
Tests for Component objects

Author: Andrew Lee
"""
import pytest

from pyomo.environ import ConcreteModel, Set

from idaes.core.components import (Component, Solute, Solvent,
                                   Ion, Anion, Cation)
from idaes.core.phases import (LiquidPhase, VaporPhase, SolidPhase, Phase,
                               PhaseType)
from idaes.core.util.exceptions import ConfigurationError


class TestComponent():
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.comp = Component()
        m.comp2 = Component()

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert m.comp.config.valid_phase_types is None
        assert m.comp.config.elemental_composition is None
        assert not m.comp.config._component_list_exists
        assert m.config.henry_components is None

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.component_list, Set)

        for j in m.component_list:
            assert j in ["comp", "comp2"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        with pytest.raises(TypeError,
                           match="comp Generic Component objects do not "
                           "support is_solute\(\) method. Use a Solvent or "
                           "Solute Component instead."):
            m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        with pytest.raises(TypeError,
                           match="comp Generic Component objects do not "
                           "support is_solvent\(\) method. Use a Solvent or "
                           "Solute Component instead."):
            m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        assert m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.comp3 = Component(default={
            "valid_phase_types": PhaseType.liquidPhase})

        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Phase = Phase()

        assert m.comp3._is_phase_valid(m.Liq)
        assert not m.comp3._is_phase_valid(m.Sol)
        assert not m.comp3._is_phase_valid(m.Vap)
        assert not m.comp3._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_vapor(self, m):
        m.comp4 = Component(default={
            "valid_phase_types": PhaseType.vaporPhase})

        assert not m.comp4._is_phase_valid(m.Liq)
        assert not m.comp4._is_phase_valid(m.Sol)
        assert m.comp4._is_phase_valid(m.Vap)
        assert not m.comp4._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_solid(self, m):
        m.comp5 = Component(default={
            "valid_phase_types": PhaseType.solidPhase})

        assert not m.comp5._is_phase_valid(m.Liq)
        assert m.comp5._is_phase_valid(m.Sol)
        assert not m.comp5._is_phase_valid(m.Vap)
        assert not m.comp5._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_LV(self, m):
        m.comp6 = Component(default={
            "valid_phase_types": [PhaseType.liquidPhase,
                                  PhaseType.vaporPhase]})

        assert m.comp6._is_phase_valid(m.Liq)
        assert not m.comp6._is_phase_valid(m.Sol)
        assert m.comp6._is_phase_valid(m.Vap)
        assert not m.comp6._is_phase_valid(m.Phase)


class TestSolute():
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.comp = Solute()

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert m.comp.config.valid_phase_types is None
        assert not m.comp.config._component_list_exists

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.component_list, Set)

        for j in m.component_list:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        assert m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.comp3 = Component(default={
            "valid_phase_types": PhaseType.liquidPhase})

        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Phase = Phase()

        assert m.comp3._is_phase_valid(m.Liq)
        assert not m.comp3._is_phase_valid(m.Sol)
        assert not m.comp3._is_phase_valid(m.Vap)
        assert not m.comp3._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_vapor(self, m):
        m.comp4 = Component(default={
            "valid_phase_types": PhaseType.vaporPhase})

        assert not m.comp4._is_phase_valid(m.Liq)
        assert not m.comp4._is_phase_valid(m.Sol)
        assert m.comp4._is_phase_valid(m.Vap)
        assert not m.comp4._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_solid(self, m):
        m.comp5 = Component(default={
            "valid_phase_types": PhaseType.solidPhase})

        assert not m.comp5._is_phase_valid(m.Liq)
        assert m.comp5._is_phase_valid(m.Sol)
        assert not m.comp5._is_phase_valid(m.Vap)
        assert not m.comp5._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_LV(self, m):
        m.comp6 = Component(default={
            "valid_phase_types": [PhaseType.liquidPhase,
                                  PhaseType.vaporPhase]})

        assert m.comp6._is_phase_valid(m.Liq)
        assert not m.comp6._is_phase_valid(m.Sol)
        assert m.comp6._is_phase_valid(m.Vap)
        assert not m.comp6._is_phase_valid(m.Phase)


class TestSovent():
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.comp = Solvent()

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert m.comp.config.valid_phase_types is None
        assert not m.comp.config._component_list_exists

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.component_list, Set)

        for j in m.component_list:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert not m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        assert m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.comp3 = Component(default={
            "valid_phase_types": PhaseType.liquidPhase})

        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Phase = Phase()

        assert m.comp3._is_phase_valid(m.Liq)
        assert not m.comp3._is_phase_valid(m.Sol)
        assert not m.comp3._is_phase_valid(m.Vap)
        assert not m.comp3._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_vapor(self, m):
        m.comp4 = Component(default={
            "valid_phase_types": PhaseType.vaporPhase})

        assert not m.comp4._is_phase_valid(m.Liq)
        assert not m.comp4._is_phase_valid(m.Sol)
        assert m.comp4._is_phase_valid(m.Vap)
        assert not m.comp4._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_solid(self, m):
        m.comp5 = Component(default={
            "valid_phase_types": PhaseType.solidPhase})

        assert not m.comp5._is_phase_valid(m.Liq)
        assert m.comp5._is_phase_valid(m.Sol)
        assert not m.comp5._is_phase_valid(m.Vap)
        assert not m.comp5._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_LV(self, m):
        m.comp6 = Component(default={
            "valid_phase_types": [PhaseType.liquidPhase,
                                  PhaseType.vaporPhase]})

        assert m.comp6._is_phase_valid(m.Liq)
        assert not m.comp6._is_phase_valid(m.Sol)
        assert m.comp6._is_phase_valid(m.Vap)
        assert not m.comp6._is_phase_valid(m.Phase)


class TestIon():
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.comp = Ion()

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert "valid_phase_types" not in m.comp.config
        assert m.comp.config.charge is None
        assert not m.comp.config._component_list_exists

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.component_list, Set)

        for j in m.component_list:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(AttributeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Phase = Phase()

        assert m.comp._is_phase_valid(m.Liq)
        assert not m.comp._is_phase_valid(m.Sol)
        assert not m.comp._is_phase_valid(m.Vap)
        assert not m.comp._is_phase_valid(m.Phase)


class TestAnion():
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.comp = Anion(default={"charge": -1})

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert "valid_phase_types" not in m.comp.config
        assert m.comp.config.charge == -1
        assert not m.comp.config._component_list_exists

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.component_list, Set)

        for j in m.component_list:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_invalid_charge(self, m):
        with pytest.raises(ConfigurationError,
                           match="an received invalid value for charge "
                           "configuration argument."
                           " Anions must have a negative charge."):
            m.an = Anion(default={"charge": +1})

    @pytest.mark.unit
    def test_no_charge(self, m):
        with pytest.raises(ConfigurationError,
                           match="an was not provided with a value "
                           "for charge."):
            m.an = Anion()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(AttributeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Phase = Phase()

        assert m.comp._is_phase_valid(m.Liq)
        assert not m.comp._is_phase_valid(m.Sol)
        assert not m.comp._is_phase_valid(m.Vap)
        assert not m.comp._is_phase_valid(m.Phase)


class TestCation():
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.comp = Cation(default={"charge": +1})

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert "valid_phase_types" not in m.comp.config
        assert m.comp.config.charge == +1
        assert not m.comp.config._component_list_exists

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.component_list, Set)

        for j in m.component_list:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_invalid_charge(self, m):
        with pytest.raises(ConfigurationError,
                           match="cat received invalid value for charge "
                           "configuration argument."
                           " Cations must have a positive charge."):
            m.cat = Cation(default={"charge": -1})

    @pytest.mark.unit
    def test_no_charge(self, m):
        with pytest.raises(ConfigurationError,
                           match="cat was not provided with a value "
                           "for charge."):
            m.cat = Cation()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(AttributeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Phase = Phase()

        assert m.comp._is_phase_valid(m.Liq)
        assert not m.comp._is_phase_valid(m.Sol)
        assert not m.comp._is_phase_valid(m.Vap)
        assert not m.comp._is_phase_valid(m.Phase)
