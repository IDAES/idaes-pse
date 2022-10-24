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
Tests for Component objects

Author: Andrew Lee
"""
import pytest
import types

from pyomo.environ import ConcreteModel, Set, Param, Var, units as pyunits

from idaes.core.base.components import (
    Component,
    Solute,
    Solvent,
    Ion,
    Anion,
    Cation,
    Apparent,
)
from idaes.core.base.phases import (
    LiquidPhase,
    VaporPhase,
    SolidPhase,
    Phase,
    PhaseType,
    AqueousPhase,
)
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
from idaes.core.base.property_meta import PropertyClassMetadata, UnitSet


class TestComponent:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        m.comp = Component()
        m.comp2 = Component()

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert m.comp.config.valid_phase_types is None
        assert m.comp.config.elemental_composition is None
        assert not m.comp.config._component_list_exists
        assert m.comp.config.henry_component is None
        assert m.comp.config.has_vapor_pressure

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.component_list, Set)

        for j in m.component_list:
            assert j in ["comp", "comp2"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        with pytest.raises(
            TypeError,
            match="comp Generic Component objects do not "
            "support is_solute\(\) method. Use a Solvent or "
            "Solute Component instead.",
        ):
            m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        with pytest.raises(
            TypeError,
            match="comp Generic Component objects do not "
            "support is_solvent\(\) method. Use a Solvent or "
            "Solute Component instead.",
        ):
            m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(TypeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.comp3 = Component(valid_phase_types=PhaseType.liquidPhase)

        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Aqu = AqueousPhase()
        m.Phase = Phase()

        assert m.comp3._is_phase_valid(m.Liq)
        assert not m.comp3._is_phase_valid(m.Sol)
        assert not m.comp3._is_phase_valid(m.Vap)
        assert not m.comp3._is_phase_valid(m.Aqu)
        assert not m.comp3._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_vapor(self, m):
        m.comp4 = Component(valid_phase_types=PhaseType.vaporPhase)

        assert not m.comp4._is_phase_valid(m.Liq)
        assert not m.comp4._is_phase_valid(m.Sol)
        assert m.comp4._is_phase_valid(m.Vap)
        assert not m.comp4._is_phase_valid(m.Aqu)
        assert not m.comp4._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_solid(self, m):
        m.comp5 = Component(valid_phase_types=PhaseType.solidPhase)

        assert not m.comp5._is_phase_valid(m.Liq)
        assert m.comp5._is_phase_valid(m.Sol)
        assert not m.comp5._is_phase_valid(m.Vap)
        assert not m.comp5._is_phase_valid(m.Aqu)
        assert not m.comp5._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_aqueous(self, m):
        m.comp6 = Component(valid_phase_types=PhaseType.aqueousPhase)

        assert not m.comp6._is_phase_valid(m.Liq)
        assert not m.comp6._is_phase_valid(m.Sol)
        assert not m.comp6._is_phase_valid(m.Vap)
        # Generic components are never valid in the aqueous phase
        assert not m.comp6._is_phase_valid(m.Aqu)
        assert not m.comp6._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_LV(self, m):
        m.comp7 = Component(
            valid_phase_types=[PhaseType.liquidPhase, PhaseType.vaporPhase]
        )

        assert m.comp7._is_phase_valid(m.Liq)
        assert not m.comp7._is_phase_valid(m.Sol)
        assert m.comp7._is_phase_valid(m.Vap)
        assert not m.comp7._is_phase_valid(m.Aqu)
        assert not m.comp7._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_create_parameters(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        m.comp = Component(
            parameter_data={
                "mw": 10,
                "pressure_crit": 1e5,
                "temperature_crit": 500,
            }
        )

        assert isinstance(m.comp.mw, Param)
        assert m.comp.mw.value == 10

        assert isinstance(m.comp.pressure_crit, Var)
        assert m.comp.pressure_crit.value == 1e5

        assert isinstance(m.comp.temperature_crit, Var)
        assert m.comp.temperature_crit.value == 500

    @pytest.mark.unit
    def test_create_parameters_convert(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        m.comp = Component(
            parameter_data={
                "mw": (10, pyunits.g / pyunits.mol),
                "pressure_crit": (1, pyunits.bar),
                "temperature_crit": (900, pyunits.degR),
            }
        )

        assert isinstance(m.comp.mw, Param)
        assert m.comp.mw.value == 1e-2

        assert isinstance(m.comp.pressure_crit, Var)
        assert m.comp.pressure_crit.value == 1e5

        assert isinstance(m.comp.temperature_crit, Var)
        assert m.comp.temperature_crit.value == 500


class TestSolute:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

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

        assert isinstance(m.solute_set, Set)

        for j in m.solute_set:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(TypeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.comp3 = Solute(valid_phase_types=PhaseType.liquidPhase)

        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Aqu = AqueousPhase()
        m.Phase = Phase()

        assert m.comp3._is_phase_valid(m.Liq)
        assert not m.comp3._is_phase_valid(m.Sol)
        assert not m.comp3._is_phase_valid(m.Vap)
        assert not m.comp3._is_phase_valid(m.Aqu)
        assert not m.comp3._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_vapor(self, m):
        m.comp4 = Solute(valid_phase_types=PhaseType.vaporPhase)

        assert not m.comp4._is_phase_valid(m.Liq)
        assert not m.comp4._is_phase_valid(m.Sol)
        assert m.comp4._is_phase_valid(m.Vap)
        assert not m.comp4._is_phase_valid(m.Aqu)
        assert not m.comp4._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_solid(self, m):
        m.comp5 = Solute(valid_phase_types=PhaseType.solidPhase)

        assert not m.comp5._is_phase_valid(m.Liq)
        assert m.comp5._is_phase_valid(m.Sol)
        assert not m.comp5._is_phase_valid(m.Vap)
        assert not m.comp5._is_phase_valid(m.Aqu)
        assert not m.comp5._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_aqueous(self, m):
        m.comp6 = Solute(valid_phase_types=PhaseType.aqueousPhase)

        assert not m.comp6._is_phase_valid(m.Liq)
        assert not m.comp6._is_phase_valid(m.Sol)
        assert not m.comp6._is_phase_valid(m.Vap)
        assert m.comp6._is_phase_valid(m.Aqu)
        assert not m.comp6._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_LV(self, m):
        m.comp7 = Solute(
            valid_phase_types=[PhaseType.liquidPhase, PhaseType.vaporPhase]
        )

        assert m.comp7._is_phase_valid(m.Liq)
        assert not m.comp7._is_phase_valid(m.Sol)
        assert m.comp7._is_phase_valid(m.Vap)
        assert not m.comp7._is_phase_valid(m.Aqu)
        assert not m.comp7._is_phase_valid(m.Phase)


class TestSovent:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

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

        assert isinstance(m.solvent_set, Set)

        for j in m.solvent_set:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert not m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(TypeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.comp3 = Solvent(valid_phase_types=PhaseType.liquidPhase)

        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Aqu = AqueousPhase()
        m.Phase = Phase()

        assert m.comp3._is_phase_valid(m.Liq)
        assert not m.comp3._is_phase_valid(m.Sol)
        assert not m.comp3._is_phase_valid(m.Vap)
        assert not m.comp3._is_phase_valid(m.Aqu)
        assert not m.comp3._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_vapor(self, m):
        m.comp4 = Solvent(valid_phase_types=PhaseType.vaporPhase)

        assert not m.comp4._is_phase_valid(m.Liq)
        assert not m.comp4._is_phase_valid(m.Sol)
        assert m.comp4._is_phase_valid(m.Vap)
        assert not m.comp4._is_phase_valid(m.Aqu)
        assert not m.comp4._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_solid(self, m):
        m.comp5 = Solvent(valid_phase_types=PhaseType.solidPhase)

        assert not m.comp5._is_phase_valid(m.Liq)
        assert m.comp5._is_phase_valid(m.Sol)
        assert not m.comp5._is_phase_valid(m.Vap)
        assert not m.comp5._is_phase_valid(m.Aqu)
        assert not m.comp5._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_LV(self, m):
        m.comp6 = Solvent(
            valid_phase_types=[PhaseType.liquidPhase, PhaseType.vaporPhase]
        )

        assert m.comp6._is_phase_valid(m.Liq)
        assert not m.comp6._is_phase_valid(m.Sol)
        assert m.comp6._is_phase_valid(m.Vap)
        assert not m.comp6._is_phase_valid(m.Aqu)
        assert not m.comp6._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_aqueous(self, m):
        m.comp7 = Solvent(valid_phase_types=PhaseType.aqueousPhase)

        assert not m.comp7._is_phase_valid(m.Liq)
        assert not m.comp7._is_phase_valid(m.Sol)
        assert not m.comp7._is_phase_valid(m.Vap)
        assert m.comp7._is_phase_valid(m.Aqu)
        assert not m.comp7._is_phase_valid(m.Phase)


class TestIon:
    @pytest.mark.unit
    def test_electrolye(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        with pytest.raises(
            NotImplementedError,
            match="comp The IonData component class is intended as a base "
            "class for the AnionData and CationData classes, and should "
            "not be used directly",
        ):
            m.comp = Ion(_electrolyte=True)

    @pytest.mark.unit
    def test_not_electrolye(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        with pytest.raises(
            PropertyPackageError,
            match="comp Ion Component types should only be used with " "Aqueous Phases",
        ):
            m.comp = Ion(_electrolyte=False)


class TestAnion:
    @pytest.mark.unit
    def test_not_electrolye(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        with pytest.raises(
            PropertyPackageError,
            match="comp Ion Component types should only be used with " "Aqueous Phases",
        ):
            m.comp = Anion(_electrolyte=False)

    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        m.anion_set = Set()

        m.comp = Anion(charge=-1, _electrolyte=True)

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert "valid_phase_types" not in m.comp.config
        assert m.comp.config.charge == -1
        assert not m.comp.config._component_list_exists
        assert not m.comp.config.has_vapor_pressure
        with pytest.raises(ValueError):
            m.comp.config.has_vapor_pressure = True

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.anion_set, Set)

        for j in m.anion_set:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_invalid_charge(self, m):
        with pytest.raises(
            ConfigurationError,
            match="an received invalid value for charge "
            "configuration argument."
            " Anions must have a negative charge.",
        ):
            m.an = Anion(charge=+1, _electrolyte=True)

    @pytest.mark.unit
    def test_no_charge(self, m):
        with pytest.raises(
            ConfigurationError, match="an was not provided with a value " "for charge."
        ):
            m.an = Anion(_electrolyte=True)

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(AttributeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Aqu = AqueousPhase()
        m.Phase = Phase()

        assert not m.comp._is_phase_valid(m.Liq)
        assert not m.comp._is_phase_valid(m.Sol)
        assert not m.comp._is_phase_valid(m.Vap)
        assert m.comp._is_phase_valid(m.Aqu)
        assert not m.comp._is_phase_valid(m.Phase)


class TestCation:
    @pytest.mark.unit
    def test_not_electrolye(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        with pytest.raises(
            PropertyPackageError,
            match="comp Ion Component types should only be used with " "Aqueous Phases",
        ):
            m.comp = Cation(_electrolyte=False)

    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        m.cation_set = Set()

        m.comp = Cation(charge=+1, _electrolyte=True)

        return m

    @pytest.mark.unit
    def test_config(self, m):
        assert "valid_phase_types" not in m.comp.config
        assert m.comp.config.charge == +1
        assert not m.comp.config._component_list_exists
        assert not m.comp.config.has_vapor_pressure
        with pytest.raises(ValueError):
            m.comp.config.has_vapor_pressure = True

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m.cation_set, Set)

        for j in m.cation_set:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_invalid_charge(self, m):
        with pytest.raises(
            ConfigurationError,
            match="cat received invalid value for charge "
            "configuration argument."
            " Cations must have a positive charge.",
        ):
            m.cat = Cation(charge=-1, _electrolyte=True)

    @pytest.mark.unit
    def test_no_charge(self, m):
        with pytest.raises(
            ConfigurationError, match="cat was not provided with a value " "for charge."
        ):
            m.cat = Cation(_electrolyte=True)

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(AttributeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Aqu = AqueousPhase()
        m.Phase = Phase()

        assert not m.comp._is_phase_valid(m.Liq)
        assert not m.comp._is_phase_valid(m.Sol)
        assert not m.comp._is_phase_valid(m.Vap)
        assert m.comp._is_phase_valid(m.Aqu)
        assert not m.comp._is_phase_valid(m.Phase)


class TestApparent:
    @pytest.fixture(scope="class")
    def m(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        m._apparent_set = Set()

        m.comp = Apparent(dissociation_species={"comp": 1}, _electrolyte=True)

        return m

    @pytest.mark.unit
    def test_not_electrolye(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        m.comp = Apparent(dissociation_species={"comp": 1}, _electrolyte=False)

        for j in m.component_list:
            assert j in ["comp"]
        for j in m.solute_set:
            assert j in ["comp"]

        assert m.comp.config.dissociation_species == {"comp": 1}

    @pytest.mark.unit
    def test_no_dissociation_species(self):
        m = ConcreteModel()

        m.meta_object = PropertyClassMetadata()
        m.meta_object._default_units = UnitSet(
            temperature=pyunits.K,
            mass=pyunits.kg,
            length=pyunits.m,
            time=pyunits.s,
            amount=pyunits.mol,
        )

        def get_metadata(self):
            return m.meta_object

        m.get_metadata = types.MethodType(get_metadata, m)

        with pytest.raises(
            ConfigurationError,
            match="dissoication_species argument was not set. "
            "Apparent components require the dissociation species to be "
            "defined.",
        ):
            m.comp = Apparent(_electrolyte=False)

    @pytest.mark.unit
    def test_config(self, m):
        assert m.comp.config.valid_phase_types is None
        assert not m.comp.config._component_list_exists

    @pytest.mark.unit
    def test_populate_component_list(self, m):
        assert isinstance(m._apparent_set, Set)

        for j in m._apparent_set:
            assert j in ["comp"]

    @pytest.mark.unit
    def test_is_solute(self, m):
        assert m.comp.is_solute()

    @pytest.mark.unit
    def test_is_solvent(self, m):
        assert not m.comp.is_solvent()

    @pytest.mark.unit
    def test_is_phase_valid_no_assignment(self, m):
        with pytest.raises(TypeError):
            m.comp._is_phase_valid("foo")

    @pytest.mark.unit
    def test_is_phase_valid_liquid(self, m):
        m.comp3 = Apparent(
            valid_phase_types=PhaseType.liquidPhase,
            dissociation_species={"comp": 1},
            _electrolyte=True,
        )

        m.Liq = LiquidPhase()
        m.Sol = SolidPhase()
        m.Vap = VaporPhase()
        m.Aqu = AqueousPhase()
        m.Phase = Phase()

        assert m.comp3._is_phase_valid(m.Liq)
        assert not m.comp3._is_phase_valid(m.Sol)
        assert not m.comp3._is_phase_valid(m.Vap)
        assert not m.comp3._is_phase_valid(m.Aqu)
        assert not m.comp3._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_vapor(self, m):
        m.comp4 = Apparent(
            valid_phase_types=PhaseType.vaporPhase,
            dissociation_species={"comp": 1},
            _electrolyte=True,
        )

        assert not m.comp4._is_phase_valid(m.Liq)
        assert not m.comp4._is_phase_valid(m.Sol)
        assert m.comp4._is_phase_valid(m.Vap)
        assert not m.comp4._is_phase_valid(m.Aqu)
        assert not m.comp4._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_solid(self, m):
        m.comp5 = Apparent(
            valid_phase_types=PhaseType.solidPhase,
            dissociation_species={"comp": 1},
            _electrolyte=True,
        )

        assert not m.comp5._is_phase_valid(m.Liq)
        assert m.comp5._is_phase_valid(m.Sol)
        assert not m.comp5._is_phase_valid(m.Vap)
        assert not m.comp5._is_phase_valid(m.Aqu)
        assert not m.comp5._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_aqueous(self, m):
        m.comp6 = Apparent(
            valid_phase_types=PhaseType.aqueousPhase,
            dissociation_species={"comp": 1},
            _electrolyte=True,
        )

        assert not m.comp6._is_phase_valid(m.Liq)
        assert not m.comp6._is_phase_valid(m.Sol)
        assert not m.comp6._is_phase_valid(m.Vap)
        assert m.comp6._is_phase_valid(m.Aqu)
        assert not m.comp6._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_phase_valid_LV(self, m):
        m.comp7 = Apparent(
            valid_phase_types=[PhaseType.liquidPhase, PhaseType.vaporPhase],
            dissociation_species={"comp": 1},
            _electrolyte=True,
        )

        assert m.comp7._is_phase_valid(m.Liq)
        assert not m.comp7._is_phase_valid(m.Sol)
        assert m.comp7._is_phase_valid(m.Vap)
        assert not m.comp7._is_phase_valid(m.Aqu)
        assert not m.comp7._is_phase_valid(m.Phase)

    @pytest.mark.unit
    def test_is_aqueous_phase_valid(self, m):
        assert m.comp._is_aqueous_phase_valid()
