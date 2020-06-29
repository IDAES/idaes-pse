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
Tests for Phase objects

Author: Andrew Lee
"""
from pyomo.environ import ConcreteModel, Set

from idaes.core.phases import (Phase, LiquidPhase, SolidPhase, VaporPhase,
                               PhaseType)
import pytest

@pytest.mark.unit
def test_PhaseType():
    assert len(PhaseType) == 4
    for i in PhaseType.__members__:
        assert i in ["undefined", "liquidPhase", "vaporPhase", "solidPhase"]


@pytest.mark.unit
def test_config():
    m = ConcreteModel()

    m.phase = Phase()

    assert len(m.phase.config) == 5
    for k, v in m.phase.config.items():
        if k == "_phase_list_exists":
            assert not v
        elif k == "parameter_data":
            assert v == {}
        else:
            assert v is None


@pytest.mark.unit
def test_populate_phase_list():
    m = ConcreteModel()

    m.phase = Phase()
    m.phase2 = Phase()

    assert isinstance(m.phase_list, Set)

    for p in m.phase_list:
        assert p in ["phase", "phase2"]


@pytest.mark.unit
def test_is_phase_generic():
    m = ConcreteModel()

    m.phase = Phase()

    assert not m.phase.is_liquid_phase()
    assert not m.phase.is_solid_phase()
    assert not m.phase.is_vapor_phase()


@pytest.mark.unit
def test_is_phase_old_style_liquid():
    m = ConcreteModel()

    m.Liq = Phase()

    assert m.Liq.is_liquid_phase()
    assert not m.Liq.is_solid_phase()
    assert not m.Liq.is_vapor_phase()


@pytest.mark.unit
def test_is_phase_old_style_solid():
    m = ConcreteModel()

    m.Sol = Phase()

    assert not m.Sol.is_liquid_phase()
    assert m.Sol.is_solid_phase()
    assert not m.Sol.is_vapor_phase()


@pytest.mark.unit
def test_is_phase_old_style_vapor():
    m = ConcreteModel()

    m.Vap = Phase()

    assert not m.Vap.is_liquid_phase()
    assert not m.Vap.is_solid_phase()
    assert m.Vap.is_vapor_phase()


@pytest.mark.unit
def test_phase_list_exists():
    m = ConcreteModel()

    m.phase = Phase(default={"_phase_list_exists": True})
    m.phase2 = Phase(default={"_phase_list_exists": True})

    assert not hasattr(m, "phase_list")


@pytest.mark.unit
def test_LiquidPhase():
    m = ConcreteModel()

    m.phase = LiquidPhase()

    assert m.phase.is_liquid_phase()
    assert not m.phase.is_solid_phase()
    assert not m.phase.is_vapor_phase()


@pytest.mark.unit
def test_SolidPhase():
    m = ConcreteModel()

    m.phase = SolidPhase()

    assert not m.phase.is_liquid_phase()
    assert m.phase.is_solid_phase()
    assert not m.phase.is_vapor_phase()


@pytest.mark.unit
def test_VaporPhase():
    m = ConcreteModel()

    m.phase = VaporPhase()

    assert not m.phase.is_liquid_phase()
    assert not m.phase.is_solid_phase()
    assert m.phase.is_vapor_phase()
