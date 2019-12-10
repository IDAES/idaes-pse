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
import pytest

from pyomo.environ import (ConcreteModel,
                           Set,
                           SolverFactory,
                           TerminationCondition,
                           value)

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic_eos.cubic_prop_pack import \
    cubic_roots_available, CubicParameterBlock, CubicStateBlock, CubicEoS


# Set module level pyest marker
pytestmark = pytest.mark.cubic_root
prop_available = cubic_roots_available()


def test_CubicEoS():
    assert len(CubicEoS) == 2
    assert CubicEoS.PR
    assert CubicEoS.SRK


class TestParameterBlock(object):
    def test_build_default(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": True})

        m.fs.params = CubicParameterBlock()

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ('Vap', 'Liq')

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 2
        for p in m.fs.params.phase_list:
            assert p in ["Vap", "Liq"]

    def test_build_VL(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": True})

        m.fs.params = CubicParameterBlock(default={
                "valid_phase": ("Vap", "Liq")})

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ('Vap', 'Liq')

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 2
        for p in m.fs.params.phase_list:
            assert p in ["Vap", "Liq"]

    def test_build_LV(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": True})

        m.fs.params = CubicParameterBlock(default={
                "valid_phase": ("Liq", "Vap")})

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ('Liq', 'Vap')

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 2
        for p in m.fs.params.phase_list:
            assert p in ["Vap", "Liq"]

    def test_build_L(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": True})

        m.fs.params = CubicParameterBlock(default={
                "valid_phase": ("Liq")})

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ('Liq')

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 1
        for p in m.fs.params.phase_list:
            assert p in ["Liq"]

    def test_build_V(self):
        m = ConcreteModel()

        m.fs = FlowsheetBlock(default={"dynamic": True})

        m.fs.params = CubicParameterBlock(default={
                "valid_phase": ("Vap")})

        assert m.fs.params.state_block_class is CubicStateBlock
        assert m.fs.params.config.valid_phase == ('Vap')

        assert isinstance(m.fs.params.phase_list, Set)
        assert len(m.fs.params.phase_list) == 1
        for p in m.fs.params.phase_list:
            assert p in ["Vap"]

#
#class TestStateBlock(object):
#    def test_config(self):
#        m = ConcreteModel()
#
#        m.fs = FlowsheetBlock(default={"dynamic": False})
#
#        m.fs.params = CubicParameterBlock()
#        m.fs.params.component_list = Set(initialize=["a", "b"])
#        m.fs.params.cubic_type = CubicEoS.PR
#
#        m.fs.props = m.fs.params.state_block_class(
#                [1], default={"parameters": m.fs.params})
