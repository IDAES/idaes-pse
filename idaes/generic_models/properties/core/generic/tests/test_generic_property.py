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
Tests for generic property package core code

Author: Andrew Lee
"""
import pytest
from sys import modules

from pyomo.environ import Block, ConcreteModel, Expression, Param, Set, Var
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.generic_models.properties.core.generic.generic_property import (
        GenericPropertyPackageError,
        # get_method,
        GenericParameterData)
from idaes.generic_models.properties.core.generic.tests import dummy_eos

from idaes.core import declare_process_block_class, Component, Phase, LiquidPhase
from idaes.core.util.exceptions import (ConfigurationError,
                                        PropertyPackageError)
from idaes.core.util.misc import add_object_reference

# -----------------------------------------------------------------------------
@declare_process_block_class("DummyParameterBlock")
class DummyParameterData(GenericParameterData):
    def configure(self):
        self.configured = True


class TestGenericParameterBlock(object):
    def test_build(self):
        m = ConcreteModel()
        m.params = DummyParameterBlock(default={
                "components": {"a": {}, "b": {}, "c": {}},
                "phases": {
                    1: {"type": LiquidPhase,
                        "component_list": ["a", "b"],
                        "equation_of_state": "foo"},
                    2: {"equation_of_state": "bar"}},
                "state_definition": "foo"})

        assert m.params.configured

        assert isinstance(m.params.component_list, Set)
        assert len(m.params.component_list) == 3
        for j in m.params.component_list:
            assert j in ["a", "b", "c"]
            assert isinstance(m.params.get_component(j), Component)

        assert isinstance(m.params.phase_list, Set)
        assert len(m.params.phase_list) == 2
        for p in m.params.phase_list:
            assert p in ["1", "2"]
        assert isinstance(m.params.get_phase("1"), LiquidPhase)
        assert isinstance(m.params.get_phase("2"), Phase)

        assert isinstance(m.params._phase_component_set, Set)
        assert len(m.params._phase_component_set) == 5
        for i in m.params._phase_component_set:
            assert i in [("1", "a"), ("1", "b"),
                         ("2", "a"), ("2", "b"), ("2", "c")]
