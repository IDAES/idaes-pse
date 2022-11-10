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
Tests for config utility methods.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel, Set
from pyomo.dae import ContinuousSet
from pyomo.network import Port
from idaes.core import (
    declare_process_block_class,
    PhysicalParameterBlock,
    StateBlock,
    StateBlockData,
    ReactionParameterBlock,
    useDefault,
)
from idaes.core.base.phases import PhaseType as PT
from idaes.core.util.config import (
    is_physical_parameter_block,
    is_reaction_parameter_block,
    is_state_block,
    is_port,
    is_time_domain,
    is_transformation_method,
    is_transformation_scheme,
    DefaultBool,
)
from idaes.core.util.exceptions import ConfigurationError


@declare_process_block_class("ParameterBlock")
class _ParameterBlock(PhysicalParameterBlock):
    def build(self):
        pass


@pytest.mark.unit
def test_is_physical_parameter_block_passes():
    # Make an instance of a Parameter Block
    p = ParameterBlock()

    # Check that is_physical_parameter_block returns the ParameterBlock
    assert p == is_physical_parameter_block(p)


@pytest.mark.unit
def test_is_physical_parameter_block_useDefault():
    assert useDefault == is_physical_parameter_block(useDefault)


@pytest.mark.unit
def test_is_physical_parameter_block_fails():
    # Test that is_physical_parameter_block returns ConfigurationError with
    # wrong input
    m = ConcreteModel()

    with pytest.raises(ConfigurationError):
        is_physical_parameter_block(m)  # Non Parameter Block Pyomo object
    with pytest.raises(ConfigurationError):
        is_physical_parameter_block("foo")  # str
    with pytest.raises(ConfigurationError):
        is_physical_parameter_block(1)  # int


@declare_process_block_class("RParameterBlock")
class _RParameterBlock(ReactionParameterBlock):
    def build(self):
        pass


@pytest.mark.unit
def test_is_reaction_parameter_block_passes():
    # Make an instance of a Parameter Block
    r = RParameterBlock()

    # Check that is_reaction_parameter_block returns the ReactionParameterBlock
    assert r == is_reaction_parameter_block(r)


@pytest.mark.unit
def test_is_reaction_parameter_block_useDefault():
    # No useDefault option for is_reaction_parameter_block
    with pytest.raises(ConfigurationError):
        is_reaction_parameter_block(useDefault)


@pytest.mark.unit
def test_is_reaction_parameter_block_fails():
    # Test that is_reaction_parameter_block returns ConfigurationError with
    # wrong input
    m = ConcreteModel()

    with pytest.raises(ConfigurationError):
        is_reaction_parameter_block(m)  # Non Parameter Block Pyomo object
    with pytest.raises(ConfigurationError):
        is_reaction_parameter_block("foo")  # str
    with pytest.raises(ConfigurationError):
        is_reaction_parameter_block(1)  # int


@declare_process_block_class("TestStateBlock", block_class=StateBlock)
class StateTestBlockData(StateBlockData):
    def build(self):
        pass


@pytest.mark.unit
def test_is_state_block_passes():
    # Make an instance of a TestStateBlock
    s = TestStateBlock()

    # Check that is_state_block returns the TestStateBlock
    assert s == is_state_block(s)


@pytest.mark.unit
def test_is_state_block_fails():
    # Test that is_state_block returns ConfigurationError with wrong input
    m = ConcreteModel()

    with pytest.raises(ConfigurationError):
        is_state_block(m)  # Non Parameter Block Pyomo object
    with pytest.raises(ConfigurationError):
        is_state_block("foo")  # str
    with pytest.raises(ConfigurationError):
        is_state_block(1)  # int


@pytest.mark.unit
def test_is_port():
    # Test that is_port passes a valid port
    m = ConcreteModel()
    m.c = Port()
    assert isinstance(is_port(m.c), Port)


@pytest.mark.unit
def test_is_port_errors():
    # Test that is_port returns errors when not given a Port
    with pytest.raises(ConfigurationError):
        is_port("foo")  # str
    with pytest.raises(ConfigurationError):
        is_port(["foo", "bar"])  # list of strs
    with pytest.raises(ConfigurationError):
        is_port({"foo": "bar"})  # dict
    with pytest.raises(ConfigurationError):
        is_port(1.0)  # float
    with pytest.raises(ConfigurationError):
        is_port(1)  # int


@pytest.mark.unit
def test_is_time_domain():
    # Test that is_time_domain accepts Sets and ContinuousSets
    m = ConcreteModel()

    m.s = Set(initialize=[1, 2, 3, 4])
    m.cs = ContinuousSet(bounds=[0, 1])

    assert isinstance(is_time_domain(m.s), Set)
    assert isinstance(is_time_domain(m.cs), ContinuousSet)


@pytest.mark.unit
def test_is_time_domain_errors():
    # Test that is_time_domain returns errors when not Set or ContinuousSet

    with pytest.raises(ConfigurationError):
        assert is_time_domain("foo")
    with pytest.raises(ConfigurationError):
        assert is_time_domain(["foo", "bar"])
    with pytest.raises(ConfigurationError):
        assert is_time_domain(("foo", "bar"))
    with pytest.raises(ConfigurationError):
        assert is_time_domain({"foo": "bar"})
    with pytest.raises(ConfigurationError):
        assert is_time_domain(1)
    with pytest.raises(ConfigurationError):
        assert is_time_domain(1.0)


@pytest.mark.unit
def test_is_transformation_method():
    assert is_transformation_method("dae.finite_difference") == "dae.finite_difference"

    assert is_transformation_method("dae.collocation") == "dae.collocation"

    with pytest.raises(ConfigurationError):
        is_transformation_method("dea.finite_difference")


@pytest.mark.unit
def test_is_transformation_scheme():
    assert is_transformation_scheme("BACKWARD") == "BACKWARD"
    assert is_transformation_scheme("FORWARD") == "FORWARD"
    assert is_transformation_scheme("LAGRANGE-RADAU") == "LAGRANGE-RADAU"
    assert is_transformation_scheme("LAGRANGE-LEGENDRE") == "LAGRANGE-LEGENDRE"

    with pytest.raises(ConfigurationError):
        is_transformation_scheme("foo")


@pytest.mark.unit
def test_DefaultBool():
    assert DefaultBool(useDefault) is useDefault
    assert DefaultBool(True)
    assert not DefaultBool(False)
    with pytest.raises(ValueError):
        DefaultBool("foo")
