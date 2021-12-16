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
Tests for common methods used by generic framework

Author: A Lee
"""
import pytest

from types import MethodType

from idaes.generic_models.properties.core.generic.utility import (
    GenericPropertyPackageError, get_method, get_phase_method,
    get_component_object, get_bounds_from_config)

from pyomo.environ import Block, units as pyunits
from pyomo.common.config import ConfigBlock, ConfigValue
from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError


@pytest.fixture
def frame():
    m = Block(concrete=True)
    m.params = Block()
    m.params.config = ConfigBlock()
    m.params.config.declare("test_arg", ConfigValue())
    m.params.config.declare("state_bounds", ConfigValue())

    def get_component(self, comp):
        return getattr(self, comp)
    m.params.get_component = MethodType(get_component, m.params)
    m.params.get_phase = MethodType(get_component, m.params)

    m.params.comp = Block()
    m.params.comp.config = ConfigBlock()
    m.params.comp.config.declare("test_arg_2", ConfigValue())

    return m


@pytest.mark.unit
def test_generic_property_package_error():
    with pytest.raises(
            PropertyPackageError,
            match="Generic Property Package instance block called for "
            "prop, but was not provided with a method "
            "for this property. Please add a method for this property "
            "in the property parameter configuration."):
        raise GenericPropertyPackageError("block", "prop")


@pytest.mark.unit
def test_get_component_object(frame):
    assert get_component_object(frame, "comp") is frame.params.comp


class TestGetMethod():
    @pytest.mark.unit
    def test_get_method_invalid_name(self, frame):
        with pytest.raises(
                AttributeError,
                match="ScalarBlock Generic Property Package called for "
                "invalid configuration option foo. Please contact the "
                "developer of the property package."):
            get_method(frame, "foo")

    @pytest.mark.unit
    def test_get_method_none(self, frame):
        with pytest.raises(
                GenericPropertyPackageError,
                match="Generic Property Package instance ScalarBlock "
                "called for test_arg, but was not provided with a "
                "method for this property. Please add a method for "
                "this property in the property parameter "
                "configuration."):
            get_method(frame, "test_arg")

    @pytest.mark.unit
    def test_get_method_not_callable(self, frame):
        frame.params.config.test_arg = "foo"
        with pytest.raises(
                ConfigurationError,
                match="ScalarBlock Generic Property Package received "
                "invalid value for argument test_arg. Value must be a "
                "method, a class with a method named expression or a "
                "module containing one of the previous."):
            get_method(frame, "test_arg")

    @pytest.mark.unit
    def test_get_method_simple(self, frame):
        def test_arg():
            return "bar"

        frame.params.config.test_arg = test_arg

        mthd = get_method(frame, "test_arg")
        assert mthd() == "bar"

    @pytest.mark.unit
    def test_get_method_class_w_method(self, frame):
        class TestClass():
            def test_arg():
                return "bar"

        frame.params.config.test_arg = TestClass

        mthd = get_method(frame, "test_arg")
        assert mthd() == "bar"

    @pytest.mark.unit
    def test_get_method_class_w_return_expression(self, frame):
        class TestClass():
            def return_expression(*args, **kwargs):
                return "bar"

        frame.params.config.test_arg = TestClass

        mthd = get_method(frame, "test_arg")
        assert mthd() == "bar"

    @pytest.mark.unit
    def test_get_method_phase(self, frame):
        def test_arg():
            return "bar"

        frame.params.config.test_arg = {"test_phase": test_arg}

        mthd = get_method(frame, "test_arg", phase="test_phase")
        assert mthd() == "bar"

    @pytest.mark.unit
    def test_get_method_comp(self, frame):
        def test_arg():
            return "bar"

        frame.params.comp.config.test_arg_2 = test_arg

        mthd = get_method(frame, "test_arg_2", comp="comp")
        assert mthd() == "bar"

    @pytest.mark.unit
    def test_get_method_phase_comp(self, frame):
        def test_arg():
            return "bar"

        frame.params.comp.config.test_arg_2 = {"test_phase": test_arg}

        mthd = get_method(frame, "test_arg_2", comp="comp", phase="test_phase")
        assert mthd() == "bar"


class TestGetPhaseMethod():
    @pytest.mark.unit
    def test_get_phase_method_invalid_name(self, frame):
        with pytest.raises(
                AttributeError,
                match="ScalarBlock Generic Property Package called for "
                "invalid configuration option foo. Please contact the "
                "developer of the property package."):
            get_phase_method(frame, "foo", "comp")

    @pytest.mark.unit
    def test_get_phase_method_none(self, frame):
        with pytest.raises(
                GenericPropertyPackageError,
                match="Generic Property Package instance ScalarBlock "
                "called for test_arg_2, but was not provided with a "
                "method for this property. Please add a method for "
                "this property in the property parameter "
                "configuration."):
            get_phase_method(frame, "test_arg_2", "comp")

    @pytest.mark.unit
    def test_get_phase_method_not_callable(self, frame):
        frame.params.comp.config.test_arg_2 = "foo"
        with pytest.raises(
                ConfigurationError,
                match="ScalarBlock Generic Property Package received "
                "invalid value for argument test_arg_2. Value must be a "
                "method, a class with a method named expression or a "
                "module containing one of the previous."):
            get_phase_method(frame, "test_arg_2", "comp")

    @pytest.mark.unit
    def test_get_phase_method_simple(self, frame):
        def test_arg():
            return "bar"

        frame.params.comp.config.test_arg_2 = test_arg

        mthd = get_phase_method(frame, "test_arg_2", "comp")
        assert mthd() == "bar"

    @pytest.mark.unit
    def test_get_phase_method_class_w_method(self, frame):
        class TestClass():
            def test_arg_2():
                return "bar"

        frame.params.comp.config.test_arg_2 = TestClass

        mthd = get_phase_method(frame, "test_arg_2", "comp")
        assert mthd() == "bar"

    @pytest.mark.unit
    def test_get_phase_method_class_w_return_expression(self, frame):
        class TestClass():
            def return_expression():
                return "bar"

        frame.params.comp.config.test_arg_2 = TestClass

        mthd = get_phase_method(frame, "test_arg_2", "comp")
        assert mthd() == "bar"


class TestGetBoundsFromConfig():
    @pytest.mark.unit
    def test_no_state_bounds(self, frame):
        bounds, value = get_bounds_from_config(frame, "foo", "bar")

        assert bounds == (None, None)
        assert value is None

    @pytest.mark.unit
    def test_no_state_bounds_3(self, frame):
        frame.params.config.state_bounds = {
            "test_state": (1, 2, 3)}
        bounds, value = get_bounds_from_config(frame, "test_state", "bar")

        assert bounds == (1, 3)
        assert value == 2

    @pytest.mark.unit
    def test_no_state_bounds_4(self, frame):
        frame.params.config.state_bounds = {
            "test_state": (1, 2, 3, pyunits.km)}
        bounds, value = get_bounds_from_config(frame, "test_state", pyunits.m)

        assert bounds == (1000, 3000)
        assert value == 2000
