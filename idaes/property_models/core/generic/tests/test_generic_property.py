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

from pyomo.environ import Block, ConcreteModel
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.property_models.core.generic.generic_property import (
        GenericPropertyPackageError, get_method)

from idaes.core.util.exceptions import ConfigurationError, PropertyPackageError
from idaes.core.util.misc import add_object_reference


def test_GenericPropertyPackageError():
    with pytest.raises(PropertyPackageError):
        raise GenericPropertyPackageError("block", "prop")


# Dummy method for testing get_method
def dummy_option(b):
    return b.dummy_response


class TestGetMethod(object):
    @pytest.fixture(scope="class")
    def frame(self):
        m = ConcreteModel()

        # Create a dummy parameter block
        m.params = Block()

        # Add necessary parameters to parameter block
        m.params.config = ConfigBlock()
        m.params.config.declare("dummy_option", ConfigValue(default=None))

        # Create a dummy state block
        m.props = Block([1])
        m.props[1].config = ConfigBlock()
        add_object_reference(m.props[1], "_params", m.params)

        m.props[1].dummy_response = "foo"

        return m

    def test_None(self, frame):
        with pytest.raises(GenericPropertyPackageError):
            get_method(frame.props[1], "dummy_option")

    def test_invlaid_option(self, frame):
        with pytest.raises(AttributeError):
            get_method(frame.props[1], "foo")

    def test_method(self, frame):
        # Test that get_method works when provided with the method directly
        frame.params.config.dummy_option = dummy_option

        assert get_method(frame.props[1], "dummy_option") is dummy_option
        assert get_method(frame.props[1], "dummy_option")(frame.props[1]) is \
            frame.props[1].dummy_response

    def test_module(self, frame):
        # Test that get_method works when pointed to a module with the method
        frame.params.config.dummy_option = modules[__name__]

        assert get_method(frame.props[1], "dummy_option") is dummy_option
        assert get_method(frame.props[1], "dummy_option")(frame.props[1]) is \
            frame.props[1].dummy_response

    def test_not_method_or_module(self, frame):
        # Test that get_method works when provided with the method directly
        frame.params.config.dummy_option = "foo"

        with pytest.raises(ConfigurationError):
            get_method(frame.props[1], "dummy_option")
