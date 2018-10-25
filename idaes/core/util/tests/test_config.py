##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
# 
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for config utility methods.

Author: Andrew Lee
"""
import pytest
from pyomo.environ import ConcreteModel
from pyomo.network import Port
from idaes.core.util.config import is_parameter_block, list_of_floats, \
                                    list_of_strings, is_port

# Import a Property file for testing
from os.path import abspath, dirname, join
from pyutilib.misc import import_file
example = join(dirname(abspath(__file__)), '..', '..', '..', '..',
               'examples', 'core', 'examples',
               'steam_properties_demo.py')
fe = import_file(example)


def test_is_parameter_block_passes():
    # Make an instance of a Parameter Block
    p = fe.PropertyParameterBlock()

    # Check that is_parameter_block returns the ParameterBlock
    assert p == is_parameter_block(p)


def test_is_parameter_block_fails():
    # Test that is_parameter_block returns TypeError with wrong input
    m = ConcreteModel()

    with pytest.raises(TypeError):
        is_parameter_block(m)  # Non Parameter Block Pyomo object
    with pytest.raises(TypeError):
        is_parameter_block("foo")  # str
    with pytest.raises(TypeError):
        is_parameter_block(1)  # int


def test_list_of_strings():
    # Test list_of_strings=returns correctly
    assert list_of_strings(1) == ['1']  # int
    assert list_of_strings([1, 2, 3]) == ['1', '2', '3']  # list of ints
    assert list_of_strings(1.0) == ['1.0']  # float
    # list of floats
    assert list_of_strings([1.0, 2.0, 3.0]) == ['1.0', '2.0', '3.0']
    assert list_of_strings("foo") == ["foo"]  # str
    assert list_of_strings(["foo", "bar"]) == ["foo", "bar"]  # list of strs


def test_list_of_strings_errors():
    # Test that list_of_strings fails correctly
    with pytest.raises(ValueError):
        list_of_floats({"foo": "bar"})  # dict


def test_list_of_floats():
    # Test list_of_floats returns correctly
    assert list_of_floats(1) == [1.0]  # int
    assert list_of_floats([1, 2, 3]) == [1.0, 2.0, 3.0]  # list of ints
    assert list_of_floats(1.0) == [1.0]  # float
    assert list_of_floats([1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]  # list of floats


def test_list_of_floats_errors():
    # Test that list_of_floats fails correctly
    with pytest.raises(ValueError):
        list_of_floats("foo")  # str
    with pytest.raises(ValueError):
        list_of_floats(["foo", "bar"])  # list of strs
    with pytest.raises(ValueError):
        list_of_floats({"foo": "bar"})  # dict


def test_is_port():
    # Test that is_port passes a valid port
    m = ConcreteModel()
    m.c = Port()
    assert isinstance(is_port(m.c), Port)


def test_is_port_errors():
    # Test that is_port returns errors when not given a Port
    with pytest.raises(TypeError):
        is_port("foo")  # str
    with pytest.raises(TypeError):
        is_port(["foo", "bar"])  # list of strs
    with pytest.raises(TypeError):
        is_port({"foo": "bar"})  # dict
    with pytest.raises(TypeError):
        is_port(1.0)  # float
    with pytest.raises(TypeError):
        is_port(1)  # int
