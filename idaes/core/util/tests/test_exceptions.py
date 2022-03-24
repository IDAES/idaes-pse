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
Tests for custom exceptions.
"""

import pytest
from idaes.core.util.exceptions import *

__author__ = "Andrew Lee"


@pytest.mark.unit
def test_BalanceTypeNotSupportedError():
    with pytest.raises(NotImplementedError):
        raise BalanceTypeNotSupportedError()
    with pytest.raises(IdaesError):
        raise BalanceTypeNotSupportedError()


@pytest.mark.unit
def test_ConfigurationError():
    with pytest.raises(ValueError):
        raise ConfigurationError()
    with pytest.raises(IdaesError):
        raise ConfigurationError()


@pytest.mark.unit
def test_DynamicError():
    with pytest.raises(ValueError):
        raise DynamicError()
    with pytest.raises(IdaesError):
        raise DynamicError()


@pytest.mark.unit
def test_BurntToast():
    with pytest.raises(IdaesError):
        raise BurntToast()


@pytest.mark.unit
def test_PropertyNotSupportedError():
    with pytest.raises(AttributeError):
        raise PropertyNotSupportedError()
    with pytest.raises(IdaesError):
        raise PropertyNotSupportedError()


@pytest.mark.unit
def test_PropertyPackageError():
    with pytest.raises(AttributeError):
        # This MUST be an AttributeError due to behaviour in Pyomo
        raise PropertyPackageError()
    with pytest.raises(IdaesError):
        raise PropertyPackageError()


@pytest.mark.unit
def test_InitializationError():
    with pytest.raises(ArithmeticError):
        raise InitializationError()
    with pytest.raises(IdaesError):
        raise InitializationError()


@pytest.mark.unit
def test_UserModelError():
    with pytest.raises(ValueError):
        raise UserModelError()
    with pytest.raises(IdaesError):
        raise UserModelError()
