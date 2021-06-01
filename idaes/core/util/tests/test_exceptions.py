###############################################################################
# ** Copyright Notice **
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
###############################################################################
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


@pytest.mark.unit
def test_ConfigurationError():
    with pytest.raises(ValueError):
        raise ConfigurationError()


@pytest.mark.unit
def test_DynamicError():
    with pytest.raises(ValueError):
        raise DynamicError()


@pytest.mark.unit
def test_BurntToast():
    with pytest.raises(Exception):
        raise BurntToast()


@pytest.mark.unit
def test_PropertyNotSupportedError():
    with pytest.raises(AttributeError):
        raise PropertyNotSupportedError()


@pytest.mark.unit
def test_PropertyPackageError():
    with pytest.raises(AttributeError):
        # This MUST be an AttributeError due to behaviour in Pyomo
        raise PropertyPackageError()
