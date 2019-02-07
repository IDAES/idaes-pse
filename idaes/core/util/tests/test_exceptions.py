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
Tests for custom exceptions.
"""

import pytest
from idaes.core.util.exceptions import *

__author__ = "Andrew Lee"


def test_BalanceTypeNotSupportedError():
    with pytest.raises(NotImplementedError):
        raise BalanceTypeNotSupportedError()


def test_ConfigurationError():
    with pytest.raises(ValueError):
        raise ConfigurationError()


def test_DynamicError():
    with pytest.raises(ValueError):
        raise DynamicError()


def test_BurntToast():
    with pytest.raises(Exception):
        raise BurntToast()


def test_PropertyNotSupportedError():
    with pytest.raises(AttributeError):
        raise PropertyNotSupportedError()


def test_PropertyPackageError():
    with pytest.raises(AttributeError):
        # This MUST be an AttributeError due to behaviour in Pyomo
        raise PropertyPackageError()
