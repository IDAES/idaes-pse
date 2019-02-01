
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
This module contains miscalaneous utility functions for use in IDAES models.
"""

__author__ = "Andrew Lee"

import pytest

from pyomo.environ import ConcreteModel, Set

from idaes.core.util.misc import add_object_reference


def test_add_object_reference():
    m = ConcreteModel()

    m.s = Set(initialize=[1, 2, 3])

    add_object_reference(m, "test_ref", m.s)

    assert hasattr(m, "test_ref")
    assert m.test_ref == m.s


def test_add_object_reference_fail():
    m = ConcreteModel()

    with pytest.raises(AttributeError):
        add_object_reference(m, "test_ref", m.s)
