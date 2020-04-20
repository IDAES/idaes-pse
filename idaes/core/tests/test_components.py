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
Tests for Component objects

Author: Andrew Lee
"""
from pyomo.environ import ConcreteModel, Set

from idaes.core.components import Component


def test_config():
    m = ConcreteModel()

    m.comp = Component()

    assert len(m.comp.config) == 1
    assert not m.comp.config._component_list_exists


def test_populate_component_list():
    m = ConcreteModel()

    m.comp = Component()
    m.comp2 = Component()

    assert isinstance(m.component_list, Set)

    for j in m.component_list:
        assert j in ["comp", "comp2"]
