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
Testing fixtures and methods for generic property packages
Authors: Andrew Lee
"""

import pytest

from pyomo.environ import ConcreteModel, Block, Var
from idaes.core.util.misc import add_object_reference


@pytest.fixture()
def generic_frame():
    m = ConcreteModel()
    m.params = Block()

    m.params.temperature_ref = Var(initialize=273.15)
    m.params.pressure_ref = Var(initialize=1e5)

    m.props = Block([1])
    add_object_reference(m.props[1], "_params", m.params)

    m.props[1].temperature = Var(initialize=298.15)
    m.props[1].pressure = Var(initialize=101325)

    return m

