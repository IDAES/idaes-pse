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
Tests for 0D heat exchanger models.

Author: John Eslick
"""

from idaes.unit_models.power_generation import TurbineInletStage
from idaes.property_models import iapws95_ph
from idaes.ui.report import degrees_of_freedom

@pytest.fixture()
def build_heater():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95_ph.Iapws95ParameterBlock()
    m.fs.heater = TurbineInletStage(default={"property_package": m.fs.properties})
    return m
