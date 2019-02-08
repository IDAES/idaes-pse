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
Tests for turbine multistage model.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory

from idaes.core import FlowsheetBlock
from idaes.unit_models.power_generation import (
    TurbineMultistage, TurbineStage, TurbineInletStage, TurbineOutletStage)
from idaes.property_models import iapws95_ph
from idaes.ui.report import degrees_of_freedom

# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None

@pytest.fixture()
def build_turbine():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95_ph.Iapws95ParameterBlock()
    # use the default number of stages
    m.fs.turb = TurbineMultistage(default={"property_package": m.fs.properties})
    return m

def test_basic_build(build_turbine):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine
    turb = m.fs.turb
    assert(isinstance(turb.inlet_stage, TurbineInletStage))
    assert(isinstance(turb.outlet_stage, TurbineOutletStage))
    assert(isinstance(turb.hp_stages, TurbineStage))
    assert(isinstance(turb.ip_stages, TurbineStage))
    assert(isinstance(turb.lp_stages, TurbineStage))
    assert(0 not in turb.hp_stages)
    assert(1 in turb.hp_stages)
    assert(2 in turb.hp_stages)
    assert(3 not in turb.hp_stages)
    assert(0 not in turb.ip_stages)
    assert(1 in turb.ip_stages)
    assert(10 in turb.ip_stages)
    assert(11 not in turb.ip_stages)
    assert(0 not in turb.lp_stages)
    assert(1 in turb.lp_stages)
    assert(5 in turb.lp_stages)
    assert(6 not in turb.lp_stages)
