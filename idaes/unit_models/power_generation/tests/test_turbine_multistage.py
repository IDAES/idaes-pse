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
from idaes.unit_models import Separator, Mixer
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
def build_turbine_for_buid_test():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95_ph.Iapws95ParameterBlock()
    # use the default number of stages
    m.fs.turb = TurbineMultistage(default={
        "property_package": m.fs.properties,
        "hp_split_locations": [0],
        "ip_split_locations": [1,2],
        "lp_split_locations": [3,4],
        "hp_mix_locations": [0],
        "ip_mix_locations": [1,2],
        "lp_mix_locations": [3,4]})
    return m

def test_basic_build(build_turbine_for_buid_test):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine_for_buid_test
    turb = m.fs.turb
    assert(isinstance(turb.inlet_stage, TurbineInletStage))
    assert(isinstance(turb.outlet_stage, TurbineOutletStage))
    assert(isinstance(turb.hp_stages, TurbineStage))
    assert(isinstance(turb.ip_stages, TurbineStage))
    assert(isinstance(turb.lp_stages, TurbineStage))
    assert(isinstance(turb.hp_split, Separator))
    assert(isinstance(turb.ip_split, Separator))
    assert(isinstance(turb.lp_split, Separator))
    assert(isinstance(turb.hp_mix, Mixer))
    assert(isinstance(turb.ip_mix, Mixer))
    assert(isinstance(turb.lp_mix, Mixer))
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
    assert(0 in turb.hp_split)
    assert(1 in turb.ip_split)
    assert(2 in turb.ip_split)
    assert(3 in turb.lp_split)
    assert(4 in turb.lp_split)
    assert(0 in turb.hp_mix)
    assert(1 in turb.ip_mix)
    assert(2 in turb.ip_mix)
    assert(3 in turb.lp_mix)
    assert(4 in turb.lp_mix)
