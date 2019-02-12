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
        "lp_split_num_outlets": {3:3}})
    return m

def test_basic_build(build_turbine_for_buid_test):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine_for_buid_test
    turb = m.fs.turb
    #assert(isinstance(turb.inlet_stage, TurbineInletStage))
    """
    assert(0 not in turb.inlet_stage)
    assert(1 in turb.inlet_stage)
    assert(2 in turb.inlet_stage)
    assert(3 in turb.inlet_stage)
    assert(4 in turb.inlet_stage)
    assert(5 not in turb.inlet_stage)
    assert(isinstance(turb.outlet_stage, TurbineOutletStage))
    assert(isinstance(turb.hp_stages, TurbineStage))
    assert(isinstance(turb.ip_stages, TurbineStage))
    assert(isinstance(turb.lp_stages, TurbineStage))
    assert(isinstance(turb.hp_split, Separator))
    assert(isinstance(turb.ip_split, Separator))
    assert(isinstance(turb.lp_split, Separator))
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
    assert((1,1) in turb.hp_stream)
    assert((1,2) not in turb.hp_stream)
    assert((2,1) not in turb.hp_stream)
    assert((2,2) not in turb.hp_stream)
    assert((1,1) in turb.ip_stream)
    assert((1,2) in turb.ip_stream)
    assert((2,1) in turb.ip_stream)
    assert((2,2) in turb.ip_stream)
    for i in [3,4,5,6,7,8,9]:
        assert((i,1) in turb.ip_stream)
        assert((i,2) not in turb.ip_stream)
    assert((10,1) not in turb.ip_stream)
    assert((10,2) not in turb.ip_stream)
    for i in [1,2]:
        assert((i,1) in turb.lp_stream)
        assert((i,2) not in turb.lp_stream)
    for i in [3,4]:
        assert((i,1) in turb.lp_stream)
        assert((i,2) in turb.lp_stream)
    assert((5,1) not in turb.lp_stream)
    assert((5,2) not in turb.lp_stream)

    # though to test everything here so to make sure connections are right
    # will fix split fractions except one for each splitter.  in this case
    # with no disconnections that should leave the number of variables in the
    # inlet port as the degrees of freedom
    """
    '''
    turb.hp_split[0].split_fraction[0,"outlet_2"].fix(0.05)
    turb.ip_split[1].split_fraction[0,"outlet_2"].fix(0.05)
    turb.ip_split[2].split_fraction[0,"outlet_2"].fix(0.05)
    turb.lp_split[3].split_fraction[0,"outlet_2"].fix(0.05)
    turb.lp_split[3].split_fraction[0,"outlet_3"].fix(0.05)
    turb.lp_split[4].split_fraction[0,"outlet_2"].fix(0.05)
    '''
    turb.inlet_split.split_fraction[0,"outlet_2"].fix(0.25)
    turb.inlet_split.split_fraction[0,"outlet_3"].fix(0.25)
    turb.inlet_split.split_fraction[0,"outlet_4"].fix(0.25)
    turb.display()
    assert(degrees_of_freedom(m)==3)
