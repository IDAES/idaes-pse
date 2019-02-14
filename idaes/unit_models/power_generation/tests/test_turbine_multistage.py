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

@pytest.fixture()
def build_turbine_for_run_test():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95_ph.Iapws95ParameterBlock()
    # roughly based on NETL baseline studies
    m.fs.turb = TurbineMultistage(default={
        "property_package": m.fs.properties,
        "num_hp": 7,
        "num_ip": 14,
        "num_lp": 11,
        "hp_split_locations": [4,7],
        "ip_split_locations": [5, 14],
        "lp_split_locations": [4,7,9,11],
        "hp_disconnect": [7],
        "ip_split_num_outlets": {14:3}})
    return m

def test_initialize(build_turbine_for_run_test):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine_for_run_test
    turb = m.fs.turb

    hin = iapws95_ph.htpx(T=880, P=2.4233e7)
    m.fs.turb.inlet_split.inlet.enth_mol[0].fix(hin)
    m.fs.turb.inlet_split.inlet.flow_mol[0].fix(26000)
    m.fs.turb.inlet_split.inlet.pressure[0].fix(2.4233e7)

    turb.hp_split[4].split_fraction[0,"outlet_2"].fix(0.05)
    turb.hp_split[7].split_fraction[0,"outlet_2"].fix(0.05)
    turb.ip_split[5].split_fraction[0,"outlet_2"].fix(0.05)
    turb.ip_split[14].split_fraction[0,"outlet_2"].fix(0.05)
    turb.ip_split[14].split_fraction[0,"outlet_3"].fix(0.10)
    turb.lp_split[4].split_fraction[0,"outlet_2"].fix(0.05)
    turb.lp_split[7].split_fraction[0,"outlet_2"].fix(0.05)
    turb.lp_split[9].split_fraction[0,"outlet_2"].fix(0.05)
    turb.lp_split[11].split_fraction[0,"outlet_2"].fix(0.05)

    turb.inlet_split.split_fraction[0,"outlet_1"].value = 0.25
    turb.inlet_split.split_fraction[0,"outlet_2"].value = 0.25
    turb.inlet_split.split_fraction[0,"outlet_3"].value = 0.25
    turb.inlet_split.split_fraction[0,"outlet_4"].value = 0.25
    #turb.display()

    turb.initialize()
    assert(degrees_of_freedom(m)==3)
