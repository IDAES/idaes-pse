#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for turbine stage model.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel

from idaes.core import FlowsheetBlock
from idaes.models_extra.power_generation.unit_models.helm import HelmTurbineStage
from idaes.models.properties import iapws95
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    activated_equalities_generator,
)
from idaes.core.solvers import get_solver

# Set up solver
solver = get_solver()


@pytest.fixture()
def build_turbine():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineStage(property_package=m.fs.properties)
    return m


@pytest.fixture()
def build_turbine_dyn():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineStage(dynamic=False, property_package=m.fs.properties)
    return m


@pytest.mark.unit
def test_basic_build(build_turbine):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine


@pytest.mark.component
def test_initialize(build_turbine):
    """Initialize a turbine model"""
    m = build_turbine
    # set inlet
    m.fs.turb.inlet.enth_mol[0].value = 70000
    m.fs.turb.inlet.flow_mol[0].value = 15000
    m.fs.turb.inlet.pressure[0].value = 8e6
    m.fs.turb.efficiency_isentropic.fix(0.8)
    m.fs.turb.ratioP.fix(0.7)
    m.fs.turb.initialize(outlvl=4)

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4
    assert degrees_of_freedom(m) == 3  # inlet was't fixed and still shouldn't be


@pytest.mark.unit
def test_report(build_turbine):
    m = build_turbine
    m.fs.turb.report()
