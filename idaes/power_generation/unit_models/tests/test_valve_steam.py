##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Tests steam/water valves

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory

from idaes.core import FlowsheetBlock
from idaes.power_generation.unit_models import SteamValve
from idaes.generic_models.properties import iapws95
from idaes.core.util.model_statistics import (
        degrees_of_freedom,
        activated_equalities_generator)

prop_available = iapws95.iapws95_available()

# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
else:
    solver = None


@pytest.fixture()
def build_valve_vapor():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.valve = SteamValve(default={"property_package": m.fs.properties})
    return m


@pytest.fixture()
def build_valve_liquid():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.valve = SteamValve(default={"property_package": m.fs.properties,
                                     "phase": "Liq"})
    return m


@pytest.mark.unit
def test_vapor_steady_state_build(build_valve_vapor):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_valve_vapor


@pytest.mark.unit
def test_liquid_steady_state_build(build_valve_liquid):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_valve_liquid


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_vapor_steady_state_initialize(build_valve_vapor):
    """Initialize a turbine model"""
    m = build_valve_vapor
    # set inlet
    hin = iapws95.htpx(T=880, P=2.4233e7)
    # set inlet
    m.fs.valve.inlet.enth_mol[0].value = hin
    m.fs.valve.inlet.flow_mol[0].value = 26000/4.0
    m.fs.valve.inlet.pressure[0].value = 2.5e7
    m.fs.valve.Cv.fix(0.01)

    m.fs.valve.initialize(outlvl=1)

    eq_cons = activated_equalities_generator(m)

    for c in eq_cons:
        assert(abs(c.body() - c.lower) < 1e-4)
    assert(degrees_of_freedom(m)==3) #inlet was't fixed and still shouldn't be


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.unit
def test_report(build_valve_vapor):
    """Initialize a turbine model"""
    m = build_valve_vapor
    m.fs.valve.report()
