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
Tests for turbine inlet model.

Author: John Eslick
"""
import pytest

from pyomo.environ import (
    ConcreteModel, SolverFactory, units as pyunits, value,
    TransformationFactory)
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock
from idaes.power_generation.unit_models import TurbineInletStage
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


@pytest.fixture(scope="module")
def build_turbine():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = TurbineInletStage(
        default={"property_package": m.fs.properties})
    return m


@pytest.fixture()
def build_turbine_dyn():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": True})
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = TurbineInletStage(default={
        "dynamic": False,
        "property_package": m.fs.properties})
    return m


@pytest.mark.integration
def test_basic_build(build_turbine):
    """Make a turbine model and make sure it doesn't throw exception"""
    assert_units_consistent(build_turbine)


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize(build_turbine):
    """Initialize a turbine model"""
    hin = value(iapws95.htpx(T=880*pyunits.K, P=2.4233e7*pyunits.Pa))
    # set inlet
    build_turbine.fs.turb.inlet.enth_mol[0].value = hin
    build_turbine.fs.turb.inlet.flow_mol[0].value = 26000/4.0
    build_turbine.fs.turb.inlet.pressure[0].value = 2.4233e7

    build_turbine.fs.turb.initialize(outlvl=1)

    eq_cons = activated_equalities_generator(build_turbine)
    for c in eq_cons:
        assert(abs(c.body() - c.lower) < 1e-4)

    assert degrees_of_freedom(build_turbine) == 3


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_dyn(build_turbine_dyn):
    """Initialize a turbine model"""
    m = build_turbine_dyn
    hin = iapws95.htpx(T=880*pyunits.K, P=2.4233e7*pyunits.Pa)
    discretizer = TransformationFactory('dae.finite_difference')
    discretizer.apply_to(m, nfe=4, wrt=m.fs.time, scheme='BACKWARD')

    # fix inlet
    m.fs.turb.inlet[:].enth_mol.fix(hin)
    m.fs.turb.inlet[:].flow_mol.fix(26000/4.0)
    m.fs.turb.inlet[:].pressure.fix(2.4233e7)
    m.fs.turb.flow_coeff[:].fix()

    assert(degrees_of_freedom(m) == 0)
    m.fs.turb.initialize()

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert(abs(c.body() - c.lower) < 1e-4)

    m.display()

    assert(degrees_of_freedom(m) == 0)


@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.unit
def test_report(build_turbine):
    build_turbine.fs.turb.report()
