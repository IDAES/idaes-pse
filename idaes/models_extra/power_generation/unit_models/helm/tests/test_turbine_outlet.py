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
Tests for turbine outlet model.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, TransformationFactory, units as pyunits

from idaes.core import FlowsheetBlock
from idaes.models_extra.power_generation.unit_models.helm import HelmTurbineOutletStage
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
    m.fs.turb = HelmTurbineOutletStage(property_package=m.fs.properties)
    return m


@pytest.fixture()
def build_turbine_dyn():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineOutletStage(dynamic=False, property_package=m.fs.properties)
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
    m.fs.turb.inlet.enth_mol[0].value = 47115
    m.fs.turb.inlet.flow_mol[0].value = 15000
    m.fs.turb.inlet.pressure[0].value = 8e4
    m.fs.turb.outlet.pressure[0].fix(4e4)

    m.fs.turb.initialize(outlvl=1)

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4
    assert degrees_of_freedom(m) == 2  # inlet was't fixed and still shouldn't be


@pytest.mark.component
def test_initialize_calc_cf(build_turbine):
    """Initialize a turbine model"""
    m = build_turbine
    # set inlet
    m.fs.turb.inlet.enth_mol[0].value = 47115
    m.fs.turb.inlet.flow_mol[0].value = 15000
    m.fs.turb.inlet.pressure[0].value = 8e4
    m.fs.turb.outlet.pressure[0].fix(4e4)

    m.fs.turb.initialize(calculate_cf=True)

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4

    m.fs.turb.inlet.enth_mol[0].fix()
    m.fs.turb.inlet.pressure[0].fix()

    solver.solve(m)
    assert m.fs.turb.inlet.flow_mol[0].value == pytest.approx(15000)
    assert degrees_of_freedom(m) == 0


@pytest.mark.component
def test_initialize_calc_cf_dyn(build_turbine_dyn):
    """Initialize a turbine model"""
    m = build_turbine_dyn
    discretizer = TransformationFactory("dae.finite_difference")
    discretizer.apply_to(m, nfe=4, wrt=m.fs.time, scheme="BACKWARD")
    # set inlet
    m.fs.turb.inlet.enth_mol.fix(47115)
    for t in m.fs.turb.inlet.flow_mol:
        m.fs.turb.inlet.flow_mol[t].value = 15000
    m.fs.turb.inlet.pressure.fix(8e4)
    m.fs.turb.outlet.pressure.fix(4e4)
    m.fs.turb.flow_coeff.fix()

    assert degrees_of_freedom(m) == 0
    m.fs.turb.initialize(calculate_cf=True)
    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4
    solver.solve(m)
    assert m.fs.turb.inlet.flow_mol[0].value == pytest.approx(15000)
    assert degrees_of_freedom(m) == 0


@pytest.mark.unit
def test_report(build_turbine):
    m = build_turbine
    m.fs.turb.report()
