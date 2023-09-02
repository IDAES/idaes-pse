#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for turbine inlet model.

Author: John Eslick
"""
import pytest

from pyomo.environ import ConcreteModel, TransformationFactory, units as pyunits, value

from idaes.core import FlowsheetBlock
from idaes.models_extra.power_generation.unit_models.helm import HelmTurbineInletStage
from idaes.models.properties import iapws95
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    activated_equalities_generator,
)
from idaes.core.solvers import get_solver
from idaes.models.properties.general_helmholtz import helmholtz_available

# Set up solver
solver = get_solver()


@pytest.fixture()
def build_turbine():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineInletStage(property_package=m.fs.properties)
    return m


@pytest.fixture()
def build_turbine_dyn():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 4], time_units=pyunits.s)
    m.fs.properties = iapws95.Iapws95ParameterBlock()
    m.fs.turb = HelmTurbineInletStage(property_package=m.fs.properties)
    return m


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.unit
def test_basic_build(build_turbine):
    """Make a turbine model and make sure it doesn't throw exception"""
    m = build_turbine


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.component
def test_initialize(build_turbine):
    """Initialize a turbine model"""
    m = build_turbine
    hin = value(iapws95.htpx(T=880 * pyunits.K, P=2.4233e7 * pyunits.Pa))
    # set inlet
    m.fs.turb.inlet.enth_mol[0].value = hin
    m.fs.turb.inlet.flow_mol[0].value = 26000 / 4.0
    m.fs.turb.inlet.pressure[0].value = 2.4233e7
    m.fs.turb.initialize(outlvl=1)
    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4
    assert degrees_of_freedom(m) == 3  # inlet was't fixed and still shouldn't be


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.component
def test_initialize_calc_cf(build_turbine):
    """Initialize a turbine model"""
    m = build_turbine
    hin = value(iapws95.htpx(T=880 * pyunits.K, P=2.4233e7 * pyunits.Pa))
    # set inlet
    m.fs.turb.inlet.enth_mol[0].value = hin
    m.fs.turb.inlet.flow_mol[0].value = 26000 / 4.0
    m.fs.turb.inlet.pressure[0].value = 2.4233e7
    m.fs.turb.ratioP[0].value = 0.6

    m.fs.turb.initialize(outlvl=1, calculate_cf=True)
    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4

    m.fs.turb.inlet.enth_mol[0].value = hin
    m.fs.turb.inlet.flow_mol[0].value = 26000 / 4.0
    m.fs.turb.inlet.pressure[0].value = 2.4233e7
    m.fs.turb.inlet.fix()
    solver.solve(m)
    assert m.fs.turb.ratioP[0].value == pytest.approx(0.6)

    assert degrees_of_freedom(m) == 0


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.component
def test_initialize_dyn(build_turbine_dyn):
    """Initialize a turbine model"""
    m = build_turbine_dyn
    hin = value(iapws95.htpx(T=880 * pyunits.K, P=2.4233e7 * pyunits.Pa))
    discretizer = TransformationFactory("dae.finite_difference")
    discretizer.apply_to(m, nfe=4, wrt=m.fs.time, scheme="BACKWARD")

    # fix inlet
    m.fs.turb.inlet.enth_mol.fix(hin)
    m.fs.turb.inlet.flow_mol.fix(26000 / 4.0)
    m.fs.turb.inlet.pressure.fix(2.4233e7)
    m.fs.turb.flow_coeff[:].fix()

    m.fs.turb.initialize()

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4

    assert degrees_of_freedom(m) == 0


@pytest.mark.skipif(not helmholtz_available(), reason="General Helmholtz not available")
@pytest.mark.component
def test_initialize_dyn2(build_turbine_dyn):
    """Initialize a turbine model"""
    m = build_turbine_dyn
    hin = iapws95.htpx(T=880 * pyunits.K, P=2.4233e7 * pyunits.Pa)
    discretizer = TransformationFactory("dae.finite_difference")
    discretizer.apply_to(m, nfe=4, wrt=m.fs.time, scheme="BACKWARD")

    # fix inlet
    m.fs.turb.inlet.enth_mol.fix(hin)
    m.fs.turb.inlet.flow_mol.fix(26000 / 4.0)
    m.fs.turb.inlet.pressure.fix(2.4233e7)
    m.fs.turb.flow_coeff[:].fix()
    m.fs.turb.ratioP[0].value = 0.6  # used for flow coeff calc

    m.fs.turb.initialize(calculate_cf=True)

    eq_cons = activated_equalities_generator(m)
    for c in eq_cons:
        assert abs(c.body() - c.lower) < 1e-4
    assert m.fs.turb.ratioP[0].value == pytest.approx(0.6)
    assert degrees_of_freedom(m) == 0


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.unit
def test_get_stream_table_contents(build_turbine):
    stable = build_turbine.fs.turb._get_stream_table_contents()

    expected = {
        "Units": {
            "Mass Flow": getattr(pyunits.pint_registry, "kg/s"),
            "Molar Flow": getattr(pyunits.pint_registry, "mol/s"),
            "Molar Enthalpy": getattr(pyunits.pint_registry, "J/mol"),
            "P": getattr(pyunits.pint_registry, "Pa"),
            "T": getattr(pyunits.pint_registry, "K"),
            "Vapor Fraction": getattr(pyunits.pint_registry, "dimensionless"),
        },
        "Inlet": {
            "Mass Flow": pytest.approx(0.01801527, rel=1e-5),
            "Molar Flow": pytest.approx(1.0, rel=1e-5),
            "Molar Enthalpy": pytest.approx(0.01102139, rel=1e-5),
            "P": pytest.approx(11032300, rel=1e-5),
            "T": pytest.approx(270.4877, rel=1e-5),
            "Vapor Fraction": pytest.approx(0.0, abs=1e-5),
        },
        "Outlet": {
            "Mass Flow": pytest.approx(0.01801527, rel=1e-5),
            "Molar Flow": pytest.approx(1.0, rel=1e-5),
            "Molar Enthalpy": pytest.approx(0.01102139, rel=1e-5),
            "P": pytest.approx(11032300, rel=1e-5),
            "T": pytest.approx(270.4877, rel=1e-5),
            "Vapor Fraction": pytest.approx(0.0, abs=1e-5),
        },
    }

    assert stable.to_dict() == expected
