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

import pytest
import pyomo.environ as pyo

from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.unit_models import FWH0DDynamic
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.testing import initialization_tester
from idaes.core.solvers import get_solver

prop_available = iapws95.iapws95_available()

# Set up solver
solver = get_solver()


@pytest.fixture(scope="module")
def build_unit():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()
    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(
        dynamic=True,
        time_set=[0, 60],
        time_units=pyo.units.s,
        default_property_package=iapws95.Iapws95ParameterBlock(),
    )
    m.fs.properties = m.fs.config.default_property_package
    m.fs.unit = FWH0DDynamic(
        has_desuperheat=True,
        has_drain_cooling=True,
        has_drain_mixer=True,
        condense={
            "cold_side": {"has_pressure_change": True},
            "hot_side": {"has_pressure_change": True},
            "has_holdup": True,
        },
        desuperheat={"dynamic": False},
        cooling={"dynamic": False, "has_holdup": False},
        property_package=m.fs.properties,
    )
    m.discretizer = pyo.TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=2, wrt=m.fs.time, scheme="BACKWARD")
    m.fs.unit.set_initial_condition()
    m.fs.unit.desuperheat.hot_side_inlet.flow_mol.fix(100)
    m.fs.unit.desuperheat.hot_side_inlet.flow_mol.unfix()
    m.fs.unit.desuperheat.hot_side_inlet.pressure.fix(201325)
    m.fs.unit.desuperheat.hot_side_inlet.enth_mol.fix(60000)
    m.fs.unit.drain_mix.drain.flow_mol.fix(1)
    m.fs.unit.drain_mix.drain.pressure.fix(201325)
    m.fs.unit.drain_mix.drain.enth_mol.fix(20000)
    m.fs.unit.cooling.cold_side_inlet.flow_mol.fix(400)
    m.fs.unit.cooling.cold_side_inlet.pressure.fix(101325)
    m.fs.unit.cooling.cold_side_inlet.enth_mol.fix(3000)
    m.fs.unit.condense.area.fix(600)
    m.fs.unit.condense.overall_heat_transfer_coefficient.fix(1010)
    m.fs.unit.desuperheat.area.fix(85)
    m.fs.unit.desuperheat.overall_heat_transfer_coefficient.fix(145)
    m.fs.unit.cooling.area.fix(100)
    m.fs.unit.cooling.overall_heat_transfer_coefficient.fix(675)
    m.fs.unit.condense.cold_side.deltaP[:].fix(0)
    m.fs.unit.condense.level.fix(0.275)
    m.fs.unit.condense.heater_diameter.fix(1.4)
    m.fs.unit.condense.vol_frac_shell.fix(0.675)
    m.fs.unit.condense.cond_sect_length.fix(6.4)
    m.fs.unit.condense.cold_side.volume.fix(1.3)

    return m


@pytest.mark.unit
def test_basic_build(build_unit):
    """Make a model and make sure it doesn't throw exception"""
    m = build_unit
    assert degrees_of_freedom(m) == 0
    # Check unit config arguments
    assert len(m.fs.unit.config) == 10
    assert m.fs.unit.config.has_desuperheat
    assert m.fs.unit.config.has_drain_cooling
    assert m.fs.unit.config.has_drain_mixer
    assert m.fs.unit.config.condense
    assert m.fs.unit.config.property_package is m.fs.properties


@pytest.mark.skipif(not iapws95.iapws95_available(), reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_unit(build_unit):
    initialization_tester(build_unit, dof=0)


# TODO
# @pytest.mark.integration
# def test_units(build_unit):
#     # ToDo: dynamic models cannot assert units consistent
#     # assert_units_consistent(build_unit)


@pytest.mark.integration
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_fwh_model(build_unit):
    m = build_unit
    # need to initialize because integration test do not run sequentially
    m.fs.unit.initialize()
    solver.options = {
        "tol": 1e-7,
        "linear_solver": "ma27",
        "max_iter": 50,
    }
    # initial drain flow rate
    drain_flow0 = m.fs.unit.cooling.hot_side_outlet.flow_mol[0].value
    assert pytest.approx(67.20607, abs=1e-2) == pyo.value(
        m.fs.unit.cooling.hot_side_outlet.flow_mol[0]
    )

    # change drain flow rate to increase water level
    for t in m.fs.time:
        if t >= 30:
            m.fs.unit.cooling.hot_side_outlet.flow_mol[t].fix(drain_flow0 * 0.95)
    m.fs.unit.condense.level.unfix()
    m.fs.unit.condense.level[0].fix()
    solver.solve(m, tee=True)
    assert degrees_of_freedom(m) == 0
    assert pytest.approx(0.2759, abs=1e-2) == pyo.value(m.fs.unit.condense.level[60])
