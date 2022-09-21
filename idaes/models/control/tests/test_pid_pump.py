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

__author__ = "Jinliang Ma"

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

from idaes.core import FlowsheetBlock
from idaes.models_extra.power_generation.unit_models.watertank import WaterTank
from idaes.models.control.controller import PIDController, ControllerType
import idaes.core.util.scaling as iscale

# TODO: Should have a test for this that does not depend on models_extra
from idaes.models_extra.power_generation.unit_models.helm import (
    HelmValve as WaterValve,
    HelmIsentropicCompressor as WaterPump,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models.properties import iapws95
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values
from idaes.core.solvers import get_solver

import pytest


def set_scaling_factors(m):
    """Set scaling factors for variables and expressions. These are used for
    variable scaling and used by the framework to scale constraints.

    Args:
        m: plant model to set scaling factors for.

    Returns:
        None
    """
    # Set scaling factors for boiler system
    iscale.set_scaling_factor(m.fs.tank.control_volume.energy_holdup, 1e-10)
    iscale.set_scaling_factor(m.fs.tank.control_volume.material_holdup, 1e-6)
    if m.dynamic:
        for t, c in m.fs.tank.control_volume.energy_accumulation_disc_eq.items():
            iscale.constraint_scaling_transform(c, 1e-6)

    # scaling factor for control valves
    for t in m.fs.time:
        iscale.set_scaling_factor(
            m.fs.valve.control_volume.properties_in[t].flow_mol, 0.001
        )

    # Calculate calculated scaling factors
    iscale.calculate_scaling_factors(m)


@pytest.fixture
def m():
    m_ss = get_model(dynamic=False)
    m_dyn = get_model(dynamic=True)
    copy_non_time_indexed_values(m_dyn.fs, m_ss.fs, copy_fixed=True)
    for t in m_dyn.fs.time:
        copy_values_at_time(m_dyn.fs, m_ss.fs, t, 0.0, copy_fixed=True)
    m_dyn.fs.controller.mv_ref.value = m_dyn.fs.valve.valve_opening[0].value
    # calculate integral error assuming error is zero
    m_dyn.fs.controller.integral_of_error[:].value = 0

    dof = degrees_of_freedom(m_dyn)
    assert dof == 0
    run_dynamic(m_dyn)
    return m_dyn


def get_model(dynamic=False):
    m = pyo.ConcreteModel(name="Testing PID controller model")
    if dynamic:
        m.dynamic = True
        m.fs = FlowsheetBlock(
            dynamic=True, time_set=[0, 50, 1000], time_units=pyo.units.s
        )
    else:
        m.dynamic = False
        m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()

    # water pump
    m.fs.pump = WaterPump(dynamic=False, property_package=m.fs.prop_water)

    # water tank
    m.fs.tank = WaterTank(
        tank_type="simple_tank", has_holdup=True, property_package=m.fs.prop_water
    )

    m.fs.valve = WaterValve(
        dynamic=False, has_holdup=False, phase="Liq", property_package=m.fs.prop_water
    )

    if dynamic:
        m.fs.controller = PIDController(
            process_var=m.fs.tank.tank_level,
            manipulated_var=m.fs.valve.valve_opening,
            type=ControllerType.PI,
        )

        m.discretizer = pyo.TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=20, wrt=m.fs.time, scheme="BACKWARD")

        m.fs.controller.gain_p.fix(-1e-1)
        m.fs.controller.gain_i.fix(-1e-2)
        m.fs.controller.setpoint.fix(5)
        m.fs.controller.mv_ref.fix(0.5)
        # m.fs.controller.integral_of_error[0].fix(0)

    m.fs.pump_to_tank = Arc(source=m.fs.pump.outlet, destination=m.fs.tank.inlet)
    m.fs.tank_to_valve = Arc(source=m.fs.tank.outlet, destination=m.fs.valve.inlet)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)

    m.fs.pump.efficiency_isentropic.fix(0.8)
    m.fs.pump.deltaP.fix(5e4)

    m.fs.tank.tank_cross_sect_area.fix(20)
    m.fs.tank.tank_level[:].fix(5)

    m.fs.valve.Cv.fix(63.4676)
    m.fs.valve.valve_opening.fix(0.5)
    m.fs.pump.inlet.flow_mol[:].fix(7000)
    m.fs.pump.inlet.enth_mol[:].fix(3924)
    m.fs.pump.inlet.pressure[:].fix(1e7)

    m.fs.valve.outlet.pressure[:].fix(10050000)

    set_scaling_factors(m)

    if dynamic:
        m.fs.tank.set_initial_condition()
        m.fs.tank.tank_level.unfix()
        m.fs.tank.tank_level[0].fix()
        m.fs.valve.valve_opening.unfix()
    else:
        m.fs.pump.initialize()

        m.fs.tank.inlet.flow_mol[:].value = m.fs.pump.outlet.flow_mol[0].value
        m.fs.tank.inlet.enth_mol[:].value = m.fs.pump.outlet.enth_mol[0].value
        m.fs.tank.inlet.pressure[:].value = m.fs.pump.outlet.pressure[0].value

        m.fs.tank.initialize()

        m.fs.valve.inlet.flow_mol[:].value = m.fs.tank.outlet.flow_mol[0].value
        m.fs.valve.inlet.enth_mol[:].value = m.fs.tank.outlet.enth_mol[0].value
        m.fs.valve.inlet.pressure[:].value = m.fs.tank.outlet.pressure[0].value

        # Solve for outlet pressure
        m.fs.valve.initialize()

        m.fs.valve.valve_opening.unfix()
        dof = degrees_of_freedom(m)
        assert dof == 0
        solver = get_solver()
        solver.solve(m, tee=True)

    return m


def run_dynamic(m):
    solver = get_solver()

    # add step change
    for t in m.fs.time:
        if t >= 50:
            m.fs.pump.inlet.flow_mol[t].fix(m.fs.pump.inlet.flow_mol[0].value * 1.2)
        else:
            m.fs.pump.inlet.flow_mol[t].fix(m.fs.pump.inlet.flow_mol[0].value)
    solver.solve(m, tee=True)

    return m


@pytest.mark.integration
def test_pid(m):
    assert 0.5000 == pytest.approx(pyo.value(m.fs.valve.valve_opening[0.0]), abs=1e-3)
    assert 0.6156 == pytest.approx(
        pyo.value(m.fs.valve.valve_opening[406.25]), abs=1e-3
    )
    assert 0.60268 == pytest.approx(pyo.value(m.fs.valve.valve_opening[1000]), abs=1e-3)
    assert 5 == pytest.approx(pyo.value(m.fs.tank.tank_level[0.0]), abs=1e-3)
    assert 4.9762 == pytest.approx(pyo.value(m.fs.tank.tank_level[406.25]), abs=1e-3)
    assert 4.99923 == pytest.approx(pyo.value(m.fs.tank.tank_level[1000]), abs=1e-3)
    assert 7000.0 == pytest.approx(pyo.value(m.fs.tank.outlet.flow_mol[0.0]), abs=1e-3)
    assert 8598.66642 == pytest.approx(
        pyo.value(m.fs.tank.outlet.flow_mol[406.25]), abs=1e-3
    )
    assert 8436.87597 == pytest.approx(
        pyo.value(m.fs.tank.outlet.flow_mol[1000]), abs=1e-3
    )
