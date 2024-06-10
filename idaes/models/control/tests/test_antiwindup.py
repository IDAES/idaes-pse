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
Tests to ensure the controller behaves as expected when antiwindup is and isn't
activated. It turns out that you need to work really hard in order to get
windup issues in an (approximately) first order system. What to do? Use an
(approximately) second order system, of course!

           valve_1   +----+                +----+   valve_3
  steam ----|><|-->--| ta |    valve_2     | ta |----|><|-->--- steam
                     | nk |-----|><|--->---| nk |
                     | _1 |                | _2 |
                     +----+                +----+


"""

__author__ = "Douglas Allan, John Eslick"

import pytest
import numpy as np
import matplotlib.pyplot as plt
import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models import Heater, Valve
from idaes.models.properties import iapws95
from idaes.core.util.initialization import propagate_state
from idaes.models.control.controller import (
    PIDController,
    ControllerType,
    ControllerMVBoundType,
    ControllerAntiwindupType,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.math import safe_sqrt
from idaes.core.solvers import get_solver, petsc
import idaes.logger as idaeslog
from idaes.core.util.plot import plot_grid_dynamic
import idaes.core.util.scaling as iscale


def set_indexed_variable_bounds(var, bounds):
    for idx, subvar in var.items():
        subvar.bounds = bounds


def _valve_pressure_flow_cb(b):
    """
    For vapor F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    umeta = b.config.property_package.get_metadata().get_derived_units

    b.Cv = pyo.Var(
        initialize=0.1,
        doc="Valve flow coefficient",
        units=umeta("amount") / umeta("time") / umeta("pressure"),
    )
    b.Cv.fix()

    b.flow_var = pyo.Reference(b.control_volume.properties_in[:].flow_mol)
    b.pressure_flow_equation_scale = lambda x: x**2

    @b.Constraint(b.flowsheet().time)
    def pressure_flow_equation(b2, t):
        Po = b2.control_volume.properties_out[t].pressure
        Pi = b2.control_volume.properties_in[t].pressure
        F = b2.control_volume.properties_in[t].flow_mol
        Cv = b2.Cv
        fun = b2.valve_function[t]
        return F**2 == Cv**2 * (Pi**2 - Po**2) * fun**2


def _add_setpoint_step(m, time=1, value=6.0e5):
    """Easy function to add an inlet pressure step change

    Args:
        m (ConcreteModel): A valve and tank model
        time (float): Time for the step to occur
        value (float): New pressure at time

    Returns:
        None
    """
    for t in m.fs.time:
        if t >= time:
            m.fs.ctrl.setpoint[t].fix(value)


def create_model(
    steady_state=True,
    time_set=None,
    time_units=pyo.units.s,
    nfe=5,
    tee=False,
    calc_integ=True,
    derivative_on_error=False,
    antiwindup=ControllerAntiwindupType.NONE,
    initial_valve1_opening=1.0,
):
    """Create a test model and solver

    Args:
        steady_state (bool): If True, create a steady state model, otherwise
            create a dynamic model
        time_set (list): The beginning and end point of the time domain
        time_units (Pyomo Unit object): Units of time domain
        nfe (int): Number of finite elements argument for the DAE
            transformation.
        calc_integ (bool): If True, calculate in the initial condition for
            the integral term, else use a fixed variable (fs.ctrl.integral_of_error[0]),
            False is the better option if you have a value from a previous
            time period

    Returns
        (tuple): (ConcreteModel, Solver)
    """
    if steady_state:
        fs_cfg = {"dynamic": False}
        model_name = "Steam Tank, Steady State"
    else:
        fs_cfg = {"dynamic": True, "time_set": time_set, "time_units": time_units}
        model_name = "Steam Tank, Dynamic"

    if time_set is None:
        time_set = [0, 3]

    m = pyo.ConcreteModel(name=model_name)
    m.fs = FlowsheetBlock(**fs_cfg)
    # Create a property parameter block
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        phase_presentation=iapws95.PhaseType.G
    )
    m.fs.prop_water.set_default_scaling("flow_mol", 1e-2)
    # Create the valve and tank models
    m.fs.valve_1 = Valve(
        dynamic=False,
        has_holdup=False,
        pressure_flow_callback=_valve_pressure_flow_cb,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.tank_1 = Heater(
        has_holdup=True,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.valve_2 = Valve(
        dynamic=False,
        has_holdup=False,
        pressure_flow_callback=_valve_pressure_flow_cb,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.tank_2 = Heater(
        has_holdup=True,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.valve_3 = Valve(
        dynamic=False,
        has_holdup=False,
        pressure_flow_callback=_valve_pressure_flow_cb,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    if not steady_state:
        # Add a controller
        m.fs.ctrl = PIDController(
            process_var=m.fs.tank_2.control_volume.properties_out[:].pressure,
            manipulated_var=m.fs.valve_1.valve_opening,
            calculate_initial_integral=calc_integ,
            mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
            controller_type=ControllerType.PID,  # rather use PI, but testing all terms
            derivative_on_error=derivative_on_error,
            antiwindup_type=antiwindup,
        )

    # Connect the models
    m.fs.v1_to_tank_1 = Arc(source=m.fs.valve_1.outlet, destination=m.fs.tank_1.inlet)
    m.fs.tank_1_to_v2 = Arc(source=m.fs.tank_1.outlet, destination=m.fs.valve_2.inlet)
    m.fs.v2_to_tank_2 = Arc(source=m.fs.valve_2.outlet, destination=m.fs.tank_2.inlet)
    m.fs.tank_2_to_v3 = Arc(source=m.fs.tank_2.outlet, destination=m.fs.valve_3.inlet)

    # Add the stream constraints and do the DAE transformation
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    if not steady_state:
        pyo.TransformationFactory("dae.finite_difference").apply_to(
            m.fs, nfe=nfe, wrt=m.fs.time, scheme="BACKWARD"
        )

    # Fix the derivative variables to zero at time 0 (steady state assumption)
    m.fs.fix_initial_conditions()

    # Fix the input variables
    m.fs.valve_1.inlet.enth_mol.fix(55000)
    m.fs.valve_1.inlet.pressure.fix(6e5)
    m.fs.valve_3.outlet.pressure.fix(101325)
    for valve in [m.fs.valve_1, m.fs.valve_2, m.fs.valve_3]:
        valve.Cv.fix(0.001)
        set_indexed_variable_bounds(valve.deltaP, [None, 0])
    m.fs.valve_1.valve_opening.fix(initial_valve1_opening)
    m.fs.valve_2.valve_opening.fix(1)
    m.fs.valve_3.valve_opening.fix(1)
    m.fs.tank_1.heat_duty.fix(0)
    m.fs.tank_1.control_volume.volume.fix(15.0)
    m.fs.tank_2.heat_duty.fix(0)
    m.fs.tank_2.control_volume.volume.fix(15.0)
    if not steady_state:
        m.fs.ctrl.gain_p.fix(5e-6)
        m.fs.ctrl.gain_i.fix(2e-6)
        m.fs.ctrl.gain_d.fix(1e-6)
        m.fs.ctrl.derivative_term[m.fs.time.first()].fix(0)
        m.fs.ctrl.setpoint.fix(1.5e5)
        m.fs.ctrl.mv_ref.fix(0)
        m.fs.ctrl.mv_lb = 0.0
        m.fs.ctrl.mv_ub = 1.0
        if antiwindup == ControllerAntiwindupType.BACK_CALCULATION:
            m.fs.ctrl.gain_b.fix(1)

    for t in m.fs.time:
        # For debugging purposes.
        for valve in [m.fs.valve_1, m.fs.valve_2, m.fs.valve_3]:
            iscale.set_scaling_factor(valve.control_volume.work[t], 1e-6)
            iscale.set_scaling_factor(valve.valve_opening[t], 1)
        for tank in [m.fs.tank_1, m.fs.tank_2]:
            iscale.set_scaling_factor(tank.control_volume.heat[t], 1e-6)
            iscale.set_scaling_factor(tank.control_volume.volume[t], 1)
    iscale.calculate_scaling_factors(m)
    iscale.scale_time_discretization_equations(m.fs, m.fs.time, 1)

    # Initialize the model

    solver = get_solver(
        options={
            "max_iter": 300,
            "nlp_scaling_method": "user-scaling",
            "linear_solver": "ma57",
        }
    )

    for t in m.fs.time:
        m.fs.valve_1.inlet.flow_mol[t] = 100  # initial guess on flow
    # simple initialize
    m.fs.valve_1.initialize(outlvl=idaeslog.DEBUG)
    propagate_state(m.fs.v1_to_tank_1)
    m.fs.tank_1.initialize()
    propagate_state(m.fs.tank_1_to_v2)
    m.fs.valve_2.initialize()
    propagate_state(m.fs.v2_to_tank_2)
    m.fs.tank_2.initialize()
    propagate_state(m.fs.tank_2_to_v3)
    # Can't specify both flow and outlet pressure so free the outlet pressure
    # for initialization and refix it after.  Inlet flow gets fixed in init
    op = {}
    for t in m.fs.time:
        op[t] = pyo.value(m.fs.valve_3.outlet.pressure[t])
        m.fs.valve_3.outlet.pressure[t].unfix()
    m.fs.valve_3.initialize()
    for t in m.fs.time:
        m.fs.valve_3.outlet.pressure[t].fix(op[t])
    if not steady_state:
        m.fs.ctrl.deactivate()
        m.fs.valve_1.valve_opening.fix()
    print("Model solve for initialization")
    solver.solve(m, tee=tee)
    if not steady_state:
        m.fs.ctrl.activate()
        m.fs.valve_1.valve_opening.unfix()
        m.fs.valve_1.valve_opening[m.fs.time.first()].fix()
        print("Model solve for controller initialization")
        solver.solve(m, tee=tee)
    # Return the model and solver
    return m, solver


def get_valve_openings(P_inlet, P_out_1, P_out_2):
    # First calculate the two steady states that should be achieved in the test
    # don't worry these steady state problems solve super fast
    m_steady, solver = create_model(tee=False)
    m_steady.fs.valve_1.inlet.pressure.fix(P_inlet)
    m_steady.fs.tank_2.control_volume.properties_out[0].pressure.fix(P_out_1)
    m_steady.fs.valve_1.valve_opening[0].unfix()
    solver.solve(m_steady, tee=False)
    s1_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])
    print(f"Steady state 1 valve opening: {s1_valve}")

    m_steady.fs.tank_2.control_volume.properties_out[0].pressure.fix(P_out_2)
    solver.solve(m_steady, tee=False)
    s2_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])
    print(f"Steady state 2 valve opening: {s2_valve}")
    return s1_valve, s2_valve


@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
@pytest.mark.component
def test_setpoint_change_windup():
    s1_valve, s2_valve = get_valve_openings(6.0e5, 1.5e5, 3.0e5)

    tsim = 40

    # Next create a model for the 0 to 5 sec time period
    m_dynamic, solver = create_model(
        steady_state=False,
        time_set=[0, 2, tsim],
        nfe=2,
        calc_integ=True,
        tee=True,
        derivative_on_error=False,
        initial_valve1_opening=s1_valve,
        antiwindup=ControllerAntiwindupType.NONE,
    )
    # return m_dynamic, solver
    # Retune controller to result in windup
    m_dynamic.fs.ctrl.gain_p.fix(5e-6)
    m_dynamic.fs.ctrl.gain_i.fix(1e-5)
    m_dynamic.fs.ctrl.gain_d.fix(0)
    # Add a step change in outlet setpoint
    _add_setpoint_step(m_dynamic, time=m_dynamic.fs.time.at(3), value=3.0e5)

    res = petsc.petsc_dae_by_time_element(
        m_dynamic,
        time=m_dynamic.fs.time,
        between=[
            m_dynamic.fs.time.first(),
            m_dynamic.fs.time.at(2),
            m_dynamic.fs.time.last(),
        ],
        ts_options={
            "--ts_type": "beuler",  # Backwards Euler
            "--ts_adapt_type": "basic",
            "--ts_dt": 1,
            "--ts_save_trajectory": 1,
            # "--snes_monitor": "",
            # "--ksp_monitor": "",
            # "--ts_monitor": "",
            "--ts_max_snes_failures": 1000,
        },
    )
    # read the trajectory data, and make it easy by interpolating a time point
    # every second
    tj = res.trajectory
    tj2 = tj.interpolate(np.linspace(0, tsim, tsim + 1))
    tf = m_dynamic.fs.time.last()

    # Test for huge value  of integral term and oscillatory control trajectory
    assert tj2.get_vec(m_dynamic.fs.ctrl.mv_integral_component[tf])[8] >= 2.95
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[13] == pytest.approx(
        1, abs=1e-4
    )
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[16] <= 0.32
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[19] >= 0.82
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[23] <= 0.56

    return m_dynamic, solver, tj2


@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
@pytest.mark.component
def test_setpoint_change_conditional_integration():
    s1_valve, s2_valve = get_valve_openings(6.0e5, 1.5e5, 3.0e5)

    tsim = 40

    # Next create a model for the 0 to 5 sec time period
    m_dynamic, solver = create_model(
        steady_state=False,
        time_set=[0, 2, tsim],
        nfe=2,
        calc_integ=True,
        tee=True,
        derivative_on_error=False,
        initial_valve1_opening=s1_valve,
        antiwindup=ControllerAntiwindupType.CONDITIONAL_INTEGRATION,
    )
    # return m_dynamic, solver
    # Retune controller to result in windup
    m_dynamic.fs.ctrl.gain_p.fix(5e-6)
    m_dynamic.fs.ctrl.gain_i.fix(1e-5)
    m_dynamic.fs.ctrl.gain_d.fix(0)
    # Add a step change in outlet setpoint
    _add_setpoint_step(m_dynamic, time=m_dynamic.fs.time.at(3), value=3.0e5)

    res = petsc.petsc_dae_by_time_element(
        m_dynamic,
        time=m_dynamic.fs.time,
        between=[
            m_dynamic.fs.time.first(),
            m_dynamic.fs.time.at(2),
            m_dynamic.fs.time.last(),
        ],
        ts_options={
            "--ts_type": "beuler",  # Backwards Euler
            "--ts_adapt_type": "basic",
            "--ts_dt": 1,
            "--ts_save_trajectory": 1,
            "--ts_max_snes_failures": 1000,
        },
    )
    # read the trajectory data, and make it easy by interpolating a time point
    # every second
    tj = res.trajectory
    tj2 = tj.interpolate(np.linspace(0, tsim, tsim + 1))
    tf = m_dynamic.fs.time.last()

    # Test that integral component of MV doesn't exceed MV bound
    assert np.all(tj2.get_vec(m_dynamic.fs.ctrl.mv_integral_component[tf]) <= 1)
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[10] >= 0.585
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[14] <= 0.72
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[18] >= 0.63
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[28] == pytest.approx(
        s2_valve, abs=1e-3
    )

    return m_dynamic, solver, tj2


@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
@pytest.mark.component
def test_setpoint_change_back_calculation():
    s1_valve, s2_valve = get_valve_openings(6.0e5, 1.5e5, 3.0e5)

    tsim = 40

    # Next create a model for the 0 to 5 sec time period
    m_dynamic, solver = create_model(
        steady_state=False,
        time_set=[0, 2, tsim],
        nfe=2,
        calc_integ=True,
        tee=True,
        derivative_on_error=False,
        initial_valve1_opening=s1_valve,
        antiwindup=ControllerAntiwindupType.BACK_CALCULATION,
    )
    # return m_dynamic, solver
    # Retune controller to result in windup
    m_dynamic.fs.ctrl.gain_p.fix(5e-6)
    m_dynamic.fs.ctrl.gain_i.fix(1e-5)
    m_dynamic.fs.ctrl.gain_d.fix(0)
    m_dynamic.fs.ctrl.gain_b.fix(1.5)
    # Add a step change in outlet setpoint
    _add_setpoint_step(m_dynamic, time=m_dynamic.fs.time.at(3), value=3.0e5)

    res = petsc.petsc_dae_by_time_element(
        m_dynamic,
        time=m_dynamic.fs.time,
        between=[
            m_dynamic.fs.time.first(),
            m_dynamic.fs.time.at(2),
            m_dynamic.fs.time.last(),
        ],
        ts_options={
            "--ts_type": "beuler",  # Backwards Euler
            "--ts_adapt_type": "basic",
            "--ts_dt": 1,
            "--ts_save_trajectory": 1,
            "--ts_max_snes_failures": 1000,
        },
    )
    # read the trajectory data, and make it easy by interpolating a time point
    # every second
    tj = res.trajectory
    tj2 = tj.interpolate(np.linspace(0, tsim, tsim + 1))
    tf = m_dynamic.fs.time.last()

    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[5] == pytest.approx(
        1, abs=1e-3
    )
    # Test that integral component of MV doesn't far exceed MV bound
    assert np.all(tj2.get_vec(m_dynamic.fs.ctrl.mv_integral_component[tf]) <= 1.2)
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[10] >= 0.565
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[14] <= 0.73
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[18] >= 0.62
    assert tj2.get_vec(m_dynamic.fs.valve_1.valve_opening[tf])[40] == pytest.approx(
        s2_valve, abs=1e-3
    )

    return m_dynamic, solver, tj2


if __name__ == "__main__":
    m, solver, tj2 = test_setpoint_change_windup()

    fig = plt.figure(figsize=(4, 5.5))
    axes = fig.subplots(nrows=3, ncols=1)
    time = tj2.get_vec("_time")
    tf = m.fs.time.last()

    axes[0].plot(time, tj2.get_vec(m.fs.valve_1.valve_opening[tf]))
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("opening (fraction open)")

    axes[1].plot(
        time, tj2.get_vec(m.fs.tank_2.control_volume.properties_out[tf].pressure) / 1000
    )
    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("tank pressure (kPa)")

    axes[2].plot(time, tj2.get_vec(m.fs.ctrl.mv_integral_component[tf]))
    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("integral component (fraction open)")

    fig.tight_layout()
