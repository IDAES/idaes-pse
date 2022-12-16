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
Basic PID control test using a tank with steam flowing through it.  The
controller is set up to maintain a tank pressure by adjusting the inlet valve
position.

           valve_1   +----+
  steam ----|><|-->--| ta |    valve_2
                     | nk |-----|><|--->--- steam
                     +----+

To test there are two basic things:
  1. That the problem with control goes to the steady state solution
  2. That the dynamic problem with control can be split smoothly across two
     different models representing different adjacent time periods.
"""

__author__ = "John Eslick"

import pytest
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
)
from idaes.core.solvers import get_solver
from idaes.core.util.plot import plot_grid_dynamic


def _valve_pressure_flow_cb(b):
    """
    For vapor F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    umeta = b.config.property_package.get_metadata().get_derived_units

    b.Cv = pyo.Var(
        initialize=0.1,
        doc="Valve flow coefficent",
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


def _add_inlet_pressure_step(m, time=1, value=6.0e5):
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
            m.fs.valve_1.inlet.pressure[t].fix(value)


def create_model(
    steady_state=True,
    time_set=None,
    time_units=pyo.units.s,
    nfe=5,
    tee=False,
    calc_integ=True,
):
    """Create a test model and solver

    Args:
        steady_state (bool): If True, create a steady state model, otherwise
            create a dynamic model
        time_set (list): The begining and end point of the time domain
        time_units (Pyomo Unit object): Units of time domain
        nfe (int): Number of finite elements argument for the DAE
            transformation.
        calc_integ (bool): If True, calculate in the initial condition for
            the integral term, else use a fixed variable (fs.ctrl.err_i0),
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
        phase_presentation=iapws95.PhaseType.LG
    )
    # Create the valve and tank models
    m.fs.valve_1 = Valve(
        dynamic=False,
        has_holdup=False,
        pressure_flow_callback=_valve_pressure_flow_cb,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.tank = Heater(
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
    if not steady_state:
        # Add a controller
        m.fs.ctrl = PIDController(
            process_var=m.fs.tank.control_volume.properties_out[:].pressure,
            manipulated_var=m.fs.valve_1.valve_opening,
            calculate_initial_integral=calc_integ,
            mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
            type=ControllerType.PID,
        )

    # The control volume block doesn't assume the two phases are in equilibrium
    # by default, so I'll make that assumption here, I don't actually expect
    # liquid to form but who knows. The phase_fraction in the control volume is
    # volumetric phase fraction hence the densities.
    @m.fs.tank.Constraint(m.fs.time)
    def vol_frac_vap(b, t):
        return (
            b.control_volume.properties_out[t].phase_frac["Vap"]
            * b.control_volume.properties_out[t].dens_mol
            / b.control_volume.properties_out[t].dens_mol_phase["Vap"]
        ) == (b.control_volume.phase_fraction[t, "Vap"])

    # Connect the models
    m.fs.v1_to_tank = Arc(source=m.fs.valve_1.outlet, destination=m.fs.tank.inlet)
    m.fs.tank_to_v2 = Arc(source=m.fs.tank.outlet, destination=m.fs.valve_2.inlet)

    # Add the stream constraints and do the DAE transformation
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    if not steady_state:
        pyo.TransformationFactory("dae.finite_difference").apply_to(
            m.fs, nfe=nfe, wrt=m.fs.time, scheme="BACKWARD"
        )

    # Fix the derivative variables to zero at time 0 (steady state assumption)
    m.fs.fix_initial_conditions()

    # Fix the input variables
    m.fs.valve_1.inlet.enth_mol.fix(50000)
    m.fs.valve_1.inlet.pressure.fix(5e5)
    m.fs.valve_2.outlet.pressure.fix(101325)
    m.fs.valve_1.Cv.fix(0.001)
    m.fs.valve_2.Cv.fix(0.001)
    m.fs.valve_1.valve_opening.fix(1)
    m.fs.valve_2.valve_opening.fix(1)
    m.fs.tank.heat_duty.fix(0)
    m.fs.tank.control_volume.volume.fix(2.0)
    if not steady_state:
        m.fs.ctrl.gain_p.fix(1e-6)
        m.fs.ctrl.gain_i.fix(1e-5)
        m.fs.ctrl.gain_d.fix(1e-6)
        m.fs.ctrl.derivative_of_error[m.fs.time.first()].fix(0)
        m.fs.ctrl.setpoint.fix(3e5)
        m.fs.ctrl.mv_ref.fix(0)
        m.fs.ctrl.mv_lb = 0.0
        m.fs.ctrl.mv_ub = 1.0

    # Initialize the model
    solver = get_solver(options={"max_iter": 50})

    for t in m.fs.time:
        m.fs.valve_1.inlet.flow_mol[t] = 100  # initial guess on flow
    # simple initialize
    m.fs.valve_1.initialize()
    propagate_state(m.fs.v1_to_tank)
    m.fs.tank.initialize()
    propagate_state(m.fs.tank_to_v2)
    # Can't specify both flow and outlet pressure so free the outlet pressure
    # for initialization and refix it after.  Inlet flow gets fixed in init
    op = {}
    for t in m.fs.time:
        op[t] = pyo.value(m.fs.valve_2.outlet.pressure[t])
        m.fs.valve_2.outlet.pressure[t].unfix()
    m.fs.valve_2.initialize()
    for t in m.fs.time:
        m.fs.valve_2.outlet.pressure[t].fix(op[t])
    if not steady_state:
        m.fs.ctrl.deactivate()
        m.fs.valve_1.valve_opening.fix()
    solver.solve(m, tee=tee)
    if not steady_state:
        m.fs.ctrl.activate()
        m.fs.valve_1.valve_opening.unfix()
        m.fs.valve_1.valve_opening[m.fs.time.first()].fix()
        solver.solve(m, tee=tee)
    # Return the model and solver
    return m, solver


@pytest.mark.integration
def test_pid():
    """This test is pretty course-grained, but it should cover everything"""

    # First calculate the two steady states that should be achieved in the test
    # don't worry these steady state problems solve super fast
    m_steady, solver = create_model(tee=False)
    m_steady.fs.valve_1.inlet.pressure.fix(5.0e5)
    m_steady.fs.tank.control_volume.properties_out[0].pressure.fix(3e5)
    m_steady.fs.valve_1.valve_opening[0].unfix()
    solver.solve(m_steady, tee=False)
    s1_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])

    m_steady.fs.valve_1.inlet.pressure.fix(5.5e5)
    solver.solve(m_steady, tee=False)
    s2_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])

    # Next create a model for the 0 to 5 sec time period
    m_dynamic, solver = create_model(
        steady_state=False, time_set=[0, 12], nfe=30, calc_integ=True, tee=False
    )

    # Add a step change right in inlet pressure
    _add_inlet_pressure_step(m_dynamic, time=6, value=5.5e5)
    solver.solve(m_dynamic, tee=False)

    # Check that we reach the expected steady state (almost) by t = 5.6 and t=12
    assert pyo.value(
        m_dynamic.fs.valve_1.valve_opening[m_dynamic.fs.time.at(15)]
    ) == pytest.approx(s1_valve, abs=0.001)
    assert pyo.value(
        m_dynamic.fs.valve_1.valve_opening[m_dynamic.fs.time.last()]
    ) == pytest.approx(s2_valve, abs=0.001)

    return m_dynamic, solver


if __name__ == "__main__":
    m, solver = test_pid()

    plot_grid_dynamic(
        x=m.fs.time,
        xlabel="time (s)",
        y=[
            m.fs.valve_1.valve_opening,
            m.fs.tank.control_volume.properties_out[:].pressure,
            m.fs.valve_1.control_volume.properties_in[:].pressure,
        ],
        ylabel=[
            "opening (fraction open)",
            "tank pressure (kPa)",
            "inlet pressure (kPa)",
        ],
        yunits=[
            None,
            pyo.units.kPa,
            pyo.units.kPa,
        ],
        cols=3,
        rows=1,
        to_file="pid_steam_tank_pressure.pdf",
    )
