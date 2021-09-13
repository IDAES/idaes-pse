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
from idaes.generic_models.unit_models import Heater, Valve
from idaes.generic_models.properties import iapws95
from idaes.core.util import copy_port_values as _set_port
from idaes.core.util.plot import stitch_dynamic
from idaes.generic_models.control import PIDBlock, PIDForm
import idaes.core.util.scaling as iscale
from idaes.core.util import get_solver


def _valve_pressure_flow_cb(b):
    """
    For vapor F = Cv*sqrt(Pi**2 - Po**2)*f(x)
    """
    umeta = b.config.property_package.get_metadata().get_derived_units

    b.Cv = pyo.Var(
        initialize=0.1,
        doc="Valve flow coefficent",
        units=umeta("amount")/umeta("time")/umeta("pressure")
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
        return F ** 2 == Cv ** 2 * (Pi ** 2 - Po ** 2) * fun ** 2


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
    calc_integ=True,
    form=PIDForm.standard,
    tee=False
):
    """ Create a test model and solver

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
        form: whether the equations are written in the standard or velocity
            form

    Returns
        (tuple): (ConcreteModel, Solver)
    """
    if steady_state:
        fs_cfg = {"dynamic": False}
        model_name = "Steam Tank, Steady State"
    else:
        fs_cfg = {
            "dynamic": True, "time_set": time_set, "time_units": time_units}
        model_name = "Steam Tank, Dynamic"

    if time_set is None:
        time_set = [0, 3]

    m = pyo.ConcreteModel(name=model_name)
    m.fs = FlowsheetBlock(default=fs_cfg)
    # Create a property parameter block
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG})
    # Create the valve and tank models
    m.fs.valve_1 = Valve(default={
        "dynamic": False,
        "has_holdup": False,
        "pressure_flow_callback": _valve_pressure_flow_cb,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "property_package": m.fs.prop_water})
    m.fs.tank = Heater(default={
        "has_holdup": True,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "property_package": m.fs.prop_water})
    m.fs.valve_2 = Valve(default={
        "dynamic": False,
        "has_holdup": False,
        "pressure_flow_callback": _valve_pressure_flow_cb,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "property_package": m.fs.prop_water})

    # Connect the models
    m.fs.v1_to_t = Arc(source=m.fs.valve_1.outlet, destination=m.fs.tank.inlet)
    m.fs.t_to_v2 = Arc(source=m.fs.tank.outlet, destination=m.fs.valve_2.inlet)

    # The control volume block doesn't assume the two phases are in equilibrium
    # by default, so I'll make that assumption here, I don't actually expect
    # liquid to form but who knows. The phase_fraction in the control volume is
    # volumetric phase fraction hence the densities.
    @m.fs.tank.Constraint(m.fs.time)
    def vol_frac_vap(b, t):
        return (b.control_volume.properties_out[t].phase_frac["Vap"] *
                b.control_volume.properties_out[t].dens_mol /
                b.control_volume.properties_out[t].dens_mol_phase["Vap"]) == (
                    b.control_volume.phase_fraction[t, "Vap"])

    # Add the stream constraints and do the DAE transformation
    pyo.TransformationFactory('network.expand_arcs').apply_to(m.fs)
    if not steady_state:
        pyo.TransformationFactory('dae.finite_difference').apply_to(
            m.fs, nfe=nfe, wrt=m.fs.time, scheme='BACKWARD')

    # Fix the derivative variables to zero at time 0 (steady state assumption)
    m.fs.fix_initial_conditions()

    # A tank pressure reference that's directly time-indexed
    m.fs.tank_pressure = pyo.Reference(
        m.fs.tank.control_volume.properties_out[:].pressure)

    # Add a controller
    m.fs.ctrl = PIDBlock(default={"pv": m.fs.tank_pressure,
                                  "output": m.fs.valve_1.valve_opening,
                                  "upper": 1.0,
                                  "lower": 0.0,
                                  "calculate_initial_integral": calc_integ,
                                  "pid_form": form})
    m.fs.ctrl.deactivate()  # Don't want controller turned on by default

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
    m.fs.ctrl.gain.fix(1e-6)
    m.fs.ctrl.time_i.fix(0.1)
    m.fs.ctrl.time_d.fix(0.1)
    m.fs.ctrl.setpoint.fix(3e5)

    # Initialize the model
    solver = get_solver(options={"max_iter":50})

    for t in m.fs.time:
        m.fs.valve_1.inlet.flow_mol[t] = 100  # initial guess on flow
    # simple initialize
    m.fs.valve_1.initialize()
    _set_port(m.fs.tank.inlet, m.fs.valve_1.outlet)
    m.fs.tank.initialize()
    _set_port(m.fs.valve_2.inlet, m.fs.tank.outlet)
    # Can't specify both flow and outlet pressure so free the outlet pressure
    # for initialization and refix it after.  Inlet flow gets fixed in init
    op = {}
    for t in m.fs.time:
        op[t] = pyo.value(m.fs.valve_2.outlet.pressure[t])
        m.fs.valve_2.outlet.pressure[t].unfix()
    m.fs.valve_2.initialize()
    for t in m.fs.time:
        m.fs.valve_2.outlet.pressure[t].fix(op[t])
    solver.solve(m, tee=tee)
    # Return the model and solver
    return m, solver


def tpid(form, tee=False):
    """This test is pretty course-grained, but it should cover everything"""

    # First calculate the two steady states that should be achieved in the test
    # don't worry these steady state problems solve super fast
    m_steady, solver = create_model(tee=tee)
    m_steady.fs.tank_pressure[0].fix(3e5)
    m_steady.fs.valve_1.valve_opening[0].unfix()
    solver.solve(m_steady, tee=tee)
    s1_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])
    solver.solve(m_steady, tee=tee)
    m_steady.fs.valve_1.inlet.pressure.fix(5.5e5)
    solver.solve(m_steady, tee=tee)
    s2_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])

    # Next create a model for the 0 to 5 sec time period
    m_dynamic, solver = create_model(
        steady_state=False,
        time_set=[0, 5],
        nfe=10,
        calc_integ=True,
        form=form,
        tee=tee
    )
    # Turn on control and solve since the setpoint is different than the
    # steady-state solution, stuff happens, also pressure step at t=5
    m_dynamic.fs.ctrl.activate()
    m_dynamic.fs.valve_1.valve_opening.unfix()
    m_dynamic.fs.valve_1.valve_opening[0].fix()

    # Add a step change right at the end of the interval, to make sure we can
    # get a continuous solution across the two models
    _add_inlet_pressure_step(m_dynamic, time=4.5, value=5.5e5)
    iscale.calculate_scaling_factors(m_dynamic)
    solver.solve(m_dynamic, tee=tee)

    # Now create a model for the 5 to 10 second interval and set the inital
    # conditions of the first model to the final (unsteady) state of the
    # previous time interval model
    m_dynamic2, solver = create_model(
        steady_state=False,
        time_set=[5, 10],
        nfe=10,
        calc_integ=False,
        form=form,
    )

    m_dynamic2.fs.valve_1.valve_opening.fix(
        m_dynamic.fs.valve_1.valve_opening[5].value)
    m_dynamic2.fs.valve_1.inlet.enth_mol.fix(
        m_dynamic.fs.valve_1.inlet.enth_mol[5].value)
    m_dynamic2.fs.valve_1.inlet.pressure.fix(
        m_dynamic.fs.valve_1.inlet.pressure[5].value)
    for i in m_dynamic2.fs.tank.control_volume.material_accumulation:
        if i[0] == 5:
            m_dynamic2.fs.tank.control_volume.material_accumulation[i].value =\
                m_dynamic.fs.tank.control_volume.material_accumulation[i].value
    for i in m_dynamic2.fs.tank.control_volume.energy_accumulation:
        if i[0] == 5:
            m_dynamic2.fs.tank.control_volume.energy_accumulation[i].value =\
                m_dynamic.fs.tank.control_volume.energy_accumulation[i].value
    m_dynamic2.fs.ctrl.err_d0.fix(pyo.value(m_dynamic.fs.ctrl.err_d[5]))
    m_dynamic2.fs.ctrl.err_i0.fix(pyo.value(m_dynamic.fs.ctrl.err_i_end))

    iscale.calculate_scaling_factors(m_dynamic2)
    # As a lazy form of initialization, solve the steady state problem before
    # turning on the controller.
    solver.solve(m_dynamic2, tee=tee)

    # Now turn on control and solve again.
    m_dynamic2.fs.ctrl.activate()
    m_dynamic2.fs.valve_1.valve_opening.unfix()
    m_dynamic2.fs.valve_1.valve_opening[5].fix()
    solver.solve(m_dynamic2, tee=tee)

    # Check that the model hit steady state about. The tolerance is loose
    # because I used a low nfe in an attempt to speed it up
    t = m_dynamic.fs.time.get_lower_element_boundary(4.45)
    assert s1_valve == pytest.approx(
        pyo.value(m_dynamic.fs.valve_1.valve_opening[t]), abs=5e-3)

    # The second should hit steady state pretty quick if the two models
    # line up smoothly.
    t = m_dynamic2.fs.time.get_lower_element_boundary(7)
    assert s2_valve == pytest.approx(
        pyo.value(m_dynamic2.fs.valve_1.valve_opening[t]), abs=5e-3)

    stitch_time = stitch_dynamic(m_dynamic.fs.time, m_dynamic2.fs.time)
    stitch_valve = stitch_dynamic(m_dynamic.fs.valve_1.valve_opening,
                                  m_dynamic2.fs.valve_1.valve_opening)

    # test the stitch function used to plot time dependent variables across
    # the two models
    for i, t in enumerate(m_dynamic.fs.time):
        assert t == stitch_time[i]
        assert m_dynamic.fs.valve_1.valve_opening[t].value == stitch_valve[i]
    for j, t in enumerate(m_dynamic2.fs.time):
        i = j + len(m_dynamic.fs.time)
        assert t == stitch_time[i]
        assert m_dynamic2.fs.valve_1.valve_opening[t].value == stitch_valve[i]
    return m_dynamic, m_dynamic2, solver



@pytest.mark.integration
def test_pid_velocity():
    """This test is pretty course-grained, but it should cover everything"""
    tpid(PIDForm.velocity)


@pytest.mark.integration
def test_pid_standard():
    """This test is pretty course-grained, but it should cover everything"""
    tpid(PIDForm.standard)

if __name__ == '__main__':
    m, m2, solver = tpid(PIDForm.velocity, tee=True)
