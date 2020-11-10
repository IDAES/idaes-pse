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
from idaes.power_generation.unit_models import SteamValve
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.generic_models.unit_models import Heater
from idaes.generic_models.properties import iapws95
from idaes.core.util import copy_port_values as _set_port
from idaes.core.util.plot import stitch_dynamic
from idaes.generic_models.control import PIDBlock, PIDForm

solver_available = pyo.SolverFactory('ipopt').available()
prop_available = iapws95.iapws95_available()


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
    time_set=[0, 3],
    time_units=pyo.units.s,
    nfe=5,
    calc_integ=True,
    form=PIDForm.standard
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

    m = pyo.ConcreteModel(name=model_name)
    m.fs = FlowsheetBlock(default=fs_cfg)
    # Create a property parameter block
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG})
    # Create the valve and tank models
    m.fs.valve_1 = SteamValve(default={
        "dynamic": False,
        "has_holdup": False,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "property_package": m.fs.prop_water})
    m.fs.tank = Heater(default={
        "has_holdup": True,
        "material_balance_type": MaterialBalanceType.componentTotal,
        "property_package": m.fs.prop_water})
    m.fs.valve_2 = SteamValve(default={
        "dynamic": False,
        "has_holdup": False,
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
    solver = pyo.SolverFactory("ipopt")
    solver.options = {'tol': 1e-6,
                      'linear_solver': "ma27",
                      'max_iter': 100}
    for t in m.fs.time:
        m.fs.valve_1.inlet.flow_mol = 100  # initial guess on flow
    # simple initialize
    m.fs.valve_1.initialize(outlvl=1)
    _set_port(m.fs.tank.inlet, m.fs.valve_1.outlet)
    m.fs.tank.initialize(outlvl=1)
    _set_port(m.fs.valve_2.inlet, m.fs.tank.outlet)
    m.fs.valve_2.initialize(outlvl=1)
    solver.solve(m, tee=True)

    # Return the model and solver
    return m, solver


def tpid(form):
    """This test is pretty course-grained, but it should cover everything"""

    # First calculate the two steady states that should be achieved in the test
    # don't worry these steady state problems solve super fast
    m_steady, solver = create_model()
    m_steady.fs.tank_pressure[0].fix(3e5)
    m_steady.fs.valve_1.valve_opening[0].unfix()
    solver.solve(m_steady, tee=True)
    s1_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])
    solver.solve(m_steady, tee=True)
    m_steady.fs.valve_1.inlet.pressure.fix(5.5e5)
    solver.solve(m_steady, tee=True)
    s2_valve = pyo.value(m_steady.fs.valve_1.valve_opening[0])

    # Next create a model for the 0 to 5 sec time period
    m_dynamic, solver = create_model(
        steady_state=False,
        time_set=[0, 5],
        nfe=10,
        calc_integ=True,
        form=form,
    )

    # Turn on control and solve since the setpoint is different than the
    # steady-state solution, stuff happens, also pressure step at t=5
    m_dynamic.fs.ctrl.activate()
    m_dynamic.fs.valve_1.valve_opening.unfix()
    m_dynamic.fs.valve_1.valve_opening[0].fix()

    # Add a step change right at the end of the interval, to make sure we can
    # get a continuous solution across the two models
    _add_inlet_pressure_step(m_dynamic, time=4.5, value=5.5e5)
    solver.solve(m_dynamic, tee=True)

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

    # As a lazy form of initialization, solve the steady state problem before
    # turning on the controller.
    solver.solve(m_dynamic2, tee=True)

    # Now turn on control and solve again.
    m_dynamic2.fs.ctrl.activate()
    m_dynamic2.fs.valve_1.valve_opening.unfix()
    m_dynamic2.fs.valve_1.valve_opening[5].fix()
    solver.solve(m_dynamic2, tee=True)

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
        assert m_dynamic.fs.valve_1.valve_opening[t] == stitch_valve[i]
    for j, t in enumerate(m_dynamic2.fs.time):
        i = j + len(m_dynamic.fs.time)
        assert t == stitch_time[i]
        assert m_dynamic2.fs.valve_1.valve_opening[t] == stitch_valve[i]


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_pid_velocity():
    """This test is pretty course-grained, but it should cover everything"""
    tpid(PIDForm.velocity)


@pytest.mark.integration
@pytest.mark.solver
@pytest.mark.skipif(not prop_available, reason="IAPWS not available")
@pytest.mark.skipif(not solver_available, reason="Solver not available")
def test_pid_standard():
    """This test is pretty course-grained, but it should cover everything"""
    tpid(PIDForm.standard)
