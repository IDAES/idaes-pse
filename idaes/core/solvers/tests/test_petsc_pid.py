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
"""This uses the PID controller model to do a little more testing of the
PETSc TS solver.  Specifically it tests constraints with an explicit time var."""

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
from idaes.core.solvers import petsc
from idaes.core.util.math import smooth_max, smooth_min
import numpy as np


def _valve_pressure_flow_cb(b):
    """
    Callback for valve pressure-flow relation F = Cv*sqrt(Pi**2 - Po**2)*f(x)
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


def create_model(
    time_set=None,
    time_units=pyo.units.s,
    nfe=5,
    tee=False,
    calc_integ=True,
):
    """Create a test model and solver

    Args:
        time_set (list): The beginning and end point of the time domain
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
    # Add a controller
    m.fs.ctrl = PIDController(
        process_var=m.fs.tank.control_volume.properties_out[:].pressure,
        manipulated_var=m.fs.valve_1.valve_opening,
        calculate_initial_integral=calc_integ,
        mv_bound_type=ControllerMVBoundType.SMOOTH_BOUND,
        type=ControllerType.PI,
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

    # Fix controller settings
    m.fs.ctrl.gain_p.fix(1e-6)
    m.fs.ctrl.gain_i.fix(1e-5)
    # m.fs.ctrl.gain_d.fix(1e-6)
    # m.fs.ctrl.derivative_of_error[m.fs.time.first()].fix(0)
    m.fs.ctrl.setpoint.fix(3e5)
    m.fs.ctrl.mv_ref.fix(0)
    m.fs.ctrl.mv_lb = 0.0
    m.fs.ctrl.mv_ub = 1.0

    for t in m.fs.time:
        m.fs.valve_1.inlet.flow_mol[t] = 100  # initial guess on flow
    # simple initialize
    m.fs.valve_1.initialize()
    propagate_state(m.fs.v1_to_tank)
    m.fs.tank.initialize()
    propagate_state(m.fs.tank_to_v2)
    # Can't specify both flow and outlet pressure so free the outlet pressure
    # for initialization and refix it after.  Inlet flow gets fixed in init, but
    # is unfixed for the final problem
    m.fs.valve_2.outlet.pressure.unfix()
    m.fs.valve_2.initialize()
    m.fs.valve_2.outlet.pressure.fix(101325)
    m.fs.valve_1.valve_opening.unfix()
    m.fs.valve_1.valve_opening[m.fs.time.first()].fix()
    # Return the model and solver
    return m


@pytest.mark.integration
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_petsc_with_pid_model():
    m = create_model(
        time_set=[0, 24],
        nfe=5,
        calc_integ=True,
    )
    # time_var will be an explicit time variable we can use in constraints.
    m.fs.time_var = pyo.Var(m.fs.time)

    # We'll add a constraint to calculate the inlet pressure based on time,
    # so we need to unfix pressure.
    m.fs.valve_1.control_volume.properties_in[:].pressure.unfix()

    # The solver will directly set the time variable for the DAE solve, but
    # solving the initial conditions is just a system of nonlinear equations,
    # so we need to fix the initial time.
    m.fs.time_var[0].fix(m.fs.time.first())

    # We could break up the time domain and solve this in pieces, but creative use
    # of min and max will let us create the ramping function we want.
    # From 10s to 12s ramp inlet pressure from 500,000 Pa to 600,000 Pa
    @m.fs.Constraint(m.fs.time)
    def inlet_pressure_eqn(b, t):
        return b.valve_1.control_volume.properties_in[t].pressure == smooth_min(
            600000, smooth_max(500000, 50000 * (b.time_var[t] - 10) + 500000)
        )

    res = petsc.petsc_dae_by_time_element(
        m,
        time=m.fs.time,
        between=[m.fs.time.first(), m.fs.time.at(3), m.fs.time.last()],
        timevar=m.fs.time_var,
        ts_options={
            "--ts_type": "beuler",
            "--ts_dt": 0.1,
            "--ts_monitor": "",  # set initial step to 0.1
            "--ts_save_trajectory": 1,
        },
    )
    assert isinstance(res.results, list)

    # read the trajectory data, and make it easy by interpolating a time point
    # every second
    tj = res.trajectory
    tj2 = tj.interpolate(np.linspace(0, 24, 25))

    # For more details about the problem behavior see the PETSc examples in the
    # examples repo.

    # make sure the inlet pressure is initially 5e5 Pa
    assert pyo.value(
        tj2.get_vec(m.fs.valve_1.control_volume.properties_in[24].pressure)[5]
    ) == pytest.approx(5e5)
    # make sure the inlet pressure ramped up to 6e5 Pa
    assert pyo.value(
        tj2.get_vec(m.fs.valve_1.control_volume.properties_in[24].pressure)[20]
    ) == pytest.approx(6e5)
    # make sure after the controller comes on the presure goes to the set point
    assert pyo.value(
        tj2.get_vec(m.fs.tank.control_volume.properties_out[24].pressure)[9]
    ) == pytest.approx(3e5)
    # make sure after ramping inlet pressure the tank pressure gets back to the
    # setpoint
    assert pyo.value(
        tj2.get_vec(m.fs.tank.control_volume.properties_out[24].pressure)[22]
    ) == pytest.approx(3e5)

    # Test derivatives.  There is no discretization equation at t=0 and
    # there is no liquid, so check the calculations of the vapor energy and
    # material in the tank
    petsc.calculate_time_derivatives(m, m.fs.time)

    der = (
        m.fs.tank.control_volume.material_holdup[m.fs.time.last(), "Vap", "H2O"]
        - m.fs.tank.control_volume.material_holdup[m.fs.time.at(5), "Vap", "H2O"]
    ) / (m.fs.time.last() - m.fs.time.at(5))
    assert pyo.value(
        m.fs.tank.control_volume.material_accumulation[m.fs.time.last(), "Vap", "H2O"]
    ) == pytest.approx(pyo.value(der))
    der = (
        m.fs.tank.control_volume.energy_holdup[m.fs.time.last(), "Vap"]
        - m.fs.tank.control_volume.energy_holdup[m.fs.time.at(5), "Vap"]
    ) / (m.fs.time.last() - m.fs.time.at(5))
    assert pyo.value(
        m.fs.tank.control_volume.energy_accumulation[m.fs.time.last(), "Vap"]
    ) == pytest.approx(pyo.value(der))
    der = (
        m.fs.ctrl.integral_of_error[m.fs.time.last()]
        - m.fs.ctrl.integral_of_error[m.fs.time.at(5)]
    ) / (m.fs.time.last() - m.fs.time.at(5))
    assert pyo.value(
        m.fs.ctrl.integral_of_error_dot[m.fs.time.last()]
    ) == pytest.approx(pyo.value(der))
