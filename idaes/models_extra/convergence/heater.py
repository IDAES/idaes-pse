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

__author__ = "John Eslick"

import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.models_extra.power_generation.unit_models.helm import HelmValve
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.models.unit_models import Heater
from idaes.models.properties import iapws95
from idaes.core.util.initialization import propagate_state as _set_port
from idaes.core.solvers import get_solver


def create_model_steady_state(f=100, p=5e5, h=5e4):
    """Create a steady state heater model."""
    m = pyo.ConcreteModel(name="Dynamic Heater Test")
    m.fs = FlowsheetBlock(dynamic=False)
    # Create a property parameter block
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        phase_presentation=iapws95.PhaseType.MIX
    )
    m.fs.heater = Heater(
        has_holdup=False,
        has_pressure_change=True,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )

    m.fs.heater.inlet.enth_mol.fix(50000)
    m.fs.heater.inlet.pressure.fix(5e5)
    m.fs.heater.inlet.flow_mol.fix(100)
    m.fs.heater.heat_duty.fix(0)
    m.fs.heater.deltaP.fix(0)

    m.fs.heater.initialize()

    solver = get_solver(options={"max_iter": 25})
    return m, solver


def create_model_dynamic(
    p_in=5e5,
    p_out=101325,
    h=5e4,
    time_set=None,
    nfe=10,
    dae_transform="dae.collocation",
    dae_scheme="LAGRANGE-RADAU",
):
    """Create a test dynamic heater model and solver.  For dynamic testing
    a valve is added to the heater outlet to set up pressure driven flow.

    Args:
        time_set (list): The begining and end point of the time domain
        nfe (int): Number of finite elements argument for the DAE transformation
    Returns:
        (tuple): (ConcreteModel, Solver)
    """
    if time_set is None:
        time_set = [0, 5]

    m = pyo.ConcreteModel(name="Dynamic Heater Test")
    m.fs = FlowsheetBlock(dynamic=True, time_set=time_set)
    # Create a property parameter block
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        phase_presentation=iapws95.PhaseType.MIX
    )
    # Create the valve and heater models
    m.fs.pipe = HelmValve(
        dynamic=False,
        has_holdup=False,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.heater = Heater(
        has_holdup=True,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )
    m.fs.valve = HelmValve(
        dynamic=False,
        has_holdup=False,
        material_balance_type=MaterialBalanceType.componentTotal,
        property_package=m.fs.prop_water,
    )

    # Connect the models
    m.fs.v1_to_t = Arc(source=m.fs.pipe.outlet, destination=m.fs.heater.inlet)
    m.fs.t_to_v2 = Arc(source=m.fs.heater.outlet, destination=m.fs.valve.inlet)

    # Add the stream constraints and do the DAE transformation
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    pyo.TransformationFactory("dae.finite_difference").apply_to(
        m.fs, nfe=nfe, wrt=m.fs.time, scheme="BACKWARD"
    )

    # Fix the derivative variables to zero at time 0 (steady state assumption)
    m.fs.fix_initial_conditions()

    # A heater pressure reference that's directly time-indexed
    m.fs.heater_pressure = pyo.Reference(
        m.fs.heater.control_volume.properties_out[:].pressure
    )

    # Fix the input variables
    m.fs.pipe.inlet.enth_mol.fix(h)
    m.fs.pipe.inlet.pressure.fix(p_in)
    m.fs.valve.outlet.pressure.fix(p_out)
    m.fs.pipe.Cv.fix(0.001)
    m.fs.valve.Cv.fix(0.001)
    m.fs.pipe.valve_opening.fix(1)
    m.fs.valve.valve_opening.fix(1)
    m.fs.heater.heat_duty.fix(0)
    m.fs.heater.control_volume.volume.fix(2.0)

    # Initialize the model
    solver = get_solver(options={"max_iter": 25})
    for t in m.fs.time:
        m.fs.pipe.inlet.flow_mol = 250  # initial guess on flow
    # simple initialize
    m.fs.pipe.initialize()
    _set_port(m.fs.heater.inlet, m.fs.pipe.outlet)
    m.fs.heater.initialize()
    _set_port(m.fs.valve.inlet, m.fs.heater.outlet)
    m.fs.valve.initialize()
    solver.solve(m, tee=True)

    # Return the model and solver
    return m, solver
