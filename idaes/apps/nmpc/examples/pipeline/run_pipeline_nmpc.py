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
import pyomo.environ as pyo
from pyomo.dae.flatten import flatten_dae_components

from workspace.pipelines.examples.pipeline_model import make_model

from idaes.apps.nmpc.dynamic_data import (
    find_nearest_index,
    interval_data_from_time_series,
    load_inputs_into_model,
)
from idaes.apps.nmpc.cost_expressions import (
    get_tracking_cost_from_constant_setpoint,
)

import matplotlib.pyplot as plt

"""
Script to run an NMPC simulation with a single natural gas pipeline model.
"""


def run_nmpc(
        simulation_horizon=20.0,
        controller_horizon=20.0,
        sample_period=2.0,
        ):
    """
    Runs a simulation where outlet (demand) flow rate is perturbed
    at a specified time.
    """
    #
    # Make steady model (for initial conditions) and dynamic model
    #
    nxfe = 4
    m_steady = make_model(dynamic=False, nxfe=nxfe)
    ntfe_per_sample = 4
    n_cycles = round(simulation_horizon/sample_period)
    m_plant = make_model(
        horizon=sample_period,
        ntfe=ntfe_per_sample,
        nxfe=nxfe,
    )
    time = m_plant.fs.time
    t0 = time.first()
    tf = time.last()

    ipopt = pyo.SolverFactory("ipopt")

    # Fix "inputs" in dynamic model: inlet pressure and outlet flow rate
    space = m_plant.fs.pipeline.control_volume.length_domain
    x0 = space.first()
    xf = space.last()
    m_plant.fs.pipeline.control_volume.flow_mass[:, xf].fix()
    m_plant.fs.pipeline.control_volume.pressure[:, x0].fix()

    #
    # Make steady model for setpoint
    #
    m_setpoint = make_model(dynamic=False, nxfe=nxfe)
    m_setpoint.fs.pipeline.control_volume.flow_mass[:, x0].fix(
        5.0e5*pyo.units.kg/pyo.units.hr
    )
    m_setpoint.fs.pipeline.control_volume.pressure[:, x0].fix(
        57.0*pyo.units.bar
    )
    ipopt.solve(m_setpoint, tee=True)

    #
    # Extract data from setpoint model
    #
    setpoint_scalar_vars, setpoint_dae_vars = flatten_dae_components(
        m_setpoint, m_setpoint.fs.time, pyo.Var
    )
    setpoint_data = {
        # HACK: Pass keys through ComponentUID twice to ensure consistency
        # of string representation.
        str(pyo.ComponentUID(var.referent)): var[t0].value
        for var in setpoint_dae_vars
    }

    #
    # Load initial inputs into steady model for initial conditions
    #
    m_steady.fs.pipeline.control_volume.flow_mass[:, x0].fix(
        3.0e5*pyo.units.kg/pyo.units.hr
    )
    m_steady.fs.pipeline.control_volume.pressure[:, x0].fix(
        57.0*pyo.units.bar
    )

    ipopt.solve(m_steady, tee=True)

    #
    # Extract data from steady state model
    #
    steady_scalar_vars, steady_dae_vars = flatten_dae_components(
        m_steady, m_steady.fs.time, pyo.Var
    )
    initial_data = {
        str(pyo.ComponentUID(var.referent)): var[t0].value
        for var in steady_dae_vars
    }
    scalar_data = {
        str(pyo.ComponentUID(var)): var.value for var in steady_scalar_vars
    }

    #
    # Load data into dynamic model
    #
    for name, val in scalar_data.items():
        m_plant.find_component(name).set_value(val)
    for name, val in initial_data.items():
        var = m_plant.find_component(name)
        for t in time:
            var[t].set_value(val)

    # Solve as a sanity check -- should be square with zero infeasibility
    ipopt.solve(m_plant, tee=True)

    #
    # Initialize data structure for simulation data
    #
    scalar_vars, dae_vars = flatten_dae_components(m_plant, time, pyo.Var)
    simulation_data = (
        [t0],
        {
            str(pyo.ComponentUID(var.referent)): [var[t0].value]
            for var in dae_vars
        },
    )

    #
    # Construct dynamic model for controller
    #
    samples_per_controller = round(controller_horizon/sample_period)
    ntfe_per_controller = ntfe_per_sample * samples_per_controller
    m_controller = make_model(
        horizon=controller_horizon,
        ntfe=ntfe_per_controller,
        nxfe=nxfe,
    )
    # Fix inputs at initial condition
    m_controller.fs.pipeline.control_volume.pressure[t0, x0].fix()
    m_controller.fs.pipeline.control_volume.flow_mass[t0, xf].fix()

    #
    # Construct tracking objective
    #
    cv = m_controller.fs.pipeline.control_volume
    tracking_variables = [
        pyo.Reference(cv.pressure[:, x0]),
        pyo.Reference(cv.pressure[:, xf]),
        pyo.Reference(cv.flow_mass[:, x0]),
        pyo.Reference(cv.flow_mass[:, xf]),
    ]
    weight_data = {
        "fs.pipeline.control_volume.flow_mass[*,%s]" % x0: 1e-10,
        "fs.pipeline.control_volume.flow_mass[*,%s]" % xf: 1e-10,
        "fs.pipeline.control_volume.pressure[*,%s]" % x0: 1e-2,
        "fs.pipeline.control_volume.pressure[*,%s]" % xf: 1e-2,
    }
    m_controller.tracking_cost = get_tracking_cost_from_constant_setpoint(
        tracking_variables,
        m_controller.fs.time,
        setpoint_data,
        weight_data=weight_data,
    )
    m_controller.tracking_objective = pyo.Objective(
        expr=sum(m_controller.tracking_cost.values())
    )

    #
    # Constrain inputs piecewise constant
    #
    piecewise_constant_vars = [
        pyo.Reference(cv.pressure[:, x0]),
        pyo.Reference(cv.flow_mass[:, xf]),
    ]
    m_controller.piecewise_constant_vars_set = pyo.Set(
        initialize=list(range(len(piecewise_constant_vars)))
    )
    sample_points = [
        t0 + sample_period*i for i in range(samples_per_controller+1)
    ]
    sample_point_set = set(sample_points)
    def piecewise_constant_rule(m, i, t):
        var = piecewise_constant_vars[i]
        if t in sample_point_set:
            return pyo.Constraint.Skip
        else:
            t_next = m_controller.fs.time.next(t)
            return var[t] == var[t_next]
    m_controller.piecewise_constant_constraint = pyo.Constraint(
        m_controller.piecewise_constant_vars_set,
        m_controller.fs.time,
        rule=piecewise_constant_rule,
    )

    #
    # Initialize dynamic model to initial steady state
    #
    for name, val in scalar_data.items():
        m_controller.find_component(name).set_value(val)
    for name, val in initial_data.items():
        var = m_controller.find_component(name)
        for t in time:
            var[t].set_value(val)

    #
    # Initialize data structure for controller inputs
    #
    input_names = [
        "fs.pipeline.control_volume.flow_mass[*,%s]" % xf,
        "fs.pipeline.control_volume.pressure[*,%s]" % x0,
    ]
    applied_inputs = (
        [t0],
        {
            name: [m_controller.find_component(name)[t0].value]
            for name in input_names
        },
    )

    # This will be necessary for controller initialization
    _, controller_dae_vars = flatten_dae_components(
        m_controller, m_controller.fs.time, pyo.Var
    )

    for i in range(n_cycles):
        # time.first() in the model corresponds to sim_t0 in "simulation time"
        # time.last() in the model corresponds to sim_tf in "simulation time"
        sim_t0 = i*sample_period
        sim_tf = (i+1)*sample_period

        #
        # Solve dynamic optimization problem to get inputs
        #
        ipopt.solve(m_controller, tee=True)

        ts = sample_points[1]
        #
        # Extract first inputs from controller
        #
        sample_interval = [t0, ts]
        input_data = {}
        for name in input_names:
            var = m_controller.find_component(name)
            input_data[name] = [var[t0].value, var[ts].value]
        extracted_inputs = (sample_interval, input_data)

        #
        # Extend data structure of applied inputs
        #
        applied_inputs[0].append(sim_t0 + ts)
        for name, values in applied_inputs[1].items():
            values.append(extracted_inputs[1][name][1])

        #
        # Load inputs from controller into plant
        #
        inputs_for_model = interval_data_from_time_series(extracted_inputs)
        load_inputs_into_model(m_plant, m_plant.fs.time, inputs_for_model)

        ipopt.solve(m_plant, tee=True)

        #
        # Extract time series data from solved model
        #
        # Initial conditions have already been accounted for.
        # Note that this is only correct because we're using an implicit
        # time discretization.
        non_initial_time = list(time)[1:]
        model_data = (
            non_initial_time,
            {
                str(pyo.ComponentUID(var.referent)): [
                    var[t].value for t in non_initial_time
                ] for var in dae_vars
            },
        )

        #
        # Apply offset to data from model
        #
        new_time_points = [t+sim_t0 for t in non_initial_time]
        new_sim_data = (new_time_points, dict(model_data[1]))

        #
        # Extend simulation data with result of new simulation
        #
        simulation_data[0].extend(new_time_points)
        for name, values in simulation_data[1].items():
            values.extend(new_sim_data[1][name])

        #
        # Re-initialize controller model
        #
        seen = set()
        tf = m_controller.fs.time.last()
        for var in controller_dae_vars:
            if id(var[t0]) in seen:
                continue
            else:
                seen.add(id(var[t0]))
            for t in m_controller.fs.time:
                ts = t + sample_period
                idx = m_controller.fs.time.find_nearest_index(ts)
                if idx is None:
                    # ts is outside the controller's horizon
                    var[t].set_value(var[tf].value)
                else:
                    ts = m_controller.fs.time.at(idx)
                    var[t].set_value(var[ts].value)

        #
        # Re-initialize model to final values.
        # This includes setting new initial conditions.
        #
        tf = m_plant.fs.time.last()
        for var in dae_vars:
            final_value = var[tf].value
            for t in m_plant.fs.time:
                var[t].set_value(final_value)
            controller_var = m_controller.find_component(var.referent)
            controller_var[t0].set_value(final_value)

    return simulation_data, applied_inputs


def plot_states_from_data(data, names, show=False):
    time, state_data = data
    plt.rcParams.update({"font.size": 16})
    for i, name in enumerate(names):
        values = state_data[name]
        fig, ax = plt.subplots()
        ax.plot(time, values, linewidth=3)
        ax.set_xlabel("Time (hr)")
        ax.set_title(name)
        fig.tight_layout()
        if show:
            fig.show()
        else:
            fname = "state%s.png" % i
            fig.savefig(fname, transparent=True)


def plot_inputs_from_data(data, names, show=False):
    time, input_data = data
    plt.rcParams.update({"font.size": 16})
    for i, name in enumerate(names):
        values = input_data[name]
        fig, ax = plt.subplots()
        ax.step(time, values, linewidth=3)
        ax.set_xlabel("Time (hr)")
        ax.set_title(name)
        fig.tight_layout()
        if show:
            fig.show()
        else:
            fname = "input%s.png" % i
            fig.savefig(fname, transparent=True)
        

if __name__ == "__main__":
    simulation_data, applied_inputs = run_nmpc(
        simulation_horizon=24,
    )
    plot_states_from_data(
        simulation_data,
        [
            "fs.pipeline.control_volume.flow_mass[*,0.0]",
            "fs.pipeline.control_volume.pressure[*,1.0]",
        ],
        show=True,
    )
    plot_inputs_from_data(
        applied_inputs,
        [
            "fs.pipeline.control_volume.flow_mass[*,1.0]",
            "fs.pipeline.control_volume.pressure[*,0.0]",
        ],
        show=True,
    )
