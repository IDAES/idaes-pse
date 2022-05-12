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

import idaes.core as idaes

from idaes.apps.nmpc.examples.pipeline.pipeline_model import (
    make_model,
    get_simulation_inputs,
)

from idaes.apps.nmpc.dynamic_data import (
    find_nearest_index,
    interval_data_from_time_series,
    load_inputs_into_model,
)

"""
Script for running "rolling-horizon" type simulation with pipeline model
"""

def run_simulation(
        simulation_horizon=20.0,
        t_ptb=4.0,
        ):
    """
    Runs a simulation where outlet (demand) flow rate is perturbed
    at a specified time.
    """
    #
    # Make steady model (for initial conditions) and dynamic model
    #
    m_steady = make_model(dynamic=False)
    model_horizon = 2.0
    ntfe = 4
    m = make_model(horizon=model_horizon, ntfe=ntfe)
    time = m.fs.time
    t0 = time.first()
    tf = time.last()

    # Fix "inputs" in dynamic model: inlet pressure and outlet flow rate
    space = m.fs.pipeline.control_volume.length_domain
    x0 = space.first()
    xf = space.last()
    m.fs.pipeline.control_volume.flow_mass[:, xf].fix()
    m.fs.pipeline.control_volume.pressure[:, x0].fix()

    input_sequence = get_simulation_inputs(
        simulation_horizon=simulation_horizon,
        t_ptb=t_ptb,
    )
    n_cycles = len(input_sequence[0])-1
    simulation_horizon = input_sequence[0][-1]

    #
    # Load initial inputs into steady model
    #
    initial_inputs = {
        name: values[0] for name, values in input_sequence[1].items()
    }
    t0_steady = m_steady.fs.time.first()
    for name, val in initial_inputs.items():
        m_steady.find_component(name)[t0_steady].fix(val)

    ipopt = pyo.SolverFactory("ipopt")
    ipopt.solve(m_steady, tee=True)

    #
    # Extract data from steady state model
    #
    steady_scalar_vars, steady_dae_vars = flatten_dae_components(
        m_steady, m_steady.fs.time, pyo.Var
    )
    initial_data = {
        str(pyo.ComponentUID(var.referent)): var[t0_steady].value
        for var in steady_dae_vars
    }
    scalar_data = {
        str(pyo.ComponentUID(var)): var.value for var in steady_scalar_vars
    }

    #
    # Load data into dynamic model
    #
    for name, val in scalar_data.items():
        m.find_component(name).set_value(val)
    for name, val in initial_data.items():
        var = m.find_component(name)
        for t in time:
            var[t].set_value(val)

    # Solve as a sanity check -- should be square with zero infeasibility
    ipopt.solve(m, tee=True)

    #
    # Initialize data structure for simulation data
    #
    scalar_vars, dae_vars = flatten_dae_components(m, time, pyo.Var)
    simulation_data = (
        [t0],
        {
            str(pyo.ComponentUID(var.referent)): [var[t0].value]
            for var in dae_vars
        },
    )

    simulation_time = input_sequence[0]
    for i in range(n_cycles):
        # time.first() in the model corresponds to sim_t0 in "simulation time"
        # time.last() in the model corresponds to sim_tf in "simulation time"
        sim_t0 = i*model_horizon
        sim_tf = (i+1)*model_horizon

        #
        # Extract inputs of sequence that are between sim_t0 and sim_tf
        #
        idx_t0 = find_nearest_index(simulation_time, sim_t0)
        idx_tf = find_nearest_index(simulation_time, sim_tf)
        extracted_inputs = (
            simulation_time[idx_t0:idx_tf + 1],
            {
                name: values[idx_t0:idx_tf + 1]
                for name, values in input_sequence[1].items()
            },
        )
        extracted_inputs = interval_data_from_time_series(extracted_inputs)
        
        #
        # Apply offset to time points so they are valid for model
        #
        offset = sim_t0
        inputs_for_model = {
            name: {
                (interval[0]-offset, interval[1]-offset): val
                for interval, val in inputs.items()
            } for name, inputs in extracted_inputs.items()
        }

        load_inputs_into_model(m, time, inputs_for_model)

        ipopt.solve(m, tee=True)

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
        new_time_points = [t+offset for t in non_initial_time]
        new_sim_data = (new_time_points, dict(model_data[1]))

        #
        # Extend simulation data with result of new simulation
        #
        simulation_data[0].extend(new_time_points)
        for name, values in simulation_data[1].items():
            values.extend(new_sim_data[1][name])

        #
        # Re-initialize model to final values.
        # This includes setting new initial conditions.
        #
        for var in dae_vars:
            for t in time:
                var[t].set_value(var[tf].value)

    return simulation_data


if __name__ == "__main__":
    run_simulation()
