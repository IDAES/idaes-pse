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
This script solves a dynamic optimization problem with a simple pipeline
network. The network has two nodes, one supply and one demand, with a
compressor and a pipeline between them.
The dynamic optimization problem penalizes deviation from target demand
"""
import pyomo.environ as pyo
from pyomo.dae.flatten import flatten_dae_components

from idaes.models_extra.gas_distribution.flowsheets.simple_network_model import (
    make_simple_model,
    add_objective_to_model,
)

from idaes.apps.nmpc.dynamic_data import (
    interval_data_from_time_series,
    load_inputs_into_model,
)

import matplotlib.pyplot as plt


def _plot_time_indexed_data(
        data,
        names,
        show=True,
        prefix=None,
        ):
    """ 
    data:
    (
        [t0, ...],
        {
            str(cuid): [value0, ...],
        },
    )
    names: list of str(cuids) that we will actually plot.
    """
    if prefix is None:
        prefix = ""
    plt.rcParams.update({"font.size": 16})
    time, name_map = data
    for i, name in enumerate(names):
        values = name_map[name]
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.plot(time, values, linewidth=3)

        # Logic for adding useful names and axis labels to plots we care
        # about using for other purposes.
        ax.set_xlabel("Time (hr)")
        ax.set_title(name)
        fig.tight_layout()
        if show:
            plt.show()
        else:
            fig.savefig(prefix + "state%s.png" % i, transparent=True)


def run_dynamic_optimization():
    nxfe = 4

    ipopt = pyo.SolverFactory("ipopt")

    #
    # Create steady state model for initial conditions and extract data
    #
    m = make_simple_model(nxfe=nxfe, dynamic=False)
    prop = m.fs.properties
    mw = prop.natural_gas.mw
    inlet_flow = 3.0e5 * pyo.units.kg/pyo.units.hr / mw
    t0 = m.fs.time.first()
    m.fs.nodes[0].supplies[0].flow_mol[t0].fix(inlet_flow)
    m.fs.nodes[0].state[t0].pressure.fix(50.0*pyo.units.bar)
    m.fs.compressor.boost_pressure[t0].fix(7.0*pyo.units.bar)
    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)
    scalar_vars, dae_vars = flatten_dae_components(m, m.fs.time, pyo.Var)
    initial_scalar_data = {
        str(pyo.ComponentUID(var)): var.value
        for var in scalar_vars
    }
    initial_data = {
        str(pyo.ComponentUID(var.referent)): var[t0].value
        for var in dae_vars
    }
    ###

    p_demand_min = 50.5
    #
    # Create steady model for target conditions and extract data
    #
    add_objective_to_model(m, dynamic=False)
    m.fs.nodes[0].supplies[0].flow_mol[t0].unfix()
    m.fs.compressor.boost_pressure[t0].unfix()
    m.fs.compressor.boost_pressure[t0].setlb(0.0)
    m.fs.compressor.boost_pressure[t0].setub(100.0)
    m.fs.nodes[1].state[t0].pressure.setlb(50.5)
    m.fs.nodes[1].demands[0].flow_mol[t0].unfix()
    m.target_demand[:].set_value(
        pyo.value(50.0 * 1e4 * 0.72 / mw)*pyo.units.kmol/pyo.units.hr
    )
    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)
    target_data = {
        str(pyo.ComponentUID(var.referent)): var[t0].value
        for var in dae_vars
    }
    ###

    #
    # Construct dynamic model
    #
    horizon = 20.0
    ntfe = 20
    time_scheme = "BACKWARD"
    m = make_simple_model(
        horizon=horizon,
        ntfe=ntfe,
        nxfe=nxfe,
        dynamic=True,
        time_scheme=time_scheme,
    )
    add_objective_to_model(m, dynamic=True)
    time = m.fs.time
    t0 = time.first()
    ###

    #
    # Set terminal costs from target steady state
    #
    space = m.fs.pipeline.control_volume.length_domain
    for x in space:
        m.terminal_pressure[x].set_value(
            # Ideally I would have a map to the variable in the dynamic model,
            # then a map from that variable to its terminal parameter.
            target_data["fs.pipeline.control_volume.pressure[*,%s]" % x]
        )
        m.terminal_flow_mass[x].set_value(
            target_data["fs.pipeline.control_volume.flow_mass[*,%s]" % x]
        )

    #
    # Fix degrees of freedom
    #
    m.fs.compressor.boost_pressure[:].fix()
    m.fs.nodes[0].state[:].pressure.fix()
    m.fs.nodes[1].demands[0].flow_mol[:].fix()
    ###

    #
    # Initialize dynamic model with steady state data
    #
    for name, val in initial_data.items():
        var = m.find_component(name)
        for t in m.fs.time:
            var[t].set_value(val)
    for name, val in initial_scalar_data.items():
        var = m.find_component(name)
        var.set_value(val)
    m.fs.pipeline.control_volume.material_accumulation[...].set_value(0.0)
    m.fs.pipeline.control_volume.flow_mass_dt[:,:].set_value(0.0)
    ###

    #
    # Sanity check solve - should have zero degrees of freedom and no
    # infeasibility
    #
    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)

    #
    # Construct target demand sequence and load into model
    #
    sample_points = [4.0, 20.0]
    val = pyo.value(50.0 * 1e4 * 0.72 / mw)
    demand_name = str(pyo.ComponentUID(m.target_demand))
    input_series_data = (
        sample_points,
        {
            demand_name: [val, val],
        },
    )
    # Note that here we are setting the value of a mutable parameter
    input_interval_data = interval_data_from_time_series(input_series_data)
    load_inputs_into_model(m, time, input_interval_data)
    ###

    #
    # Make sure proper degrees of freedom are unfixed
    #
    m.fs.nodes[0].state[:].pressure.fix()
    m.fs.compressor.boost_pressure[:].setub(100.0)
    for t in time:
        if t != t0:
            m.fs.compressor.boost_pressure[t].unfix()
            m.fs.nodes[1].demands[0].flow_mol[t].unfix()

            m.fs.nodes[1].state[t].pressure.setlb(p_demand_min)

    m.fs.pipeline.control_volume.area.fix()
    m.fs.pipeline.diameter_eqn.deactivate()

    res = ipopt.solve(m, tee=True)
    pyo.assert_optimal_termination(res)

    #
    # Extract values we care about and plot
    #
    scalar_vars, dae_vars = flatten_dae_components(m, time, pyo.Var)
    sim_data = (
        list(time),
        {
            str(pyo.ComponentUID(var.referent)): [var[t].value for t in time]
            for var in dae_vars
        },
    )
    return sim_data


if __name__ == "__main__":
    sim_data = run_dynamic_optimization()
    names_to_plot = [
        "fs.nodes[0].state[*].flow_mol",
        "fs.nodes[1].state[*].flow_mol",
        "fs.compressor.boost_pressure[*]",
        "fs.nodes[1].state[*].pressure",
    ]
    _plot_time_indexed_data(sim_data, names_to_plot, show=True)
