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
This file contains utilities for constructing a flowsheet for the
simplest gas pipeline network, which contains a supply node, a
compressor, a pipeline, and a demand node.
"""
import pyomo.environ as pyo
import idaes.core as idaes
import pyomo.network as network
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline
from idaes.models_extra.gas_distribution.unit_models.compressor import (
    IsothermalCompressor as Compressor,
)
from idaes.models_extra.gas_distribution.unit_models.node import PipelineNode


def make_simple_model(
        dynamic=True,
        nxfe=2,
        space_method="dae.finite_difference",
        space_scheme="FORWARD",
        ntfe=40,
        horizon=20.0,
        time_method="dae.finite_difference",
        time_scheme="BACKWARD",
        ):
    m = pyo.ConcreteModel()
    default = {"dynamic": dynamic}
    if dynamic:
        default["time_set"] = [0.0, horizon]
        default["time_units"] = pyo.units.hr
    m.fs = idaes.FlowsheetBlock(default=default)
    m.fs.properties = NaturalGasParameterBlock()
    node_configs = [
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 0,
            "n_outlet_pipelines": 1,
            "n_supplies": 1,
            "n_demands": 0,
        },
        {
            "property_package": m.fs.properties,
            "n_inlet_pipelines": 1,
            "n_outlet_pipelines": 0,
            "n_supplies": 0,
            "n_demands": 1,
        },
    ]
    m.fs.node_set = pyo.Set(initialize=list(range(len(node_configs))))
    node_configs = {i: config for i, config in enumerate(node_configs)}
    m.fs.nodes = PipelineNode(m.fs.node_set, initialize=node_configs)
    nodes = m.fs.nodes
    pipeline_config = {
        "property_package": m.fs.properties,
        "finite_elements": nxfe,
        "transformation_method": space_method,
        "transformation_scheme": space_scheme,
        "has_holdup": True,
    }
    m.fs.pipeline = GasPipeline(default=pipeline_config)
    pipeline = m.fs.pipeline
    compressor_config = {"property_package": m.fs.properties}
    m.fs.compressor = Compressor(default=compressor_config)
    compressor = m.fs.compressor
    m._compressor_to_pipeline = network.Arc(
        ports=(compressor.outlet_port, pipeline.inlet_port)
    )

    m.fs.nodes[0].add_pipeline_to_outlet(compressor)
    m.fs.nodes[1].add_pipeline_to_inlet(pipeline)
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m)

    cv = m.fs.pipeline.control_volume
    m.fs.pipeline.diameter.fix(0.92*pyo.units.m)
    cv.length.fix(300*pyo.units.km)
    x0 = cv.length_domain.first()
    xf = cv.length_domain.last()
    j = next(iter(m.fs.properties.component_list))
    if dynamic:
        # Fix initial conditions
        t0 = m.fs.time.first()
        for x in cv.length_domain:
            if x != x0:
                cv.pressure[t0, x].fix()
            if x != xf:
                cv.flow_mass[t0, x].fix()
        # Apply transformation
        disc = pyo.TransformationFactory(time_method)
        disc.apply_to(m, nfe=ntfe, wrt=m.fs.time, scheme=time_scheme)

        # Deactivate constraints
        # I want to deactivate differential equations at (t0, xf).
        # Material balance already doesn't exist here.
        # TODO: This constraint deactivation should go in the unit model,
        # and which constraint we deactivate (x0 vs xf) depends on the
        # discretization.
        cv.momentum_balance[t0, xf].deactivate()

    # Fix "dynamic inputs." This needs to be done after a potential
    # discretization transformation.
    nodes[0].state[:].mole_frac_comp[j].fix()
    # Fixing node temperature specifies the temperature of the supply
    # and the pipeline inlet.
    nodes[0].state[:].temperature.fix()
    #nodes[0].state[:].pressure.fix()
    #nodes[1].state[:].flow_mol.fix()

    return m


def add_objective_to_model(
        m,
        dynamic=True,
        add_terminal_penalty=None,
        ):
    if add_terminal_penalty is None:
        add_terminal_penalty = dynamic
    m.supply_cost_coef = pyo.Param(initialize=0.0, mutable=True)
    m.electricity_cost_coef = pyo.Param(initialize=0.1, mutable=True)
    #m.demand_cost_coef = pyo.Param(initialize=1e6, mutable=True)
    # Demand cost coefficient needs to be scaled down for the units
    # of my model. The coefficient in the DAE code is 1e6 ((1e4 SCM)/hr)^-2.
    #
    # This needs some information from the property package
    demand_conv_factor = (
        1e4
        * (0.72*pyo.units.kg / pyo.units.m**3)
        / (18.0*pyo.units.kg / pyo.units.kmol)
    )
    m.demand_cost_coef = pyo.Param(
        initialize=pyo.value(1e6*demand_conv_factor**-2),
        mutable=True,
    )
    space = m.fs.pipeline.control_volume.length_domain
    x0 = space.first()
    xf = space.last()
    if dynamic:
        m.terminal_pressure_coef = pyo.Param(initialize=1e6, mutable=True)
        # Terminal flow penalty uses units of kg/hr instead of kmol/hr
        terminal_flow_conv_factor = 1e4 * (0.72*pyo.units.kg / pyo.units.m**3)
        m.terminal_flow_coef = pyo.Param(
            initialize=pyo.value(1e6*terminal_flow_conv_factor**-2),
            mutable=True,
        )
        m.terminal_pressure = pyo.Param(
            space,
            initialize=57.0,
            units=pyo.units.bar,
            mutable=True,
        )
        m.terminal_flow_mass = pyo.Param(
            space,
            initialize=16666.7,
            units=pyo.units.kmol/pyo.units.hr,
            mutable=True,
        )

    time = m.fs.time
    t0 = time.first()
    tf = time.last()
    supply = m.fs.nodes[0].supplies[0]
    demand = m.fs.nodes[1].demands[0]
    compressor = m.fs.compressor

    m.target_demand = pyo.Param(
        time,
        mutable=True,
        units=pyo.units.kmol/pyo.units.hr,
        initialize=16667.0,
    )

    if dynamic:
        m.supply_cost = pyo.Expression(expr=sum(
            m.supply_cost_coef * supply.flow_mol[t] * (t - time.prev(t))
            for t in time if t != t0
        ))

        m.boost_cost = pyo.Expression(expr=sum(
            m.electricity_cost_coef * compressor.power[t] * (t - time.prev(t))
            for t in time if t != t0
        ))

        m.demand_penalty = pyo.Expression(expr=sum(
            m.demand_cost_coef
            * (demand.flow_mol[t] - m.target_demand[t])**2
            * (t - time.prev(t))
            for t in time if t != t0
        ))

        m.terminal_pressure_cost = pyo.Expression(expr=sum(
            m.terminal_pressure_coef * (
                m.fs.pipeline.control_volume.pressure[tf, x]
                - m.terminal_pressure[x]
            )**2
            for x in space
        ))

        m.terminal_flow_mass_cost = pyo.Expression(expr=sum(
            m.terminal_flow_coef * (
                m.fs.pipeline.control_volume.flow_mass[tf, x]
                - m.terminal_flow_mass[x]
            )**2
            for x in space
        ))
    else:
        m.supply_cost = pyo.Expression(expr=sum(
            m.supply_cost_coef * supply.flow_mol[t] for t in time
        ))

        m.boost_cost = pyo.Expression(expr=sum(
            m.electricity_cost_coef * compressor.power[t] for t in time
        ))

        m.demand_penalty = pyo.Expression(expr=sum(
            m.demand_cost_coef
            * (demand.flow_mol[t] - m.target_demand[t])**2
            for t in time
        ))

    cost_expr = m.supply_cost + m.boost_cost + m.demand_penalty
    if dynamic and add_terminal_penalty:
        cost_expr += m.terminal_pressure_cost
        cost_expr += m.terminal_flow_mass_cost
    cost_expr *= 1e-6

    m.obj = pyo.Objective(expr=cost_expr)
