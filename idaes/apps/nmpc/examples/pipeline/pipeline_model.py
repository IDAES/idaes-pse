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
import idaes.core as idaes
from idaes.models_extra.gas_distribution.properties.natural_gas import (
    NaturalGasParameterBlock,
)
from idaes.models_extra.gas_distribution.unit_models.pipeline import GasPipeline

"""
Utilities for constructing and simulating a gas pipeline model
"""

def get_simulation_inputs(
        simulation_horizon=20.0,
        model_horizon=2.0,
        initial_flow=3.0e5,
        perturbed_flow=3.6e5,
        initial_pressure=57.0,
        perturbed_pressure=57.0,
        t_ptb=4.0,
        ):
    """
    Flow in kg/hr. Pressure in bar.
    """
    n_cycles = round(simulation_horizon/model_horizon)
    simulation_time = [model_horizon*i for i in range(n_cycles+1)]
    input_sequence = (
        simulation_time,
        {
            "fs.pipeline.control_volume.flow_mass[*,1.0]": [
                initial_flow if t <= t_ptb else perturbed_flow
                for t in simulation_time
                ],
            "fs.pipeline.control_volume.pressure[*,0.0]": [
                initial_pressure if t <= t_ptb else perturbed_pressure
                for t in simulation_time
                ],
        },
    )
    return input_sequence


def make_model(
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
    pipeline_config = {
        "property_package": m.fs.properties,
        "finite_elements": nxfe,
        "transformation_method": space_method,
        "transformation_scheme": space_scheme,
        "has_holdup": True,
    }
    m.fs.pipeline = GasPipeline(default=pipeline_config)
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
    cv.properties[:, x0].mole_frac_comp[j].fix(1.0)
    cv.properties[:, x0].temperature.fix(293.15*pyo.units.K)
    return m
