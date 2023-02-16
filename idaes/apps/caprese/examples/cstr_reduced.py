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
Example for Caprese's module for NMPC.
"""
import random
from idaes.apps.caprese.nmpc import NMPCSim
from idaes.apps.caprese.dynamic_block import DynamicBlock
from idaes.apps.caprese.controller import ControllerBlock
from idaes.apps.caprese.util import apply_noise_with_bounds
from idaes.apps.caprese.categorize import (
    categorize_dae_variables_and_constraints,
    VariableCategory,
    ConstraintCategory,
)

VC = VariableCategory
CC = ConstraintCategory

import pyomo.environ as pyo
from pyomo.dae.flatten import flatten_dae_components
from pyomo.dae.initialization import solve_consistent_initial_conditions
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.incidence_analysis.interface import IncidenceGraphInterface
from pyomo.core.expr.calculus.derivatives import reverse_ad

import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model

import numpy as np
import scipy.sparse as sps
import pandas as pd
import matplotlib.pyplot as plt

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if pyo.SolverFactory("ipopt").available():
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
        "tol": 1e-6,
        "bound_push": 1e-8,
        "halt_on_ampl_error": "yes",
        "linear_solver": "ma57",
    }
else:
    solver = None


class PlotData(object):
    def __init__(self, group, location, name=None, t_switch=None):
        # Would really like a PlotData class that is constructed based on an
        # NMPCVar object that contains necessary setpoint/reference
        # information, instead of having to access that in the NMPCVarGroup
        time = group.index_set
        if t_switch == None:
            t_switch = group.t0

        self.name = name

        var = group.varlist[location]
        initial = group.reference[location]
        setpoint = group.setpoint[location]
        self.data_series = pd.Series(
            [var[t].value for t in time], index=[t for t in time]
        )
        self.setpoint_series = pd.Series(
            [initial if t < t_switch else setpoint for t in time]
        )

    def plot(self):
        # fig, ax can be formatted to the user's liking
        fig, ax = plt.subplots()
        if self.name is not None:
            self.data_series.plot(label=self.name)
        else:
            self.data_series.plot()
        return fig, ax


def main(plot_switch=False):

    # This tests the same model constructed in the test_nmpc_constructor_1 file
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
    sample_time = 0.5
    m_plant = make_model(horizon=sample_time, ntfe=5, ntcp=2)
    time_plant = m_plant.fs.time

    solve_consistent_initial_conditions(m_plant, time_plant, solver)

    #####
    # Flatten and categorize controller model
    #####
    model = m_controller
    time = model.fs.time
    t0 = time.first()
    t1 = time[2]
    scalar_vars, dae_vars = flatten_dae_components(
        model,
        time,
        pyo.Var,
    )
    scalar_cons, dae_cons = flatten_dae_components(
        model,
        time,
        pyo.Constraint,
    )
    inputs = [
        model.fs.mixer.S_inlet.flow_vol,
        model.fs.mixer.E_inlet.flow_vol,
    ]
    measurements = [
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "C"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "E"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "S"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "P"]),
        model.fs.cstr.outlet.temperature,
    ]
    model.fs.cstr.control_volume.material_holdup[:, "aq", "Solvent"].fix()
    model.fs.cstr.total_flow_balance.deactivate()
    var_partition, con_partition = categorize_dae_variables_and_constraints(
        model,
        dae_vars,
        dae_cons,
        time,
        input_vars=inputs,
    )
    controller = ControllerBlock(
        model=model,
        time=time,
        measurements=measurements,
        category_dict={None: var_partition},
    )
    controller.construct()

    solve_consistent_initial_conditions(m_controller, time, solver)
    controller.initialize_to_initial_conditions()

    m_controller._dummy_obj = pyo.Objective(expr=0)
    nlp = PyomoNLP(m_controller)
    igraph = IncidenceGraphInterface(nlp)
    m_controller.del_component(m_controller._dummy_obj)
    diff_vars = [var[t1] for var in var_partition[VC.DIFFERENTIAL]]
    alg_vars = [var[t1] for var in var_partition[VC.ALGEBRAIC]]
    deriv_vars = [var[t1] for var in var_partition[VC.DERIVATIVE]]
    diff_eqns = [con[t1] for con in con_partition[CC.DIFFERENTIAL]]
    alg_eqns = [con[t1] for con in con_partition[CC.ALGEBRAIC]]

    # Assemble and factorize "derivative Jacobian"
    dfdz = nlp.extract_submatrix_jacobian(diff_vars, diff_eqns)
    dfdy = nlp.extract_submatrix_jacobian(alg_vars, diff_eqns)
    dgdz = nlp.extract_submatrix_jacobian(diff_vars, alg_eqns)
    dgdy = nlp.extract_submatrix_jacobian(alg_vars, alg_eqns)
    dfdzdot = nlp.extract_submatrix_jacobian(deriv_vars, diff_eqns)
    fact = sps.linalg.splu(dgdy.tocsc())
    dydz = fact.solve(dgdz.toarray())
    deriv_jac = dfdz - dfdy.dot(dydz)
    fact = sps.linalg.splu(dfdzdot.tocsc())
    dzdotdz = -fact.solve(deriv_jac)

    # Use some heuristic on the eigenvalues of the derivative Jacobian
    # to identify fast states.
    w, V = np.linalg.eig(dzdotdz)
    w_max = np.max(np.abs(w))
    (fast_modes,) = np.where(np.abs(w) > w_max / 2)
    fast_states = []
    for idx in fast_modes:
        evec = V[:, idx]
        _fast_states, _ = np.where(np.abs(evec) > 0.5)
        fast_states.extend(_fast_states)
    fast_states = set(fast_states)

    # Store components necessary for model reduction in a model-
    # independent form.
    fast_state_derivs = [
        pyo.ComponentUID(var_partition[VC.DERIVATIVE][idx].referent, context=model)
        for idx in fast_states
    ]
    fast_state_diffs = [
        pyo.ComponentUID(var_partition[VC.DIFFERENTIAL][idx].referent, context=model)
        for idx in fast_states
    ]
    fast_state_discs = [
        pyo.ComponentUID(con_partition[CC.DISCRETIZATION][idx].referent, context=model)
        for idx in fast_states
    ]

    # Perform pseudo-steady state model reduction on the fast states
    # and re-categorize
    for cuid in fast_state_derivs:
        var = cuid.find_component_on(m_controller)
        var.fix(0.0)
    for cuid in fast_state_diffs:
        var = cuid.find_component_on(m_controller)
        var[t0].unfix()
    for cuid in fast_state_discs:
        con = cuid.find_component_on(m_controller)
        con.deactivate()

    var_partition, con_partition = categorize_dae_variables_and_constraints(
        model,
        dae_vars,
        dae_cons,
        time,
        input_vars=inputs,
    )
    controller.del_component(model)

    # Re-construct controller block with new categorization
    measurements = [
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "C"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "E"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "S"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "P"]),
    ]
    controller = ControllerBlock(
        model=model,
        time=time,
        measurements=measurements,
        category_dict={None: var_partition},
    )
    controller.construct()

    #####
    # Construct dynamic block for plant
    #####
    model = m_plant
    time = model.fs.time
    t0 = time.first()
    t1 = time[2]
    scalar_vars, dae_vars = flatten_dae_components(
        model,
        time,
        pyo.Var,
    )
    scalar_cons, dae_cons = flatten_dae_components(
        model,
        time,
        pyo.Constraint,
    )
    inputs = [
        model.fs.mixer.S_inlet.flow_vol,
        model.fs.mixer.E_inlet.flow_vol,
    ]
    measurements = [
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "C"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "E"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "S"]),
        pyo.Reference(model.fs.cstr.outlet.conc_mol[:, "P"]),
    ]
    model.fs.cstr.control_volume.material_holdup[:, "aq", "Solvent"].fix()
    model.fs.cstr.total_flow_balance.deactivate()

    var_partition, con_partition = categorize_dae_variables_and_constraints(
        model,
        dae_vars,
        dae_cons,
        time,
        input_vars=inputs,
    )
    plant = DynamicBlock(
        model=model,
        time=time,
        measurements=measurements,
        category_dict={None: var_partition},
    )
    plant.construct()

    p_t0 = plant.time.first()
    c_t0 = controller.time.first()
    p_ts = plant.sample_points[1]
    c_ts = controller.sample_points[1]

    controller.set_sample_time(sample_time)
    plant.set_sample_time(sample_time)

    # We now perform the "RTO" calculation: Find the optimal steady state
    # to achieve the following setpoint
    setpoint = [
        (controller.mod.fs.cstr.outlet.conc_mol[0, "P"], 0.4),
        # (controller.mod.fs.cstr.outlet.conc_mol[0, 'S'], 0.01),
        (controller.mod.fs.cstr.outlet.conc_mol[0, "S"], 0.1),
        (controller.mod.fs.cstr.control_volume.energy_holdup[0, "aq"], 300),
        (controller.mod.fs.mixer.E_inlet.flow_vol[0], 0.1),
        (controller.mod.fs.mixer.S_inlet.flow_vol[0], 2.0),
        (controller.mod.fs.cstr.volume[0], 1.0),
    ]
    setpoint_weights = [
        (controller.mod.fs.cstr.outlet.conc_mol[0, "P"], 1.0),
        (controller.mod.fs.cstr.outlet.conc_mol[0, "S"], 1.0),
        (controller.mod.fs.cstr.control_volume.energy_holdup[0, "aq"], 1.0),
        (controller.mod.fs.mixer.E_inlet.flow_vol[0], 1.0),
        (controller.mod.fs.mixer.S_inlet.flow_vol[0], 1.0),
        (controller.mod.fs.cstr.volume[0], 1.0),
    ]

    # Some of the "differential variables" that have been fixed in the
    # model file are different from the measurements listed above. We
    # unfix them here so the RTO solve is not overconstrained.
    # (The RTO solve will only automatically unfix inputs and measurements.)
    controller.mod.fs.cstr.control_volume.material_holdup[0, ...].unfix()
    controller.mod.fs.cstr.control_volume.energy_holdup[0, ...].unfix()
    # controller.mod.fs.cstr.volume[0].unfix()
    controller.mod.fs.cstr.control_volume.material_holdup[0, "aq", "Solvent"].fix()

    controller.add_setpoint_objective(setpoint, setpoint_weights)
    controller.solve_setpoint(solver)

    # Now we are ready to construct the tracking NMPC problem
    tracking_weights = [
        *((v, 1.0) for v in controller.vectors.differential[:, 0]),
        *((v, 1.0) for v in controller.vectors.input[:, 0]),
    ]

    controller.add_tracking_objective(tracking_weights)

    controller.constrain_control_inputs_piecewise_constant()

    controller.initialize_to_initial_conditions()

    # Solve the first control problem
    controller.vectors.input[...].unfix()
    controller.vectors.input[:, 0].fix()
    solver.solve(controller, tee=True)

    # For a proper NMPC simulation, we must have noise.
    # We do this by treating inputs and measurements as Gaussian random
    # variables with the following variances (and bounds).
    cstr = controller.mod.fs.cstr
    variance = [
        (cstr.outlet.conc_mol[0.0, "S"], 0.01),
        (cstr.outlet.conc_mol[0.0, "E"], 0.005),
        (cstr.outlet.conc_mol[0.0, "C"], 0.01),
        (cstr.outlet.conc_mol[0.0, "P"], 0.005),
        (cstr.outlet.temperature[0.0], 1.0),
        (cstr.volume[0.0], 0.05),
    ]
    controller.set_variance(variance)
    measurement_variance = [v.variance for v in controller.MEASUREMENT_BLOCK[:].var]
    measurement_noise_bounds = [
        (0.0, var[c_t0].ub) for var in controller.MEASUREMENT_BLOCK[:].var
    ]

    mx = plant.mod.fs.mixer
    variance = [
        (mx.S_inlet_state[0.0].flow_vol, 0.02),
        (mx.E_inlet_state[0.0].flow_vol, 0.001),
    ]
    plant.set_variance(variance)
    input_variance = [v.variance for v in plant.INPUT_BLOCK[:].var]
    input_noise_bounds = [(0.0, var[p_t0].ub) for var in plant.INPUT_BLOCK[:].var]

    random.seed(100)

    # Extract inputs from controller and inject them into plant
    inputs = controller.generate_inputs_at_time(c_ts)
    plant.inject_inputs(inputs)

    # This "initialization" really simulates the plant with the new inputs.
    plant.vectors.input[:, :].fix()
    plant.initialize_by_solving_elements(solver)
    plant.vectors.input[:, :].fix()
    solver.solve(plant, tee=True)

    for i in range(1, 11):
        print("\nENTERING NMPC LOOP ITERATION %s\n" % i)
        measured = plant.generate_measurements_at_time(p_ts)
        plant.advance_one_sample()
        plant.initialize_to_initial_conditions()
        measured = apply_noise_with_bounds(
            measured,
            measurement_variance,
            random.gauss,
            measurement_noise_bounds,
        )

        controller.advance_one_sample()
        controller.load_measurements(measured)

        solver.solve(controller, tee=True)

        inputs = controller.generate_inputs_at_time(c_ts)
        inputs = apply_noise_with_bounds(
            inputs,
            input_variance,
            random.gauss,
            input_noise_bounds,
        )
        plant.inject_inputs(inputs)

        plant.initialize_by_solving_elements(solver)
        solver.solve(plant)

    import pdb

    pdb.set_trace()


if __name__ == "__main__":
    main()
