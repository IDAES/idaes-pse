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
from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory, Reference
from pyomo.dae.initialization import solve_consistent_initial_conditions
import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model
import pandas as pd
import matplotlib.pyplot as plt
from idaes.apps.caprese.data_manager import (
                    ControllerDataManager,)

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {
            'tol': 1e-6,
            'bound_push': 1e-8,
            'halt_on_ampl_error': 'yes',
            'linear_solver': 'ma57',
            }
else:
    solver = None

def main():
    # This tests the same model constructed in the test_nmpc_constructor_1 file
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
    sample_time = 0.5
    m_plant = make_model(horizon=sample_time, ntfe=5, ntcp=2)
    time_plant = m_plant.fs.time

    simulation_horizon = 60
    n_samples_to_simulate = round(simulation_horizon/sample_time)

    samples_to_simulate = [time_plant.first() + i*sample_time
                           for i in range(1, n_samples_to_simulate)]

    # We must identify for the controller which variables are our
    # inputs and measurements.
    inputs = [
            m_plant.fs.mixer.S_inlet.flow_vol[0],
            m_plant.fs.mixer.E_inlet.flow_vol[0],
            ]
    measurements = [
            m_controller.fs.cstr.outlet.conc_mol[0, 'C'],
            m_controller.fs.cstr.outlet.conc_mol[0, 'E'],
            m_controller.fs.cstr.outlet.conc_mol[0, 'S'],
            m_controller.fs.cstr.outlet.conc_mol[0, 'P'],
            m_controller.fs.cstr.outlet.temperature[0],
            m_controller.fs.cstr.volume[0],
            ]
    
    # Construct the "NMPC simulator" object
    nmpc = NMPCSim(
            plant_model=m_plant,
            plant_time_set=m_plant.fs.time,
            controller_model=m_controller, 
            controller_time_set=m_controller.fs.time,
            inputs_at_t0=inputs,
            measurements=measurements,
            sample_time=sample_time,
            )

    plant = nmpc.plant
    controller = nmpc.controller

    p_t0 = nmpc.plant.time.first()
    c_t0 = nmpc.controller.time.first()
    p_ts = nmpc.plant.sample_points[1]
    c_ts = nmpc.controller.sample_points[1]

    #--------------------------------------------------------------------------
    # Declare variables of interest for plotting.
    # It's ok not declaring anything. The data manager will still save some 
    # important data, but the user should use the default string of CUID for plotting afterward.
    states_of_interest = [nmpc.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','S'],
                          nmpc.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','E'],
                          nmpc.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','C'],
                          nmpc.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','P'],
                          nmpc.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','Solvent'],
                          nmpc.plant.mod.fs.cstr.control_volume.energy_holdup[:,'aq'],
                          ]
    ref_soi = [Reference(var) for var in states_of_interest]
    
    inputs_of_interest = [Reference(nmpc.plant.mod.fs.mixer.S_inlet_state[:].flow_vol),
                          Reference(nmpc.plant.mod.fs.mixer.E_inlet_state[:].flow_vol),]
    
    data_manager = ControllerDataManager(plant, 
                                        controller,
                                        ref_soi,
                                        inputs_of_interest,)
    #--------------------------------------------------------------------------
    
    solve_consistent_initial_conditions(plant, plant.time, solver)
    solve_consistent_initial_conditions(controller, controller.time, solver)

    # We now perform the "RTO" calculation: Find the optimal steady state
    # to achieve the following setpoint
    setpoint = [
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'P'], 0.4),
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'S'], 0.0),
            (controller.mod.fs.cstr.control_volume.energy_holdup[0, 'aq'], 300),
            (controller.mod.fs.mixer.E_inlet.flow_vol[0], 0.1),
            (controller.mod.fs.mixer.S_inlet.flow_vol[0], 2.0),
            (controller.mod.fs.cstr.volume[0], 1.0),
            ]
    setpoint_weights = [
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'P'], 1.),
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'S'], 1.),
            (controller.mod.fs.cstr.control_volume.energy_holdup[0, 'aq'], 1.),
            (controller.mod.fs.mixer.E_inlet.flow_vol[0], 1.),
            (controller.mod.fs.mixer.S_inlet.flow_vol[0], 1.),
            (controller.mod.fs.cstr.volume[0], 1.),
            ]

    # Some of the "differential variables" that have been fixed in the
    # model file are different from the measurements listed above. We
    # unfix them here so the RTO solve is not overconstrained.
    # (The RTO solve will only automatically unfix inputs and measurements.)
    nmpc.controller.mod.fs.cstr.control_volume.material_holdup[0,...].unfix()
    nmpc.controller.mod.fs.cstr.control_volume.energy_holdup[0,...].unfix()
    nmpc.controller.mod.fs.cstr.volume[0].unfix()

    nmpc.controller.add_setpoint_objective(setpoint, setpoint_weights)
    nmpc.controller.solve_setpoint(solver)

    # Now we are ready to construct the tracking NMPC problem
    tracking_weights = [
            *((v, 1.) for v in nmpc.controller.vectors.differential[:,0]),
            *((v, 1.) for v in nmpc.controller.vectors.input[:,0]),
            ]

    nmpc.controller.add_tracking_objective(tracking_weights)

    nmpc.controller.constrain_control_inputs_piecewise_constant()
    
    nmpc.controller.initialize_to_initial_conditions()
    
    
    # Solve the first control problem
    nmpc.controller.vectors.input[...].unfix()
    nmpc.controller.vectors.input[:,0].fix()
    solver.solve(nmpc.controller, tee=True)
    data_manager.save_controller_data(iteration = 0)


    # For a proper NMPC simulation, we must have noise.
    # We do this by treating inputs and measurements as Gaussian random
    # variables with the following variances (and bounds).
    cstr = nmpc.controller.mod.fs.cstr
    variance = [
            (cstr.outlet.conc_mol[0.0, 'S'], 0.2),
            (cstr.outlet.conc_mol[0.0, 'E'], 0.05),
            (cstr.outlet.conc_mol[0.0, 'C'], 0.1),
            (cstr.outlet.conc_mol[0.0, 'P'], 0.05),
            (cstr.outlet.temperature[0.0], 5.),
            (cstr.volume[0.0], 0.05),
            ]
    nmpc.controller.set_variance(variance)
    measurement_variance = [v.variance for v in controller.measurement_vars]
    measurement_noise_bounds = [
            (0.0, var[c_t0].ub) for var in controller.measurement_vars
            ]

    mx = plant.mod.fs.mixer
    variance = [
            (mx.S_inlet_state[0.0].flow_vol, 0.02),
            (mx.E_inlet_state[0.0].flow_vol, 0.001),
            ]
    nmpc.plant.set_variance(variance)
    input_variance = [v.variance for v in plant.input_vars]
    input_noise_bounds = [(0.0, var[p_t0].ub) for var in plant.input_vars]

    random.seed(246)

    # Extract inputs from controller and inject them into plant
    inputs = controller.generate_inputs_at_time(c_ts)
    plant.inject_inputs(inputs)
    
    data_manager.save_initial_plant_data()

    # This "initialization" really simulates the plant with the new inputs.
    nmpc.plant.initialize_by_solving_elements(solver)
    solver.solve(nmpc.plant)
    data_manager.save_plant_data(iteration = 0)

    for i in range(1, 11):
        print('\nENTERING NMPC LOOP ITERATION %s\n' % i)
        measured = nmpc.plant.generate_measurements_at_time(p_ts)
        nmpc.plant.advance_one_sample()
        nmpc.plant.initialize_to_initial_conditions()
        measured = apply_noise_with_bounds(
                measured,
                measurement_variance,
                random.gauss,
                measurement_noise_bounds,
                )

        nmpc.controller.advance_one_sample()
        nmpc.controller.load_measurements(measured, 
                                          target = "measurement",
                                          timepoint = controller.time.first())

        solver.solve(nmpc.controller, tee=True)
        data_manager.save_controller_data(iteration = i)

        inputs = controller.generate_inputs_at_time(c_ts)
        inputs = apply_noise_with_bounds(
                inputs,
                input_variance,
                random.gauss,
                input_noise_bounds,
                )
        plant.inject_inputs(inputs)
        
        nmpc.plant.initialize_by_solving_elements(solver)
        solver.solve(nmpc.plant)
        data_manager.save_plant_data(iteration = i)

    data_manager.plot_setpoint_tracking_results(ref_soi)
    data_manager.plot_control_input(inputs_of_interest)

    return nmpc, data_manager

if __name__ == '__main__':
    nmpc, data_manager = main()