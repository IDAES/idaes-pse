##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
Example for Caprese's module for NMPC.
Main script for running the example.
"""
import random
from idaes.apps.caprese.dynamic_builder import DynamicSim
from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory, Reference
from pyomo.dae.initialization import solve_consistent_initial_conditions
import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_rodrigo.cstr_rodrigo_model import make_model
from idaes.apps.caprese.data_manager import PlantDataManager
from idaes.apps.caprese.data_manager import ControllerDataManager
from idaes.apps.caprese.plotlibrary import (
        plot_setpoint_tracking_results,
        plot_control_input)

__author__ = "Kuan-Han Lin"



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
    m_controller = make_model(horizon=10, ntfe=5, ntcp=2, bounds=True)
    sample_time = 2.
    m_plant = make_model(horizon=sample_time, ntfe=2, ntcp=2, bounds = True)
    time_plant = m_plant.t

    simulation_horizon = 20
    n_samples_to_simulate = round(simulation_horizon/sample_time)

    samples_to_simulate = [time_plant.first() + i*sample_time
                           for i in range(1, n_samples_to_simulate)]

    # We must identify for the controller which variables are our
    # inputs and measurements.
    inputs = [
            m_plant.Tjinb[0],
            ]
    measurements = [
            m_plant.Tall[0, "T"],
            m_plant.Tall[0, "Tj"],
            m_plant.Ca[0],
            ]

    # Construct the "NMPC simulator" object
    nmpc = DynamicSim(
            plant_model=m_plant,
            plant_time_set=m_plant.t,
            controller_model=m_controller,
            controller_time_set=m_controller.t,
            inputs_at_t0=inputs,
            measurements_at_t0=measurements,
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
    # important data.
    states_of_interest = [Reference(nmpc.plant.mod.Ca[:]),
                          Reference(nmpc.plant.mod.Tall[:, "T"])]

    plant_data= PlantDataManager(plant, states_of_interest)
    controller_data= ControllerDataManager(controller, states_of_interest)
    #--------------------------------------------------------------------------

    solve_consistent_initial_conditions(plant, plant.time, solver)
    solve_consistent_initial_conditions(controller, controller.time, solver)

    # We now perform the "RTO" calculation: Find the optimal steady state
    # to achieve the following setpoint
    setpoint = [(controller.mod.Ca[0], 0.018)]
    setpoint_weights = [(controller.mod.Ca[0], 1.)]

    # nmpc.controller.add_setpoint_objective(setpoint, setpoint_weights)
    # nmpc.controller.solve_setpoint(solver)
    nmpc.controller.add_single_time_optimization_objective(setpoint,
                                                           setpoint_weights)
    nmpc.controller.solve_single_time_optimization(solver,
                                                   ic_type = "measurement_var",
                                                   require_steady = True,
                                                   load_setpoints = True)

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
    controller_data.save_controller_data(iteration = 0)

    #-------------------------------------------------------------------------
    #noise for measurements
    variance = [
        (nmpc.controller.mod.Tall[0, "T"], 0.05),
        (nmpc.controller.mod.Tall[0, "Tj"], 0.02),
        (nmpc.controller.mod.Ca[0], 1.0E-5),
        ]
    nmpc.controller.set_variance(variance)
    measurement_variance = [
            v.variance for v in controller.MEASUREMENT_BLOCK[:].var
            ]
    measurement_noise_bounds = [
            (var[c_t0].lb, var[c_t0].ub)
            for var in controller.MEASUREMENT_BLOCK[:].var
            ]

    # noise for inputs
    variance = [
        (plant.mod.Tjinb[0], 0.01),
        ]
    nmpc.plant.set_variance(variance)
    input_variance = [v.variance for v in plant.INPUT_BLOCK[:].var]
    input_noise_bounds = [
        (var[p_t0].lb, var[p_t0].ub) for var in plant.INPUT_BLOCK[:].var
        ]

    random.seed(246)
    #-------------------------------------------------------------------------

    # Extract inputs from controller and inject them into plant
    inputs = controller.generate_inputs_at_time(c_ts)
    plant.inject_inputs(inputs)

    plant_data.save_initial_plant_data()

    # This "initialization" really simulates the plant with the new inputs.
    nmpc.plant.initialize_by_solving_elements(solver)
    nmpc.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(nmpc.plant, tee = True)
    plant_data.save_plant_data(iteration = 0)

    for i in range(1, n_samples_to_simulate +1):
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
        nmpc.controller.load_initial_conditions(measured)

        solver.solve(nmpc.controller, tee=True)
        controller_data.save_controller_data(iteration = i)

        inputs = controller.generate_inputs_at_time(c_ts)
        inputs = apply_noise_with_bounds(
                inputs,
                input_variance,
                random.gauss,
                input_noise_bounds,
                )
        plant.inject_inputs(inputs)

        nmpc.plant.initialize_by_solving_elements(solver)
        nmpc.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(nmpc.plant, tee = True)
        plant_data.save_plant_data(iteration = i)

    plot_setpoint_tracking_results(states_of_interest,
                                   plant_data.plant_df,
                                   controller_data.setpoint_df)

    inputs_to_plot = [Reference(nmpc.plant.mod.Tjinb[:])]
    plot_control_input(inputs_to_plot, plant_data.plant_df)

    return nmpc, plant_data, controller_data

if __name__ == '__main__':
    nmpc, plant_data, controller_data = main()
