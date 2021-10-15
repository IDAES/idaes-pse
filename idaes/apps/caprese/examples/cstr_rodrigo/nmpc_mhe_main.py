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
Example for Caprese's module for NMPC/MHE.
"""
import random
from idaes.apps.caprese.dynamic_builder import DynamicSim
from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory, Reference
from pyomo.dae.initialization import solve_consistent_initial_conditions
# import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_rodrigo.cstr_rodrigo_model import make_model
from idaes.apps.caprese.data_manager import (
        PlantDataManager,
        ControllerDataManager,
        EstimatorDataManager
        )
from idaes.apps.caprese.plotlibrary import (
        plot_setpoint_tracking_results,
        plot_control_input,
        plot_estimation_results,)

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
    m_estimator = make_model(horizon=10., ntfe=10, ntcp=2, bounds=True)
    m_controller = make_model(horizon=10, ntfe=5, ntcp=2, bounds=True)
    sample_time = 2.
    m_plant = make_model(horizon=sample_time, ntfe=2, ntcp=2, bounds = True)
    time_plant = m_plant.t

    simulation_horizon = 20
    n_samples_to_simulate = round(simulation_horizon/sample_time)

    samples_to_simulate = [time_plant.first() + i*sample_time
                           for i in range(1, n_samples_to_simulate)]

    # We must identify for the dynamic system which variables are our
    # inputs and measurements.
    inputs = [
            m_plant.Tjinb[0],
            ]
    measurements = [
            m_plant.Tall[0, "T"],
            # m_plant.Tall[0, "Tj"],
            m_plant.Ca[0],
            ]

    # Construct the "Dynamic simulator" object
    dyna = DynamicSim(
            plant_model = m_plant,
            plant_time_set = m_plant.t,
            estimator_model = m_estimator, 
            estimator_time_set = m_estimator.t,
            controller_model = m_controller, 
            controller_time_set = m_controller.t,
            inputs_at_t0 = inputs,
            measurements_at_t0 = measurements,
            sample_time = sample_time,
            )

    plant = dyna.plant
    estimator = dyna.estimator
    controller = dyna.controller

    p_t0 = dyna.plant.time.first()
    e_t0 = dyna.estimator.time.first()
    c_t0 = dyna.controller.time.first()
    p_ts = dyna.plant.sample_points[1]
    e_ts = dyna.estimator.sample_points[1]
    c_ts = dyna.controller.sample_points[1]

    #--------------------------------------------------------------------------
    # Declare variables of interest for plotting.
    # It's ok not declaring anything. The data manager will still save some 
    # important data, but the user should use the default string of CUID for plotting afterward.
    states_of_interest = [Reference(dyna.plant.mod.Ca[:]),
                          Reference(dyna.plant.mod.Tall[:, "T"])]
    inputs_of_interest = [Reference(dyna.plant.mod.Tjinb[...])]

    plant_data = PlantDataManager(plant,
                                  states_of_interest,
                                  inputs_of_interest)

    controller_data = ControllerDataManager(controller,
                                            states_of_interest,
                                            inputs_of_interest)

    estimator_data = EstimatorDataManager(estimator,
                                          states_of_interest)

    #--------------------------------------------------------------------------
    # Plant setup
    solve_consistent_initial_conditions(plant, plant.time, solver)


    # Controller setup
    solve_consistent_initial_conditions(controller, controller.time, solver)

    # We now perform the "RTO" calculation: Find the optimal steady state
    # to achieve the following setpoint
    setpoint = [(controller.mod.Ca[0], 0.018)]
    setpoint_weights = [(controller.mod.Ca[0], 1.)]

    dyna.controller.add_single_time_optimization_objective(setpoint,
                                                           setpoint_weights)
    dyna.controller.solve_single_time_optimization(solver,
                                                   ic_type = "differential_var",
                                                   require_steady = True,
                                                   load_setpoints = True)

    # Now we are ready to construct the tracking NMPC problem
    tracking_weights = [
            *((v, 1.) for v in dyna.controller.vectors.differential[:,0]),
            *((v, 1.) for v in dyna.controller.vectors.input[:,0]),
            ]

    dyna.controller.add_tracking_objective(tracking_weights)

    dyna.controller.constrain_control_inputs_piecewise_constant()

    dyna.controller.initialize_to_initial_conditions()

    # Set up input noises that will be applied to control inputs
    variance = [
        (plant.mod.Tjinb[0], 0.01),
        ]
    dyna.plant.set_variance(variance)
    input_variance = [v.variance for v in plant.INPUT_BLOCK[:].var]
    input_noise_bounds = [(var[p_t0].lb, var[p_t0].ub) for var in plant.INPUT_BLOCK[:].var]
    #--------------------------------------------------------------------------


    # Estimator setup
    # Here we solve for a steady state and use it to fill in past measurements
    desired_ss = [(estimator.mod.Ca[0], 0.021)]
    ss_weights = [(estimator.mod.Ca[0], 1.)]
    dyna.estimator.initialize_past_info_with_steady_state(desired_ss, ss_weights, solver)

    # Now we are ready to construct the objective function for MHE
    model_disturbance_weights = [
            (estimator.mod.Ca[0], 1.),
            (estimator.mod.Tall[0, "T"], 1.),
            (estimator.mod.Tall[0, "Tj"], 1.),
            ]

    measurement_noise_weights = [
            (estimator.mod.Ca[0], 100.),
            (estimator.mod.Tall[0, "T"], 20.),
            ]   

    dyna.estimator.add_noise_minimize_objective(model_disturbance_weights,
                                                measurement_noise_weights)

    # Set up measurement noises that will be applied to measurements
    variance = [
        (dyna.estimator.mod.Tall[0, "T"], 0.05),
        (dyna.estimator.mod.Ca[0], 1.0E-2),
        ]
    dyna.estimator.set_variance(variance)
    measurement_variance = [
        v.variance for v in estimator.MEASUREMENT_BLOCK[:].var
        ]
    measurement_noise_bounds = [
        (var[e_t0].lb, var[e_t0].ub) for var in estimator.MEASUREMENT_BLOCK[:].var
        ]
    random.seed(246)
    #--------------------------------------------------------------------------

    plant_data.save_initial_plant_data()    

    # Solve the first control problem
    dyna.controller.vectors.input[...].unfix()
    dyna.controller.vectors.input[:,0].fix()
    solver.solve(dyna.controller, tee=True)
    controller_data.save_controller_data(iteration = 0)

    # Extract inputs from controller and inject them into plant
    inputs = controller.generate_inputs_at_time(c_ts)
    plant.inject_inputs(inputs)

    # This "initialization" really simulates the plant with the new inputs.
    dyna.plant.initialize_by_solving_elements(solver)
    dyna.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(dyna.plant, tee = True)
    plant_data.save_plant_data(iteration = 0)

    # Extract measurements from the plant and inject them into MHE
    measurements = dyna.plant.generate_measurements_at_time(p_ts)
    dyna.estimator.load_measurements(measurements,
                                     target = "actualmeasurement",
                                     timepoint = estimator.time.last())
    dyna.estimator.load_inputs_into_last_sample(inputs)

    # Solve the first estimation problem
    dyna.estimator.check_var_con_dof(skip_dof_check = False)
    solver.solve(dyna.estimator, tee=True)
    estimator_data.save_estimator_data(iteration = 0)

    for i in range(1, 11):
        print('\nENTERING MHE LOOP ITERATION %s\n' % i)

        estimates = dyna.estimator.generate_estimates_at_time(estimator.time.last())

        dyna.controller.advance_one_sample()
        dyna.controller.load_initial_conditions(estimates)    

        solver.solve(dyna.controller, tee = True)
        controller_data.save_controller_data(iteration = i)

        dyna.plant.advance_one_sample()
        dyna.plant.initialize_to_initial_conditions()
        inputs = controller.generate_inputs_at_time(c_ts)
        inputs = apply_noise_with_bounds(
                inputs,
                input_variance,
                random.gauss,
                input_noise_bounds,
                )
        plant.inject_inputs(inputs)

        dyna.plant.initialize_by_solving_elements(solver)
        dyna.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(dyna.plant, tee = True)    
        plant_data.save_plant_data(iteration = i)

        measurements = dyna.plant.generate_measurements_at_time(p_ts)
        measurements = apply_noise_with_bounds(
                    measurements,
                    measurement_variance,
                    random.gauss,
                    measurement_noise_bounds,
                    )
        dyna.estimator.advance_one_sample()
        dyna.estimator.load_measurements(measurements,
                                         target = "actualmeasurement",
                                         timepoint = estimator.time.last())
        dyna.estimator.load_inputs_into_last_sample(inputs)

        dyna.estimator.check_var_con_dof(skip_dof_check = False)
        solver.solve(dyna.estimator, tee = True)
        estimator_data.save_estimator_data(iteration = i)

    plot_setpoint_tracking_results(states_of_interest,
                                   plant_data.plant_df,
                                   controller_data.setpoint_df)

    plot_control_input(inputs_of_interest,
                       plant_data.plant_df)

    plot_estimation_results(states_of_interest,
                            plant_data.plant_df,
                            estimator_data.estimator_df)

    return dyna, plant_data, controller_data, estimator_data

if __name__ == '__main__':
    dyna, plant_data, controller_data, estimator_data = main()
