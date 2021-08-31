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
Example for Caprese's module for MHE.
"""
import random
from idaes.apps.caprese.dynamic_builder import DynamicSim
from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory, Reference
from pyomo.dae.initialization import solve_consistent_initial_conditions
# import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_rodrigo.cstr_rodrigo_model import make_model
from idaes.apps.caprese.data_manager import EstimatorDataManager

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
    sample_time = 2.
    m_plant = make_model(horizon=sample_time, ntfe=4, ntcp=2, bounds = True)
    time_plant = m_plant.t

    simulation_horizon = 20
    n_samples_to_simulate = round(simulation_horizon/sample_time)

    samples_to_simulate = [time_plant.first() + i*sample_time
                           for i in range(1, n_samples_to_simulate)]

    # We must identify for the estimator which variables are our
    # inputs and measurements.
    inputs = [
            m_plant.Tjinb[0],
            ]
    measurements = [
            m_plant.Tall[0, "T"],
            # m_plant.Tall[0, "Tj"],
            m_plant.Ca[0],
            ]
    
    # Construct the "MHE simulator" object
    mhe = DynamicSim(
            plant_model=m_plant,
            plant_time_set=m_plant.t,
            estimator_model=m_estimator, 
            estimator_time_set=m_estimator.t,
            inputs_at_t0=inputs,
            measurements_at_t0=measurements,
            sample_time=sample_time,
            )

    plant = mhe.plant
    estimator = mhe.estimator
    
    p_t0 = mhe.plant.time.first()
    e_t0 = mhe.estimator.time.first()
    p_ts = mhe.plant.sample_points[1]
    e_ts = mhe.estimator.sample_points[1]
    #--------------------------------------------------------------------------
    # Declare variables of interest for plotting.
    # It's ok not declaring anything. The data manager will still save some 
    # important data, but the user should use the default string of CUID for plotting afterward.
    states_of_interest = [Reference(mhe.plant.mod.Ca[:]),
                          Reference(mhe.plant.mod.Tall[:, "T"])]

    # Set up data manager to save estimation data
    data_manager = EstimatorDataManager(plant, 
                                       estimator,
                                       states_of_interest,)
    #--------------------------------------------------------------------------
    solve_consistent_initial_conditions(plant, plant.time, solver)
    
    # Here we solve for a steady state and use it to fill in past measurements
    desired_ss = [(estimator.mod.Ca[0], 0.021)]
    ss_weights = [(estimator.mod.Ca[0], 1.)]
    mhe.estimator.initialize_past_info_with_steady_state(desired_ss, ss_weights, solver)
        
    # Now we are ready to construct the objective function for MHE
    model_disturbance_variances = [
            (estimator.mod.Ca[0], 1.),
            (estimator.mod.Tall[0, "T"], 1.),
            (estimator.mod.Tall[0, "Tj"], 1.),
            ]

    measurement_noise_variances = [
            (mhe.estimator.mod.Tall[0, "T"], 0.05),
            (mhe.estimator.mod.Ca[0], 1.0E-2),
            ]  
    
    mhe.estimator.add_noise_minimize_objective(model_disturbance_variances,
                                               measurement_noise_variances,
                                               givenform = "variance")
    
    #-------------------------------------------------------------------------
    # Set up measurement noises that will be applied to measurements
    mhe.estimator.set_variance(measurement_noise_variances)
    measurement_variance = [v.variance for v in estimator.measurement_vars]
    measurement_noise_bounds = [
            (var[e_t0].lb, var[e_t0].ub) for var in estimator.measurement_vars
            ]
    #-------------------------------------------------------------------------
    
    data_manager.save_initial_plant_data()
    
    # This "initialization" really simulates the plant with the new inputs.
    mhe.plant.initialize_by_solving_elements(solver)
    mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(mhe.plant, tee = True)
    data_manager.save_plant_data(iteration = 0)
    
    # Extract measurements from the plant and inject them into MHE
    measurements = mhe.plant.generate_measurements_at_time(p_ts)
    mhe.estimator.load_measurements(measurements,
                                    target = "actualmeasurement",
                                    timepoint = estimator.time.last())
    mhe.estimator.load_inputs_for_MHE([mhe.plant.mod.Tjinb[p_ts].value])
    
    # Solve the first estimation problem
    mhe.estimator.check_var_con_dof(skip_dof_check = False)
    solver.solve(mhe.estimator, tee=True)
    data_manager.save_estimator_data(iteration = 0)
    
    cinput = {ind: 250.+ind*5 if ind<=5 else 260.-ind*5 for ind in range(1, 11)}
    
    
    for i in range(1,11):
        print('\nENTERING MHE LOOP ITERATION %s\n' % i)
        
        mhe.plant.advance_one_sample()
        mhe.plant.initialize_to_initial_conditions()
        inputs = [cinput[i]]
        mhe.plant.inject_inputs(inputs)
        
        mhe.plant.initialize_by_solving_elements(solver)
        mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(mhe.plant, tee = True)
        data_manager.save_plant_data(iteration = i)
        
        measurements = mhe.plant.generate_measurements_at_time(p_ts)
        measurements = apply_noise_with_bounds(
                    measurements,
                    measurement_variance,
                    random.gauss,
                    measurement_noise_bounds,
                    )
        
        mhe.estimator.advance_one_sample()
        mhe.estimator.load_measurements(measurements,
                                        target = "actualmeasurement",
                                        timepoint = estimator.time.last())
        mhe.estimator.load_inputs_for_MHE(inputs)
        
        mhe.estimator.check_var_con_dof(skip_dof_check = False)
        # mhe.estimator.vectors.modeldisturbance[...].fix(0.0)
        solver.solve(mhe.estimator, tee=True)
        data_manager.save_estimator_data(iteration = i)
        
    data_manager.plot_estimation_results(states_of_interest)
    return mhe, data_manager

if __name__ == '__main__':
    mhe, data_manager = main()