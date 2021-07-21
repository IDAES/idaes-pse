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
from idaes.apps.caprese.mhe import MHESim
from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory
from pyomo.dae.initialization import solve_consistent_initial_conditions
# import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model
import pandas as pd

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
    m_estimator = make_model(horizon=8, ntfe=16, ntcp=2, bounds=True)
    sample_time = 0.5
    m_plant = make_model(horizon=sample_time, ntfe=5, ntcp=2)
    time_plant = m_plant.fs.time

    simulation_horizon = 20
    n_samples_to_simulate = round(simulation_horizon/sample_time)

    samples_to_simulate = [time_plant.first() + i*sample_time
                           for i in range(1, n_samples_to_simulate)]

    # We must identify for the estimator which variables are our
    # inputs and measurements.
    inputs = [
            m_plant.fs.mixer.S_inlet.flow_vol[0],
            m_plant.fs.mixer.E_inlet.flow_vol[0],
            ]
    measurements = [
            m_estimator.fs.cstr.outlet.conc_mol[0, 'C'],
            m_estimator.fs.cstr.outlet.conc_mol[0, 'E'],
            m_estimator.fs.cstr.outlet.conc_mol[0, 'S'],
            m_estimator.fs.cstr.outlet.conc_mol[0, 'P'],
            m_estimator.fs.cstr.outlet.temperature[0],
            m_estimator.fs.cstr.volume[0],
            ]
    
    # Construct the "MHE simulator" object
    mhe = MHESim(
            plant_model=m_plant,
            plant_time_set=m_plant.fs.time,
            estimator_model=m_estimator, 
            estimator_time_set=m_estimator.fs.time,
            inputs_at_t0=inputs,
            measurements=measurements,
            sample_time=sample_time,
            )
    
    plant = mhe.plant
    estimator = mhe.estimator
    
    p_t0 = mhe.plant.time.first()
    c_t0 = mhe.estimator.time.first()
    p_ts = mhe.plant.sample_points[1]
    c_ts = mhe.estimator.sample_points[1]

    solve_consistent_initial_conditions(plant, plant.time, solver)
    
    estimator.mod.fs.cstr.volume[0].unfix()
    
    # Here we solve for a steady state and use it to fill in past measurements
    desired_ss = [
        (estimator.mod.fs.cstr.outlet.conc_mol[0, 'P'], 0.4),
        (estimator.mod.fs.cstr.outlet.conc_mol[0, 'S'], 0.0),
        (estimator.mod.fs.cstr.control_volume.energy_holdup[0, 'aq'], 300),
        (estimator.mod.fs.mixer.E_inlet.flow_vol[0], 0.1),
        (estimator.mod.fs.mixer.S_inlet.flow_vol[0], 2.0),
        (estimator.mod.fs.cstr.volume[0], 1.0),
        ]
    ss_weights = [
        (estimator.mod.fs.cstr.outlet.conc_mol[0, 'P'], 1.),
        (estimator.mod.fs.cstr.outlet.conc_mol[0, 'S'], 1.),
        (estimator.mod.fs.cstr.control_volume.energy_holdup[0, 'aq'], 1.),
        (estimator.mod.fs.mixer.E_inlet.flow_vol[0], 1.),
        (estimator.mod.fs.mixer.S_inlet.flow_vol[0], 1.),
        (estimator.mod.fs.cstr.volume[0], 1.),
        ]
    mhe.estimator.initialize_past_info_with_steady_state(desired_ss, ss_weights, solver)
    
    # Now we are ready to construct the objective function for MHE
    model_disturbance_weights = [
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','S'], 1.),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','E'], 1.),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','C'], 1.),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','P'], 1.),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','Solvent'], 1.),
            (estimator.mod.fs.cstr.control_volume.energy_holdup[0,'aq'], 1.),
            ]

    measurement_noise_weights = [
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'C'], 100.),
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'E'], 20.),
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'S'], 50.),
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'P'], 20.),
            (estimator.mod.fs.cstr.outlet.temperature[0], 10.),
            (estimator.mod.fs.cstr.volume[0], 20.),
            ]   
    
    mhe.estimator.add_noise_minimize_objective(model_disturbance_weights,
                                               measurement_noise_weights)
    
    #-------------------------------------------------------------------------
    # Set up measurement noises that will be applied to measurements
    cstr = mhe.estimator.mod.fs.cstr
    variance = [
            (cstr.outlet.conc_mol[0.0, 'C'], 0.01),
            (cstr.outlet.conc_mol[0.0, 'E'], 0.05),
            (cstr.outlet.conc_mol[0.0, 'S'], 0.02),
            (cstr.outlet.conc_mol[0.0, 'P'], 0.05),
            (cstr.outlet.temperature[0.0], 0.1),
            (cstr.volume[0.0], 0.05),
            ]
    mhe.estimator.set_variance(variance)
    measurement_variance = [v.variance for v in estimator.measurement_vars]
    measurement_noise_bounds = [
            (0.0, var[c_t0].ub) for var in estimator.measurement_vars
            ]
    #-------------------------------------------------------------------------
    
    # Set up pandas dataframe to save plant data
    mhe.plant.initialize_plant_dataframe()
    
    # This "initialization" really simulates the plant with the new inputs.
    mhe.plant.initialize_by_solving_elements(solver)
    mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(mhe.plant, tee = True)
    mhe.plant.record_plant_data()
    
    # Extract measurements from the plant and inject them into MHE
    measurements = mhe.plant.generate_measurements_at_time(p_ts)
    mhe.estimator.load_measurements(measurements,
                                    target = "actualmeasurement",
                                    timepoint = estimator.time.last())
    mhe.estimator.load_inputs_for_MHE([mhe.plant.mod.fs.mixer.S_inlet.flow_vol[p_ts].value,
                                       mhe.plant.mod.fs.mixer.E_inlet.flow_vol[p_ts].value])
    
    # Set up pandas dataframe to save estimation results
    mhe.estimator.initialize_estimator_dataframe()
    
    # Solve the first estimation problem
    mhe.estimator.check_var_con_dof(skip_dof_check = False)
    solver.solve(mhe.estimator, tee=True)
    mhe.estimator.record_estimator_data()
        
    cinput1 = [0.5608456705408656, 3.4818166997491384, 5.0, 0.9629431563506397, 2.0623866186035156, 
               4.9999999797327686, 2.285805028476981, 3.913753219840146, 3.4585265451075538, 5.0]
    cinput2 = [0.28666548361218924, 0.01, 0.01, 0.01, 0.12654063510571273,
               0.01, 0.9999996329001195, 0.242203179025321, 0.7110096123027149, 0.01]
    
    for i in range(1,11):
        print('\nENTERING MHE LOOP ITERATION %s\n' % i)
        
        mhe.plant.advance_one_sample()
        mhe.plant.initialize_to_initial_conditions()
        cinput = [cinput1[i-1], cinput2[i-1]]
        mhe.plant.inject_inputs(cinput)
        
        mhe.plant.initialize_by_solving_elements(solver)
        mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(mhe.plant, tee = True)
        mhe.plant.record_plant_data()
            
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
        mhe.estimator.load_inputs_for_MHE(cinput)
        
        mhe.estimator.check_var_con_dof(skip_dof_check = False)
        # mhe.estimator.vectors.modeldisturbance[...].fix(0.0)
        solver.solve(mhe.estimator, tee=True)
        mhe.estimator.record_estimator_data()
        
        
    return mhe

if __name__ == '__main__':
    mhe = main() 

    # Plot the estiamtion result (estimated states vs. real states)
    vars_soi = (mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','S'],
                mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','E'],
                mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','C'],
                mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','P'],
                mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','Solvent'],
                mhe.plant.mod.fs.cstr.control_volume.energy_holdup[:,'aq'],
                )
    mhe.plot_estimation_result(vars_soi, xlabel = "time",)