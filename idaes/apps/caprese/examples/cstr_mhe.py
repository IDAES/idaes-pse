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
from idaes.apps.caprese.examples.cstr_model import make_model
from idaes.apps.caprese.data_manager import PlantDataManager
from idaes.apps.caprese.data_manager import EstimatorDataManager
from idaes.apps.caprese.plotlibrary import plot_estimation_results

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
            m_plant.fs.cstr.outlet.conc_mol[0, 'C'],
            m_plant.fs.cstr.outlet.conc_mol[0, 'E'],
            m_plant.fs.cstr.outlet.conc_mol[0, 'S'],
            m_plant.fs.cstr.outlet.conc_mol[0, 'P'],
            m_plant.fs.cstr.outlet.temperature[0],
            m_plant.fs.cstr.volume[0],
            ]

    # Construct the "MHE simulator" object
    mhe = DynamicSim(
            plant_model=m_plant,
            plant_time_set=m_plant.fs.time,
            estimator_model=m_estimator,
            estimator_time_set=m_estimator.fs.time,
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
    states_of_interest = (Reference(mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','S']),
                          Reference(mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','E']),
                          # Reference(mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','C']),
                          # Reference(mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','P']),
                          # Reference(mhe.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','Solvent']),
                          Reference(mhe.plant.mod.fs.cstr.control_volume.energy_holdup[:,'aq']),
                          )

    # Set up data manager to save estimation data
    plant_data = PlantDataManager(plant, states_of_interest)
    estimator_data = EstimatorDataManager(estimator, states_of_interest)
    #--------------------------------------------------------------------------
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
    measurement_variance = [
            v.variance for v in estimator.MEASUREMENT_BLOCK[:].var
            ]
    measurement_noise_bounds = [
            (0.0, var[e_t0].ub) for var in estimator.MEASUREMENT_BLOCK[:].var
            ]
    #-------------------------------------------------------------------------

    plant_data.save_initial_plant_data()

    # This "initialization" really simulates the plant with the new inputs.
    mhe.plant.initialize_by_solving_elements(solver)
    mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(mhe.plant, tee = True)
    plant_data.save_plant_data(iteration = 0)

    # Extract measurements from the plant and inject them into MHE
    measurements = mhe.plant.generate_measurements_at_time(p_ts)
    mhe.estimator.load_measurements(measurements,
                                    timepoint = estimator.time.last())
    mhe.estimator.load_inputs_into_last_sample(
        mhe.plant.generate_inputs_at_time(p_ts))

    # Solve the first estimation problem
    mhe.estimator.check_var_con_dof(skip_dof_check = False)
    solver.solve(mhe.estimator, tee=True)
    estimator_data.save_estimator_data(iteration = 0)

    cinput1 = [0.56, 3.48, 5.00, 0.96, 2.06,
               5.00, 2.29, 3.91, 3.46, 5.0]
    cinput2 = [0.29, 0.01, 0.01, 0.01, 0.13,
               0.01, 1.00, 0.24, 0.71, 0.01]

    for i in range(1,11):
        print('\nENTERING MHE LOOP ITERATION %s\n' % i)

        mhe.plant.advance_one_sample()
        mhe.plant.initialize_to_initial_conditions()
        cinput = [cinput1[i-1], cinput2[i-1]]
        mhe.plant.inject_inputs(cinput)

        mhe.plant.initialize_by_solving_elements(solver)
        mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(mhe.plant, tee = True)
        plant_data.save_plant_data(iteration = i)

        measurements = mhe.plant.generate_measurements_at_time(p_ts)
        measurements = apply_noise_with_bounds(
                   measurements,
                   measurement_variance,
                   random.gauss,
                   measurement_noise_bounds,
                   )

        mhe.estimator.advance_one_sample()
        mhe.estimator.load_measurements(measurements,
                                        timepoint = estimator.time.last())
        mhe.estimator.load_inputs_into_last_sample(cinput)

        mhe.estimator.check_var_con_dof(skip_dof_check = False)
        # mhe.estimator.vectors.modeldisturbance[...].fix(0.0)
        solver.solve(mhe.estimator, tee=True)
        estimator_data.save_estimator_data(iteration = i)

    plot_estimation_results(states_of_interest,
                            plant_data.plant_df,
                            estimator_data.estimator_df)

    return mhe, plant_data, estimator_data

if __name__ == '__main__':
    mhe, plant_data, estimator_data = main()
