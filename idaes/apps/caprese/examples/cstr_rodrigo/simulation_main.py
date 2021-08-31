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
Example for Caprese's module for simulation of a plant.
"""
import random
from idaes.apps.caprese.dynamic_builder import DynamicSim
# from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory, Reference
from pyomo.dae.initialization import solve_consistent_initial_conditions
# import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_rodrigo.cstr_rodrigo_model import make_model
from idaes.apps.caprese.data_manager import PlantDataManager

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
    sample_time = 2.
    m_plant = make_model(horizon=sample_time, ntfe=4, ntcp=2, bounds = True)
    time_plant = m_plant.t
    
    # We must identify for the plant which variables are our
    # inputs and measurements.
    inputs = [
            m_plant.Tjinb[0],
            ]
    measurements = [
            m_plant.Tall[0, "T"],
            # m_plant.Tall[0, "Tj"],
            m_plant.Ca[0],
            ]
    
    # Construct the "plant simulator" object
    simulator = DynamicSim(
                    plant_model=m_plant,
                    plant_time_set=m_plant.t,
                    inputs_at_t0=inputs,
                    measurements_at_t0=measurements,
                    sample_time=sample_time,
                    )
    
    plant = simulator.plant
    
    p_t0 = simulator.plant.time.first()
    p_ts = simulator.plant.sample_points[1]
    #--------------------------------------------------------------------------
    # Declare variables of interest for plotting.
    # It's ok not declaring anything. The data manager will still save some 
    # important data, but the user should use the default string of CUID for plotting afterward.
    states_of_interest = [Reference(simulator.plant.mod.Ca[:]),
                          Reference(simulator.plant.mod.Tall[:, "T"])]
    inputs_of_interest = [Reference(simulator.plant.mod.Tjinb[...])]
    
    # Set up data manager to save plant data
    data_manager = PlantDataManager(plant, 
                                    states_of_interest,
                                    inputs_of_interest)
    #--------------------------------------------------------------------------
    solve_consistent_initial_conditions(plant, plant.time, solver)

    input_list = {ind: 250.+ind*5 if ind<=5 else 260.-ind*5 for ind in range(0, 11)}
    
    data_manager.save_initial_plant_data()
    
    plant.inject_inputs([input_list[0]])

    # This "initialization" really simulates the plant with the new inputs.
    simulator.plant.initialize_by_solving_elements(solver)
    simulator.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(simulator.plant, tee = True)
    data_manager.save_plant_data(iteration = 0)

    for i in range(1,11):
        print('\nENTERING SIMULATOR LOOP ITERATION %s\n' % i)
        
        simulator.plant.advance_one_sample()
        simulator.plant.initialize_to_initial_conditions()
        simulator.plant.inject_inputs([input_list[i]])
        
        simulator.plant.initialize_by_solving_elements(solver)
        simulator.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(simulator.plant, tee = True)
        data_manager.save_plant_data(iteration = i)
        
    data_manager.plot_plant_state_evolution(states_of_interest)
    data_manager.plot_control_input(inputs_of_interest)
    
    return simulator, data_manager

if __name__ == '__main__':
    simulator, data_manager = main()