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
from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory, Reference
from pyomo.dae.initialization import solve_consistent_initial_conditions
# import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model
from idaes.apps.caprese.data_manager import PlantDataManager
from idaes.apps.caprese.plotlibrary import (
        plot_plant_state_evolution,
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
    sample_time = 0.5
    m_plant = make_model(horizon=sample_time, ntfe=5, ntcp=2)
    time_plant = m_plant.fs.time

    # We must identify for the plant which variables are our
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

    # Construct the "plant simulator" object
    simulator = DynamicSim(
                    plant_model=m_plant,
                    plant_time_set=m_plant.fs.time,
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
    states_of_interest = [Reference(simulator.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','S']),
                          Reference(simulator.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','E']),
                          # Reference(simulator.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','C']),
                          # Reference(simulator.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','P']),
                          # Reference(simulator.plant.mod.fs.cstr.control_volume.material_holdup[:,'aq','Solvent']),
                          Reference(simulator.plant.mod.fs.cstr.control_volume.energy_holdup[:,'aq']),
                          ]   

    inputs_of_interest = [Reference(simulator.plant.mod.fs.mixer.S_inlet_state[:].flow_vol),
                          # Reference(simulator.plant.mod.fs.mixer.E_inlet_state[:].flow_vol),
                          ]

    # Set up data manager to save plant data
    data_manager = PlantDataManager(plant, 
                                    states_of_interest,
                                    inputs_of_interest)
    #--------------------------------------------------------------------------
    solve_consistent_initial_conditions(plant, plant.time, solver)

    cinput1 = [0.5608456705408656, 3.4818166997491384, 5.0, 0.9629431563506397, 2.0623866186035156, 
               4.9999999797327686, 2.285805028476981, 3.913753219840146, 3.4585265451075538, 5.0]
    cinput2 = [0.28666548361218924, 0.01, 0.01, 0.01, 0.12654063510571273,
               0.01, 0.9999996329001195, 0.242203179025321, 0.7110096123027149, 0.01]

    data_manager.save_initial_plant_data()

    cinput = [cinput1[0], cinput2[0]]
    plant.inject_inputs(cinput)

    # This "initialization" really simulates the plant with the new inputs.
    simulator.plant.initialize_by_solving_elements(solver)
    simulator.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(simulator.plant, tee = True)
    data_manager.save_plant_data(iteration = 0)

    for i in range(1, 10):
        print('\nENTERING SIMULATOR LOOP ITERATION %s\n' % i)

        simulator.plant.advance_one_sample()
        simulator.plant.initialize_to_initial_conditions()
        cinput = [cinput1[i], cinput2[i]]
        simulator.plant.inject_inputs(cinput)

        simulator.plant.initialize_by_solving_elements(solver)
        simulator.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(simulator.plant, tee = True)
        data_manager.save_plant_data(iteration = i)

    plot_plant_state_evolution(states_of_interest, data_manager.plant_df)
    plot_control_input(inputs_of_interest, data_manager.plant_df)

    return simulator, data_manager

if __name__ == '__main__':
    simulator, data_manager = main()  
