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
# from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory
from pyomo.dae.initialization import solve_consistent_initial_conditions
# import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model
import pandas as pd
import matplotlib.pyplot as plt

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
                [var[t].value for t in time],
                index=[t for t in time])
        self.setpoint_series = pd.Series(
                [initial if t < t_switch else setpoint for t in time])

    def plot(self):
        # fig, ax can be formatted to the user's liking
        fig, ax = plt.subplots()
        if self.name is not None:
            self.data_series.plot(label=self.name)
        else:
            self.data_series.plot()
        return fig, ax

def main(plot_switch=False):
    m_estimator = make_model(horizon=1, ntfe=4, ntcp=2, bounds=True)
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
    #measurement_error_constraints are only defined at sample points
    solve_consistent_initial_conditions(estimator, estimator.time, solver, suppress_warnings = True)
    
    estimator.mod.fs.cstr.volume[0].unfix()
    
    # Now we are ready to construct the objective function for MHE
    model_disturbance_weights = [
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','S'], 0.1),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','E'], 0.1),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','C'], 0.1),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','P'], 0.1),
            (estimator.mod.fs.cstr.control_volume.material_holdup[0,'aq','Solvent'], 0.1),
            (estimator.mod.fs.cstr.control_volume.energy_holdup[0,'aq'], 0.1),
            ]

    measurement_noise_weights = [
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'C'], 10.),
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'E'], 10.),
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'S'], 10.),
            (estimator.mod.fs.cstr.outlet.conc_mol[0, 'P'], 10.),
            (estimator.mod.fs.cstr.outlet.temperature[0], 0.1),
            (estimator.mod.fs.cstr.volume[0], 1.),
            ]   
    
    mhe.estimator.add_noise_minimize_objective(model_disturbance_weights,
                                               measurement_noise_weights)
    
    
    mhe.estimator.MHE_initialize_to_initial_conditions()
    
    # This "initialization" really simulates the plant with the new inputs.
    mhe.plant.initialize_by_solving_elements(solver)
    mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
    solver.solve(mhe.plant, tee = True)
    
    measurements = mhe.plant.generate_measurements_at_time(p_ts)
    #apply measurement error here
    mhe.estimator.load_measurements_for_MHE(measurements)
    
    # Solve the first estimation problem
    mhe.estimator.check_var_con_dof(skip_dof_check = False)
    solver.solve(mhe.estimator, tee=True)
    
    soi_string = ["S", "E", "C", "P", "Solvent", "Energy"]
    plant_rec = {}
    MHE_rec = {}
    for ss in soi_string:
        plant_rec[ss] = []
        MHE_rec[ss] = []
        
    cinput1 = [0.5608456705408656, 3.4818166997491384, 5.0, 0.9629431563506397, 2.0623866186035156, 
               4.9999999797327686, 2.285805028476981, 3.913753219840146, 3.4585265451075538, 5.0]
    cinput2 = [0.28666548361218924, 0.01, 0.01, 0.01, 0.12654063510571273,
               0.01, 0.9999996329001195, 0.242203179025321, 0.7110096123027149, 0.01]
    
    for i in range(1,11):
        print('\nENTERING NMPC LOOP ITERATION %s\n' % i)
        
        mhe.plant.advance_one_sample()
        mhe.plant.initialize_to_initial_conditions()
        #inject inputs here if it's updated
        cinput = [cinput1[i-1], cinput2[i-1]]
        mhe.plant.inject_inputs(cinput)
        
        mhe.plant.initialize_by_solving_elements(solver)
        mhe.plant.vectors.input[...].fix() #Fix the input to solve the plant
        solver.solve(mhe.plant)
        
        # def save_plant_data(plant, rec_dict):
        plant_rec["S"].append(plant.mod.fs.cstr.control_volume.material_holdup[p_ts,'aq','S'].value)
        plant_rec["E"].append(plant.mod.fs.cstr.control_volume.material_holdup[p_ts,'aq','E'].value)
        plant_rec["C"].append(plant.mod.fs.cstr.control_volume.material_holdup[p_ts,'aq','C'].value)
        plant_rec["P"].append(plant.mod.fs.cstr.control_volume.material_holdup[p_ts,'aq','P'].value)
        plant_rec["Solvent"].append(plant.mod.fs.cstr.control_volume.material_holdup[p_ts,'aq','Solvent'].value)
        plant_rec["Energy"].append(plant.mod.fs.cstr.control_volume.energy_holdup[p_ts,'aq'].value)
            
        measurements = mhe.plant.generate_measurements_at_time(p_ts)
        #apply measurement error here
        
        mhe.estimator.MHE_advance_one_sample()
        mhe.estimator.load_measurements_for_MHE(measurements)
        mhe.estimator.load_inputs_for_MHE(cinput)
        # mhe.estimator.vectors.input[...].fix(cinput[i])
        
        mhe.estimator.check_var_con_dof(skip_dof_check = False)
        # mhe.estimator.vectors.modeldisturbance[...].fix(0.0)
        solver.solve(mhe.estimator, tee=True)
        
        tl = mhe.estimator.time.last()
        MHE_rec["S"].append(estimator.mod.fs.cstr.control_volume.material_holdup[tl,'aq','S'].value)
        MHE_rec["E"].append(estimator.mod.fs.cstr.control_volume.material_holdup[tl,'aq','E'].value)
        MHE_rec["C"].append(estimator.mod.fs.cstr.control_volume.material_holdup[tl,'aq','C'].value)
        MHE_rec["P"].append(estimator.mod.fs.cstr.control_volume.material_holdup[tl,'aq','P'].value)
        MHE_rec["Solvent"].append(estimator.mod.fs.cstr.control_volume.material_holdup[tl,'aq','Solvent'].value)
        MHE_rec["Energy"].append(estimator.mod.fs.cstr.control_volume.energy_holdup[tl,'aq'].value)
    
    return mhe, plant_rec, MHE_rec

if __name__ == '__main__':
    mhe, plant_rec, MHE_rec = main()
    
    import matplotlib.pyplot as plt
    soi_string = ["S", "E", "C", "P", "Solvent", "Energy"]
    for ind, item in enumerate(soi_string):
        plt.figure(ind)
        plt.title(item)
        plt.plot(plant_rec[item], "r")
        plt.plot(MHE_rec[item], "b")
        
    plt.show()
