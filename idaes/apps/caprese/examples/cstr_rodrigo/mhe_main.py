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
from cstr_rodrigo_model2 import make_model
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
    m_estimator = make_model(horizon=4., ntfe=4, ntcp=2, bounds=True)
    sample_time = 2.
    m_plant = make_model(horizon=sample_time, ntfe=2, ntcp=2, bounds = True)
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
            m_estimator.Tall[0, "T"],
            # m_estimator.Tall[0, "Tj"],
            m_estimator.Ca[0],
            ]
    
    # Construct the "MHE simulator" object
    mhe = MHESim(
            plant_model=m_plant,
            plant_time_set=m_plant.t,
            estimator_model=m_estimator, 
            estimator_time_set=m_estimator.t,
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
        
    # Now we are ready to construct the objective function for MHE
    model_disturbance_weights = [
            (estimator.mod.Ca[0], 0.1),
            (estimator.mod.Tall[0, "T"], 0.2),
            (estimator.mod.Tall[0, "Tj"], 0.3),
            ]

    measurement_noise_weights = [
            (estimator.mod.Ca[0], 10.),
            (estimator.mod.Tall[0, "T"], 20.),
            ]   
    
    mhe.estimator.add_noise_minimize_objective(model_disturbance_weights,
                                               measurement_noise_weights)
    
    mhe.estimator.MHE_initialize_to_initial_conditions()
    
    # This "initialization" really simulates the plant with the new inputs.
    mhe.plant.initialize_by_solving_elements(solver)
    mhe.plant.mod.Tjinb.fix() #Fix the input to solve the plant
    solver.solve(mhe.plant, tee = True)
    
    measurements = mhe.plant.generate_measurements_at_time(p_ts)
    #apply measurement error here
    mhe.estimator.load_measurements_for_MHE(measurements)
    
    # Solve the first estimation problem
    mhe.estimator.check_var_con_dof()
    solver.solve(mhe.estimator, tee=True)
    
    
    for i in range(1,11):
        print('\nENTERING NMPC LOOP ITERATION %s\n' % i)
        
        mhe.plant.advance_one_sample()
        mhe.plant.initialize_to_initial_conditions()
        #inject inputs here if it's updated
        
        mhe.plant.initialize_by_solving_elements(solver)
        mhe.plant.mod.Tjinb.fix() #Fix the input to solve the plant
        solver.solve(mhe.plant)
        
        measurements = mhe.plant.generate_measurements_at_time(p_ts)
        #apply measurement error here
        
        mhe.estimator.MHE_advance_one_sample()
        mhe.estimator.load_measurements_for_MHE(measurements)
        
        mhe.estimator.check_var_con_dof()
        solver.solve(mhe.estimator, tee=True)
        
    
    return mhe

if __name__ == '__main__':
    mhe = main()