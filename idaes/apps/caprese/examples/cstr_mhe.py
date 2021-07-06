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
# from pyomo.dae.initialization import solve_consistent_initial_conditions
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
    
    return mhe

    plant = nmpc.plant
    estimator = nmpc.estimator
    
    p_t0 = nmpc.plant.time.first()
    c_t0 = nmpc.estimator.time.first()
    p_ts = nmpc.plant.sample_points[1]
    c_ts = nmpc.estimator.sample_points[1]

if __name__ == '__main__':
    mhe = main()
