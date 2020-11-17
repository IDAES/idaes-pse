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
"""
import random
from idaes.apps.caprese.nmpc import NMPCSim
from idaes.apps.caprese.util import apply_noise_with_bounds
from pyomo.environ import SolverFactory
from pyomo.dae.initialization import solve_consistent_initial_conditions
import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model
import pandas as pd
import matplotlib.pyplot as plt

__author__ = "Robert Parker"


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

    # This tests the same model constructed in the test_nmpc_constructor_1 file
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
    sample_time = 0.5
    m_plant = make_model(horizon=sample_time, ntfe=5, ntcp=2)
    time_plant = m_plant.fs.time

    simulation_horizon = 60
    n_samples_to_simulate = round(simulation_horizon/sample_time)

    samples_to_simulate = [time_plant.first() + i*sample_time
                           for i in range(1, n_samples_to_simulate)]

    inputs = [
            m_plant.fs.mixer.S_inlet.flow_vol[0],
            m_plant.fs.mixer.E_inlet.flow_vol[0],
            ]
    measurements = [
            m_controller.fs.cstr.outlet.conc_mol[0, 'C'],
            m_controller.fs.cstr.outlet.conc_mol[0, 'E'],
            m_controller.fs.cstr.outlet.conc_mol[0, 'S'],
            m_controller.fs.cstr.outlet.conc_mol[0, 'P'],
            m_controller.fs.cstr.outlet.temperature[0],
            ]
    
    nmpc = NMPCSim(
            plant_model=m_plant,
            plant_time_set=m_plant.fs.time,
            controller_model=m_controller, 
            controller_time_set=m_controller.fs.time,
            inputs_at_t0=inputs,
            measurements=measurements,
            sample_time=sample_time,
            )

    plant = nmpc.plant
    controller = nmpc.controller

    solve_consistent_initial_conditions(plant, plant.time, solver)
    solve_consistent_initial_conditions(controller, controller.time, solver)

    setpoint = [
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'P'], 0.4),
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'S'], 0.0),
            (controller.mod.fs.cstr.control_volume.energy_holdup[0, 'aq'], 300),
            (controller.mod.fs.mixer.E_inlet.flow_vol[0], 0.1),
            (controller.mod.fs.mixer.S_inlet.flow_vol[0], 2.0),
            ]
    setpoint_weights = [
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'P'], 1.),
            (controller.mod.fs.cstr.outlet.conc_mol[0, 'S'], 1.),
            (controller.mod.fs.cstr.control_volume.energy_holdup[0, 'aq'], 1.),
            (controller.mod.fs.mixer.E_inlet.flow_vol[0], 1.),
            (controller.mod.fs.mixer.S_inlet.flow_vol[0], 1.),
            ]

    nmpc.controller.add_setpoint_objective(setpoint, setpoint_weights)
    nmpc.controller.solve_setpoint(solver)

    tracking_weights = [
            *((v, 1.) for v in nmpc.controller.vectors.differential[:,0]),
            *((v, 1.) for v in nmpc.controller.vectors.input[:,0]),
            ]

    nmpc.controller.add_tracking_objective(tracking_weights)

    nmpc.controller.constrain_control_inputs_piecewise_constant()
    
    nmpc.controller.initialize_to_initial_conditions()
    
    nmpc.controller.vectors.input[...].unfix()
    nmpc.controller.vectors.input[:,0].fix()
    solver.solve(nmpc.controller, tee=True)

    cv = controller.mod.fs.cstr.control_volume
    variance = [
            (cv.material_holdup[0.0,'aq','S'], 0.2),
            (cv.material_holdup[0.0,'aq','E'], 0.05),
            (cv.material_holdup[0.0,'aq','C'], 0.1),
            (cv.material_holdup[0.0,'aq','P'], 0.05),
            (cv.energy_holdup[0.0,'aq'], 5.),
            (cv.volume[0.0], 0.05),
            ]
    nmpc.controller.set_variance(variance)
    measurement_variance = [v.variance for v in controller.measurement_vars]
    t0 = nmpc.controller.time.first()
    measurement_noise_bounds = [
            (0.0, var[t0].ub) for var in controller.measurement_vars
            ]

    mx = plant.mod.fs.mixer
    variance = [
            (mx.S_inlet_state[0.0].flow_vol, 0.02),
            (mx.E_inlet_state[0.0].flow_vol, 0.001),
            ]
    nmpc.plant.set_variance(variance)
    input_variance = [v.variance for v in plant.input_vars]
    t0 = nmpc.plant.time.first()
    input_noise_bounds = [(0.0, var[t0].ub) for var in plant.input_vars]

    c_ts = nmpc.controller.sample_points[1]
    p_ts = nmpc.plant.sample_points[1]
    inputs = controller.generate_inputs_at_time(c_ts)
    plant.inject_inputs(inputs)

    nmpc.plant.initialize_by_solving_elements(solver)

    for i in range(0,10):
        print('\nENTERING NMPC LOOP ITERATION %s\n' % i)
        measured = nmpc.plant.generate_measurements_at_time(p_ts)
        nmpc.plant.shift_back_one_sample()
        nmpc.plant.initialize_to_initial_conditions()
        measured = apply_noise_with_bounds(
                measured,
                measurement_variance,
                random.gauss,
                measurement_noise_bounds,
                )

        nmpc.controller.shift_back_one_sample()
        nmpc.controller.load_measurements(measured)

        solver.solve(nmpc.controller, tee=True)

        inputs = controller.generate_inputs_at_time(c_ts)
        inputs = apply_noise_with_bounds(
                inputs,
                input_variance,
                random.gauss,
                input_noise_bounds,
                )
        plant.inject_inputs(inputs)
        
        nmpc.plant.initialize_by_solving_elements(solver)


if __name__ == '__main__':
    main()

