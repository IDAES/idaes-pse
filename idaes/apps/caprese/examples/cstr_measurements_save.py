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
from idaes.apps.caprese.nmpc import *
from idaes.apps.caprese.mhe import *
import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model
from idaes.apps.caprese.util import (InputHistory, InputRecord, search_history,
        histories_from_json)
import pandas as pd
import matplotlib.pyplot as plt
import json

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'linear_solver': 'ma57',
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes'}
else:
    solver = None

def main(plot_switch=False):

    # This tests the same model constructed in the test_nmpc_constructor_1 file
    m_plant = make_model(horizon=6, ntfe=60, ntcp=2)
    m_estimator = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
    sample_time = 0.5
    time_plant = m_plant.fs.time
    time_estimator = m_estimator.fs.time

    plant_horizon = time_plant.last() - time_plant.first()
    estimator_horizon = time_estimator.last() - time_estimator.first()
    n_plant_samples = plant_horizon/sample_time
    assert n_plant_samples == int(n_plant_samples)
    n_plant_samples = int(n_plant_samples)

    plant_sample_points = [time_plant.first() + i*sample_time
                           for i in range(1, n_plant_samples+1)]
    for t in plant_sample_points:
        assert t in time_plant
    # Six samples per horizon, five elements per sample

    initial_plant_inputs = [m_plant.fs.mixer.S_inlet.flow_vol[0],
                            m_plant.fs.mixer.E_inlet.flow_vol[0]]

    measurements_at_t0 = [
            m_plant.fs.cstr.control_volume.properties_out[0].conc_mol['S'],
            m_plant.fs.cstr.control_volume.properties_out[0].conc_mol['E'],
            m_plant.fs.cstr.control_volume.properties_out[0].conc_mol['P'],
            m_plant.fs.cstr.control_volume.properties_out[0].conc_mol['C'],
            m_plant.fs.cstr.control_volume.properties_out[0].temperature,
            ]

    mhe = MHESim(m_plant, time_plant, m_estimator, time_estimator,
            measurements_at_t0, inputs_at_t0=initial_plant_inputs)

    mhe.set_sample_time(sample_time)

    input_history_list = []
    for i in range(mhe.plant._MHE_NAMESPACE.n_input_vars):
        input_history_list.append(InputHistory(sample_time=sample_time))
    
    input_history = input_history_list[0]
    for i, t in enumerate(plant_sample_points):
        if i == 0:
            prev = time_plant.first()
        else:
            prev = plant_sample_points[i-1]
        if i % 2 == 0:
            value = 1.8
        else:
            value = 2.3
        start = prev
        end = t
        input_history.append(InputRecord(start=start, end=end, value=value))

    input_history = input_history_list[1]
    for i, t in enumerate(plant_sample_points):
        if i == 0:
            prev = time_plant.first()
        else:
            prev = plant_sample_points[i-1]
        if i % 2 == 0:
            value = 0.13
        else:
            value = 0.08
        start = prev
        end = t
        input_history.append(InputRecord(start=start, end=end, value=value))

    mhe.set_inputs_from_history(mhe.plant, input_history_list)

    print('Initializing plant...')
    mhe.initialize_plant(solve_initial_conditions=True)

    sigma_0 = 0.03
    variance_list = []
    for var in mhe.plant._MHE_NAMESPACE.measurement_vars:
        variance = sigma_0 * var[estimator_horizon].value
        variance_list.append((var[0], variance))

    mhe.set_measurement_variance(mhe.plant, variance_list)
    
    measurement_history = mhe.generate_measurement_history_from_plant(
            t_start=0, t_stop=estimator_horizon)

    filename = 'measurements.json'
    with open(filename, 'w') as f:
        f.write(json.dumps(measurement_history))

    new_measurement_history = histories_from_json(filename)

    print('loading measurements into estimator')
    mhe.generate_initial_measurements_from_history(new_measurement_history)

    import pdb; pdb.set_trace()


if __name__ == '__main__':
    plot_switch = False
    main(plot_switch)

