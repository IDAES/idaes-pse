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
Example for Cappresse's module for NMPC.
"""
from idaes.dynamic.cappresse import nmpc
from idaes.dynamic.cappresse.nmpc import *
import idaes.logger as idaeslog
from idaes.dynamic.cappresse.tests.cstr_for_testing import make_model

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes'}
else:
    solver = None

def main():

    # This tests the same model constructed in the test_nmpc_constructor_1 file
    m_plant = make_model(horizon=6, ntfe=60, ntcp=2)
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
    m_steady = make_model(steady=True)
    sample_time = 0.5
    time_plant = m_plant.fs.time

    n_plant_samples = (time_plant.last()-time_plant.first())/sample_time
    assert n_plant_samples == int(n_plant_samples)
    n_plant_samples = int(n_plant_samples)

    plant_sample_points = [time_plant.first() + i*sample_time
                           for i in range(1, n_plant_samples)]
    for t in plant_sample_points:
        assert t in time_plant
    # Six samples per horizon, five elements per sample

    initial_plant_inputs = [m_plant.fs.mixer.S_inlet.flow_rate[0],
                            m_plant.fs.mixer.E_inlet.flow_rate[0]]
    
    nmpc = NMPCSim(m_plant.fs, m_controller.fs, initial_plant_inputs,
            solver=solver, outlvl=idaeslog.DEBUG, 
            sample_time=sample_time)
    
    p_mod = nmpc.p_mod
    c_mod = nmpc.c_mod
    
    set_point = [(c_mod.cstr.outlet.conc_mol[0, 'P'], 0.4),
                 (c_mod.cstr.outlet.conc_mol[0, 'S'], 0.0),
                 (c_mod.cstr.control_volume.energy_holdup[0, 'aq'], 300),
                 (c_mod.mixer.E_inlet.flow_rate[0], 0.1),
                 (c_mod.mixer.S_inlet.flow_rate[0], 2.0)]
    # Interestingly, this (st.st. set point) solve converges infeasible
    # if energy_holdup set point is not 300. (Needs higher weight?)

    weight_tolerance = 5e-7
    
    # Weight overwrite expects a list of (VarData, value) tuples
    # in the STEADY MODEL
    weight_overwrite = [(m_steady.fs.mixer.E_inlet.flow_rate[0], 20.0)]
    
    nmpc.add_setpoint(set_point,
            steady_model=m_steady.fs,
            outlvl=idaeslog.DEBUG,
            steady_weight_overwrite=weight_overwrite,
            steady_weight_tol=weight_tolerance)
    
    nmpc.add_pwc_constraints()
    
    nmpc.initialize_control_problem(strategy='initial_conditions')
    
    nmpc.solve_control_problem()

    nmpc.inject_inputs_into_plant(time_plant.first(),
                                  add_noise=True)

    nmpc.simulate_plant(time_plant.first())

    for t in plant_sample_points:
        nmpc.load_initial_conditions_from_plant(t, 
                                                add_noise=True)

        nmpc.initialize_control_problem(strategy='from_previous')

        nmpc.solve_control_problem()

        nmpc.inject_inputs_into_plant(t,
                                      add_noise=True)
        
        nmpc.simulate_plant(t)

if __name__ == '__main__':
    main()

