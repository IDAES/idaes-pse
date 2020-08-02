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
from idaes.apps.caprese import NMPCSim, ControlInitOption
from pyomo.environ import SolverFactory
import idaes.logger as idaeslog
from idaes.apps.caprese.examples.cstr_model import make_model
import pandas as pd
import matplotlib.pyplot as plt

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes',
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
    m_plant = make_model(horizon=6, ntfe=60, ntcp=2)
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
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

    initial_plant_inputs = [m_plant.fs.mixer.S_inlet.flow_vol[0],
                            m_plant.fs.mixer.E_inlet.flow_vol[0]]
    
    nmpc = NMPCSim(plant_model=m_plant.fs, 
                   plant_time_set=m_plant.fs.time,
                   controller_model=m_controller.fs, 
                   controller_time_set=m_controller.fs.time,
                   inputs_at_t0=initial_plant_inputs,
                   solver=solver, outlvl=idaeslog.DEBUG,
                   sample_time=sample_time)

    plant = nmpc.plant
    controller = nmpc.controller

    nmpc.solve_consistent_initial_conditions(plant)
    nmpc.solve_consistent_initial_conditions(controller)
    
    set_point = [(controller.cstr.outlet.conc_mol[0, 'P'], 0.4),
                 (controller.cstr.outlet.conc_mol[0, 'S'], 0.0),
                 (controller.cstr.control_volume.energy_holdup[0, 'aq'], 300),
                 (controller.mixer.E_inlet.flow_vol[0], 0.1),
                 (controller.mixer.S_inlet.flow_vol[0], 2.0)]
    # Interestingly, this (st.st. set point) solve converges infeasible
    # if energy_holdup set point is not 300. (Needs higher weight?)

    weight_tolerance = 5e-7
    
    # Weight overwrite expects a list of (VarData, value) tuples
    # in the STEADY MODEL
    weight_override = [(controller.mixer.E_inlet.flow_vol[0], 20.0)]

    nmpc.calculate_full_state_setpoint(set_point,
            objective_weight_override=weight_override,
            objective_weight_tolerance=weight_tolerance,
            outlvl=idaeslog.DEBUG,
            allow_inconsistent=False,
            tolerance=1e-6)

    nmpc.add_setpoint_to_controller()
    
    nmpc.constrain_control_inputs_piecewise_constant()
    
    nmpc.initialize_control_problem(
            control_init_option=ControlInitOption.FROM_INITIAL_CONDITIONS)
    
    nmpc.solve_control_problem()

    nmpc.inject_control_inputs_into_plant(time_plant.first(),
                                  add_input_noise=True)

    nmpc.simulate_plant(time_plant.first())

    for t in plant_sample_points:
        nmpc.transfer_current_plant_state_to_controller(t,
                                                add_plant_noise=True)

        nmpc.initialize_control_problem(
                control_init_option=ControlInitOption.FROM_PREVIOUS)

        nmpc.solve_control_problem()

        nmpc.inject_control_inputs_into_plant(t,
                                      add_input_noise=True)
        
        nmpc.simulate_plant(t)

    # TODO: add option for specifying "user-interest variables"

    if plot_switch:
        temp_info = plant._NMPC_NAMESPACE.var_locator[
                plant.cstr.outlet.temperature[0.]]
        temp_location = temp_info.location
        temp_group = temp_info.group
        temperature_data = PlotData(temp_group, temp_location, name='Temperature')
        fig, ax = temperature_data.plot()
        fig.savefig(temperature_data.name)
    
        P_info = plant._NMPC_NAMESPACE.var_locator[
                plant.cstr.outlet.conc_mol[0.,'P']]
        P_location = P_info.location
        P_group = P_info.group
        P_data = PlotData(P_group, P_location, name='P_conc')
        fig, ax = P_data.plot()
        fig.savefig(P_data.name)
    
        S_info = plant._NMPC_NAMESPACE.var_locator[
                plant.cstr.outlet.conc_mol[0.,'S']]
        S_location = S_info.location
        S_group = S_info.group
        S_data = PlotData(S_group, S_location, name='S_conc')
        fig, ax = S_data.plot()
        fig.savefig(S_data.name)


if __name__ == '__main__':
    plot_switch = False
    main(plot_switch)

