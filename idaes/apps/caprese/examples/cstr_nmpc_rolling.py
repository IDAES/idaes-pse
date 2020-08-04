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
from idaes.apps.caprese import NMPCSim, ControlInitOption, VariableCategory
from idaes.apps.caprese.common.config import PlantHorizonType
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
                      'halt_on_ampl_error': 'yes'}
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

    m_plant = make_model(horizon=0.8, ntfe=8, ntcp=2)
    # Here I am intentionally using a longer plant horizon than my sample
    # time. This behavior is supported, although not recommended unless
    # the plant is being used for multiple applications that require different
    # sampling times.
    m_controller = make_model(horizon=3, ntfe=30, ntcp=2, bounds=True)
    sample_time = 0.5
    time_plant = m_plant.fs.time

    simulation_length = 6
    samples_to_simulate = int(simulation_length/sample_time)

    initial_plant_inputs = [m_plant.fs.mixer.S_inlet.flow_vol[0],
                            m_plant.fs.mixer.E_inlet.flow_vol[0]]
    
    nmpc = NMPCSim(plant_model=m_plant.fs, 
                   plant_time_set=m_plant.fs.time,
                   controller_model=m_controller.fs, 
                   controller_time_set=m_controller.fs.time,
                   inputs_at_t0=initial_plant_inputs,
                   solver=solver, 
                   outlvl=idaeslog.DEBUG,
                   sample_time=sample_time,
                   plant_horizon_type=PlantHorizonType.ROLLING)

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
            outlvl=idaeslog.DEBUG)

    nmpc.add_setpoint_to_controller()
    
    nmpc.constrain_control_inputs_piecewise_constant()
    
    nmpc.initialize_control_problem(
            control_init_option=ControlInitOption.FROM_INITIAL_CONDITIONS)
    
    nmpc.solve_control_problem()

    planned_input_history = nmpc.initialize_input_history_from(controller, 
            t_real=0)

    nmpc.inject_control_inputs_into_plant(add_input_noise=True)
    # FIXME: Input noise is added to controller variables because of
    #        the silly way I generate noisy values. This is fixed in the
    #        MHE branch. Should port over to this one.

    nmpc.simulate_plant(time_plant.first())

    variables_of_interest = [
            VariableCategory.DIFFERENTIAL,
            VariableCategory.INPUT,
            plant.cstr.outlet.temperature[0],
            ]

    plant_history = nmpc.initialize_history_from_plant(0, 0.5,
            variables=variables_of_interest)

    actual_input_history = nmpc.initialize_input_history_from(plant, t_real=0)
    # FIXME: Addition of input noise appears to be broken.
    #        Not true, my implementation is just wonky.

    for i in range(samples_to_simulate):

        nmpc.transfer_current_plant_state_to_controller(time_plant.first(),
                                                add_plant_noise=True)

        nmpc.cycle_plant()
        
        nmpc.initialize_control_problem(
                control_init_option=ControlInitOption.FROM_PREVIOUS)

        nmpc.solve_control_problem()

        nmpc.extend_input_history_from(planned_input_history, controller)
        # FIXME: This breaks due to inconsistent inputs. Do I enforce a 
        #        pwc constraint at t == 0? No. Breaks due to my bad noise
        #        implementation. I apply noise, then back-shift.

        nmpc.inject_control_inputs_into_plant(time_plant.first(),
                                      add_input_noise=True)
        
        nmpc.simulate_plant(time_plant.first())

        nmpc.extend_input_history_from(actual_input_history, plant)

        nmpc.extend_history_from_plant(plant_history)
    import pdb; pdb.set_trace()

    # get outlet temperature series from plant history
    # utility function: 
    #   temperature_series = data_series_from_var(
    #               plant_history, 
    #               time,
    #               plant.cstr.outlet.temperature[0],
    #               )
    #   nmpc.get_series(history, temperature[0])
    #   history.from_var(plant.cstr.outlet.temperature[0], time)
    # plot it.


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

