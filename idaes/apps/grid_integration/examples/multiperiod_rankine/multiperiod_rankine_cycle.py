#This script provides functions to create a multiperiod rankine cycle model
import pyomo.environ as pyo
import numpy as np
from random import random
from idaes.apps.grid_integration.multiperiod.multiperiod import MultiPeriodModel
from idaes.apps.rankine.simple_rankine_cycle import create_model, set_inputs, initialize_model, close_flowsheet_loop, add_operating_cost

#Create a steady-state ranking cycle model, not yet setup for multi-period
def create_ss_rankine_model():
    p_lower_bound = 175
    p_upper_bound = 450

    m = pyo.ConcreteModel()
    m.rankine = create_model(heat_recovery=True)
    m.rankine = set_inputs(m.rankine)
    m.rankine = initialize_model(m.rankine)
    m.rankine = close_flowsheet_loop(m.rankine)
    m.rankine = add_operating_cost(m.rankine)

    # Setting bounds for net cycle power output for a capex design
    m.rankine.fs.eq_min_power = pyo.Constraint(
        expr=m.rankine.fs.net_cycle_power_output >= p_lower_bound*1e6)

    m.rankine.fs.eq_max_power = pyo.Constraint(
        expr=m.rankine.fs.net_cycle_power_output <= p_upper_bound*1e6)

    m.rankine.fs.boiler.inlet.flow_mol[0].unfix()
    m.rankine.fs.boiler.inlet.flow_mol[0].setlb(1)

    return m

#Create a multiperiod capable steady-state rankine cycle model. This is a
#user-provided function to a MultiPeriod class
def create_mp_rankine_block():
    m = create_ss_rankine_model()
    turbine_ramp_rate = 100
    battery_ramp_rate = 50
    b1 = m.rankine

    #Add coupling variable (next_power_output) to represent power output in next time period
    b1.power_output = pyo.Expression(expr = b1.fs.net_cycle_power_output*1e-6)  #MW
    b1.next_power_output = pyo.Var(within=pyo.NonNegativeReals, initialize=1.5) #MW

    #Use coupling variable to add ramping constraint
    b1.ramp1 = pyo.Constraint(expr=b1.power_output - b1.next_power_output <= turbine_ramp_rate)
    b1.ramp2 = pyo.Constraint(expr=b1.next_power_output - b1.power_output <= turbine_ramp_rate)

    #Add battery integration to rankine cycle
    b1.P_to_battery = pyo.Var(within=pyo.NonNegativeReals,initialize = 0.0)
    b1.P_to_grid = pyo.Var(within=pyo.NonNegativeReals,initialize = 0.0)
    b1.P_total = pyo.Constraint(expr = b1.power_output == b1.P_to_battery + b1.P_to_grid)

    #Simple battery model (soc = state of charge). We create a coupling variable called next_soc
    m.battery = pyo.Block()
    b2=m.battery

    #soc = state of charge
    b2.soc = pyo.Var(within=pyo.NonNegativeReals,initialize=0.0, bounds=(0,100))
    b2.next_soc = pyo.Var(within=pyo.NonNegativeReals,initialize=0.0, bounds=(0,100))

    #Amount discharged to grid this time period (assume discharge is positive)
    b2.efficiency = np.sqrt(0.88)
    b2.discharge = pyo.Var(initialize = 0.0)
    b2.energy_change = pyo.Constraint(expr = b2.next_soc == b2.soc - b2.discharge/b2.efficiency + b1.P_to_battery*b2.efficiency)
    b2.energy_down_ramp = pyo.Constraint(expr = b2.soc - b2.next_soc <= battery_ramp_rate)
    b2.energy_up_ramp = pyo.Constraint(expr = b2.next_soc - b2.soc <= battery_ramp_rate)

    return m

#the power output and battery state are linked between time periods
def get_rankine_link_variable_pairs(t1,t2):
    """
        t1: current time block
        t2: next time block
    """
    return [(t1.rankine.next_power_output,t2.rankine.power_output),
            (t1.battery.next_soc,t2.battery.soc)]

#the final power output and battery state must be the same as the intial power output and battery state
def get_rankine_periodic_variable_pairs(t1,t2):
    """
        t1: final time block
        t2: first time block
    """
    return [(t1.battery.next_soc,t2.battery.soc)]


def create_multiperiod_rankine_model(n_time_points=4):
    """
        Create a multi-period rankine cycle object. This object contains a pyomo model
        with a block for each time instance.

        n_time_points: Number of time blocks to create
    """
    n_time_points = n_time_points
    mp_rankine = MultiPeriodModel(n_time_points, create_mp_rankine_block, get_rankine_link_variable_pairs)

    #create the multiperiod object
    mp_rankine.build_multi_period_model()
    return mp_rankine
