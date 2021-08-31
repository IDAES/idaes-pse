from simple_rankine_cycle import create_model, set_inputs, initialize_model, close_flowsheet_loop, add_operating_cost
from idaes.apps.multiperiod.multiperiod import MultiPeriodModel
import pyomo.environ as pyo
import numpy as np
from random import random

"""
    Example driver for multiperiod.py using a simple rankine cycle
"""

#Create a steady-state ranking cycle model, note yet multi-period capable
def create_ss_rankine_model():
    p_lower_bound = 300
    p_upper_bound = 450

    m = pyo.ConcreteModel()
    m.rankine = create_model(heat_recovery=True)
    m.rankine = set_inputs(m.rankine)
    m.rankine = initialize_model(m.rankine)
    m.rankine = close_flowsheet_loop(m.rankine)
    m.rankine = add_operating_cost(m.rankine)

    # Setting bounds for net cycle power output for the capex plant
    m.rankine.fs.eq_min_power = pyo.Constraint(
        expr=m.rankine.fs.net_cycle_power_output >= p_lower_bound*1e6)

    m.rankine.fs.eq_max_power = pyo.Constraint(
        expr=m.rankine.fs.net_cycle_power_output <= p_upper_bound*1e6)

    m.rankine.fs.boiler.inlet.flow_mol[0].unfix()
    m.rankine.fs.boiler.inlet.flow_mol[0].setlb(1)

    #Now add coupling information
    return m

#Create a multiperiod capable steady-state rankine cycle model. This is a
#user-provided function to a MultiPeriod class
def create_mp_rankine_block(lmp_signal):
    m = create_ss_rankine_model()
    b1 = m.rankine

    #Add new coupling variable for next time
    b1.power_output = pyo.Expression(expr = b1.fs.net_cycle_power_output*1e-6)
    b1.next_power_output = pyo.Var(within=pyo.NonNegativeReals,bounds = (300,450), initialize=1.5)

    #Add ramping constraint
    b1.ramp1 = pyo.Constraint(expr=b1.power_output - b1.next_power_output <= 100)
    b1.ramp2 = pyo.Constraint(expr=b1.next_power_output - b1.power_output <= 100)

    #Add battery integration
    b1.P_to_battery = pyo.Var(within=pyo.NonNegativeReals,initialize = 0.0)
    b1.P_to_grid = pyo.Var(within=pyo.NonNegativeReals,initialize = 0.0)
    b1.P_total = pyo.Constraint(expr = b1.power_output == b1.P_to_battery + b1.P_to_grid)

    #Simple battery model
    m.battery = pyo.Block()
    b2=m.battery

    b2.soc = pyo.Var(within=pyo.NonNegativeReals,initialize=0.0, bounds=(0,100))
    b2.next_soc = pyo.Var(within=pyo.NonNegativeReals,initialize=0.0, bounds=(0,100))

    #Amount discharged to grid this time period (assume discharge is positive)
    b2.discharge = pyo.Var(initialize = 0.0)
    b2.energy_change = pyo.Constraint(expr = b2.next_soc == b2.soc - b2.discharge + b1.P_to_battery)
    b2.energy_down_ramp = pyo.Constraint(expr = b2.soc - b2.next_soc <= 50)
    b2.energy_up_ramp = pyo.Constraint(expr = b2.next_soc - b2.soc <= 50)

    #Objective function
    m.lmp_signal = lmp_signal
    m.revenue = lmp_signal*(b1.P_to_grid + b2.discharge)
    m.cost = pyo.Expression(expr=-(m.revenue - b1.fs.operating_cost))
    return m

#Function to retrieve variable pairs linked between time steps. This is a
#user-provided function to a MultiPeriod class
def get_rankine_link_variable_pairs(b1,b2):
    return [(b1.rankine.next_power_output,b2.rankine.power_output),
            (b1.battery.next_soc,b2.battery.soc)]

#########################################
#check that a steady-state block solves
#########################################
b_ss = create_ss_rankine_model()
opt = pyo.SolverFactory('ipopt')
opt.solve(b_ss)
opt.solve(b_ss, tee=True)

#########################################
# Check that a multiperiod-capable block solves
#########################################
b_mp = create_mp_rankine_block(lmp_signal=30.0)
b_mp.obj = pyo.Objective(expr=b_mp.cost)
b_mp.battery.soc.fix(0)
opt = pyo.SolverFactory('ipopt')
opt.solve(b_mp)
opt.solve(b_mp, tee=True)

#########################################
#create and solve a multiperiod model
#########################################
n_time_points = 4
mp_rankine = MultiPeriodModel(n_time_points, create_mp_rankine_block, get_rankine_link_variable_pairs)#, state_variable_func=get_rankine_link_variable_pairs)

#data passed to each block time period
time_points = np.arange(0,n_time_points)
data_points = [{"lmp_signal":random()*100} for _ in range(n_time_points)]
data_kwargs = dict(zip(time_points,data_points))

#create the multiperiod object
mp_rankine.build_multi_period_model(data_kwargs)

#retrieve pyomo model and active process blocks (i.e. time blocks)
m = mp_rankine.pyomo_model
blks = mp_rankine.get_active_process_blocks()

#set objective function, fix initial state, and solve
m.obj = pyo.Objective(expr = sum([blk.cost for blk in blks]))
blks[0].battery.soc.fix(0)
opt = pyo.SolverFactory('ipopt')
opt.solve(m,tee=True)

#check results
for block in blks:
    print("Dispatch [MW]: ", pyo.value(block.rankine.power_output), "Charge [MWh]: ",pyo.value(block.battery.soc))

#########################################
#advance time steps
#########################################
mp_rankine.advance_time(lmp_signal=0.0)
blks = mp_rankine.get_active_process_blocks()

#Update data (objective price and battery storage level)
m.obj.expr = sum([blk.cost for blk in blks])
blks[0].battery.soc.fix(pyo.value(blks[0].battery.soc))
blks[0].fix_power = pyo.Constraint(expr = blks[0].rankine.power_output == pyo.value(blks[0].rankine.power_output))
opt.solve(m,tee=True)

#Check solution.
for block in blks:
    print("Dispatch [MW]: ", pyo.value(block.rankine.power_output), "Charge [MWh]: ",pyo.value(block.battery.soc))

#advance time again
mp_rankine.advance_time(lmp_signal=200.0)
blks = mp_rankine.get_active_process_blocks()
m.obj.expr = sum([blk.cost for blk in blks])
blks[0].battery.soc.fix(pyo.value(blks[0].battery.soc))
blks[0].fix_power = pyo.Constraint(expr = blks[0].rankine.power_output == pyo.value(blks[0].rankine.power_output))
opt.solve(m,tee=True)
for block in blks:
    print("Dispatch [MW]: ", pyo.value(block.rankine.power_output), "Charge [MWh]: ",pyo.value(block.battery.soc))
