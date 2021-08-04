##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
Cost evaluation using IDAES Costing Framework for EMRE NGCC Model
Base Case : B31B - NETL Baseline Report Rev 4
Author: A. Deshpande, Alex Noring, M. Zamarripa
"""

# NGCC cost evaluation using IDAES Costing Framework

import pytest
import idaes
from idaes.power_generation.costing.power_plant_costing import \
     (get_PP_costing,
      build_flowsheet_cost_constraint,
      costing_initialization,
      get_fixed_OM_costs,
      get_variable_OM_costs,
      initialize_fixed_OM_costs,
      initialize_variable_OM_costs)
from idaes.core.util.model_statistics import (degrees_of_freedom)
import pyomo.environ as pyo
from pyomo.environ import *
from idaes.core import FlowsheetBlock
from pyomo.environ import units as pyunits
# **Important Correlations for Indirect Process Variables**

# Assume linear scaling of cooling tower duty based on condenser duty values of the current case and B31B
# Cooling tower duty = Cooling tower duty B31B * (Condenser duty)/(Condenser duty B31B)

# Assume linear scaling of circulating water flowrate based on LP turbine exhaust steam - current case and B31B
# Circulating water flowrate = Circulating water flowrate B31B * (LP turbine exhaust steam )/(LP turbine exhaust steam  B31B)


# Assume linear scaling of auxilliary load based on feedwater and circulating water flowrate - current case and B31B
# feedwater_circwater_flow = feedwater_flowrate + cir_water_flowrate
# Auxilliary Load = Auxilliary Load B31B * (feedwater_circwater_flow)/(feedwater_circwater_flow  B31B)


# Assume linear scaling of raw water withdrawal based on condenser duty values of the current case and B31B
# raw_water_withdrawal = raw_water_withdrawal B31B * (Condenser duty)/(Condenser duty B31B)


# Assume linear scaling of process water discharge based on condenser duty values of the current case and B31B
# process_water_discharge = process_water_discharge B31B * (Condenser duty)/(Condenser duty B31B)


def get_ngcc_costing(m, evaluate_cost=False):

    m.fs.get_costing(year='2018')

    # Accounts corresponding to each component of the NGCC plant
    # Reference - Exhibit 5-17, NETL Baseline Report Rev 4
    NGCC_accounts = ['3.1', '3.2', '3.3', '3.4', '3.5', '3.6', '3.7', '3.9',
                     '6.1', '6.3', '6.4', '6.5', '7.1', '7.2', '7.3', '7.4',
                     '7.5', '7.6', '8.1', '8.2', '8.3', '8.4', '8.5', '9.1',
                     '9.2', '9.3', '9.4', '9.5', '9.6', '9.7', '11.1', '11.2',
                     '11.3', '11.4', '11.5', '11.6', '11.7', '11.8', '11.9',
                     '12.1', '12.2', '12.3', '12.4', '12.5', '12.6', '12.7',
                     '12.8', '12.9', '13.1', '13.2', '13.3', '14.1', '14.3',
                     '14.4', '14.5', '14.6', '14.7', '14.8', '14.9', '14.10']

    # Accounts with Feedwater Flow to HP section of HRSG, as the
    # reference/scaling parameter

    FW_accounts = ['3.1', '3.3', '8.4']

    m.fs.b1 = pyo.Block()
    feedwater_flowrate = m.fs.HP_ECON1.side_1.properties_in[0].flow_mol *\
        18.01528*0.00220462*3600  # mol/s to lb/hr conversion

    if evaluate_cost:
        m.fs.b1.feedwaterflowrate = \
            pyo.Var(initialize=pyo.value(feedwater_flowrate))
        m.fs.b1.feedwaterflowrate.fix()
        get_PP_costing(m.fs.b1, FW_accounts,
                       m.fs.b1.feedwaterflowrate, 'lb/hr', 6)
    else:
        get_PP_costing(m.fs.b1, FW_accounts,
                       feedwater_flowrate, 'lb/hr', 6)

    # Accounts with Raw water withdrawal as the reference/scaling parameter
    RW_withdraw_accounts = ['3.2', '3.4', '3.5', '9.5', '14.6']

    m.fs.b2 = pyo.Block()

    condenser_duty = m.fs.main_condenser.heat_duty[0]*0.000947817*3600*1e-6
    # W converted to MMBtu/hr

    raw_water_withdrawal = 4773*condenser_duty/788  # gpm

    if evaluate_cost:
        m.fs.b2.rawwaterwithdrawal = \
            pyo.Var(initialize=pyo.value(raw_water_withdrawal))
        m.fs.b2.rawwaterwithdrawal.fix()
        get_PP_costing(m.fs.b2, RW_withdraw_accounts,
                       m.fs.b2.rawwaterwithdrawal, 'gpm', 6)
    else:
        get_PP_costing(m.fs.b2, RW_withdraw_accounts,
                       raw_water_withdrawal, 'gpm', 6)

    # Accounts with fuel gas flowrate as the reference/scaling parameter
    FuelG_accounts = ['3.6', '3.9', '6.1', '6.3', '6.4']

    m.fs.b3 = pyo.Block()

    # Obtain Fuel gas flowrate in lb/hr
    fuelgas_flowrate = m.fs.feed_fuel1.properties[0.0].flow_mol  # mol/s
    fuelgas_flowrate = \
        fuelgas_flowrate*pyo.value(m.fs.feed_fuel1.properties[0.0].mw)*1000\
        * 0.00220462*3600  # mol/s to lb/hr conversion

    if evaluate_cost:
        m.fs.b3.fuelgasflowrate = \
            pyo.Var(initialize=pyo.value(fuelgas_flowrate))
        m.fs.b3.fuelgasflowrate.fix()
        get_PP_costing(m.fs.b3, FuelG_accounts,
                       m.fs.b3.fuelgasflowrate, 'lb/hr', 6)
    else:
        get_PP_costing(m.fs.b3, FuelG_accounts,
                       fuelgas_flowrate, 'lb/hr', 6)

    # Accounts with process water discharge as the reference/scaling parameter
    PW_discharge_accounts = ['3.7']

    m.fs.b4 = pyo.Block()

    process_water_discharge = 1670*condenser_duty/788  # gpm

    if evaluate_cost:
        m.fs.b4.processwaterdischarge = \
            pyo.Var(initialize=pyo.value(process_water_discharge))
        m.fs.b4.processwaterdischarge.fix()
        get_PP_costing(m.fs.b4, PW_discharge_accounts,
                       m.fs.b4.processwaterdischarge, 'gpm', 6)
    else:
        get_PP_costing(m.fs.b4, PW_discharge_accounts,
                       process_water_discharge, 'gpm', 6)

    # Accounts with flue gas flowrate as the reference/scaling parameter
    FG_accounts = ['7.6']

    m.fs.b5 = pyo.Block()

    # Obtain Flue gas flowrate in acm

    fluegas_flowrate = m.fs.exhaust_1.properties[0.0].flow_vol*35.3147*60
    # m3/s to ft3/min conversion

    if evaluate_cost:
        m.fs.b5.fluegasflowrate = \
            pyo.Var(initialize=pyo.value(fluegas_flowrate))
        m.fs.b5.fluegasflowrate.fix()
        get_PP_costing(m.fs.b5, FG_accounts,
                       m.fs.b5.fluegasflowrate, 'acfm', 6)
    else:
        get_PP_costing(m.fs.b5, FG_accounts,
                       fluegas_flowrate, 'acfm', 6)

    # Accounts with combustion turbine gross power as the reference/scaling
    # parameter
    CT_grosspower_accounts = ['6.5']

    m.fs.b6 = pyo.Block()

    # Obtain combustion turbine gross power in kW

    CT_gross_power = -m.fs.gt_power[0]/1e3  # W to kW conversion

    if evaluate_cost:
        m.fs.b6.CTgrosspower = pyo.Var(initialize=pyo.value(CT_gross_power))
        m.fs.b6.CTgrosspower.fix()
        get_PP_costing(m.fs.b6, CT_grosspower_accounts,
                       m.fs.b6.CTgrosspower, 'kW', 6)
    else:
        get_PP_costing(m.fs.b6, CT_grosspower_accounts,
                       CT_gross_power, 'kW', 6)

    # Accounts with HRSG duty as the reference/scaling parameter
    HRSG_duty_accounts = ['7.1', '7.2']

    m.fs.b7 = pyo.Block()

    # Obtain HRSG duty in MMBtu/hr

    HRSG_duty = sum(m.fs.HP_SH4.side_2.properties_in[0.0].flow_mol_comp[k] for k in m.fs.prop_gas.component_list)*\
    m.fs.HP_SH4.side_2.properties_in[0.0].enth_mol - \
        sum(m.fs.LP_ECON.side_2.properties_out[0.0].flow_mol_comp[k] for k in m.fs.prop_gas.component_list)*\
            m.fs.LP_ECON.side_2.properties_in[0.0].enth_mol

    HRSG_duty = HRSG_duty*0.000947817*3600*1e-6  # W converted to MMBtu/hr

    if evaluate_cost:
        m.fs.b7.HRSGduty = pyo.Var(initialize=pyo.value(HRSG_duty))
        m.fs.b7.HRSGduty.fix()
        get_PP_costing(m.fs.b7, HRSG_duty_accounts,
                       m.fs.b7.HRSGduty, 'MMBtu/hr', 6)
    else:
        get_PP_costing(m.fs.b7, HRSG_duty_accounts,
                       HRSG_duty, 'MMBtu/hr', 6)

    # Accounts with gas flow to stack as the reference/scaling parameter
    Stack_flow_gas_accounts = ['7.3', '7.4', '7.5']

    m.fs.b8 = pyo.Block()

    # Obtain gas flowrate to stack in ft3/min

    stack_flow_gas = m.fs.LP_ECON.side_2.properties_out[0].flow_vol*35.3147*60
    # m3/s to ft3/min conversion

    if evaluate_cost:
        m.fs.b8.stackflowgas = pyo.Var(initialize=pyo.value(stack_flow_gas))
        m.fs.b8.stackflowgas.fix()
        get_PP_costing(m.fs.b8, Stack_flow_gas_accounts,
                       m.fs.b8.stackflowgas, 'acfm', 6)
    else:
        get_PP_costing(m.fs.b8, Stack_flow_gas_accounts,
                       stack_flow_gas, 'acfm', 6)

    # Accounts with steam turbine gross power as the reference/scaling
    # parameter
    Steam_turbine_gross_power_accounts = ['8.1', '8.2', '8.5', '14.3']

    m.fs.b9 = pyo.Block()

    # Obtain steam turbine gross power in kW

    ST_gross_power = -m.fs.steam_turbine.power[0]*1e-3  # W to kW conversion

    if evaluate_cost:
        m.fs.b9.STgrosspower = pyo.Var(initialize=pyo.value(ST_gross_power))
        m.fs.b9.STgrosspower.fix()
        get_PP_costing(m.fs.b9, Steam_turbine_gross_power_accounts,
                       m.fs.b9.STgrosspower, 'kW', 6)
    else:
        get_PP_costing(m.fs.b9, Steam_turbine_gross_power_accounts,
                       ST_gross_power, 'kW', 6)

    # Accounts with condenser duty as the reference/scaling parameter
    Condenser_duty_accounts = ['8.3']

    m.fs.b10 = pyo.Block()

    if evaluate_cost:
        m.fs.b10.condenserduty = pyo.Var(initialize=pyo.value(condenser_duty))
        m.fs.b10.condenserduty.fix()
        get_PP_costing(m.fs.b10, Condenser_duty_accounts,
                       m.fs.b10.condenserduty, 'MMBtu/hr', 6)
    else:
        get_PP_costing(m.fs.b10, Condenser_duty_accounts,
                       condenser_duty, 'MMBtu/hr', 6)

    # Accounts with cooling tower duty as the reference/scaling parameter
    Cooling_tower_accounts = ['9.1']

    m.fs.b11 = pyo.Block()

    # Obtain cooling tower duty in MMBtu/hr (includes condenser, Acid gas
    # removal, and other cooling loads)

    cooling_tower_duty = 2208*condenser_duty/788  # MMBtu/hr

    if evaluate_cost:
        m.fs.b11.coolingtowerduty = \
            pyo.Var(initialize=pyo.value(cooling_tower_duty))
        m.fs.b11.coolingtowerduty.fix()
        get_PP_costing(m.fs.b11, Cooling_tower_accounts,
                       m.fs.b11.coolingtowerduty, 'MMBtu/hr', 6)
    else:
        get_PP_costing(m.fs.b11, Cooling_tower_accounts,
                       cooling_tower_duty, 'MMBtu/hr', 6)

    # Accounts with circulating water flowrate as the reference/scaling
    # parameter

    Circ_water_accounts = ['9.2', '9.3', '9.4', '9.6', '9.7', '14.5']

    m.fs.b12 = pyo.Block()

    # Obtain circulating water flowrate in gpm

    exhaust_steam_LP_turbine = m.fs.steam_turbine.outlet_stage.control_volume.properties_out[0.0].flow_mol  # mol/s
    cir_water_flowrate = 217555*exhaust_steam_LP_turbine/((821097*453.592/18.01528)/3600)  # gpm

    if evaluate_cost:
        m.fs.b12.cirwaterflowrate = \
            pyo.Var(initialize=pyo.value(cir_water_flowrate))
        m.fs.b12.cirwaterflowrate.fix()
        get_PP_costing(m.fs.b12, Circ_water_accounts,
                       m.fs.b12.cirwaterflowrate, 'gpm', 6)
    else:
        get_PP_costing(m.fs.b12, Circ_water_accounts,
                       cir_water_flowrate, 'gpm', 6)

    # Accounts with total plant gross power as the reference/scaling parameter

    plant_gross_power_accounts = ['11.1', '11.7', '11.9', '13.1', '13.2',
                                  '13.3', '14.4', '14.7', '14.8', '14.9',
                                  '14.10']

    m.fs.b13 = pyo.Block()

    # Obtain total plant gross power in kW

    plant_gross_power = CT_gross_power + ST_gross_power  # kW

    if evaluate_cost:
        m.fs.b13.plantgrosspower = \
            pyo.Var(initialize=pyo.value(plant_gross_power))
        m.fs.b13.plantgrosspower.fix()
        get_PP_costing(m.fs.b13, plant_gross_power_accounts,
                       m.fs.b13.plantgrosspower, 'kW', 6)
    else:
        get_PP_costing(m.fs.b13, plant_gross_power_accounts,
                       plant_gross_power, 'kW', 6)

    # Accounts with auxilliary load as the reference/scaling parameter

    auxilliary_load_accounts = ['11.2', '11.3', '11.4', '11.5', '11.6', '12.1',
                                '12.2', '12.3', '12.4', '12.5', '12.6', '12.7',
                                '12.8', '12.9']

    m.fs.b14 = pyo.Block()

    # Obtain auxilliary load in kW

    feedwater_circwater_flow = feedwater_flowrate + cir_water_flowrate*8.33*60
    aux_load = 44*feedwater_circwater_flow/(1086734 + (217555*8.33*60))

    if evaluate_cost:
        m.fs.b14.auxload = pyo.Var(initialize=pyo.value(aux_load))
        m.fs.b14.auxload.fix()
        get_PP_costing(m.fs.b14, auxilliary_load_accounts,
                       m.fs.b14.auxload, 'kW', 6)
    else:
        get_PP_costing(m.fs.b14, auxilliary_load_accounts,
                       aux_load, 'kW', 6)

    # Accounts with STG,CTG output as the reference/scaling parameter

    stg_ctg_accounts = ['11.8']

    m.fs.b15 = pyo.Block()

    # Obtain STG,CTG output in kW

    stg_ctg_op = 689800  # kW

    m.fs.b15.stgctgop = pyo.Var(initialize=stg_ctg_op)
    m.fs.b15.stgctgop.fix()
    get_PP_costing(m.fs.b15, stg_ctg_accounts,
                   m.fs.b15.stgctgop, 'kW', 6)

    # Accounts with gas turbine power as the reference/scaling parameter

    gasturbine_accounts = ['14.1']

    m.fs.b16 = pyo.Block()

    # Obtain gas turbine power in kW

    gt_power = CT_gross_power  # kW

    if evaluate_cost:
        m.fs.b16.gtpower = pyo.Var(initialize=pyo.value(gt_power))
        m.fs.b16.gtpower.fix()
        get_PP_costing(m.fs.b16, gasturbine_accounts,
                       m.fs.b16.gtpower, 'kW', 6)
    else:
        get_PP_costing(m.fs.b16, gasturbine_accounts,
                       gt_power, 'kW', 6)

    # Build cost constraints
    build_flowsheet_cost_constraint(m)

    # Initialize costing
    costing_initialization(m.fs)

    # if not evaluate_cost:
        # Solve the model
    # solver = pyo.SolverFactory('ipopt')
    # solver.solve(m, tee=True)
    # else:
    #     pass

    # Obtain the total plant cost for the entire process
    # TPC = pyo.value(m.fs.flowsheet_cost)
    # print('\n\nThe total plant cost is {0} Million $'.format(TPC))
    # return TPC


def ngcc_costing_validation(m):
    #  Testing cost component values for some accounts against the
    #  NETL baseline report revision 4, case B31B

    RW_withdraw_accounts = ['3.2', '3.4', '3.5', '9.5', '14.6']

    FuelG_accounts = ['3.6', '3.9', '6.1', '6.3', '6.4']

    PW_discharge_accounts = ['3.7']

    HRSG_duty_accounts = ['7.1', '7.2']

    Condenser_duty_accounts = ['8.3']

    Cooling_tower_accounts = ['9.1']

    # Accounts with raw water withdrawal as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b2.costing.total_plant_cost[ac])
                             for ac in RW_withdraw_accounts), rel=0.35) \
        == 37.094

    # Accounts with fuel gas as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b3.costing.total_plant_cost[ac])
                             for ac in FuelG_accounts), rel=0.35) \
        == 158.415

    # Accounts with process water discharge as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b4.costing.total_plant_cost[ac])
                             for ac in PW_discharge_accounts), rel=0.35) \
        == 22.512

    # Accounts with HRSG duty as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b7.costing.total_plant_cost[ac])
                             for ac in HRSG_duty_accounts), rel=0.35) \
        == 79.724

    # Accounts with condenser duty as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b10.costing.total_plant_cost[ac])
                             for ac in Condenser_duty_accounts), rel=0.35) \
        == 9.376

    # Accounts with cooling tower duty as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b11.costing.total_plant_cost[ac])
                             for ac in Cooling_tower_accounts), rel=0.35) \
        == 20.809


def build_ngcc_OM_costs(m):
    # testing fixed OM costs
    get_fixed_OM_costs(m, 650,  # MW
                       operators_per_shift=5,
                       tech=6,
                       fixed_TPC=555.35)  # MW
    initialize_fixed_OM_costs(m)

    # fixed cost is required to estimate variable O&M costs
    # testing variable OM costs
    # m.fs.power = pyo.Var(m.fs.time,
    #                      initialize={0: 650},
    #                      units=pyunits.MW)
    # m.fs.power.fix(650)
    # power = m.fs.net_power[0]/1e6*(-1) # MW
    m.fs.net_power_cost = pyo.Var(m.fs.time,
                            initialize=650,
                            units=pyunits.MW)
    m.fs.net_power_cost.fix()
    m.fs.natural_gas = pyo.Var(m.fs.time,
                              initialize={0: 7239},
                              units=pyunits.MBtu/pyunits.day)
    m.fs.natural_gas.fix(110955)
    # m.fs.fuel_lhv = pyo.Var() # J/kg
    # m.fs.fuel_lhv.fix(47.2e6)
    # m.fs.ng_preheater.tube.properties_in[0].flow_mass
    # lhv = pyo.value(m.fs.fuel_lhv/1e6) # MJ/kg
    # fuel = pyo.value(m.tags['fuel01_F']) #kg/s
    # m.fs.inject1.gas_state[t].flow_mass*b.fuel_lhv  J/kg * kg/s
    # nat gas needs to be in pyunits.MBtu
    # 1 J/s = 0.00094781712031332 BTU/s
    m.fs.water_use = pyo.Var(m.fs.time,
                             initialize={0: 2090*1000},
                             units=pyunits.gallon/pyunits.day)
    m.fs.water_use.fix()
    #
    get_variable_OM_costs(m, m.fs.net_power_cost,
                          ["natural gas", "water"],
                          [m.fs.natural_gas, m.fs.water_use])
    # get_variable_OM_costs(m, m.fs.net_power,
    #                       ["natural gas"],
    #                       [m.fs.natural_gas])
                          # [m.fs.flow_mass]

    initialize_variable_OM_costs(m)

if __name__ == "__main__":
    import pyomo.environ as pyo
    from idaes.core import FlowsheetBlock
    m = pyo.ConcreteModel("NGCC")
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.flow_mass = pyo.Var(m.fs.time,initialize=27.3,
                             units=pyunits.kg/pyunits.s)
    m.fs.lhv = pyo.Var(m.fs.time,initialize=47.2e6,
                    units=pyunits.J/pyunits.kg)
    # m.fs.net_power_cost = pyo.Var(m.fs.time,
    #                         initialize=650,
    #                         units=pyunits.MW)
    m.fs.flow_mass.fix()
    m.fs.lhv.fix()

    build_ngcc_OM_costs(m)
    @m.fs.Constraint(m.fs.time)
    def eq1(c, t):
        return m.fs.natural_gas[t] == pyunits.convert(m.fs.flow_mass[t]*m.fs.lhv[t], pyunits.MBtu/pyunits.day)
    m.fs.natural_gas.unfix()
    solver = pyo.SolverFactory('ipopt')
    results = solver.solve(m, tee=True)
    # m.fs.fuel_lhv*pyunits.J/pyunits*s, pyunits.MBtu)
