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
Case B31A - Natural Gas Combined Cycle (NGCC) plant_gross_power
Reference: NETL Baseline Report Rev 4
Author: A. Deshpande and M. Zamarripa
"""
import pytest
from idaes.power_generation.costing.power_plant_costing import \
     (get_PP_costing,
      build_flowsheet_cost_constraint,
      costing_initialization)
from idaes.core.util.model_statistics import degrees_of_freedom
import pyomo.environ as pyo
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver

# Get default solver for testing
solver = get_solver()


@pytest.mark.unit
def test_get_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.get_costing(year='2018')

    # Accounts with Feedwater Flow to HP section of HRSG, as the
    # reference/scaling parameter - Exhibit 5-15
    FW_accounts = ['3.1', '3.3', '8.4']

    m.fs.unit = pyo.Block()
    m.fs.unit.feedwater_flowrate = pyo.Var(initialize=1085751)  # lb/hr
    m.fs.unit.feedwater_flowrate.fix()
    get_PP_costing(m.fs.unit, FW_accounts,
                   m.fs.unit.feedwater_flowrate, 'lb/hr', 6)
    assert isinstance(m.fs.unit.costing.total_plant_cost, pyo.Var)
    assert hasattr(m.fs.unit.costing, "bare_erected_cost")
    assert isinstance(m.fs.unit.costing.total_plant_cost_eq, pyo.Constraint)
    assert degrees_of_freedom(m) == 0


@pytest.mark.solver
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_get_ngcc_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})
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
    # reference/scaling parameter - Exhibit 5-15
    FW_accounts = ['3.1', '3.3', '8.4']
    m.fs.b1 = pyo.Block()

    m.fs.b1.feedwater_flowrate = pyo.Var(initialize=1085751)  # lb/hr
    m.fs.b1.feedwater_flowrate.fix()
    get_PP_costing(m.fs.b1, FW_accounts,
                   m.fs.b1.feedwater_flowrate, 'lb/hr', 6)

    # Accounts with Raw water withdrawal as the reference/scaling parameter
    # Exhibit 5-14
    RW_withdraw_accounts = ['3.2', '3.4', '3.5', '9.5', '14.6']
    m.fs.b2 = pyo.Block()

    m.fs.b2.raw_water_withdrawal = pyo.Var(initialize=2902)  # gpm
    m.fs.b2.raw_water_withdrawal.fix()
    get_PP_costing(m.fs.b2, RW_withdraw_accounts,
                   m.fs.b2.raw_water_withdrawal, 'gpm', 6)

    # Accounts with fuel gas flowrate as the reference/scaling parameter
    # Exhibit 5-15 stream 2, Exhibit 5-8
    FuelG_accounts = ['3.6', '3.9', '6.1', '6.3', '6.4']
    m.fs.b3 = pyo.Block()
    # Obtain Fuel gas flowrate in acm
    fuelgas_value = 205630  # lb/hr

    m.fs.b3.fg_flowrate = pyo.Var(initialize=fuelgas_value)  # lb/hr
    m.fs.b3.fg_flowrate.fix()
    get_PP_costing(m.fs.b3, FuelG_accounts,
                   m.fs.b3.fg_flowrate, 'lb/hr', 6)

    # Accounts with process water discharge as the reference/scaling parameter
    # Exhibit 5-14
    PW_discharge_accounts = ['3.7']
    m.fs.b4 = pyo.Block()

    m.fs.b4.process_water_discharge = pyo.Var(initialize=657)  # gpm
    m.fs.b4.process_water_discharge.fix()
    get_PP_costing(m.fs.b4, PW_discharge_accounts,
                   m.fs.b4.process_water_discharge, 'gpm', 6)

    # Accounts with flue gas flowrate as the reference/scaling parameter
    # Exhibit 5-15 stream 3, Exhibit 5-8
    FG_accounts = ['7.6']
    m.fs.b5 = pyo.Block()
    # Obtain Flue gas flowrate in acm
    fluegas_value = (8658430/60)/0.025  # ft3/min

    m.fs.b5.fg_flowrate = pyo.Var(initialize=fluegas_value)  # ft3/min
    m.fs.b5.fg_flowrate.fix()
    get_PP_costing(m.fs.b5, FG_accounts,
                   m.fs.b5.fg_flowrate, 'acfm', 6)

    # Accounts with combustion turbine gross power as the reference/scaling
    # parameter
    # Exhibit 5-9
    CT_grosspower_accounts = ['6.5']
    m.fs.b6 = pyo.Block()
    # Obtain combustion turbine gross power in kW
    CT_gross_power = 477*1000  # kW

    m.fs.b6.ct_gross_power = pyo.Var(initialize=CT_gross_power)  # kW
    m.fs.b6.ct_gross_power.fix()
    get_PP_costing(m.fs.b6, CT_grosspower_accounts,
                   m.fs.b6.ct_gross_power, 'kW', 6)

    # Accounts with HRSG duty as the reference/scaling parameter
    # Exhibit 5-8, streams 3 and 4
    HRSG_duty_accounts = ['7.1', '7.2']
    m.fs.b7 = pyo.Block()
    # Obtain HRSG duty in MMBtu/hr, overall energy balance
    HRSG_duty = -(-538.1+277.1)*8658430/(10**6)  # MMBtu/hr

    m.fs.b7.hrsg_duty = pyo.Var(initialize=HRSG_duty)  # MMBtu/hr
    m.fs.b7.hrsg_duty.fix()
    get_PP_costing(m.fs.b7, HRSG_duty_accounts,
                   m.fs.b7.hrsg_duty, 'MMBtu/hr', 6)

    # Accounts with gas flow to stack as the reference/scaling parameter
    # Exhibit 5-8, stream 4
    Stack_flow_gas_accounts = ['7.3', '7.4', '7.5']
    m.fs.b8 = pyo.Block()
    # Obtain gas flowrate to stack in ft3/min
    stack_flow_gas = (8658430/60)/0.061  # ft3/min

    m.fs.b8.stack_flow_gas = pyo.Var(initialize=stack_flow_gas)  # ft3/min
    m.fs.b8.stack_flow_gas.fix()
    get_PP_costing(m.fs.b8, Stack_flow_gas_accounts,
                   m.fs.b8.stack_flow_gas, 'acfm', 6)

    # Accounts with steam turbine gross power as the reference/scaling
    # parameter
    # Exhibit 5-9
    Steam_turbine_gross_power_accounts = ['8.1', '8.2', '8.5', '14.3']
    m.fs.b9 = pyo.Block()
    # Obtain steam turbine gross power in kW
    ST_gross_power = 263*1000  # kW

    m.fs.b9.st_gross_power = pyo.Var(initialize=ST_gross_power)  # kW
    m.fs.b9.st_gross_power.fix()
    get_PP_costing(m.fs.b9, Steam_turbine_gross_power_accounts,
                   m.fs.b9.st_gross_power, 'kW', 6)

    # Accounts with condenser duty as the reference/scaling parameter
    # Exhibit 5-9
    Condenser_duty_accounts = ['8.3']
    m.fs.b10 = pyo.Block()
    # Obtain condenser duty in MMBtu/hr
    condenser_duty = 1332  # MMBtu/hr

    m.fs.b10.cond_duty = pyo.Var(initialize=condenser_duty)  # MMBtu/hr
    m.fs.b10.cond_duty.fix()
    get_PP_costing(m.fs.b10, Condenser_duty_accounts,
                   m.fs.b10.cond_duty, 'MMBtu/hr', 6)

    # Accounts with cooling tower duty as the reference/scaling parameter
    # Exhibit 5-16
    Cooling_tower_accounts = ['9.1']
    m.fs.b11 = pyo.Block()
    # Obtain cooling tower duty in MMBtu/hr (includes condenser, Acid gas
    # removal, and other cooling loads)
    cooling_tower_duty = 1357  # MMBtu/hr

    m.fs.b11.cool_tower_duty = pyo.Var(initialize=cooling_tower_duty)
    # MMBtu/hr
    m.fs.b11.cool_tower_duty.fix()
    get_PP_costing(m.fs.b11, Cooling_tower_accounts,
                   m.fs.b11.cool_tower_duty, 'MMBtu/hr', 6)

    # Accounts with circulating water flowrate as the reference/scaling
    # parameter
    Circ_water_accounts = ['9.2', '9.3', '9.4', '9.6', '9.7', '14.5']
    m.fs.b12 = pyo.Block()
    # Obtain circulating water flowrate in gpm
    cir_water_flowrate = 217555  # gpm

    m.fs.b12.circ_water_flow = pyo.Var(initialize=cir_water_flowrate)  # gpm
    m.fs.b12.circ_water_flow.fix()
    get_PP_costing(m.fs.b12, Circ_water_accounts,
                   m.fs.b12.circ_water_flow, 'gpm', 6)

    # Accounts with total plant gross power as the reference/scaling parameter
    # Exhibit 5-9
    plant_gross_power_accounts = ['11.1', '11.7', '11.9', '13.1', '13.2',
                                  '13.3', '14.4', '14.7', '14.8', '14.9',
                                  '14.10']
    m.fs.b13 = pyo.Block()
    # Obtain total plant gross power in kW
    plant_gross_power = 740000  # kW

    m.fs.b13.gross_power = pyo.Var(initialize=plant_gross_power)  # kW
    m.fs.b13.gross_power.fix()
    get_PP_costing(m.fs.b13, plant_gross_power_accounts,
                   m.fs.b13.gross_power, 'kW', 6)

    # Accounts with auxilliary load as the reference/scaling parameter
    # Exhibit 5-9
    auxilliary_load_accounts = ['11.2', '11.3', '11.4', '11.5', '11.6', '12.1',
                                '12.2', '12.3', '12.4', '12.5', '12.6', '12.7',
                                '12.8', '12.9']
    m.fs.b14 = pyo.Block()
    # Obtain auxilliary load in kW
    aux_load = 14*1000  # kW

    m.fs.b14.auxilliary_load = pyo.Var(initialize=aux_load)  # kW
    m.fs.b14.auxilliary_load.fix()
    get_PP_costing(m.fs.b14, auxilliary_load_accounts,
                   m.fs.b14.auxilliary_load, 'kW', 6)

    # Accounts with STG,CTG output as the reference/scaling parameter
    # Case B31A Account 9 - Pg 501 rev 4 baseline report
    stg_ctg_accounts = ['11.8']
    m.fs.b15 = pyo.Block()
    # Obtain STG,CTG output in kW
    stg_ctg_op = 689800  # kW

    m.fs.b15.stg_ctg_output = pyo.Var(initialize=stg_ctg_op)  # kW
    m.fs.b15.stg_ctg_output.fix()
    get_PP_costing(m.fs.b15, stg_ctg_accounts,
                   m.fs.b15.stg_ctg_output, 'kW', 6)

    # Accounts with gas turbine power as the reference/scaling parameter
    # Exhibit 5-9
    gasturbine_accounts = ['14.1']
    m.fs.b16 = pyo.Block()
    # Obtain gas turbine power in kW
    gt_power = 477*1000  # kW

    m.fs.b16.gas_turbine_power = pyo.Var(initialize=gt_power)  # kW
    m.fs.b16.gas_turbine_power.fix()
    get_PP_costing(m.fs.b16, gasturbine_accounts,
                   m.fs.b16.gas_turbine_power, 'kW', 6)

    # Build cost constraints
    build_flowsheet_cost_constraint(m)

    # Initialize costing
    costing_initialization(m.fs)
    assert degrees_of_freedom(m) == 0

    # Solve the model
    results = solver.solve(m, tee=True)
    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok

    # Obtain the total plant costs for each account
    TPC_components = []
    TPC_components.append(sum(pyo.value(m.fs.b1.costing.total_plant_cost[ac])
                              for ac in FW_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b2.costing.total_plant_cost[ac])
                              for ac in RW_withdraw_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b3.costing.total_plant_cost[ac])
                              for ac in FuelG_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b4.costing.total_plant_cost[ac])
                              for ac in PW_discharge_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b5.costing.total_plant_cost[ac])
                              for ac in FG_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b6.costing.total_plant_cost[ac])
                              for ac in CT_grosspower_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b7.costing.total_plant_cost[ac])
                              for ac in HRSG_duty_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b8.costing.total_plant_cost[ac])
                              for ac in Stack_flow_gas_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b9.costing.total_plant_cost[ac])
                              for ac in Steam_turbine_gross_power_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b10.costing.total_plant_cost[ac])
                              for ac in Condenser_duty_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b11.costing.total_plant_cost[ac])
                              for ac in Cooling_tower_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b12.costing.total_plant_cost[ac])
                              for ac in Circ_water_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b13.costing.total_plant_cost[ac])
                              for ac in plant_gross_power_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b14.costing.total_plant_cost[ac])
                              for ac in auxilliary_load_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b15.costing.total_plant_cost[ac])
                              for ac in stg_ctg_accounts))
    TPC_components.append(sum(pyo.value(m.fs.b16.costing.total_plant_cost[ac])
                              for ac in gasturbine_accounts))

    assert pytest.approx(sum(TPC_components), abs=0.5) == 574.85

    #  Testing cost component values for some accounts against the
    #  NETL baseline report revision 4

    # Accounts with raw water withdrawal as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b2.costing.total_plant_cost[ac])
                             for ac in RW_withdraw_accounts), abs=0.5) \
        == 26.435

    # Accounts with fuel gas as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b3.costing.total_plant_cost[ac])
                             for ac in FuelG_accounts), abs=0.5) \
        == 158.415

    # Accounts with process water discharge as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b4.costing.total_plant_cost[ac])
                             for ac in PW_discharge_accounts), abs=0.5) \
        == 11.608

    # Accounts with HRSG duty as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b7.costing.total_plant_cost[ac])
                             for ac in HRSG_duty_accounts), abs=0.5) \
        == 90.794

    # Accounts with condenser duty as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b10.costing.total_plant_cost[ac])
                             for ac in Condenser_duty_accounts), abs=0.5) \
        == 14.27

    # Accounts with cooling tower duty as reference parameter
    assert pytest.approx(sum(pyo.value(m.fs.b11.costing.total_plant_cost[ac])
                             for ac in Cooling_tower_accounts), abs=0.5) \
        == 14.73
