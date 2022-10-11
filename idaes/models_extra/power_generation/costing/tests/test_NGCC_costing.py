#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Case B31A - Natural Gas Combined Cycle (NGCC) plant_gross_power
Reference: NETL-PUB-22638
Cost and Performance Baseline for Fossil Energy Plants Volume 1:
Bituminous Coal and Natural Gas to Electricity; https://doi.org/10.2172/1569246
Author: A. Deshpande and M. Zamarripa
"""
import pytest
from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)
from idaes.core.util.model_statistics import degrees_of_freedom
import pyomo.environ as pyo
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent

# Get default solver for testing
solver = get_solver()


@pytest.fixture(scope="module")
def build_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()
    CE_index_year = "2018"

    # Accounts with Feedwater Flow to HP section of HRSG, as the
    # reference/scaling parameter - Exhibit 5-15
    FW_accounts = ["3.1", "3.3", "8.4"]

    m.fs.b1 = UnitModelBlock()
    m.fs.b1.feedwater_flowrate = pyo.Var(
        initialize=1085751, units=pyunits.lb / pyunits.hr
    )
    m.fs.b1.feedwater_flowrate.fix()

    m.fs.b1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": FW_accounts,
            "scaled_param": m.fs.b1.feedwater_flowrate,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    return m


@pytest.mark.unit
def test_get_costing(build_costing):
    m = build_costing

    assert isinstance(m.fs.b1.costing.total_plant_cost, pyo.Var)
    assert hasattr(m.fs.b1.costing, "bare_erected_cost")
    assert isinstance(m.fs.b1.costing.total_plant_cost_eq, pyo.Constraint)
    assert degrees_of_freedom(m) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_units1_costing(build_costing):
    m = build_costing
    CE_index_year = "2018"

    # Accounts with Feedwater Flow to HP section of HRSG, as the
    # reference/scaling parameter - Exhibit 5-15
    FW_accounts = ["3.1", "3.3", "8.4"]

    # Accounts with Raw water withdrawal as the reference/scaling parameter
    # Exhibit 5-14
    RW_withdraw_accounts = ["3.2", "3.4", "3.5", "9.5", "14.6"]
    m.fs.b2 = UnitModelBlock()

    m.fs.b2.raw_water_withdrawal = pyo.Var(
        initialize=2902, units=pyunits.gal / pyunits.min
    )
    m.fs.b2.raw_water_withdrawal.fix()
    m.fs.b2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": RW_withdraw_accounts,
            "scaled_param": m.fs.b2.raw_water_withdrawal,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with fuel gas flowrate as the reference/scaling parameter
    # Exhibit 5-15 stream 2, Exhibit 5-8
    FuelG_accounts = ["3.6", "3.9", "6.1", "6.3", "6.4"]
    m.fs.b3 = UnitModelBlock()
    # Obtain Fuel gas flowrate in acm
    fuelgas_value = 205630  # lb/hr

    m.fs.b3.fg_flowrate = pyo.Var(
        initialize=fuelgas_value, units=pyunits.lb / pyunits.hr
    )
    m.fs.b3.fg_flowrate.fix()
    m.fs.b3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": FuelG_accounts,
            "scaled_param": m.fs.b3.fg_flowrate,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with process water discharge as the reference/scaling parameter
    # Exhibit 5-14
    PW_discharge_accounts = ["3.7"]
    m.fs.b4 = UnitModelBlock()

    m.fs.b4.process_water_discharge = pyo.Var(
        initialize=657, units=pyunits.gal / pyunits.min
    )
    m.fs.b4.process_water_discharge.fix()
    m.fs.b4.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": PW_discharge_accounts,
            "scaled_param": m.fs.b4.process_water_discharge,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Initialize costing
    QGESSCostingData.costing_initialization(m.fs.costing)
    assert degrees_of_freedom(m) == 0

    # Solve the model
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok

    # Accounts with raw water withdrawal as reference parameter
    assert pytest.approx(26.435, abs=0.5) == sum(
        pyo.value(m.fs.b2.costing.total_plant_cost[ac]) for ac in RW_withdraw_accounts
    )
    # Accounts with fuel gas as reference parameter
    assert pytest.approx(158.415, abs=0.5) == sum(
        pyo.value(m.fs.b3.costing.total_plant_cost[ac]) for ac in FuelG_accounts
    )

    # Accounts with process water discharge as reference parameter
    assert pytest.approx(11.608, abs=0.5) == sum(
        pyo.value(m.fs.b4.costing.total_plant_cost[ac]) for ac in PW_discharge_accounts
    )


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_units2_costing(build_costing):
    m = build_costing
    CE_index_year = "2018"
    # Accounts with flue gas flowrate as the reference/scaling parameter
    # Exhibit 5-15 stream 3, Exhibit 5-8
    FG_accounts = ["7.6"]
    m.fs.b5 = UnitModelBlock()
    # Obtain Flue gas flowrate in acm
    fluegas_value = (8658430 / 60) / 0.025  # ft3/min

    m.fs.b5.fg_flowrate = pyo.Var(
        initialize=fluegas_value, units=pyunits.ft**3 / pyunits.min
    )
    m.fs.b5.fg_flowrate.fix()
    m.fs.b5.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": FG_accounts,
            "scaled_param": m.fs.b5.fg_flowrate,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with combustion turbine gross power as the reference/scaling
    # parameter
    # Exhibit 5-9
    CT_grosspower_accounts = ["6.5"]
    m.fs.b6 = UnitModelBlock()
    # Obtain combustion turbine gross power in kW
    CT_gross_power = 477 * 1000  # kW

    m.fs.b6.ct_gross_power = pyo.Var(initialize=CT_gross_power, units=pyunits.kW)
    m.fs.b6.ct_gross_power.fix()
    m.fs.b6.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": CT_grosspower_accounts,
            "scaled_param": m.fs.b6.ct_gross_power,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with HRSG duty as the reference/scaling parameter
    # Exhibit 5-8, streams 3 and 4
    HRSG_duty_accounts = ["7.1", "7.2"]
    m.fs.b7 = UnitModelBlock()
    # Obtain HRSG duty in MMBtu/hr, overall energy balance
    HRSG_duty = -(-538.1 + 277.1) * 8658430 / (10**6)  # MMBtu/hr

    m.fs.b7.hrsg_duty = pyo.Var(initialize=HRSG_duty, units=pyunits.MBtu / pyunits.hr)
    m.fs.b7.hrsg_duty.fix()
    m.fs.b7.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": HRSG_duty_accounts,
            "scaled_param": m.fs.b7.hrsg_duty,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with gas flow to stack as the reference/scaling parameter
    # Exhibit 5-8, stream 4
    Stack_flow_gas_accounts = ["7.3", "7.4", "7.5"]
    m.fs.b8 = UnitModelBlock()
    # Obtain gas flowrate to stack in ft3/min
    stack_flow_gas = (8658430 / 60) / 0.061  # ft3/min

    m.fs.b8.stack_flow_gas = pyo.Var(
        initialize=stack_flow_gas, units=pyunits.ft**3 / pyunits.min
    )
    m.fs.b8.stack_flow_gas.fix()
    m.fs.b8.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": Stack_flow_gas_accounts,
            "scaled_param": m.fs.b8.stack_flow_gas,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Initialize costing
    QGESSCostingData.costing_initialization(m.fs.costing)
    assert degrees_of_freedom(m) == 0

    # Solve the model
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok
    # Accounts with HRSG duty as reference parameter
    assert pytest.approx(90.794, abs=0.1) == sum(
        pyo.value(m.fs.b7.costing.total_plant_cost[ac]) for ac in HRSG_duty_accounts
    )


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_units3_costing(build_costing):
    m = build_costing
    CE_index_year = "2018"
    # Accounts with steam turbine gross power as the reference/scaling
    # parameter
    # Exhibit 5-9
    Steam_turbine_gross_power_accounts = ["8.1", "8.2", "8.5", "14.3"]
    m.fs.b9 = UnitModelBlock()
    # Obtain steam turbine gross power in kW
    ST_gross_power = 263 * 1000  # kW

    m.fs.b9.st_gross_power = pyo.Var(initialize=ST_gross_power, units=pyunits.kW)
    m.fs.b9.st_gross_power.fix()
    m.fs.b9.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": Steam_turbine_gross_power_accounts,
            "scaled_param": m.fs.b9.st_gross_power,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with condenser duty as the reference/scaling parameter
    # Exhibit 5-9
    Condenser_duty_accounts = ["8.3"]
    m.fs.b10 = UnitModelBlock()
    # Obtain condenser duty in MMBtu/hr
    condenser_duty = 1332  # MMBtu/hr

    m.fs.b10.cond_duty = pyo.Var(
        initialize=condenser_duty, units=pyunits.MBtu / pyunits.hr
    )
    m.fs.b10.cond_duty.fix()
    m.fs.b10.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": Condenser_duty_accounts,
            "scaled_param": m.fs.b10.cond_duty,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with cooling tower duty as the reference/scaling parameter
    # Exhibit 5-16
    Cooling_tower_accounts = ["9.1"]
    m.fs.b11 = UnitModelBlock()
    # Obtain cooling tower duty in MMBtu/hr (includes condenser, Acid gas
    # removal, and other cooling loads)
    cooling_tower_duty = 1357  # MMBtu/hr

    m.fs.b11.cool_tower_duty = pyo.Var(
        initialize=cooling_tower_duty, units=pyunits.MBtu / pyunits.hr
    )
    m.fs.b11.cool_tower_duty.fix()
    m.fs.b11.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": Cooling_tower_accounts,
            "scaled_param": m.fs.b11.cool_tower_duty,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with circulating water flowrate as the reference/scaling
    # parameter
    Circ_water_accounts = ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"]
    m.fs.b12 = UnitModelBlock()
    # Obtain circulating water flowrate in gal/min
    cir_water_flowrate = 217555  # gal/min

    m.fs.b12.circ_water_flow = pyo.Var(
        initialize=cir_water_flowrate, units=pyunits.gal / pyunits.min
    )
    m.fs.b12.circ_water_flow.fix()
    m.fs.b12.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": Circ_water_accounts,
            "scaled_param": m.fs.b12.circ_water_flow,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with total plant gross power as the reference/scaling parameter
    # Exhibit 5-9
    plant_gross_power_accounts = [
        "11.1",
        "11.7",
        "11.9",
        "13.1",
        "13.2",
        "13.3",
        "14.4",
        "14.7",
        "14.8",
        "14.9",
        "14.10",
    ]
    m.fs.b13 = UnitModelBlock()
    # Obtain total plant gross power in kW
    plant_gross_power = 740000  # kW

    m.fs.b13.gross_power = pyo.Var(initialize=plant_gross_power, units=pyunits.kW)
    m.fs.b13.gross_power.fix()
    m.fs.b13.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": plant_gross_power_accounts,
            "scaled_param": m.fs.b13.gross_power,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with auxilliary load as the reference/scaling parameter
    # Exhibit 5-9
    auxilliary_load_accounts = [
        "11.2",
        "11.3",
        "11.4",
        "11.5",
        "11.6",
        "12.1",
        "12.2",
        "12.3",
        "12.4",
        "12.5",
        "12.6",
        "12.7",
        "12.8",
        "12.9",
    ]
    m.fs.b14 = UnitModelBlock()
    # Obtain auxilliary load in kW
    aux_load = 14 * 1000  # kW

    m.fs.b14.auxilliary_load = pyo.Var(initialize=aux_load, units=pyunits.kW)
    m.fs.b14.auxilliary_load.fix()
    m.fs.b14.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": auxilliary_load_accounts,
            "scaled_param": m.fs.b14.auxilliary_load,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with STG,CTG output as the reference/scaling parameter
    # Case B31A Account 9 - Pg 501 rev 4 baseline report
    stg_ctg_accounts = ["11.8"]
    m.fs.b15 = UnitModelBlock()
    # Obtain STG,CTG output in kW
    stg_ctg_op = 689800  # kW

    m.fs.b15.stg_ctg_output = pyo.Var(initialize=stg_ctg_op, units=pyunits.kW)
    m.fs.b15.stg_ctg_output.fix()
    m.fs.b15.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": stg_ctg_accounts,
            "scaled_param": m.fs.b15.stg_ctg_output,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Accounts with gas turbine power as the reference/scaling parameter
    # Exhibit 5-9
    gasturbine_accounts = ["14.1"]
    m.fs.b16 = UnitModelBlock()
    # Obtain gas turbine power in kW
    gt_power = 477 * 1000  # kW

    m.fs.b16.gas_turbine_power = pyo.Var(initialize=gt_power, units=pyunits.kW)
    m.fs.b16.gas_turbine_power.fix()
    m.fs.b16.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": gasturbine_accounts,
            "scaled_param": m.fs.b16.gas_turbine_power,
            "tech": 6,
            "ccs": "B",
            "CE_index_year": CE_index_year,
        },
    )

    # Initialize costing
    QGESSCostingData.costing_initialization(m.fs.costing)
    assert degrees_of_freedom(m) == 0

    # Solve the model
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok

    # Accounts with condenser duty as reference parameter
    assert pytest.approx(14.27, abs=0.1) == sum(
        pyo.value(m.fs.b10.costing.total_plant_cost[ac])
        for ac in Condenser_duty_accounts
    )

    # Accounts with cooling tower duty as reference parameter
    assert pytest.approx(14.73, abs=0.2) == sum(
        pyo.value(m.fs.b11.costing.total_plant_cost[ac])
        for ac in Cooling_tower_accounts
    )


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_flowsheet_costing(build_costing):
    m = build_costing
    # Build cost constraints
    m.fs.costing.build_process_costs(
        net_power=None,
        fixed_OM=False,
        variable_OM=False,
        resources=None,
        rates=None,
        prices=None,
        fuel="natural_gas",
    )

    # Initialize costing
    QGESSCostingData.costing_initialization(m.fs.costing)
    assert degrees_of_freedom(m) == 0

    # Solve the model
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal
    assert results.solver.status == pyo.SolverStatus.ok

    # Verify total plant costs
    assert pytest.approx(574.85, abs=0.1) == pyo.value(m.fs.costing.total_TPC)
