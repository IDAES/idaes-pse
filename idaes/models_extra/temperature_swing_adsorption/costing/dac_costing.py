#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Costing models for a direct air capture (DAC) plant. The fixed-bed temperature
swing adsorption model used to represent the DAC process should be passed as an
argument to the costing functions.
"""

__author__ = "Daison Yancy Caballero"

import os
import json
import textwrap
from sys import stdout
from pandas import DataFrame

from pyomo.environ import units, Var, value
from pyomo.common.fileutils import this_file_dir
from pyomo.util.calc_var_value import calculate_variable_from_constraint

from idaes.core import UnitModelBlock, UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.tables import stream_table_dataframe_to_string

from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)

directory = this_file_dir()


def get_dac_costing(tsa, costing_case="electric_boiler"):
    if costing_case == "electric_boiler":
        _get_costing_electric_boiler(tsa)
    elif costing_case == "retrofit_ngcc":
        _get_costing_retrofit_ngcc(tsa)
    else:
        raise ConfigurationError("costing case not defined.")


def _get_costing_electric_boiler(tsa):

    # load custom costing parameters
    with open(
        os.path.join(directory, "costing_params_dac_electric_boiler.json"), "r"
    ) as f:
        costing_params = json.load(f)

    # get flowsheet
    fs = tsa.flowsheet()

    # create costing block
    fs.costing = QGESSCosting()

    # EPAT accounts for:
    # 1) Sorbent Handling, 2) Sorbent Preparation & Feed,
    # 3) Feedwater & Miscellaneous Bop Systems, 9) Cooling Water System,
    # 10) Ash/Spent Sorbent Handling System, 11) Accessory Electric Plant
    # 12) Instrumentation & Control, 13) Improvements to Site,
    # 14) Buildings & Structures, 15) DAC system

    # grouping accounts
    sorbent_makeup_accounts = [
        "1.5",
        "1.6",
        "1.7",
        "1.8",
        "1.9",
        "2.5",
        "2.6",
        "2.9",
        "10.6",
        "10.7",
        "10.9",
    ]
    raw_water_withdrawal_accounts = ["3.2", "3.4", "9.5", "14.6"]
    steam_accounts = ["3.1", "3.3", "3.5"]
    process_water_discharge_accounts = ["3.7"]
    cooling_tower_accounts = ["9.1"]
    circulating_water_accounts = ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"]
    total_auxiliary_load_accounts = [
        "11.1",
        "11.2",
        "11.3",
        "11.4",
        "11.5",
        "11.6",
        "11.7",
        "11.8",
        "11.9",
        "12.4",
        "12.5",
        "12.6",
        "12.7",
        "12.8",
        "12.9",
        "13.1",
        "13.2",
        "13.3",
        "14.4",
        "14.7",
        "14.8",
        "14.9",
        "14.10",
    ]
    dac_system_accounts = [
        "15.1",
        "15.2",
        "15.3",
        "15.4",
        "15.5",
        "15.6",
        "15.7",
        "15.8",
        "15.9",
    ]

    # reference parameters for accounts
    # sorbent makeup rate - from model
    sorbent_makeup_rate = (
        units.convert(tsa.bed_volume, to_units=units.ft**3)
        * (1 - tsa.bed_voidage)
        * tsa.number_beds
        / 0.5
        / 365
        / units.d
    )  # [ft^3/d]
    # CO2 product mass flow rate - from model
    CO2_product_mass_flow = units.convert(
        tsa.mw["CO2"] * tsa.flow_mol_co2_rich_stream[0, "CO2"],
        to_units=units.lb / units.hr,
    )  # [lb/hr]
    # raw water withdrawal flow rate - from surrogates
    raw_water_withdrawal = (
        (1.0496e-03 * CO2_product_mass_flow * units.hr / units.lb + 5.3359e01)
        * units.gal
        / units.min
    )  # [gpm]
    # process water discharge flow rate - from surrogates
    process_water_discharge = (
        (5.4007e-04 * CO2_product_mass_flow * units.hr / units.lb + 1.2000e01)
        * units.gal
        / units.min
    )  # [gpm]
    # cooling tower heat duty - from surrogates
    cooling_tower_duty = (
        (3.0797e-04 * CO2_product_mass_flow * units.hr / units.lb + 2.5000e01)
        * units.MBtu
        / units.hr
    )  # [MMBtu/hr]
    # circulating water flow_rate - from surrogates
    circulating_water_flow_rate = (
        (3.0797e-02 * CO2_product_mass_flow * units.hr / units.lb + 2.5000e03)
        * units.gal
        / units.min
    )  # [gpm]
    # mole flow rate of dac feed - from model (exhaust or air)
    flow_mol_inlet = units.convert(
        tsa.flow_mol_in_total, to_units=units.kmol / units.hr
    )  # [kmol/hr]
    # compressor auxiliary load - from surrogates
    product_compressor_auxiliary_load = (
        0.0012 * flow_mol_inlet * units.hr / units.kmol - 2.2798
    ) * units.kW  # [kW]
    # boiler auxiliary load - from surrogates
    boiler_auxiliary_load = (
        0.000335999
        * units.convert(tsa.flow_mass_steam, to_units=units.lb / units.hr)
        * units.hr
        / units.lb
        * 1e3
    ) * units.kW  # covert MW to kW
    # total auxiliary load
    # calculate with fans work + compressor for CO2 pure + boiler aux load
    total_auxiliary_load = (
        product_compressor_auxiliary_load
        + boiler_auxiliary_load
        + units.convert(tsa.compressor.unit.work_mechanical[0], to_units=units.kW)
    )
    # compressor aftercooler heat exchanger duty - from surrogates
    compressor_aftercooler_heat_duty = (
        (2e-6 * flow_mol_inlet * units.hr / units.kmol - 7e-8) * units.MBtu / units.hr
    )  # [MMBtu/hr]
    # auxiliary load for 2 beds (pressure changer from TSA model estimates
    # the power required to move air or exhaust gas to all beds)
    auxiliary_load_2_beds = (
        2
        * units.convert(tsa.compressor.unit.work_mechanical[0], to_units=units.kW)
        / tsa.number_beds
    )  # [kW]

    # sorbent makeup accounts
    fs.sorbent_makeup = UnitModelBlock()
    fs.sorbent_makeup.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": sorbent_makeup_accounts,
            "scaled_param": sorbent_makeup_rate,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # raw water withdrawal accounts
    fs.raw_water_withdrawal = UnitModelBlock()
    fs.raw_water_withdrawal.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": raw_water_withdrawal_accounts,
            "scaled_param": raw_water_withdrawal,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # steam accounts
    fs.steam = UnitModelBlock()
    fs.steam.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": steam_accounts,
            "scaled_param": units.convert(
                tsa.flow_mass_steam, to_units=units.lb / units.hr
            ),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # process water discharge accounts
    fs.process_water_discharge = UnitModelBlock()
    fs.process_water_discharge.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": process_water_discharge_accounts,
            "scaled_param": process_water_discharge,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # cooling tower accounts
    fs.cooling_tower = UnitModelBlock()
    fs.cooling_tower.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": cooling_tower_accounts,
            "scaled_param": cooling_tower_duty,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # circulating water accounts
    fs.circulating_water = UnitModelBlock()
    fs.circulating_water.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": circulating_water_accounts,
            "scaled_param": circulating_water_flow_rate,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # total plant auxiliary load accounts
    fs.total_auxiliary_load = UnitModelBlock()
    fs.total_auxiliary_load.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": total_auxiliary_load_accounts,
            "scaled_param": total_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # get EPAT costing accounts for 15) DAC system
    accounts_info = {}
    for k, v in costing_params["8"]["B"].items():
        if k in dac_system_accounts:
            accounts_info[k] = v["Account Name"].replace("DAC ", "")
    accounts = list(accounts_info.keys())

    # 15.1 - DAC adsorption/desorption vessels (cost of 120 vessels)
    vessels_account = [accounts[0]]
    fs.vessels = UnitModelBlock()
    fs.vessels.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": vessels_account,
            "scaled_param": units.convert(tsa.bed_volume, to_units=units.ft**3),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.2 - DAC CO2 Compression & Drying
    product_compression_account = [accounts[1]]
    fs.product_compression = UnitModelBlock()
    fs.product_compression.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": product_compression_account,
            "scaled_param": product_compressor_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.3 - DAC CO2 Compressor Aftercooler
    compressor_aftercooler_account = [accounts[2]]
    fs.compressor_aftercooler = UnitModelBlock()
    fs.compressor_aftercooler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": compressor_aftercooler_account,
            "scaled_param": compressor_aftercooler_heat_duty,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.4 - DAC System Air Handling Duct and Dampers
    duct_dampers_account = [accounts[3]]
    fs.duct_dampers = UnitModelBlock()
    fs.duct_dampers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": duct_dampers_account,
            "scaled_param": 2
            * units.convert(tsa.flow_mass_in_total_bed, to_units=units.lb / units.hr),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.5 - DAC System Air Handling Fans
    feed_fans_account = [accounts[4]]
    fs.feed_fans = UnitModelBlock()
    fs.feed_fans.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": feed_fans_account,
            "scaled_param": auxiliary_load_2_beds,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.6 - DAC Desorption Process Gas Handling System
    desorption_gas_handling_account = [accounts[5]]
    fs.desorption_gas_handling = UnitModelBlock()
    fs.desorption_gas_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": desorption_gas_handling_account,
            "scaled_param": CO2_product_mass_flow,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.7 - DAC Steam Distribution System
    steam_distribution_account = [accounts[6]]
    fs.steam_distribution = UnitModelBlock()
    fs.steam_distribution.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": steam_distribution_account,
            "scaled_param": units.convert(
                tsa.flow_mass_steam, to_units=units.lb / units.hr
            ),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.8 - DAC System Controls Equipment
    controls_equipment_account = [accounts[7]]
    fs.controls_equipment = UnitModelBlock()
    fs.controls_equipment.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": controls_equipment_account,
            "scaled_param": total_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.9 - Electric Boiler
    electric_boiler_account = [accounts[8]]
    fs.electric_boiler = UnitModelBlock()
    fs.electric_boiler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": electric_boiler_account,
            "scaled_param": units.convert(
                tsa.flow_mass_steam, to_units=units.lb / units.hr
            ),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # total plant cost
    TPC_list = {}
    for o in fs.component_objects(descend_into=True):
        # look for costing blocks
        if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
            for k in o.costing.total_plant_cost.keys():
                if k not in ["15.1", "15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k]
                if k in ["15.1"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] / 120 * tsa.number_beds
                if k in ["15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] * tsa.number_beds / 2

    # Total plant cost of dac unit
    @fs.costing.Expression(doc="total TPC for TSA system in $MM")
    def total_TPC(b):
        return sum(TPC_list.values())

    # build variable costs components

    # gallon/day of water required
    water_rate = (
        units.convert(tsa.flow_mol_in_total, to_units=units.kmol / units.hr)
        * units.hr
        / units.kmol
        * 60792.407
        / 1629629
        * units.gallon
        / units.day
    )
    # ton/day of water_chems required
    water_chems_rate = (
        units.convert(tsa.flow_mol_in_total, to_units=units.kmol / units.hr)
        * units.hr
        / units.kmol
        * 0.1811
        / 1629629
        * units.ton
        / units.day
    )
    # ft^3/day of sorbent required
    sorbent_rate = (
        units.convert(tsa.bed_volume, to_units=units.ft**3)
        * (1 - tsa.bed_voidage)
        * tsa.number_beds
        / 0.5
        / 365
        / units.day
    )
    # kW/day of aux_power required
    aux_power_rate = total_auxiliary_load * 24 / units.day

    # kg/day of steam required
    steam_rate = units.convert(tsa.flow_mass_steam, to_units=units.kg / units.day)

    # variables for rates
    fs.costing.net_power = Var(fs.time, initialize=690, units=units.MW)
    fs.costing.water = Var(fs.time, initialize=1.0, units=units.gallon / units.day)
    fs.costing.water_chems = Var(fs.time, initialize=1.0, units=units.ton / units.day)
    fs.costing.sorbent = Var(fs.time, initialize=1.0, units=units.ft**3 / units.day)
    fs.costing.aux_power = Var(fs.time, initialize=1.0, units=units.kW / units.day)
    fs.costing.steam = Var(fs.time, initialize=1.0, units=units.kg / units.day)

    # constraints for rates
    @fs.costing.Constraint(fs.time, doc="Equation for cost of water")
    def water_eq(b, t):
        return b.water[t] == water_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of water_chems")
    def water_chems_eq(b, t):
        return b.water_chems[t] == water_chems_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of sorbent")
    def sorbent_eq(b, t):
        return b.sorbent[t] == sorbent_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of aux_power")
    def aux_power_eq(b, t):
        return b.aux_power[t] == aux_power_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of steam")
    def steam_eq(b, t):
        return b.steam[t] == steam_rate

    fs.costing.steam_eq.deactivate()
    fs.costing.steam.fix(0.0)

    fs.costing.net_power.fix()

    # resources to be costed
    resources = [
        "water",
        "water_treatment_chemicals",
        "sorbent",
        "aux_power",
        "waste_sorbent",
        "IP_steam",
    ]

    # vars for resource consumption rates
    rates = [
        fs.costing.water,
        fs.costing.water_chems,
        fs.costing.sorbent,
        fs.costing.aux_power,
        fs.costing.sorbent,
        fs.costing.steam,
    ]

    # resource prices
    prices = {
        "sorbent": 201 * units.USD_2018 / units.ft**3,
        "aux_power": 0.06 * units.USD_2018 / units.kW,
        "waste_sorbent": 0.86 * units.USD_2018 / units.ft**3,
        "IP_steam": 0.00733 * units.USD_2018 / units.kg,
    }

    # land cost for total overnight cost
    @fs.costing.Expression(doc="Land Cost [$MM]")
    def land_cost_exp(b):
        return (156000 * (tsa.number_beds / 120) ** (0.78)) * 1e-6  # scaled to Millions

    # build fixed and variables O&M cost
    fs.costing.build_process_costs(
        net_power=fs.costing.net_power,
        # arguments related to fixed OM costs
        total_plant_cost=True,
        labor_rate=38.50,
        labor_burden=30,
        operators_per_shift=8,
        tech=6,
        fixed_OM=True,
        # arguments related owners costs
        variable_OM=True,
        land_cost=fs.costing.land_cost_exp,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
        tonne_CO2_capture=tsa.total_CO2_captured_year,
    )

    # electric boiler emissions
    @fs.Expression(doc="Electric Boiler Emissions [mol/s]")
    def emissions_electric_boiler(b):
        # eq. emissions for a US average electricity grid - kgCO2e/MWh
        equi_emissions = 425.021
        return total_auxiliary_load * 1e-3 * equi_emissions * 1000 / 3600 / 44.01

    @fs.Expression(doc="Electric Boiler Emissions, PV electricity grid [mol/s]")
    def emissions_electric_boiler_pv(b):
        # eq. emissions for solar PV electricity grid - kgCO2e/MWh
        equi_emissions = 48.472
        return total_auxiliary_load * 1e-3 * equi_emissions * 1000 / 3600 / 44.01

    # initialization
    for c in fs.costing._registered_unit_costing:
        for key in c.bare_erected_cost.keys():
            calculate_variable_from_constraint(
                c.bare_erected_cost[key],
                c.bare_erected_cost_eq[key],
            )
            calculate_variable_from_constraint(
                c.total_plant_cost[key],
                c.total_plant_cost_eq[key],
            )
    calculate_variable_from_constraint(fs.costing.water[0], fs.costing.water_eq[0])
    calculate_variable_from_constraint(
        fs.costing.water_chems[0], fs.costing.water_chems_eq[0]
    )
    calculate_variable_from_constraint(fs.costing.sorbent[0], fs.costing.sorbent_eq[0])
    calculate_variable_from_constraint(
        fs.costing.aux_power[0], fs.costing.aux_power_eq[0]
    )
    fs.costing.initialize_fixed_OM_costs()
    fs.costing.initialize_variable_OM_costs()


def _get_costing_retrofit_ngcc(tsa):

    # load custom costing parameters
    with open(
        os.path.join(directory, "costing_params_dac_retrofit_ngcc.json"), "r"
    ) as f:
        costing_params = json.load(f)

    # get flowsheet
    fs = tsa.flowsheet()

    # create costing block
    fs.costing = QGESSCosting()

    # EPAT accounts for: 1) Sorbent Handling,
    # 2) Sorbent Preparation & Feed,
    # 7) HRSG, Ductwork & Stack, 8) Steam Turbine & Accessories
    # 10) Ash/Spent Sorbent Handling System, 11) Accessory Electric Plant
    # 12) Instrumentation & Control, 15) DAC system

    # grouping accounts
    sorbent_makeup_accounts = [
        "1.5",
        "1.6",
        "1.7",
        "1.8",
        "1.9",
        "2.5",
        "2.6",
        "2.9",
        "10.6",
        "10.7",
        "10.9",
    ]
    gas_flow_to_dac_accounts = ["7.3"]
    steam_accounts = ["8.4"]
    total_auxiliary_load_accounts = [
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
    dac_system_accounts = [
        "15.1",
        "15.2",
        "15.3",
        "15.4",
        "15.5",
        "15.6",
        "15.7",
        "15.8",
    ]

    # reference parameters for accounts
    # sorbent makeup rate - from model
    sorbent_makeup_rate = (
        units.convert(tsa.bed_volume, to_units=units.ft**3)
        * (1 - tsa.bed_voidage)
        * tsa.number_beds
        / 0.5
        / 365
        / units.d
    )  # [ft^3/d]
    # CO2 product mass flow rate - from model
    CO2_product_mass_flow = units.convert(
        tsa.mw["CO2"] * tsa.flow_mol_co2_rich_stream[0, "CO2"],
        to_units=units.lb / units.hr,
    )  # [lb/hr]
    # mole flow rate of dac feed - from model (exhaust or air)
    flow_mol_inlet = units.convert(
        tsa.flow_mol_in_total, to_units=units.kmol / units.hr
    )  # [kmol/hr]
    # mass flow rate of dac feed - from model (exhaust or air)
    flow_mass_inlet = units.convert(
        tsa.flow_mass_in_total, to_units=units.lb / units.hr
    )  # [lb/hr]
    # compressor auxiliary load - from surrogates
    product_compressor_auxiliary_load = (
        0.0012 * flow_mol_inlet * units.hr / units.kmol - 2.2798
    ) * units.kW  # [kW]
    # total auxiliary load
    # calculate with fans work + compressor for CO2 pure
    total_auxiliary_load = product_compressor_auxiliary_load + units.convert(
        tsa.compressor.unit.work_mechanical[0], to_units=units.kW
    )
    # compressor aftercooler heat exchanger duty - from surrogates
    compressor_aftercooler_heat_duty = (
        (2e-6 * flow_mol_inlet * units.hr / units.kmol - 7e-8) * units.MBtu / units.hr
    )  # [MMBtu/hr]
    # auxiliary load for 2 beds (pressure changer from TSA model estimates
    # the power required to move air or exhaust gas to all beds)
    auxiliary_load_2_beds = (
        2
        * units.convert(tsa.compressor.unit.work_mechanical[0], to_units=units.kW)
        / tsa.number_beds
    )  # [kW]

    # sorbent makeup accounts
    fs.sorbent_makeup = UnitModelBlock()
    fs.sorbent_makeup.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": sorbent_makeup_accounts,
            "scaled_param": sorbent_makeup_rate,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # gas flow to dac accounts
    fs.gas_flow_to_dac = UnitModelBlock()
    fs.gas_flow_to_dac.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": gas_flow_to_dac_accounts,
            "scaled_param": flow_mass_inlet,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # steam accounts
    fs.steam = UnitModelBlock()
    fs.steam.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": steam_accounts,
            "scaled_param": units.convert(
                tsa.flow_mass_steam, to_units=units.lb / units.hr
            ),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # total plant auxiliary load accounts
    fs.total_auxiliary_load = UnitModelBlock()
    fs.total_auxiliary_load.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": total_auxiliary_load_accounts,
            "scaled_param": total_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # get EPAT costing accounts for 15) DAC system
    accounts_info = {}
    for k, v in costing_params["8"]["B"].items():
        if k in dac_system_accounts:
            accounts_info[k] = v["Account Name"].replace("DAC ", "")
    accounts = list(accounts_info.keys())

    # 15.1 - DAC adsorption/desorption vessels (cost of 120 vessels)
    vessels_account = [accounts[0]]
    fs.vessels = UnitModelBlock()
    fs.vessels.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": vessels_account,
            "scaled_param": units.convert(tsa.bed_volume, to_units=units.ft**3),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.2 - DAC CO2 Compression & Drying
    product_compression_account = [accounts[1]]
    fs.product_compression = UnitModelBlock()
    fs.product_compression.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": product_compression_account,
            "scaled_param": product_compressor_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.3 - DAC CO2 Compressor Aftercooler
    compressor_aftercooler_account = [accounts[2]]
    fs.compressor_aftercooler = UnitModelBlock()
    fs.compressor_aftercooler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": compressor_aftercooler_account,
            "scaled_param": compressor_aftercooler_heat_duty,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.4 - DAC System Air Handling Duct and Dampers
    duct_dampers_account = [accounts[3]]
    fs.duct_dampers = UnitModelBlock()
    fs.duct_dampers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": duct_dampers_account,
            "scaled_param": 2
            * units.convert(tsa.flow_mass_in_total_bed, to_units=units.lb / units.hr),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.5 - DAC System Air Handling Fans
    feed_fans_account = [accounts[4]]
    fs.feed_fans = UnitModelBlock()
    fs.feed_fans.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": feed_fans_account,
            "scaled_param": auxiliary_load_2_beds,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.6 - DAC Desorption Process Gas Handling System
    desorption_gas_handling_account = [accounts[5]]
    fs.desorption_gas_handling = UnitModelBlock()
    fs.desorption_gas_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": desorption_gas_handling_account,
            "scaled_param": CO2_product_mass_flow,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.7 - DAC Steam Distribution System
    steam_distribution_account = [accounts[6]]
    fs.steam_distribution = UnitModelBlock()
    fs.steam_distribution.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": steam_distribution_account,
            "scaled_param": units.convert(
                tsa.flow_mass_steam, to_units=units.lb / units.hr
            ),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.8 - DAC System Controls Equipment
    controls_equipment_account = [accounts[7]]
    fs.controls_equipment = UnitModelBlock()
    fs.controls_equipment.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": controls_equipment_account,
            "scaled_param": total_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # total plant cost
    TPC_list = {}
    for o in fs.component_objects(descend_into=True):
        # look for costing blocks
        if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
            for k in o.costing.total_plant_cost.keys():
                if k not in ["15.1", "15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k]
                if k in ["15.1"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] / 120 * tsa.number_beds
                if k in ["15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] * tsa.number_beds / 2

    # Total plant cost of dac unit
    @fs.costing.Expression(doc="total TPC for TSA system in $MM")
    def total_TPC(b):
        return sum(TPC_list.values())

    # build variable costs components

    # gallon/day of water required
    water_rate = (
        units.convert(tsa.flow_mol_in_total, to_units=units.kmol / units.hr)
        * units.hr
        / units.kmol
        * 60792.407
        / 1629629
        * units.gallon
        / units.day
    )
    # ton/day of water_chems required
    water_chems_rate = (
        units.convert(tsa.flow_mol_in_total, to_units=units.kmol / units.hr)
        * units.hr
        / units.kmol
        * 0.1811
        / 1629629
        * units.ton
        / units.day
    )
    # ft^3/day of sorbent required
    sorbent_rate = (
        units.convert(tsa.bed_volume, to_units=units.ft**3)
        * (1 - tsa.bed_voidage)
        * tsa.number_beds
        / 0.5
        / 365
        / units.day
    )
    # kW/day of aux_power required
    aux_power_rate = total_auxiliary_load * 24 / units.day

    # kg/day of steam required
    steam_rate = units.convert(tsa.flow_mass_steam, to_units=units.kg / units.day)

    # variables for rates
    fs.costing.net_power = Var(fs.time, initialize=690, units=units.MW)
    fs.costing.water = Var(fs.time, initialize=1.0, units=units.gallon / units.day)
    fs.costing.water_chems = Var(fs.time, initialize=1.0, units=units.ton / units.day)
    fs.costing.sorbent = Var(fs.time, initialize=1.0, units=units.ft**3 / units.day)
    fs.costing.aux_power = Var(fs.time, initialize=1.0, units=units.kW / units.day)
    fs.costing.steam = Var(fs.time, initialize=1.0, units=units.kg / units.day)

    # constraints for rates
    @fs.costing.Constraint(fs.time, doc="Equation for cost of water")
    def water_eq(b, t):
        return b.water[t] == water_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of water_chems")
    def water_chems_eq(b, t):
        return b.water_chems[t] == water_chems_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of sorbent")
    def sorbent_eq(b, t):
        return b.sorbent[t] == sorbent_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of aux_power")
    def aux_power_eq(b, t):
        return b.aux_power[t] == aux_power_rate

    @fs.costing.Constraint(fs.time, doc="Equation for cost of steam")
    def steam_eq(b, t):
        return b.steam[t] == steam_rate

    fs.costing.steam_eq.deactivate()
    fs.costing.steam.fix(0.0)

    fs.costing.aux_power_eq.deactivate()
    fs.costing.aux_power.fix(0.0)

    fs.costing.net_power.fix()

    # resources to be costed
    resources = [
        "water",
        "water_treatment_chemicals",
        "sorbent",
        "aux_power",
        "waste_sorbent",
        "IP_steam",
    ]

    # vars for resource consumption rates
    rates = [
        fs.costing.water,
        fs.costing.water_chems,
        fs.costing.sorbent,
        fs.costing.aux_power,
        fs.costing.sorbent,
        fs.costing.steam,
    ]

    # resource prices
    prices = {
        "sorbent": 201 * units.USD_2018 / units.ft**3,
        "aux_power": 0.06 * units.USD_2018 / units.kW,
        "waste_sorbent": 0.86 * units.USD_2018 / units.ft**3,
        "IP_steam": 0.00733 * units.USD_2018 / units.kg,
    }

    # land cost for total overnight cost
    @fs.costing.Expression(doc="Land Cost [$MM]")
    def land_cost_exp(b):
        return (156000 * (tsa.number_beds / 120) ** (0.78)) * 1e-6  # scaled to Millions

    # build fixed and variables O&M cost
    fs.costing.build_process_costs(
        net_power=fs.costing.net_power,
        # arguments related to fixed OM costs
        total_plant_cost=True,
        labor_rate=38.50,
        labor_burden=30,
        operators_per_shift=8,
        tech=6,
        fixed_OM=True,
        # arguments related owners costs
        variable_OM=True,
        land_cost=fs.costing.land_cost_exp,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
        tonne_CO2_capture=tsa.total_CO2_captured_year,
    )


def print_dac_costing(tsa):
    fs = tsa.flowsheet()

    TPC_list = {}
    for o in fs.component_objects(descend_into=True):
        # look for costing blocks
        if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
            for k in o.costing.total_plant_cost.keys():
                if k not in ["15.1", "15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k]
                if k in ["15.1"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] / 120 * tsa.number_beds
                if k in ["15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] * tsa.number_beds / 2

    for i, k in TPC_list.items():
        print(i, value(k))


def _var_dict_costing(tsa):

    # get flowsheet
    fs = tsa.flowsheet()

    # create dir with costing summary
    var_dict = {}

    var_dict["Annualized capital cost of dac unit [$MM/year]"] = value(
        fs.costing.annualized_cost
    )
    var_dict["Fixed O&M cost of dac unit [$MM/year]"] = value(
        fs.costing.total_fixed_OM_cost
    )
    var_dict["Variable O&M cost of dac unit [$MM/year]"] = value(
        fs.costing.total_variable_OM_cost[0]
    )
    var_dict["Total annualized cost of dac unit [$MM/year]"] = value(
        fs.costing.annualized_cost
        + fs.costing.total_fixed_OM_cost
        + fs.costing.total_variable_OM_cost[0] * fs.costing.capacity_factor
    )
    var_dict["Capture cost [$/tonne CO2]"] = value(fs.costing.cost_of_capture * 1e6)

    if hasattr(fs, "emissions_electric_boiler"):
        var_dict["Electric Boiler Emissions [mol/s]"] = value(
            fs.emissions_electric_boiler
        )

    if hasattr(fs, "emissions_electric_boiler_pv"):
        var_dict["Electric Boiler Emissions, PV electricity grid [mol/s]"] = value(
            fs.emissions_electric_boiler_pv
        )

    return var_dict


def dac_costing_summary(tsa, export=False):

    fs = tsa.flowsheet()

    if not hasattr(fs, "vessels"):
        raise ConfigurationError(f"{tsa.name} does not have any costing block.")

    var_dict = _var_dict_costing(tsa)

    summary_dir = {}
    summary_dir["Value"] = {}
    summary_dir["pos"] = {}

    count = 1
    for k, v in var_dict.items():
        summary_dir["Value"][k] = value(v)
        summary_dir["pos"][k] = count
        count += 1

    df = DataFrame.from_dict(summary_dir, orient="columns")
    del df["pos"]
    if export:
        df.to_csv(f"{tsa.local_name}_summary_costing.csv")

    print("\n" + "=" * 84)
    print(f"summary costing {tsa.local_name}")
    print("-" * 84)
    stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
    print("\n" + "=" * 84 + "\n")
