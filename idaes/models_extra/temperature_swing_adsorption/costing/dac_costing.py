#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
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

__author__ = "Daison Yancy Caballero, Ryan Hughes"

import os
import json
import textwrap
from sys import stdout
from pandas import DataFrame

from pyomo.environ import (
    units,
    value,
    Var,
)
from pyomo.common.fileutils import this_file_dir

from idaes.core import UnitModelBlock, UnitModelCostingBlock
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
from idaes.core.util.tables import stream_table_dataframe_to_string
from idaes.core.util.math import smooth_max

from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)

directory = this_file_dir()
_log = idaeslog.getLogger(__name__)


def get_dac_costing_data(case):
    if case == "electric_boiler":
        fname = "costing_params_dac_electric_boiler.json"
    elif case == "retrofit_ngcc":
        fname = "costing_params_dac_retrofit_ngcc.json"
    else:
        raise ConfigurationError("costing case not defined.")

    # load custom costing parameters
    with open(os.path.join(directory, fname), "r") as f:
        costing_params = json.load(f)

    return costing_params


def get_dac_costing(unit, costing_case):
    # capital costs
    costing_params = get_dac_costing_data(costing_case)

    # get flowsheet
    fs = unit.flowsheet()

    fs.costing = QGESSCosting()

    # reference parameters for accounts ==============================
    # === accounts present for both costing cases ===
    # compressor auxiliary load - from surrogates
    _pcal_dimless = (
        0.0012
        * units.hr
        / units.kmol
        * units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr)
        - 2.2798
    )
    product_compressor_auxiliary_load = smooth_max(0, _pcal_dimless) * units.kW  # [kW]
    # compressor aftercooler heat exchanger duty - from surrogates
    _cahd_dimless = (
        2e-6
        * units.hr
        / units.kmol
        * units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr)
        - 7e-8
    )
    compressor_aftercooler_heat_duty = (
        smooth_max(0, _cahd_dimless) * units.MBtu / units.hr
    )  # [MMBtu/hr]

    auxiliary_load_2_beds = (
        2
        * units.convert(unit.compressor.unit.work_mechanical[0], to_units=units.kW)
        / unit.number_beds
    )  # [kW]

    # CO2 product mass flow rate - from model
    CO2_product_mass_flow = units.convert(
        unit.mw["CO2"] * unit.flow_mol_co2_rich_stream[0, "CO2"],
        to_units=units.lb / units.hr,
    )  # [lb/hr]

    # assume most water was knocked out in vacuum (3% enters storage)
    CO2_mole_frac = 1 - 0.03
    kmol_p_hr_total = (
        units.convert(CO2_product_mass_flow, to_units=units.kg / units.hr)
        * units.kmol
        / (44.01 * units.kg)
        / CO2_mole_frac
    )
    dens = 1.8 * units.kg / units.m**3  # from Sorbent report stream table
    MW = 43.146 * units.kg / units.kmol  # from Sorbent report stream table
    CO2_storage_throughput = kmol_p_hr_total * MW / dens

    if costing_case == "electric_boiler":
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
        # boiler auxiliary load - from surrogates
        boiler_auxiliary_load = (
            0.000335999
            * units.convert(unit.flow_mass_steam, to_units=units.lb / units.hr)
            * units.hr
            / units.lb
            * 1e3
            * units.kW
        )  # convert MW to kW
        # total auxiliary load
        # calculate with fans work + compressor for CO2 pure + boiler aux load
        total_auxiliary_load = (
            product_compressor_auxiliary_load
            + boiler_auxiliary_load
            + units.convert(unit.compressor.unit.work_mechanical[0], to_units=units.kW)
        )

        fs.raw_water_system = UnitModelBlock()
        fs.raw_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.2", "3.4", "9.5", "14.6"],
                "scaled_param": raw_water_withdrawal,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        fs.steam_system = UnitModelBlock()
        fs.steam_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.1", "3.3", "3.5"],
                "scaled_param": units.convert(
                    unit.flow_mass_steam, to_units=units.lb / units.hr
                ),
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        fs.cooling_tower = UnitModelBlock()
        fs.cooling_tower.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.1"],
                "scaled_param": cooling_tower_duty,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        fs.water_discharge_system = UnitModelBlock()
        fs.water_discharge_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["3.7"],
                "scaled_param": process_water_discharge,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        fs.cooling_water_system = UnitModelBlock()
        fs.cooling_water_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"],
                "scaled_param": circulating_water_flow_rate,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        fs.electric_systems = UnitModelBlock()
        fs.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": [
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
                ],
                "scaled_param": total_auxiliary_load,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        # Electric Boiler 15.9
        fs.electric_boiler = UnitModelBlock()
        fs.electric_boiler.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["15.9"],
                "scaled_param": units.convert(
                    unit.flow_mass_steam, to_units=units.lb / units.hr
                ),
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

    elif costing_case == "retrofit_ngcc":
        # calculate with fans work + compressor for CO2 pure
        total_auxiliary_load = product_compressor_auxiliary_load + units.convert(
            unit.compressor.unit.work_mechanical[0], to_units=units.kW
        )

        # gas flow to dac accounts
        fs.gas_flow_to_dac = UnitModelBlock()
        fs.gas_flow_to_dac.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["7.3"],
                "scaled_param": units.convert(
                    unit.flow_mass_in_total, to_units=units.lb / units.hr
                ),
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        fs.steam_flow_system = UnitModelBlock()
        fs.steam_flow_system.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": ["8.4"],
                "scaled_param": units.convert(
                    unit.flow_mass_steam, to_units=units.lb / units.hr
                ),
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )

        fs.electric_systems = UnitModelBlock()
        fs.electric_systems.costing = UnitModelCostingBlock(
            flowsheet_costing_block=fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": [
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
                ],
                "scaled_param": total_auxiliary_load,
                "tech": 8,
                "ccs": "B",
                "additional_costing_params": costing_params,
            },
        )
    else:
        raise ConfigurationError("costing case not defined.")

    # sorbent makeup accounts
    sorbent_lifespan = 1

    sorbent_makeup_rate = (
        units.convert(unit.bed_volume, to_units=units.ft**3)
        * (1 - unit.bed_voidage)
        * unit.number_beds
        / sorbent_lifespan
        / 365
        / units.day
    )

    fs.sorbent_makeup = UnitModelBlock()
    fs.sorbent_makeup.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": [
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
            ],
            "scaled_param": sorbent_makeup_rate,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    fs.vessels = UnitModelBlock()
    fs.vessels.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.1"],
            "scaled_param": units.convert(unit.bed_volume, to_units=units.ft**3),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.2 - DAC CO2 Compression & Drying
    fs.product_compression = UnitModelBlock()
    fs.product_compression.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.2"],
            "scaled_param": product_compressor_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.3 - DAC CO2 Compressor Aftercooler
    fs.compressor_aftercooler = UnitModelBlock()
    fs.compressor_aftercooler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.3"],
            "scaled_param": compressor_aftercooler_heat_duty,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.4 - DAC System Air Handling Duct and Dampers (1 system per 2 beds)
    fs.duct_dampers = UnitModelBlock()
    fs.duct_dampers.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.4"],
            "scaled_param": 2
            * units.convert(unit.flow_mass_in_total_bed, to_units=units.lb / units.hr),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.5 - DAC System Air Handling Fans (1 system per 2 beds)
    fs.feed_fans = UnitModelBlock()
    fs.feed_fans.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.5"],
            "scaled_param": auxiliary_load_2_beds,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.6 - DAC Desorption Process Gas Handling System (pure CO2 gas)
    fs.desorption_gas_handling = UnitModelBlock()
    fs.desorption_gas_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.6"],
            "scaled_param": CO2_product_mass_flow,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.7 - DAC Steam Distribution System
    fs.steam_distribution = UnitModelBlock()
    fs.steam_distribution.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.7"],
            "scaled_param": units.convert(
                unit.flow_mass_steam, to_units=units.lb / units.hr
            ),
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.8 - DAC System Controls Equipment
    fs.controls_equipment = UnitModelBlock()
    fs.controls_equipment.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.8"],
            "scaled_param": total_auxiliary_load,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.10 - CO2 interim storage vessel
    fs.CO2_storage_vessel = UnitModelBlock()
    fs.CO2_storage_vessel.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.10"],
            "scaled_param": CO2_storage_throughput,
            "tech": 8,
            "ccs": "B",
            "additional_costing_params": costing_params,
        },
    )

    # 15.11 - DAC CO2 Dryer
    fs.CO2_dryer = UnitModelBlock()
    fs.CO2_dryer.costing = UnitModelCostingBlock(
        flowsheet_costing_block=fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ["15.11"],
            "scaled_param": CO2_product_mass_flow,
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
                    TPC_list[k] = o.costing.total_plant_cost[k] * unit.number_beds
                if k in ["15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] * unit.number_beds / 2

    # Total plant cost of dac unit
    @fs.costing.Expression(doc="total TPC for TSA system in $MM")
    def total_TPC(b):
        return sum(TPC_list.values())

    # resource consumption rates for variable costing
    capacity_factor = 0.85

    @fs.costing.Expression(fs.time)
    def sorbent_rate(b):
        return sorbent_makeup_rate

    @fs.costing.Expression(fs.time)
    def water_use(b, t):
        # water use is lineraly scaled from NETL reference
        ref_water = 60792.407 * units.gal / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr) * (
            ref_water / ref_air
        )

    @fs.costing.Expression(fs.time)
    def water_treatment_chems(b, t):
        # treatment chemical use is lineraly scaled from NETL reference
        ref_chem = 0.1811 * units.ton / units.day
        ref_air = 1629629 * units.kmol / units.hr
        return units.convert(unit.flow_mol_in_total, to_units=units.kmol / units.hr) * (
            ref_chem / ref_air
        )

    fs.costing.energy_purchased = Var(
        fs.time, initialize=1.0, units=units.kW * units.hr / units.day
    )
    fs.costing.steam_rate = Var(fs.time, initialize=1.0, units=units.kg / units.day)

    @fs.costing.Constraint(fs.time)
    def energy_purchased_eq(b, t):
        return b.energy_purchased[t] == units.convert(
            total_auxiliary_load, to_units=units.kW * units.hr / units.day
        )

    @fs.costing.Constraint(fs.time, doc="Equation for cost of steam")
    def steam_eq(b, t):
        return b.steam_rate[t] == units.convert(
            unit.flow_mass_steam, to_units=units.kg / units.day
        )

    if costing_case == "retrofit_ngcc":
        fs.costing.energy_purchased_eq.deactivate()
        fs.costing.energy_purchased.fix(0.0)

    fs.costing.steam_eq.deactivate()
    fs.costing.steam_rate.fix(0.0)

    fs.costing.net_power = Var(fs.time, initialize=690, units=units.MW)
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
        fs.costing.water_use,
        fs.costing.water_treatment_chems,
        fs.costing.sorbent_rate,
        fs.costing.energy_purchased,
        fs.costing.sorbent_rate,
        fs.costing.steam_rate,
    ]

    # resource prices
    prices = {
        "sorbent": 100 * units.USD_2018 / units.ft**3,
        "aux_power": 0.06 * units.USD_2018 / units.kWh,
        "waste_sorbent": 0.86 * units.USD_2018 / units.ft**3,
        "IP_steam": 0.00733 * units.USD_2018 / units.kg,
    }

    @fs.costing.Expression()
    def land_cost_exp(b):
        return (
            156000 * (unit.number_beds / 120) ** (0.78)
        ) * 1e-6  # scaled to Millions

    fs.costing.build_process_costs(
        # net_power=None,
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
        capacity_factor=capacity_factor,
        land_cost=fs.costing.land_cost_exp,
        resources=resources,
        rates=rates,
        prices=prices,
        tonne_CO2_capture=unit.total_CO2_captured_year,
    )

    fs.costing.costing_initialization()


def print_dac_costing(unit):
    fs = unit.flowsheet()

    TPC_list = {}
    for o in fs.component_objects(descend_into=True):
        # look for costing blocks
        if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
            for k in o.costing.total_plant_cost.keys():
                if k not in ["15.1", "15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k]
                if k in ["15.1"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] / 120 * unit.number_beds
                if k in ["15.4", "15.5"]:
                    TPC_list[k] = o.costing.total_plant_cost[k] * unit.number_beds / 2

    for i, k in TPC_list.items():
        print(i, value(k))


def _var_dict_costing(unit):

    # get flowsheet
    fs = unit.flowsheet()

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


def dac_costing_summary(unit, export=False):

    fs = unit.flowsheet()

    if not hasattr(fs, "vessels"):
        raise ConfigurationError(f"{unit.name} does not have any costing block.")

    var_dict = _var_dict_costing(unit)

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
        df.to_csv(f"{unit.local_name}_summary_costing.csv")

    print("\n" + "=" * 84)
    print(f"summary costing {unit.local_name}")
    print("-" * 84)
    stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
    print("\n" + "=" * 84 + "\n")
