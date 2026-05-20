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
Python script to read costing components
This script reads the library of costing components (scaled cost, reference
parameters, costing exponents, etc.) from the json files.
First, open json file, then create a python dictionary that gets imported into
power_plant_capcost.py

Two python dictionaries that are loaded:
* BB_costing_data
* sCO2_costing_params

"""

# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

__author__ = "Costing Team (A. Noring, B. Paul, D. Caballero, and M. Zamarripa)"
__version__ = "1.0.0"

import os
import json
from sys import stdout
import textwrap

from pyomo.common.fileutils import this_file_dir
from pyomo.environ import value, units as pyunits

from idaes.core import register_idaes_currency_units
from idaes.core.util.tables import stream_table_dataframe_to_string
import idaes.logger as idaeslog

from pandas import DataFrame

directory = this_file_dir()

_log = idaeslog.getLogger(__name__)

def register_power_plant_currency_units():
    """
    Define conversion rates for US Dollars based on CE Index.
    """
    register_idaes_currency_units()
    if (
        "USD_2008_Nov" in pyunits._pint_registry  # pylint: disable=protected-access
        and "USD_2019_Sep" in pyunits._pint_registry  # pylint: disable=protected-access
        and "USD_2018_Dec" in pyunits._pint_registry  # pylint: disable=protected-access
    ):
        # Assume that custom plant units have already been registered
        # Log a message and end
        _log.info(
            "Custom plant currency units (USD_2008_Nov, USD_2019_Sep, USD_2018_Dec) "
            "already appear in Pyomo unit registry. Assuming repeated "
            "call of custom_power_plant_currency_units."
        )
    else:
        pyunits.load_definitions_from_strings(
            [
                # power plant cost account units
                # from https://toweringskills.com/financial-analysis/cost-indices/
                "USD_2008_Nov = 500/566.2 * USD_CE500",
                "USD_2019_Sep = 500/599.3 * USD_CE500",
                "USD_2018_Dec = 500/615.7 * USD_CE500",
            ]
        )

def load_BB_costing_dictionary():
    """
    The costing data dictionary contains information from the BBR4 COE
    spreadsheet and from the QGESS on capital cost scaling methodology
    (DOE/NETL-2019/1784). Specifically it includes scaling exponents,
    valid ranges for the scaled parameter, and units for those ranges.
    It is important to note the units only apply to the ranges and are
    not necessarily the units that the reference parameter value will
    be given in.. It includes the total plant cost (TPC), reference
    parameter value, and units for that value.

    Some accounts are costed using two different reference parameters, these
    accounts have been divided into two separate accounts following the naming
    convention x.x.a and x.x.b.

    This dictionary is nested with the following structure:
    tech type --> CCS --> account --> property name --> property values
    """

    # assuming the dictionary exists, load it so it is importable when called
    with open(os.path.join(directory, "BB_costing_data.json"), "r") as file:
        BB_costing_params = json.load(file)
    return BB_costing_params


def load_sCO2_costing_dictionary():
    with open(os.path.join(directory, "sCO2_costing_parameters.json"), "r") as file:
        sCO2_costing_params = json.load(file)
    return sCO2_costing_params


def define_preloaded_accounts():
    PC_preloaded_accounts = {
        "Coal Handling": ["1.1", "1.2", "1.3", "1.4", "1.9a"],
        "Sorbent Handling": ["1.5", "1.6", "1.7", "1.8", "1.9b"],
        "Coal Feed": ["2.1", "2.2", "2.9a"],
        "Sorbent Feed": ["2.5", "2.6", "2.9b"],
        "Feedwater System": ["3.1", "3.3"],
        "PC Boiler": ["4.9"],
        "Steam Turbine": ["8.1"],
        "Condenser": ["8.3"],
        "Cooling Tower": ["9.1"],
        "Circulating Water System": ["9.2", "9.3", "9.4", "9.6", "9.7"],
        "Ash Handling": ["10.6", "10.7", "10.9"],
    }

    IGCC_preloaded_accounts = {
        "Coal Handling": ["1.1", "1.2", "1.3", "1.4", "1.9"],
        "Coal Feed": ["2.1", "2.2", "2.3", "2.4", "2.9"],
        "Feedwater System": ["3.1", "3.3"],
        "Gasifier": ["4.1"],
        "Syngas Cooler": ["4.2"],
        "ASU": ["4.3a"],
        "ASU Oxidant Compression": ["4.3b"],
        "Combustion Turbine": ["6.1", "6.3"],
        "Syngas Expander": ["6.2"],
        "HRSG": ["7.1", "7.2"],
        "Steam Turbine": ["8.1"],
        "Condenser": ["8.3"],
        "Cooling Tower": ["9.1"],
        "Circulating Water System": ["9.2", "9.3", "9.4", "9.6", "9.7"],
        "Slag Handling": ["10.1", "10.2", "10.3", "10.6", "10.7", "10.8", "10.9"],
    }

    NGCC_preloaded_accounts = {
        "Feedwater System": ["3.1", "3.3"],
        "Combustion Turbine": ["6.1", "6.3"],
        "HRSG": ["7.1", "7.2"],
        "Steam Turbine": ["8.1"],
        "Condenser": ["8.3"],
        "Cooling Tower": ["9.1"],
        "Circulating Water System": ["9.2", "9.3", "9.4", "9.6", "9.7"],
    }

    AUSC_preloaded_accounts = {
        "PC Boiler": ["4.9"],
        "Steam Turbine": ["8.1"],
        "Steam Piping": ["8.4"],
    }

    return (
        PC_preloaded_accounts,
        IGCC_preloaded_accounts,
        NGCC_preloaded_accounts,
        AUSC_preloaded_accounts,
    )


def load_default_resource_prices():
    """
    Dictionary of default prices
    MUSD: the currency units are millions of USD, so its price need a 1e-6 multiplier to get USD
    """

    default_resource_prices = {
        # from NETL QGESS Rev. 4 report
        "natural_gas": 4.42 * 1e-6 * pyunits.USD_2018 / pyunits.MBtu,  # $/MMbtu
        "coal": 51.96 * 1e-6 * pyunits.USD_2018 / pyunits.ton,
        "water": 1.90e-3 * 1e-6 * pyunits.USD_2018 / pyunits.gallon,
        "water_treatment_chemicals": 550 * 1e-6 * pyunits.USD_2018 / pyunits.ton,
        "ammonia": 300 * 1e-6 * pyunits.USD_2018 / pyunits.ton,
        "SCR_catalyst": 150 * 1e-6 * pyunits.USD_2018 / pyunits.ft**3,
        "triethylene_glycol": 6.80 * 1e-6 * pyunits.USD_2018 / pyunits.gallon,
        "SCR_catalyst_waste": 2.50 * 1e-6 * pyunits.USD_2018 / pyunits.ft**3,
        "triethylene_glycol_waste": 0.35 * 1e-6 * pyunits.USD_2018 / pyunits.gallon,
        "amine_purification_unit_waste": 38 * 1e-6 * pyunits.USD_2018 / pyunits.ton,
        "thermal_reclaimer_unit_waste": 38 * 1e-6 * pyunits.USD_2018 / pyunits.ton,
        # from S. McNaul, "Screening Techno-economic Analysis of NETL Reactive
        # Capture Technology," National Energy Technology Laboratory, Pittsburgh,
        # September 30, 2022. Exhibit B-8
        "PSA_adsorbent": 150 * pyunits.USD_2018_Dec / pyunits.ft**3,
        "amine_entrainer": 861.1 * pyunits.USD_2018_Dec / pyunits.kg,
        "electricity": 71.7 * pyunits.USD_2018_Dec / pyunits.MWh,
        "reactive_membrane_replacement": 15.1 * pyunits.USD_2018_Dec / pyunits.m**2,
        "PSA_adsorbent_waste_disposal": 1.5 * pyunits.USD_2018_Dec / pyunits.ft**3,
        "nonharzardous_waste_disposal": 38 * pyunits.USD_2018_Dec / pyunits.ton,
        # from X. Ge, R. Zhang, P. Liu, B. Liu, B. Liu, Optimization and control of
        # extractive distillation for formic acid-water separation with maximum-boiling
        # azeotrope, Computers & Chemical Engineering, Vol. 169, 2023. Table 2.
        "cooling_water": 0.354 * pyunits.USD_2016 / pyunits.GJ,
        "chilled_water": 4.42 * pyunits.USD_2016 / pyunits.GJ,
        "lp_heating_steam": 7.78 * pyunits.USD_2016 / pyunits.GJ,
        "mp_heating_steam": 8.22 * pyunits.USD_2016 / pyunits.GJ,
        "hp_heating_steam": 9.83 * pyunits.USD_2016 / pyunits.GJ,
        # from Zaiz, Toufik & Lanez, Hafnaoui. (2013). ASPEN HYSYS SIMULATION AND
        # COMPARISON BETWEEN ORGANIC SOLVENTS (SULFOLANE AND DMSO) USED FOR BENZENE
        # EXTRACTION. International Journal of Chemical and Petroleum Sciences. 2.
        # 10-19. Section 7.1b on page 18.
        "sulfolane_entrainer": 3700 * pyunits.USD_2013 / pyunits.ton,
        # from S. McNaul, "Comparison of Commercial State-of-the-Art,
        # Fossil-Based Hydrogen Production Technologies," National Energy
        # Technology Laboratory, Pittsburgh, April 12, 2022. Exhibit 3-15
        "ZnO_sulfur_guard_catalyst": 600 * pyunits.USD_2018_Dec / pyunits.ft**3,
        "prereformer_catalyst": 1250 * pyunits.USD_2018_Dec / pyunits.ft**3,
        "reformer_catalyst": 525 * pyunits.USD_2018_Dec / pyunits.ft**3,
        "ZnO_sulfur_guard_catalyst_waste_disposal": 40
        * pyunits.USD_2018_Dec
        / pyunits.ft**3,
        "prereformer_catalyst_waste_disposal": 0.0
        * pyunits.USD_2018_Dec
        / pyunits.ft**3,
        "reformer_catalyst_waste_disposal": 0.0 * pyunits.USD_2018_Dec / pyunits.ft**3,
    }
    return default_resource_prices

def load_fixed_OM_data():
    # set labor types, rates, operators per shift
    labor_types = [
        "skilled",
        "unskilled",
        "supervisor",
        "maintenance",
        "technician",
        "engineer",
        "operator",  # per power plant costing
    ]

    labor_rates = [27.90, 23.26, 30.29, 24.06, 23.43, 46.82, 38.50]

    # define maintenance percentages
    # each entry has elements [0], [1], and [2]
    # [0] is how much of the maintenance cost is the labor
    # [1] is the % of TPC used to estimate the maintenance labor cost
    # [2] is the % of TPC used to estimate the maintenance material cost
    # so we get
    # maintenance labor cost = TPC * [0] * [1]
    # maintenance material cost = TPC * (1 - [0]) * [2]
    # [1] = [2] for all technologies except 8, which is specific to that reference

    maintenance_percentages = {
        1: [0.4, 0.016, 0.016],
        2: [0.4, 0.016, 0.016],
        3: [0.35, 0.03, 0.03],
        4: [0.35, 0.03, 0.03],
        5: [0.35, 0.03, 0.03],
        6: [0.4, 0.019, 0.019],
        7: [0.4, 0.016, 0.016],
        8: [0.4, 0.019, 0.016],
        9: [0.4, 0.019, 0.019],
    }

    return labor_types, labor_rates, maintenance_percentages

    def report(b, export=False):
        """
        b: flowsheet-level costing block
        export: generate csv results file of costing results
        """
        var_dict = {}

        if hasattr(b, "total_TPC"):
            var_dict["Total TPC [$MM]"] = value(b.total_TPC)

        if hasattr(b, "total_overnight_capital"):
            var_dict["Total Overnight Cost [$MM]"] = value(b.total_overnight_capital)

        if hasattr(b, "total_as_spent_cost"):
            var_dict["Total As Spent Cost [$MM]"] = value(b.total_as_spent_cost)

        if hasattr(b, "annualized_cost"):
            var_dict["Total Annualized Capital Cost [$MM/year]"] = value(
                b.annualized_cost
            )

        if hasattr(b, "total_fixed_OM_cost"):
            var_dict["Total fixed O&M cost [$MM/year]"] = value(
                b.total_fixed_OM_cost
            )

        if hasattr(b, "total_variable_OM_cost"):
            var_dict["Total variable O&M cost full capacity [$MM/year]"] = value(
                b.total_variable_OM_cost[0]
            )
            var_dict["Total variable O&M cost operating capacity [$MM/year]"] = value(
                b.total_variable_OM_cost[0] * b.capacity_factor
            )

        if (
            hasattr(b, "annualized_cost")
            and hasattr(b, "total_fixed_OM_cost")
            and hasattr(b, "total_variable_OM_cost")
        ):
            var_dict["Total Annualized Cost [$MM/year]"] = (
                value(b.annualized_cost)
                + value(b.total_fixed_OM_cost)
                + b.capacity_factor * value(b.total_variable_OM_cost[0])
            )

        if hasattr(b, "cost_of_electricity"):
            var_dict["Cost of Electricity [$/MWh]"] = value(
                b.cost_of_electricity * 1e6
            )

        if hasattr(b, "cost_of_production"):
            var_dict["Cost of Production [$/kg product]"] = value(
                b.cost_of_production * 1e6
            )

        if hasattr(b, "transport_cost_of_product"):
            var_dict["Total Transport Cost Of Product [$MM]"] = value(b.transport_cost_of_product)

        if hasattr(b, "cost_of_capture"):
            var_dict["Cost of Capture [$/tonne CO2]"] = value(b.cost_of_capture * 1e6)

        if hasattr(b, "CO2_transport_cost"):
            var_dict["Total Transport Cost of Feedstock & Captured CO2 [$MM]"] = value(b.CO2_transport_cost)

        report_dir = {}
        report_dir["Value"] = {}
        report_dir["pos"] = {}

        count = 1
        for k, v in var_dict.items():
            report_dir["Value"][k] = value(v)
            report_dir["pos"][k] = count
            count += 1

        df = DataFrame.from_dict(report_dir, orient="columns")
        del df["pos"]
        if export:
            df.to_csv(f"{b.local_name}_report.csv")

        print("\n" + "=" * 84)
        print(f"{b.local_name}")
        print("-" * 84)
        stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
        print("\n" + "=" * 84 + "\n")