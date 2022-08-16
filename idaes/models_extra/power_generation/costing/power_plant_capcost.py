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
Power Plant costing library
This method leverages NETL costing capabilities. Two main methods have been
developed to calculate the capital cost of power generation plants:

1. Fossil fueled power plants (from SCPC to IGCC) (get_PP_costing)
2. Supercritical CO2 power cycles (direct and indirect) (get_sCO2_unit_cost)

Other methods:

    - get_fixed_OM_costs() to cost fixed O&M costs
    - get_variable_OM_costs() to cost variable O&M costs
    - get_ASU_cost() to cost air separation units
    - costing_initialization() to initialize costing blocks
    - display_total_plant_costs() to display total plant cost (TPC)
    - display_bare_erected_costs() to display BEC costs
    - get_total_TPC() to display the total TPC of the entire flowsheet
    - display_flowsheet_cost() to display flowsheet cost
    - check_sCO2_costing_bounds() to display a warnning if costing model have
      been used outside the range that where designed for
"""
__author__ = "Costing Team (A. Noring and M. Zamarripa)"
__version__ = "1.0.0"

from pandas import DataFrame
from sys import stdout
import textwrap

from pyomo.environ import (
    Param,
    Var,
    Block,
    Constraint,
    Expression,
    value,
    Expr_if,
    units as pyunits,
)
from pyomo.core.base.units_container import InconsistentUnitsError
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
from idaes.core import register_idaes_currency_units
from idaes.models_extra.power_generation.costing.costing_dictionaries import (
    BB_costing_exponents,
    BB_costing_params,
    sCO2_costing_params,
)

from idaes.core import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)

from idaes.core.util.tables import stream_table_dataframe_to_string

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

# Register standard currency units
register_idaes_currency_units()

# -----------------------------------------------------------------------------
# Power Plant Costing Library
# -----------------------------------------------------------------------------

import pyomo.environ as pyo
from idaes.core import declare_process_block_class


@declare_process_block_class("QGESSCosting")
class QGESSCostingData(FlowsheetCostingBlockData):
    # Register currency and conversion rates based on CE Index
    register_idaes_currency_units()

    def build_global_params(self):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        # Set the base year for all costs
        self.base_currency = pyo.units.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyo.units.year

    def build_process_costs(
        self,
        # arguments related to Fixed OM costs
        total_plant_cost=None,
        nameplate_capacity=650,
        labor_rate=38.50,
        labor_burden=30,
        operators_per_shift=6,
        tech=1,
        # arguments related to total owners costs
        land_cost=None,
        net_power=None,
        resources=None,
        rates=None,
        prices=None,
        fixed_OM=True,
        variable_OM=False,
        fuel=None,
        chemicals=None,
        waste=None,
        tonne_CO2_capture=None,
    ):
        """
        This is where you do all your process wide costing.
        This is completely up to you, but you will have access to the
        following aggregate costs:

        1. self.aggregate_capital_cost
        2. self.aggregate_fixed_operating_cost
        3. self.aggregate_variable_operating_cost
        4. self.aggregate_flow_costs (indexed by flow type)
        """

        if total_plant_cost is None:
            self.get_total_TPC()
        # else:
        #     self.total_TPC = Var(initialize=0, bounds=(0, 1e4), doc="total TPC in $MM")

        #     @self.Constraint()
        #     def total_TPC_eq(c):
        #         return c.total_TPC == total_plant_cost

        if fixed_OM is True:
            self.get_fixed_OM_costs(
                nameplate_capacity=nameplate_capacity,
                labor_rate=labor_rate,
                labor_burden=labor_burden,
                operators_per_shift=operators_per_shift,
                tech=tech
            )

        if variable_OM is True:
            # get_variable_OM_costs(fs, production_rate, resources, rates, prices={}):
            self.get_variable_OM_costs(net_power, resources, rates, prices=prices)

        # build system costs (owner's, total overnight costs, annualized costs, COE, and cost of capture)
        if fixed_OM and variable_OM:
            # total overnight cost requires fixed costs
            self.pct_TPC = Param(
                initialize=20.2/100, doc="Fixed percentaje for other owners cost"
            )
            # self.land_cost = Param(initialize=land_cost, doc='percentaje of total plant cost')
            self.six_month_labor = Expression(
                expr=(
                    self.annual_operating_labor_cost
                    + self.maintenance_labor_cost
                    + self.admin_and_support_labor_cost
                )
                / 2
            )
            non_fuel_resources = resources  # duplicate resources list
            if fuel is not None:
                self.fuel_cost_OC = Expression(
                    expr=self.variable_operating_costs[0, fuel]
                    / (85 / 100)
                    / 12
                    * 2.25,
                    doc="Owner's costs - 2.25 months of fuel costs",
                )
                non_fuel_resources.remove(fuel)  # remove fuel from the list

            if waste is not None:
                self.waste_costs_OC = Expression(
                    expr= (sum(self.variable_operating_costs[0, i] for i in waste) / (85 / 100) / 12
                )
            )
            if chemicals is not None:
                self.chemical_costs_OC = Expression(
                    expr= (sum(self.variable_operating_costs[0, i] for i in chemicals)/ 2  # six months of chemicals
                )
            )

            self.non_fuel_and_waste_OC = Expression(
                expr=(
                    sum(self.variable_operating_costs[0, i] for i in non_fuel_resources)
                )
                / (85 / 100)
                / 12
            )

            self.total_overnight_capital = Expression(
                expr=
                self.total_TPC
                # pre production costs
                + self.six_month_labor
                + self.maintenance_material_cost / 12  # 1 month maintenance materials
                + self.non_fuel_and_waste_OC  # 1 month non fuel consumables
                + (self.waste_costs_OC if waste is not None else 0)  # 1 month waste
                # inventory capital costs
                + (self.fuel_cost_OC if fuel is not None else 0)  # 60 day fuel supply
                # Other costs
                + (self.chemical_costs_OC if waste is not None else 0) # Initial Cost for Catalyst and Chemicals
                + (land_cost if land_cost is not None else 0)
                + self.total_TPC * self.pct_TPC  # other owners costs (other + spare parts + financing + other pre-production)
            )

            self.tasc_toc_factor = Param(
                initialize=1.093,
                mutable=True,
                doc="TASC/TOC factor from Exhibit 3-7 reference 1, real for three years 1.093, real for five years 1.154, nominal for three years 1.242, nominal for five years 1.289",
            )

            self.total_as_spent_cost = Expression(
                expr=self.total_overnight_capital * self.tasc_toc_factor
            )

            self.fixed_charge_factor = Param(
                initialize=0.0707,
                mutable=True,
                doc="Fixed charge rate from Exhibit 3-5 based on CRF values, real for three/five years 0.0707, nominal for three/five years 0.0886",
            )
            self.capacity_factor = Param(
                initialize=0.85, doc="capacity factor of the plant"
            )
            self.annualized_cost = Expression(
                expr=self.fixed_charge_factor * self.total_as_spent_cost
            )

            # only build COE for power plants
            if net_power is not None:
                self.cost_of_electricity = Expression(
                    expr=(
                        self.annualized_cost
                        + self.total_fixed_OM_cost
                        + self.total_variable_OM_cost[0] * self.capacity_factor
                    )
                    / (self.capacity_factor * net_power[0])
                )
            if tonne_CO2_capture is not None:
                self.cost_of_capture = Expression(
                    expr=(
                        self.annualized_cost
                        + self.total_fixed_OM_cost
                        + self.total_variable_OM_cost[0] * self.capacity_factor
                    ) * 1e6
                    / tonne_CO2_capture
                )

    @staticmethod
    def initialize_build(self):
        """
        Here we can add intialization steps for the things we built in
        build_process_costs.

        Note that the aggregate costs will be initialized by the framework.
        """
        # TODO: For now,  no additional process level costs to initialize
        pass

    def report(self, export=False):

        var_dict = {}

        if hasattr(self, "total_TPC"):
            var_dict["Total TPC [$MM]"] = value(self.total_TPC)

        if hasattr(self, "total_overnight_capital"):
            var_dict["Total Overnigth Cost [$MM]"] = value(self.total_overnight_capital)

        if hasattr(self, "total_as_spent_cost"):
            var_dict["Total As Spent Cost [$MM]"] = value(self.total_as_spent_cost)

        if hasattr(self, "annualized_cost"):
            var_dict["Total Annualized Capital Cost [$MM/year]"] = value(
                self.annualized_cost
            )

        if hasattr(self, "total_fixed_OM_cost"):
            var_dict["Total fixed O&M cost [$MM/year]"] = value(
                self.total_fixed_OM_cost
            )

        if hasattr(self, "total_variable_OM_cost"):
            var_dict["Total variable O&M cost [$MM/year]"] = value(
                self.total_variable_OM_cost[0]
            )

        if (
            hasattr(self, "annualized_cost")
            and hasattr(self, "total_fixed_OM_cost")
            and hasattr(self, "total_variable_OM_cost")
        ):
            var_dict["Total Annualized Cost [$MM/year]"] = (
                value(self.annualized_cost)
                + value(self.total_fixed_OM_cost)
                + value(self.total_variable_OM_cost[0])
            )

        if hasattr(self, "cost_of_electricity"):
            var_dict["Cost of Electricity [$/MWh]"] = value(self.cost_of_electricity)

        if hasattr(self, "cost_of_capture"):
            var_dict["Cost of Capture [$/tonne]"] = value(self.cost_of_capture)

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
            df.to_csv(f"{self.local_name}_report.csv")

        print("\n" + "=" * 84)
        print(f"{self.local_name}")
        print("-" * 84)
        stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
        print("\n" + "=" * 84 + "\n")

    def get_PP_costing(
        blk,
        cost_accounts,
        scaled_param,
        units,
        tech,
        ccs="B",
        CE_index_base=671.1,
        additional_costing_exponents=None,
        additional_costing_params=None,
    ):
        """
        Power Plant Costing Method
        This method relies on the capital cost scaling methodologies developed
        by NETL. Report #DOE/NETL-341/013113
        Multiple vendors quotes have been used to determine the cost of several
        plant equipments (i.e. boiler, pumps, heat exchangers, etc.), other cost
        incurred during the plant operation (i.e. solids handling, etc.)

        The scaling approach uses a main equation:
            SC = RC*(SP/RP)^Exp

        where:
            SC is the scaled cost
            RC is the reference cost
            SP is the scaled operational/design parameter
            RP is the reference operational/design parameter
            Exp is the scaling exponent

        This equation is implemented when one operational/design parameter is
        used to scale the cost

        When two operational/design parameters are used to scale the cost,

        the following equation is implemented:

            SC = RC*w1*(SP1/RP1)^Exp + RC*w2*(SP2/RP2)^Exp

        where:
            SC is the scaled cost
            RC is the reference cost
            SP1 is the scaled operational/design parameter 1
            RP1 is the reference operational/design parameter 1
            w1 is the fraction of the cost scaled using
            operational/design parameter 1
            SP2 is the scaled operational/design parameter 2
            RP2 is the reference operational/design parameter 2
            w2 is the fraction of the cost scaled using
            operational/design parameter 2
            Exp is the scaling exponent

        The scaled cost is computed using ref values for different technologies.
        Categories:
        1. Supercritical PC, air-fired, with and without CO2 capture,
        Illinois No. 6 coal
        2. Subcritical PC, air-fired, with and without CO2 capture,
        Illinois No. 6 coal
        3. Two-stage, slurry-feed, oxygen-blown gasifier with and without
        CO2 capture, Illinois No. 6 coal
        4. Single-stage, slurry-feed, oxygen-blown gasifier with and without
        CO2 capture, Illinois No. 6 coal
        5. Single-stage, dry-feed, oxygen-blown, up-flow gasifier with
        and without CO2 capture, Illinois No. 6 coal
        6. Natural gas, air-fired, with and without CO2 capture
        7. Advanced Ultrasupercritical PC

        This method computes the capital cost of units and main components of the
        power plant, and requires a few arguments to build a constraint as part of
        your main model.

        Args:
            self: A block or unit model where costing constraints can be added to
            accounts: A list of accounts to be included in the total cost,
                they should all use the same reference parameter
            scaled_param: the process parameter for the system(s) being costed
            units: the units of the scaled_param, used for verification
            tech: int 1-7 representing the above catagories
            ccs: 'A' or 'B' representing no CCS or CCS

        The appropriate scaling parameters for various cost accounts can be found
        in the QGESS on capital cost scaling (Report #DOE/NETL-341/013113).
        The correct units for the reference parameters are found in the BBR4 COE
        spreadsheet.

        """
        from idaes.models_extra.power_generation.costing.power_plant_capcost_custom_dict import (
            custom_costing_exponents,
            custom_costing_params,
        )

        # ------------------------ Power Plant Cost ------------------------

        # check to see if a costing block already exists
        # if hasattr(self, "costing"):
        #     raise AttributeError(
        #         "{} already has an attribute costing. "
        #         "Check that you are not calling get_costing"
        #         " twice on the same model".format(self.name)
        #     )

        # # create a costing Block
        # self.costing = Block()
        # self.costing.library = "PP"

        # # find flowsheet block to create global costing parameters
        # try:
        #     fs = self.flowsheet()
        # except AttributeError:
        #     fs = self.parent_block()

        # # build flowsheet level parameters CE_index = year
        # if not hasattr(fs, "costing"):
        #     fs.get_costing(year="2018")

        CE_index = CE_index_base
        if additional_costing_params is not None:
            # adding custom exponents and params to the dictionary
            for key, value in additional_costing_exponents.items():
                custom_costing_exponents[key] = value
            for key, value in additional_costing_params.items():
                custom_costing_params[key] = value

        # define preloaded accounts
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

        # preloaded account handling
        if type(cost_accounts) == str:
            if tech in [1, 2]:
                cost_accounts = PC_preloaded_accounts[cost_accounts]
            elif tech in [3, 4, 5]:
                cost_accounts = IGCC_preloaded_accounts[cost_accounts]
            elif tech == 6:
                cost_accounts = NGCC_preloaded_accounts[cost_accounts]
            elif tech == 7:
                cost_accounts = AUSC_preloaded_accounts[cost_accounts]
            else:
                AttributeError("{} technology not supported".format(blk.name))

        # pull data for each account into dictionaries
        process_params = {}
        reference_units = {}
        account_names = {}
        exponents = {}
        reference_costs = {}
        reference_costs_init = {}
        reference_params = {}
        cost_scaling_fractions = {}
        engineering_fees = {}
        process_contingencies = {}
        project_contingencies = {}

        for account in cost_accounts:
            try:  # first look for data in json file info
                process_params[account] = BB_costing_exponents[str(tech)][account][
                    "Process Parameter"
                ]
                reference_units[account] = BB_costing_params[str(tech)][ccs][
                    cost_accounts[0]
                ]["Units"]
                account_names[account] = BB_costing_exponents[str(tech)][account][
                    "Account Name"
                ]
                exponents[account] = float(
                    BB_costing_exponents[str(tech)][account]["Exponent"]
                )
                reference_costs[account] = BB_costing_params[str(tech)][ccs][account][
                    "BEC"
                ]
                reference_costs_init[account] = (
                    BB_costing_params[str(tech)][ccs][account]["BEC"] * 1e-3
                )

                if type(process_params[account]) == list:
                    for i, processparam in enumerate(process_params[account]):
                        reference_params[account, processparam] = BB_costing_params[
                            str(tech)
                        ][ccs][account]["RP Value"][i]
                        cost_scaling_fractions[
                            account, processparam
                        ] = BB_costing_params[str(tech)][ccs][account][
                            "Cost scaling fraction"
                        ][
                            i
                        ]

                elif type(process_params[account]) == str:
                    reference_params[account] = BB_costing_params[str(tech)][ccs][
                        account
                    ]["RP Value"]

                engineering_fees[account] = BB_costing_params[str(tech)][ccs][account][
                    "Eng Fee"
                ]
                process_contingencies[account] = BB_costing_params[str(tech)][ccs][
                    account
                ]["Process Contingency"]
                project_contingencies[account] = BB_costing_params[str(tech)][ccs][
                    account
                ]["Project Contingency"]
            except KeyError:
                try:  # next look for data in custom dictionaries
                    process_params[account] = custom_costing_exponents[str(tech)][
                        account
                    ]["Process Parameter"]
                    reference_units[account] = custom_costing_params[str(tech)][ccs][
                        cost_accounts[0]
                    ]["Units"]
                    account_names[account] = custom_costing_exponents[str(tech)][
                        account
                    ]["Account Name"]
                    exponents[account] = float(
                        custom_costing_exponents[str(tech)][account]["Exponent"]
                    )
                    reference_costs[account] = custom_costing_params[str(tech)][ccs][
                        account
                    ]["BEC"]
                    reference_costs_init[account] = (
                        custom_costing_params[str(tech)][ccs][account]["BEC"] * 1e-3
                    )

                    if type(process_params[account]) == list:
                        for i, processparam in enumerate(process_params[account]):
                            reference_params[
                                account, processparam
                            ] = custom_costing_params[str(tech)][ccs][account][
                                "RP Value"
                            ][
                                i
                            ]

                            cost_scaling_fractions[
                                account, processparam
                            ] = custom_costing_params[str(tech)][ccs][account][
                                "Cost scaling fraction"
                            ][
                                i
                            ]

                    elif type(process_params[account]) == str:
                        reference_params[account] = custom_costing_params[str(tech)][
                            ccs
                        ][account]["RP Value"]

                    engineering_fees[account] = custom_costing_params[str(tech)][ccs][
                        account
                    ]["Eng Fee"]
                    process_contingencies[account] = custom_costing_params[str(tech)][
                        ccs
                    ][account]["Process Contingency"]
                    project_contingencies[account] = custom_costing_params[str(tech)][
                        ccs
                    ][account]["Project Contingency"]
                except KeyError:
                    print(
                        "KeyError: Account {} could not be found in the "
                        "dictionary for technology {} with CCS {}".format(
                            account, str(tech), ccs
                        )
                    )

        # check that all accounts use the same process parameter
        param_check = None
        for account in cost_accounts:
            param = process_params[account]
            if param_check is None:
                param_check = param
            elif param != param_check:
                raise ValueError(
                    "{} cost accounts selected do not use "
                    "the same process parameter".format(blk.name)
                )

        # check that the user passed the correct units
        for account in cost_accounts:
            ref_units = reference_units[account]
            if units != ref_units:
                raise ValueError(
                    "Account %s uses units of %s. "
                    "Units of %s were passed." % (cost_accounts[0], ref_units, units)
                )

        # Used by other functions for reporting results
        blk.account_names = account_names

        # define parameters
        blk.exp = Param(
            cost_accounts,
            mutable=True,
            initialize=exponents,
            doc="exponential parameter for account",
        )

        blk.ref_cost = Param(
            cost_accounts,
            mutable=True,
            initialize=reference_costs,
            doc="reference cost for account",
        )

        if type(process_params[cost_accounts[0]]) == list:
            if len(process_params[cost_accounts[0]]) > 1:
                blk.ref_param = Param(
                    cost_accounts,
                    process_params[cost_accounts[0]],
                    mutable=True,
                    initialize=reference_params,
                    doc="reference parameter for account",
                )

                blk.cost_scaling_fracs = Param(
                    cost_accounts,
                    process_params[cost_accounts[0]],
                    mutable=True,
                    initialize=cost_scaling_fractions,
                    doc="reference parameter for account",
                )
        elif type(process_params[cost_accounts[0]]) == str:
            blk.ref_param = Param(
                cost_accounts,
                mutable=True,
                initialize=reference_params,
                doc="reference parameter for account",
            )

        blk.eng_fee = Param(
            cost_accounts,
            mutable=True,
            initialize=engineering_fees,
            doc="engineering fee percentage",
        )

        blk.process_conting = Param(
            cost_accounts,
            mutable=True,
            initialize=process_contingencies,
            doc="process contingency percentage",
        )

        blk.project_conting = Param(
            cost_accounts,
            mutable=True,
            initialize=project_contingencies,
            doc="project contingency percentage",
        )

        # define variables
        blk.bare_erected_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, 1e4),
            doc="scaled bare erected cost in $MM",
        )

        blk.total_plant_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, 1e4),
            doc="total plant cost in $MM",
        )

        # rule for scaling BEC
        # reference cost is in 2018 dollars, 671.1 is CE index for 2018
        def bare_erected_cost_rule(costing, i):
            if type(process_params[i]) == list:
                if len(process_params[i]) > 1:
                    return costing.bare_erected_cost[i] * 1e3 == (
                        CE_index / CE_index_base
                    ) * costing.ref_cost[i] * sum(
                        costing.cost_scaling_fracs[i, p]
                        * (scaled_param[j] / costing.ref_param[i, p]) ** costing.exp[i]
                        for j, p in enumerate(process_params[i])
                    )
            elif type(process_params[i]) == str:
                return (
                    costing.bare_erected_cost[i] * 1e3
                    == (CE_index / CE_index_base)
                    * costing.ref_cost[i]
                    * (scaled_param / costing.ref_param[i]) ** costing.exp[i]
                )

        blk.bare_erected_cost_eq = Constraint(
            cost_accounts, rule=bare_erected_cost_rule
        )

        # rule for calculating TPC
        def total_plant_cost_rule(blk, i):
            return blk.total_plant_cost[i] == blk.bare_erected_cost[i] * (
                1 + blk.eng_fee[i] + blk.process_conting[i]
            ) * (1 + blk.project_conting[i])

        blk.total_plant_cost_eq = Constraint(cost_accounts, rule=total_plant_cost_rule)

        # rule for sum of BEC
        def BEC_sum_rule(blk):
            return sum(blk.bare_erected_cost[i] for i in cost_accounts)

        blk.bare_erected_cost_sum = Expression(rule=BEC_sum_rule)

        # rule for sum of TPC
        def TPC_sum_rule(blk):
            return sum(blk.total_plant_cost[i] for i in cost_accounts)

        blk.total_plant_cost_sum = Expression(rule=TPC_sum_rule)

        # add variable and constraint scaling
        for i in cost_accounts:
            iscale.set_scaling_factor(blk.bare_erected_cost[i], 1)
            iscale.set_scaling_factor(blk.total_plant_cost[i], 1)
            iscale.constraint_scaling_transform(
                blk.bare_erected_cost_eq[i], 1e-3, overwrite=False
            )
            iscale.constraint_scaling_transform(
                blk.total_plant_cost_eq[i], 1, overwrite=False
            )

    # -----------------------------------------------------------------------------
    # Supercritical CO2 Costing Library
    # -----------------------------------------------------------------------------
    def get_sCO2_unit_cost(
        self, equipment, scaled_param, temp_C=None, n_equip=1, custom_accounts=None
    ):
        """
        Args:
            self: pyomo Block where constraints will be made
            unit_name: the name of the SCO2 equipment to cost
            scaling_param: the scaling parameter (in appropriate units) for the
                elected equipment
            temp: the maximum temperature of the equipment. Not all types of
                equipment use a temperature correction factor, so it is optional
            n_equip: the number of pieces of equipment to cost

        Cost is in M$
        """
        # check to see if a costing block already exists
        if hasattr(self, "costing"):
            raise AttributeError(
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(self.name)
            )

        # create a costing Block
        self.costing = Block()
        self.costing.library = "sCO2"
        self.costing.equipment = equipment

        # find flowsheet block to create global costing parameters
        try:
            fs = self.flowsheet()
        except AttributeError:
            fs = self.parent_block()

        # build flowsheet level parameters CE_index = year
        if not hasattr(fs, "costing"):
            fs.get_costing(year="2017")

        CE_index = fs.costing.CE_index

        param_dict = sCO2_costing_params[equipment]

        # define parameters
        self.costing.ref_cost = Param(
            mutable=True, initialize=param_dict["a"], doc="Reference cost"
        )

        self.costing.exp = Param(
            mutable=True, initialize=param_dict["b"], doc="Scaling exponent"
        )

        self.costing.c = Param(
            mutable=True,
            initialize=param_dict["c"],
            doc="coefficient for temperature correction",
        )

        self.costing.d = Param(
            mutable=True,
            initialize=param_dict["d"],
            doc="coefficient for temperature correction",
        )

        self.costing.material_cost = Param(
            mutable=True,
            doc="material installation cost",
            initialize=param_dict["Material Cost"],
        )

        self.costing.labor_cost = Param(
            mutable=True,
            initialize=param_dict["Labor Cost"],
            doc="labor installation cost",
        )

        # estimates for the percentages of TPC will be added later

        self.costing.eng_fee = Param(
            mutable=True, initialize=0, doc="engineering fee percentage"
        )

        self.costing.process_conting = Param(
            mutable=True, initialize=0, doc="process contingency percentage"
        )

        self.costing.project_conting = Param(
            mutable=True, initialize=0, doc="project contingency percentage"
        )

        # define variables
        # n_equip is left as a fixed variable to support MINLP optimization
        self.costing.n_equip = Var(
            initialize=n_equip, doc="number of pieces of equipment"
        )

        self.costing.n_equip.fix(n_equip)

        self.costing.scaled_param = Var(
            initialize=scaled_param, bounds=(0, 1e12), doc="scaled parameter"
        )

        self.costing.temp_factor = Var(
            initialize=1, bounds=(0.9, 100), doc="temperature correction factor"
        )

        self.costing.equipment_cost = Var(
            initialize=self.costing.ref_cost,
            bounds=(0, 1e4),
            doc="equipment cost of sCO2 unit in $MM",
        )

        self.costing.bare_erected_cost = Var(
            initialize=self.costing.ref_cost,
            bounds=(0, 1e4),
            doc="bare erected cost of sCO2 unit" "in $MM",
        )

        self.costing.total_plant_cost = Var(
            initialize=self.costing.ref_cost,
            bounds=(0, 1e4),
            doc="total plant cost of sCO2 unit" "in $MM",
        )

        # divides the scaled parameter by the number of pieces of equipment
        def scaled_param_rule(costing):
            return costing.scaled_param * costing.n_equip == scaled_param

        self.costing.scaled_param_eq = Constraint(rule=scaled_param_rule)

        # check if equipment requires a temperature correction factor
        if equipment in [
            "Axial turbine",
            "Radial turbine",
            "Coal-fired heater",
            "Natural gas-fired heater",
            "Recuperator",
        ]:

            if temp_C is None:
                raise ValueError(
                    "Temperature argument is "
                    "required to cost %s equipment" % equipment
                )

            else:
                self.costing.temperature = Var(
                    initialize=500, bounds=(0, 1e6), doc="dummy var for temperature"
                )

                self.costing.temp_eq = Constraint(
                    expr=(self.costing.temperature == temp_C)
                )

                def temp_correction_rule(costing):  # rule for temp correction
                    return (
                        Expr_if(
                            costing.temperature < 550,
                            1e-6 * costing.temperature + 1,
                            1
                            + costing.c * (costing.temperature - 550)
                            + costing.d * (costing.temperature - 550) ** 2,
                        )
                        == costing.temp_factor
                    )

                self.costing.temp_correction_eq = Constraint(rule=temp_correction_rule)
        else:
            self.costing.temp_factor.fix(1)

        # rule for equipment cost
        def equipment_cost_rule(costing):
            return (
                costing.equipment_cost * 1e6
                == (CE_index / 567.5)
                * costing.n_equip
                * costing.ref_cost
                * (costing.scaled_param**costing.exp)
                * costing.temp_factor
            )

        self.costing.equipment_cost_eq = Constraint(rule=equipment_cost_rule)

        # rule for bare erected cost
        def bare_erected_cost_rule(costing):
            return costing.bare_erected_cost == costing.equipment_cost * (
                1 + costing.material_cost + costing.labor_cost
            )

        self.costing.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

        # rule for calculating total plant cost
        def total_plant_cost_rule(costing):
            return costing.total_plant_cost == costing.bare_erected_cost * (
                1 + costing.eng_fee + costing.process_conting + costing.project_conting
            )

        self.costing.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # add variable and constraint scaling
        if equipment in ["Recuperator", "Direct air cooler"]:
            iscale.set_scaling_factor(self.costing.scaled_param, 1e-5)
        else:
            iscale.set_scaling_factor(self.costing.scaled_param, 1)

        iscale.set_scaling_factor(self.costing.equipment_cost, 1e3)
        iscale.set_scaling_factor(self.costing.bare_erected_cost, 1e3)
        iscale.set_scaling_factor(self.costing.total_plant_cost, 1e3)
        iscale.constraint_scaling_transform(
            self.costing.equipment_cost_eq, 1e-6, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.costing.bare_erected_cost_eq, 1e3, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.costing.bare_erected_cost_eq, 1e3, overwrite=False
        )

    # -----------------------------------------------------------------------------
    # Air Separation Unit Costing Library
    # -----------------------------------------------------------------------------
    def get_ASU_cost(self, scaled_param):
        # scaled parameter is O2 flowrate in TPD

        params = {
            "Reference Cost": 3.26e6,
            "Reference Parameter": 13078,
            "Exponent": 0.7,
            "Eng Fee": 0.097,
            "Process": 0,
            "Project": 0.110,
        }

        # check to see if a costing block already exists
        if hasattr(self, "costing"):
            raise AttributeError(
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(self.name)
            )

        # create a costing Block
        self.costing = Block()
        self.costing.library = "ASU"

        # find flowsheet block to create global costing parameters
        try:
            fs = self.flowsheet()
        except AttributeError:
            fs = self.parent_block()

        # build flowsheet level parameters CE_index = year
        if not hasattr(fs, "costing"):
            fs.get_costing(year="2017")

        CE_index = fs.costing.CE_index

        # define parameters
        self.costing.ref_cost = Param(
            initialize=params["Reference Cost"], mutable=True, doc="ASU reference cost"
        )

        self.costing.ref_param = Param(
            initialize=params["Reference Parameter"],
            mutable=True,
            doc="ASU reference parameter value",
        )

        self.costing.exp = Param(
            initialize=params["Exponent"], mutable=True, doc="ASU scaling exponent"
        )

        self.costing.eng_fee = Param(
            mutable=True, initialize=params["Eng Fee"], doc="engineering fee percentage"
        )

        self.costing.process_conting = Param(
            mutable=True,
            initialize=params["Process"],
            doc="process contingency percentage",
        )

        self.costing.project_conting = Param(
            mutable=True,
            initialize=params["Project"],
            doc="project contingency percentage",
        )

        # define variables
        self.costing.bare_erected_cost = Var(
            initialize=params["Reference Cost"],
            bounds=(0, 1e4),
            doc="scaled bare erected cost in $MM",
        )

        self.costing.total_plant_cost = Var(
            initialize=params["Reference Cost"],
            bounds=(0, 1e4),
            doc="total plant cost in $MM",
        )

        # rule for scaling BEC
        # reference cost is in 2008 dollars, 566.2 is CE index for Nov 2008

        def bare_erected_cost_rule(costing):
            return (
                costing.bare_erected_cost * 1e3
                == (CE_index / 566.2)
                * costing.ref_cost
                * (scaled_param / costing.ref_param) ** costing.exp
            )

        self.costing.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

        # rule for calculating TPC
        def total_plant_cost_rule(costing):
            return costing.total_plant_cost == costing.bare_erected_cost * (
                1 + costing.eng_fee + costing.process_conting + costing.project_conting
            )

        self.costing.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # add variable and constraint scaling
        iscale.set_scaling_factor(self.costing.bare_erected_cost, 1)
        iscale.set_scaling_factor(self.costing.total_plant_cost, 1)

        iscale.constraint_scaling_transform(
            self.costing.bare_erected_cost_eq, 1e-3, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.costing.total_plant_cost_eq, 1, overwrite=False
        )

    # -----------------------------------------------------------------------------
    # Operation & Maintenance Costing Library
    # -----------------------------------------------------------------------------

    def get_fixed_OM_costs(
        b,
        nameplate_capacity=650,
        labor_rate=38.50,
        labor_burden=30,
        operators_per_shift=6,
        tech=1,
        fixed_TPC=None,
    ):
        """
        Creates constraints for the following fixed O&M costs in $MM/yr:
            1. Annual operating labor
            2. Maintenance labor
            3. Admin and support labor
            4. Property taxes and insurance
            5. Other fixed costs
            6. Total fixed O&M cost
            7. Maintenance materials (actually a variable cost, but scales off TPC)

        These costs apply to the project as a whole and are scaled based on the
        total TPC.

        Args:
            b: pyomo concrete model or flowsheet block
            nameplate_capacity: rated plant output in MW
            labor_rate: hourly rate of plant operators in project dollar year
            labor_burden: a percentage multiplier used to estimate non-salary
                labor expenses
            operators_per_shift: average number of operators per shift
            tech: int 1-7 representing the catagories in get_PP_costing, used to
                determine maintenance costs
            TPC_value: The TPC in $MM that will be used to determine fixed O&M
            costs. If the value is None, the function will try to use the TPC
                calculated from the individual units.

        Returns:
            None
        """
        # # check if costing block exists
        # if not hasattr(b, "costing"):
        #     b.get_costing(year="2018")

        # create and fix total_TPC if it does not exist yet
        if not hasattr(b, "total_TPC"):
            b.total_TPC = Var(initialize=0, bounds=(0, 1e4), doc="total TPC in $MM")
            if fixed_TPC is None:
                b.total_TPC.fix(100)
                _log.warning(
                    "b.costing.total_TPC does not exist and a value "
                    "for fixed_TPC was not specified, total_TPC will be "
                    "fixed to 100 MM$"
                )
            else:
                b.total_TPC.fix(fixed_TPC)
        else:
            if fixed_TPC is not None:
                _log.warning(
                    "b.costing.total_TPC already exists, the value "
                    "passed for fixed_TPC will be ignored."
                )

        # make params
        b.labor_rate = Param(initialize=labor_rate, mutable=True)
        b.labor_burden = Param(initialize=labor_burden, mutable=True)
        b.operators_per_shift = Param(initialize=operators_per_shift, mutable=True)
        b.nameplate_capacity = Param(initialize=nameplate_capacity, mutable=True)

        maintenance_percentages = {
            1: [0.4, 0.016],
            2: [0.4, 0.016],
            3: [0.35, 0.03],
            4: [0.35, 0.03],
            5: [0.35, 0.03],
            6: [0.4, 0.019],
            7: [0.4, 0.016],
        }

        b.maintenance_labor_TPC_split = Param(
            initialize=maintenance_percentages[tech][0], mutable=True
        )
        b.maintenance_labor_percent = Param(
            initialize=maintenance_percentages[tech][1], mutable=True
        )
        b.maintenance_material_TPC_split = Param(
            initialize=(1 - maintenance_percentages[tech][0]), mutable=True
        )
        b.maintenance_material_percent = Param(
            initialize=maintenance_percentages[tech][1], mutable=True
        )

        # make vars
        b.annual_operating_labor_cost = Var(
            initialize=1, bounds=(0, 1e4), doc="annual labor cost in $MM/yr"
        )
        b.maintenance_labor_cost = Var(
            initialize=1, bounds=(0, 1e4), doc="maintenance labor cost in $MM/yr"
        )
        b.admin_and_support_labor_cost = Var(
            initialize=1, bounds=(0, 1e4), doc="admin and support labor cost in $MM/yr"
        )
        b.property_taxes_and_insurance = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="property taxes and insurance cost in $MM/yr",
        )
        b.total_fixed_OM_cost = Var(
            initialize=4, bounds=(0, 1e4), doc="total fixed O&M costs in $MM/yr"
        )

        # variable for user to assign other fixed costs to, fixed to 0 by default
        b.other_fixed_costs = Var(
            initialize=0, bounds=(0, 1e4), doc="other fixed costs in $MM/yr"
        )
        b.other_fixed_costs.fix(0)

        # maintenance material cost is technically a variable cost, but it makes
        # more sense to include with the fixed costs becuase it uses TPC
        b.maintenance_material_cost = Var(
            initialize=5, bounds=(0, 1e4), doc="cost of maintenance materials in $/MWh"
        )

        # create constraints
        TPC = b.total_TPC  # quick reference to total_TPC

        # calculated from labor rate, labor burden, and operators per shift
        @b.Constraint()
        def annual_labor_cost_rule(c):
            return c.annual_operating_labor_cost * 1e6 == (
                c.operators_per_shift * c.labor_rate * (1 + c.labor_burden / 100) * 8760
            )

        # technology specific percentage of TPC
        @b.Constraint()
        def maintenance_labor_cost_rule(c):
            return c.maintenance_labor_cost == (
                TPC * c.maintenance_labor_TPC_split * c.maintenance_labor_percent
            )

        # 25% of the sum of annual operating labor and maintenance labor
        @b.Constraint()
        def admin_and_support_labor_cost_rule(c):
            return c.admin_and_support_labor_cost == (
                0.25 * (c.annual_operating_labor_cost + c.maintenance_labor_cost)
            )

        # 2% of TPC
        @b.Constraint()
        def taxes_and_insurance_cost_rule(c):
            return c.property_taxes_and_insurance == 0.02 * TPC

        # sum of fixed O&M costs
        @b.Constraint()
        def total_fixed_OM_cost_rule(c):
            return c.total_fixed_OM_cost == (
                c.annual_operating_labor_cost
                + c.maintenance_labor_cost
                + c.admin_and_support_labor_cost
                + c.property_taxes_and_insurance
                + c.other_fixed_costs
            )
        power_plant=False
        # technology specific percentage of TPC
        @b.Constraint()
        def maintenance_material_cost_rule(c):
            if power_plant is True:
                return c.maintenance_material_cost == (
                    TPC
                    * c.maintenance_material_TPC_split
                    * c.maintenance_material_percent
                    / 0.85
                    / c.nameplate_capacity
                    / 8760
                )
            else:
                return c.maintenance_material_cost == (
                    TPC
                    * c.maintenance_material_TPC_split
                    * c.maintenance_material_percent
                    / 0.85
                    * (85/100)
                )

    def get_variable_OM_costs(self, production_rate, resources, rates, prices={}):
        """
        This function is used to calculate the variable cost of producing either
        electricity in $/MWh or hydrogen in $/kg. The function is structured to
        allow the user to generate all fuel, consumable, and waste disposal costs
        at once. A total variable cost is created for each point in fs.time.

        Args:
            fs: pyomo flowsheet block
            production_rate: pyomo var indexed by fs.time representing the net
                system power or the hydrogen production rate
            resources: a list of strings for the resorces to be costed
            rates: a list of pyomo vars for resource consumption rates
            prices: a dict of resource prices to be added to the premade dictionary

        Returns:
            None.

        """

        # check the units on production_rate
        if isinstance(production_rate, Expression) and production_rate.is_indexed:
            production_units = pyunits.get_units(
                production_rate[production_rate.index_set().first()]
            )
        else:
            production_units = pyunits.get_units(production_rate)

        try:
            pyunits.convert(production_units, pyunits.MW)
            mode = "power"
            cost_units = pyunits.USD_2018 / pyunits.MWh
            unit_tag = "$/MWh"
        except InconsistentUnitsError:
            try:
                pyunits.convert(production_units, pyunits.kg / pyunits.s)
                mode = "hydrogen"
                cost_units = pyunits.USD_2018 / pyunits.kg
                unit_tag = "$/kg"
            except InconsistentUnitsError:
                raise Exception(
                    "units not compatable, make sure production rate"
                    "has dimensions of power or mass flowrate"
                )

        # assert arguments are correct types
        if type(resources) is not list:
            raise TypeError("resources argument must be a list")
        if type(rates) is not list:
            raise TypeError("rates argument must be a list")
        if type(prices) is not dict:
            raise TypeError("prices argument must be a dictionary")

        # assert lists are the same length
        if len(resources) != len(rates):
            raise Exception("resources and rates must be lists of the same length")

        # # check if flowsheet level costing block exists
        # if not hasattr(fs, "costing"):
        #     self.get_costing(year="2018")

        # dictionary of default prices
        default_prices = {
            "natural_gas": 4.42 * pyunits.USD_2018 / pyunits.MBtu,  # $/MMbtu
            "coal": 51.96 * pyunits.USD_2018 / pyunits.ton,
            "water": 1.90e-3 * pyunits.USD_2018 / pyunits.gallon,
            "water_treatment_chemicals": 550 * pyunits.USD_2018 / pyunits.ton,
            "ammonia": 300 * pyunits.USD_2018 / pyunits.ton,
            "SCR_catalyst": 150 * pyunits.USD_2018 / pyunits.ft**3,
            "triethylene_glycol": 6.80 * pyunits.USD_2018 / pyunits.gallon,
            "SCR_catalyst_waste": 2.50 * pyunits.USD_2018 / pyunits.ft**3,
            "triethylene_glycol_waste": 0.35 * pyunits.USD_2018 / pyunits.gallon,
            "amine_purification_unit waste": 38 * pyunits.USD_2018 / pyunits.ton,
            "thermal_reclaimer_unit_waste": 38 * pyunits.USD_2018 / pyunits.ton,
        }

        # add entrys from prices to defualt_prices
        for key in prices.keys():
            default_prices[key] = prices[key]

        # raise error if the user included a resource not in default_prices
        if not set(resources).issubset(default_prices.keys()):
            raise Exception(
                "A resource was included that does not contain a "
                "price. Prices exist for the following resources: "
                "{}".format(list(default_prices.keys()))
            )
        # create list of prices
        prices = [default_prices[r] for r in resources]

        # zip rates and prices into a dict accessible by resource
        resource_rates = dict(zip(resources, rates))
        resource_prices = dict(zip(resources, prices))

        # if costing power make vars in fs.costing, if costing hydrogen make vars
        # in fs.H2_costing
        # if mode == "power":
        #     costing = self.costing
        # elif mode == "hydrogen":
        #     self.H2_costing = Block()
        #     costing = self.H2_costing

        # make vars
        self.variable_operating_costs = Var(
            self.parent_block().time,
            resources,
            initialize=1,
            doc="variable operating costs",
            units=cost_units,
        )

        self.other_variable_costs = Var(
            self.parent_block().time,
            initialize=0,
            doc="a variable to include non-standard O&M costs",
            units=cost_units,
        )
        # assume the user is not using this
        self.other_variable_costs.fix(0)

        self.total_variable_OM_cost = Var(
            self.parent_block().time,
            initialize=1,
            doc="total variable operating and maintenance costs",
            units=cost_units,
        )

        # # make constraints
        # if mode == "power":

        @self.Constraint(self.parent_block().time, resources)
        def variable_cost_rule_power(c, t, r):
            # return costing.variable_operating_costs[t, r] == (
            #     pyunits.convert(
            #         resource_prices[r] * resource_rates[r][t] / production_rate[t],
            #         cost_units,
            #     )
            # )

            return (
                c.variable_operating_costs[t, r]
                == (
                    pyunits.convert(
                        resource_prices[r] * resource_rates[r][t] * 365 * 85 / 100,
                        pyunits.USD_2018,
                    )
                )
                / 1e6
            )

        # elif mode == "hydrogen":

        # @self.Constraint(self.time, resources)
        # def variable_cost_rule_hydrogen(c, t, r):
        #     # return costing.variable_operating_costs[t, r] == (
        #     #     pyunits.convert(
        #     #         resource_prices[r] * resource_rates[r][t] / production_rate[t],
        #     #         cost_units,
        #     #     )
        #     # )

        #     return (
        #         costing.variable_operating_costs[t, r]
        #         == (
        #             pyunits.convert(
        #                 resource_prices[r] * resource_rates[r][t] * 365 * 85 / 100,
        #                 pyunits.USD_2018 / pyunits.s,
        #             )
        #         )
        #         / 1e6
        #     )

        # if mode == "power" and hasattr(self, "maintenance_material_cost"):

        @self.Constraint(self.parent_block().time)
        def total_variable_cost_rule_power(c, t):
            return (
                c.total_variable_OM_cost[t]
                == sum(c.variable_operating_costs[t, r] for r in resources)
                + c.maintenance_material_cost
                + c.other_variable_costs[t]
            )

        # else:

        # @self.Constraint(self.time)
        # def total_variable_cost_rule_hydrogen(c, t):
        #     return (
        #         c.H2_costing.total_variable_OM_cost[t]
        #         == sum(
        #             c.H2_costing.variable_operating_costs[t, r] for r in resources
        #         )
        #         + c.H2_costing.other_variable_costs[t]
        #     )

        # if mode == "power":
        #     _log.warning(
        #         "The variable fs.costing.maintenance_material_cost "
        #         "could not be found, total_variable_cost_rule was "
        #         "constructed without in. get_fixed_OM_costs should be"
        #         " called before get_variable_OM_costs"
        #     )

    def initialize_fixed_OM_costs(b):
        # This method accepts either a concrete model or a flowsheet block
        if hasattr(b, "costing") and hasattr(b, "total_fixed_OM_cost"):

            calculate_variable_from_constraint(
                b.annual_operating_labor_cost, b.annual_labor_cost_rule
            )

            calculate_variable_from_constraint(
                b.maintenance_labor_cost, b.maintenance_labor_cost_rule
            )

            calculate_variable_from_constraint(
                b.admin_and_support_labor_cost,
                b.admin_and_support_labor_cost_rule,
            )

            calculate_variable_from_constraint(
                b.property_taxes_and_insurance,
                b.taxes_and_insurance_cost_rule,
            )

            calculate_variable_from_constraint(
                b.total_fixed_OM_cost, b.total_fixed_OM_cost_rule
            )

            calculate_variable_from_constraint(
                b.maintenance_material_cost,
                b.maintenance_material_cost_rule,
            )

    def initialize_variable_OM_costs(fs):
        # This method accepts only a flowsheet block
        # initialization for power generation costs
        if hasattr(fs, "costing") and hasattr(fs.costing, "variable_operating_costs"):

            for i in fs.costing.variable_operating_costs.keys():
                if hasattr(fs.costing, "variable_cost_rule_power"):
                    calculate_variable_from_constraint(
                        fs.costing.variable_operating_costs[i],
                        fs.costing.variable_cost_rule_power[i],
                    )

            for i in fs.costing.total_variable_OM_cost.keys():
                calculate_variable_from_constraint(
                    fs.costing.total_variable_OM_cost[i],
                    fs.total_variable_cost_rule_power[i],
                )

        # initialization for H2 production costs
        if hasattr(fs, "H2_costing") and hasattr(
            fs.H2_costing, "variable_operating_costs"
        ):

            if hasattr(fs.costing, "variable_cost_rule_hydrogen"):
                for i in fs.H2_costing.variable_operating_costs.keys():
                    calculate_variable_from_constraint(
                        fs.H2_costing.variable_operating_costs[i],
                        fs.H2_costing.variable_cost_rule_hydrogen[i],
                    )

            for i in fs.H2_costing.total_variable_OM_cost.keys():
                calculate_variable_from_constraint(
                    fs.H2_costing.total_variable_OM_cost[i],
                    fs.total_variable_cost_rule_hydrogen[i],
                )

    # -----------------------------------------------------------------------------
    # Costing Library Utility Functions
    # -----------------------------------------------------------------------------

    def costing_initialization(fs):
        for o in fs.component_objects(descend_into=False):
            # look for costing blocks
            if hasattr(o, "costing") and hasattr(o.costing, "library"):
                if o.costing.library == "sCO2":

                    if o.costing.equipment in [
                        "Axial turbine",
                        "Radial turbine",
                        "Coal-fired heater",
                        "Natural gas-fired heater",
                        "Recouperator",
                    ]:
                        calculate_variable_from_constraint(
                            o.costing.temperature, o.costing.temp_eq
                        )
                        calculate_variable_from_constraint(
                            o.costing.temp_factor, o.costing.temp_correction_eq
                        )

                    calculate_variable_from_constraint(
                        o.costing.scaled_param, o.costing.scaled_param_eq
                    )
                    calculate_variable_from_constraint(
                        o.costing.equipment_cost, o.costing.equipment_cost_eq
                    )
                    calculate_variable_from_constraint(
                        o.costing.bare_erected_cost, o.costing.bare_erected_cost_eq
                    )
                    calculate_variable_from_constraint(
                        o.costing.total_plant_cost, o.costing.total_plant_cost_eq
                    )

                elif o.costing.library in ["PP", "ASU"]:
                    for key in o.costing.bare_erected_cost.keys():
                        calculate_variable_from_constraint(
                            o.costing.bare_erected_cost[key],
                            o.costing.bare_erected_cost_eq[key],
                        )
                        calculate_variable_from_constraint(
                            o.costing.total_plant_cost[key],
                            o.costing.total_plant_cost_eq[key],
                        )

    def display_total_plant_costs(fs):
        print("-----Total Plant Costs-----")
        for o in fs.component_objects(descend_into=False):
            # look for costing blocks
            if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
                print(
                    "%s: $%.2f Million" % (value(o.name), value(o.costing.total_cost))
                )

    def display_bare_erected_costs(fs):
        print("-----Bare Erected Costs-----")
        for o in fs.component_objects(descend_into=False):
            # look for costing blocks
            if hasattr(o, "costing") and hasattr(o.costing, "bare_erected_cost"):
                print(
                    "%s: $%.2f Million"
                    % (value(o.name), value(o.costing.bare_erected_cost))
                )

    def display_equipment_costs(fs):
        print("-----Equipment Costs-----")
        for o in fs.component_objects(descend_into=False):
            # look for costing blocks
            if hasattr(o, "costing") and hasattr(o.costing, "equipment_cost"):
                print(
                    "%s: $%.2f Million"
                    % (value(o.name), value(o.costing.equipment_cost))
                )

    def get_total_TPC(b):
        # This method accepts either a concrete model or a flowsheet block
        # isinstance( m.fs.boiler.costing, QGESSCostingData)
        TPC_list = []
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if hasattr(o, "costing") and hasattr(o.costing, "total_plant_cost"):
                for key in o.costing.total_plant_cost.keys():
                    TPC_list.append(o.costing.total_plant_cost[key])

        b.total_TPC = Var(initialize=0, bounds=(0, 1e4), doc="total TPC in $MM")

        @b.Constraint()
        def total_TPC_eq(c):
            return c.total_TPC == sum(TPC_list)

    def display_flowsheet_cost(b):
        print("\n")
        print("Total flowsheet cost: $%.3f Million" % value(b.flowsheet_cost))

    def check_sCO2_costing_bounds(fs):
        # interate through the children of the flowsheet
        for o in fs.component_objects(descend_into=False):
            # look for costing blocks
            if hasattr(o, "costing"):
                costing = o.costing
                if costing.library == "sCO2":
                    equipment = costing.equipment
                    lower_bound = sCO2_costing_params[equipment]["Lower Bound"]
                    upper_bound = sCO2_costing_params[equipment]["Upper Bound"]
                    if value(costing.scaled_param) < lower_bound:
                        print(
                            """%s: The scaled parameter (%f) is below the lower
                            bound (%f)."""
                            % (value(o.name), value(costing.scaled_param), lower_bound)
                        )
                    elif value(costing.scaled_param) > upper_bound:
                        print(
                            """%s: The scaled parameter (%f) is above the upper
                            bound (%f)."""
                            % (value(o.name), value(costing.scaled_param), upper_bound)
                        )
                    else:
                        print(
                            """%s: The scaled parameter is within the bounds."""
                            % value(o.name)
                        )
