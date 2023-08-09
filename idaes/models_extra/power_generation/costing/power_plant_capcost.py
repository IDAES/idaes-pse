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
    - check_sCO2_costing_bounds() to display a warning if costing model have
      been used outside the range that where designed for
"""
# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

__author__ = (
    "Costing Team (A. Noring, A. Deshpande, B. Paul, D. Caballero, and M. Zamarripa)"
)
__version__ = "1.0.0"

from sys import stdout
import textwrap

from pandas import DataFrame

from pyomo.environ import (
    Param,
    Var,
    Constraint,
    Expression,
    value,
    Expr_if,
    units as pyunits,
)
from pyomo.core.base.expression import ScalarExpression
from pyomo.core.base.units_container import InconsistentUnitsError, UnitsError
from pyomo.util.calc_var_value import calculate_variable_from_constraint

import idaes.core.util.scaling as iscale
from idaes.core import (
    FlowsheetCostingBlockData,
    register_idaes_currency_units,
)
from idaes.models_extra.power_generation.costing.costing_dictionaries import (
    load_BB_costing_dictionary,
    load_sCO2_costing_dictionary,
    load_REE_costing_dictionary,
)
from idaes.models_extra.power_generation.costing.generic_ccs_capcost_custom_dict import (
    load_generic_ccs_costing_dictionary,
)

from idaes.core.util.tables import stream_table_dataframe_to_string
import idaes.logger as idaeslog
from idaes.core import declare_process_block_class

_log = idaeslog.getLogger(__name__)

# -----------------------------------------------------------------------------
# Power Plant Costing Library
# -----------------------------------------------------------------------------


def custom_power_plant_currency_units():
    """
    Define conversion rates for US Dollars based on CE Index.
    """
    register_idaes_currency_units()
    if (
        "USD_2008_Nov" in pyunits._pint_registry  # pylint: disable=protected-access
        and "USD_2019_Sep" in pyunits._pint_registry  # pylint: disable=protected-access
    ):
        # Assume that custom power plant units have already been registered
        # Log a message and end
        _log.debug(
            "Custom power plant currency units (USD_2008_Nov, USD_2019_Sep) "
            "already appear in Pyomo unit registry. Assuming repeated "
            "call of custom_power_plant_currency_units."
        )
    else:
        pyunits.load_definitions_from_strings(
            [
                "USD_2008_Nov = 500/566.2 * USD_CE500",
                "USD_2019_Sep = 500/599.3 * USD_CE500",
            ]
        )


@declare_process_block_class("QGESSCosting")
class QGESSCostingData(FlowsheetCostingBlockData):
    # Register currency and conversion rates based on CE Index
    # register_idaes_currency_units()
    custom_power_plant_currency_units()

    def build_global_params(self):
        """
        This is where we can declare any global parameters we need, such as
        Lang factors, or coefficients for costing methods that should be
        shared across the process.

        You can do what you want here, so you could have e.g. sub-Blocks
        for each costing method to separate the parameters for each method.
        """
        # Set the base year for all costs
        self.base_currency = pyunits.USD_2018
        # Set a base period for all operating costs
        self.base_period = pyunits.year

    def build_process_costs(
        self,
        # arguments related to Fixed OM costs
        total_plant_cost=None,
        nameplate_capacity=650,
        capacity_factor=0.85,
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
        chemicals_inventory=None,
        waste=None,
        transport_cost=None,
        tonne_CO2_capture=None,
        CE_index_year="2018",
    ):
        """
        This method builds process-wide costing, including fixed and variable
        operating & maintenance costs, costs of production, cost of
        electricity and cost of capture.

        Args:
            total_plant_cost: The TPC in $MM that will be used to determine fixed O&M,
                costs. If the value is None, the function will try to use the
                TPC calculated from the individual units. This quantity should
                be a Pyomo Var or Param that will contain the TPC value.
            nameplate_capacity: rated plant output in MW
            capacity_factor: multiplicative factor for normal operating
                capacity
            labor_rate: hourly rate of plant operators in project dollar year
            labor_burden: a percentage multiplier used to estimate non-salary
                labor expenses
            operators_per_shift: average number of operators per shift
            tech: int 1-7 representing the categories in get_PP_costing, used
                to determine maintenance costs
            land_cost: Expression, Var or Param to calculate land costs
            net_power: actual plant output in MW, only required if calculating
                variable costs
            resources: list setting resources to cost
            rates: list setting flow rates of resources
            prices: list setting prices of resources
            fixed_OM: True/False flag for calculating fixed O&M costs
            variable_OM: True/False flag for calculating variable O&M costs
            fuel: string setting fuel type for fuel costs
            chemicals: string setting chemicals type for chemicals costs
            chemicals_inventory: string setting chemicals type for inventory costs
            waste: string setting waste type for waste costs
            tonne_CO2_capture: Var or value to use for tonnes of CO2 capture
                in one year
            transport_cost: Expression, Var or Param to use for transport costs
                per ton of CO2 captured (note, this is not part of the TOC)
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars
        """

        try:
            CE_index_units = getattr(
                pyunits, "MUSD_" + CE_index_year
            )  # millions of USD, for base year
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        self.capacity_factor = Param(
            initialize=capacity_factor,
            mutable=True,
            doc="capacity factor of the plant",
            units=pyunits.dimensionless,
        )
        if net_power is not None:
            if not hasattr(self, "net_power"):
                self.net_power = Param(
                    self.parent_block().time,
                    initialize=net_power,
                    mutable=True,
                    units=net_power.get_units(),
                )
        else:
            self.net_power = net_power  # None

        if total_plant_cost is None:
            self.get_total_TPC(CE_index_year)

        if fixed_OM is True:
            self.get_fixed_OM_costs(
                net_power=self.net_power,
                nameplate_capacity=nameplate_capacity,
                capacity_factor=self.capacity_factor,
                labor_rate=labor_rate,
                labor_burden=labor_burden,
                operators_per_shift=operators_per_shift,
                tech=tech,
                fixed_TPC=total_plant_cost,
                CE_index_year=CE_index_year,
            )

        if variable_OM is True:
            self.get_variable_OM_costs(
                resources=resources,
                rates=rates,
                prices=prices,
                capacity_factor=self.capacity_factor,
                CE_index_year=CE_index_year,
            )

        # build system costs (owner's, total overnight costs, annualized costs,
        # COE, and cost of capture)
        if fixed_OM and variable_OM:
            # total overnight cost requires fixed costs
            self.pct_TPC = Param(
                initialize=20.2 / 100, doc="Fixed percentage for other owners cost"
            )
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
                    expr=self.variable_operating_costs[0, fuel] / 12 * 2.25,
                    doc="Owner's costs - 2.25 months of fuel costs",
                )
                non_fuel_resources.remove(fuel)  # remove fuel from the list

            if waste is not None:
                self.waste_costs_OC = Expression(
                    expr=(sum(self.variable_operating_costs[0, i] for i in waste) / 12)
                )
            if chemicals is not None:
                self.chemical_costs_OC = Expression(
                    expr=(
                        sum(self.variable_operating_costs[0, i] for i in chemicals)
                        / 2  # six months of chemicals
                    )
                )

            if chemicals_inventory is not None:
                if total_plant_cost is not None:
                    inventory_plant_basis = total_plant_cost
                else:
                    inventory_plant_basis = self.total_TPC
                self.chemical_inventory_costs_OC = Expression(
                    expr=(
                        (
                            sum(
                                self.variable_operating_costs[0, i]
                                for i in chemicals_inventory
                            )
                            / 6
                        )
                        * pyunits.year  # two months of chemicals inventory
                        + 0.005 * inventory_plant_basis
                    )
                )

            self.non_fuel_and_waste_OC = Expression(
                expr=(
                    sum(self.variable_operating_costs[0, i] for i in non_fuel_resources)
                    / 12
                )
            )

            if land_cost is not None:
                if isinstance(land_cost, (Expression, ScalarExpression)) or (
                    isinstance(land_cost, (Param, Var)) and land_cost.get_units is None
                ):
                    self.land_cost = land_cost * CE_index_units
                else:
                    self.land_cost = land_cost

            self.total_overnight_capital = Expression(
                expr=self.total_TPC
                # pre production costs
                + self.six_month_labor
                + (
                    self.chemical_inventory_costs_OC
                    if chemicals_inventory is not None
                    else 0 * CE_index_units
                )  # Initial Cost for Catalyst and Chemicals Inventory
                + (
                    self.maintenance_material_cost / 12  # 1 month materials
                    + self.non_fuel_and_waste_OC  # 1 month nonfuel consumables
                    + (
                        self.waste_costs_OC if waste is not None else 0 * CE_index_units
                    )  # 1 month waste
                    # inventory capital costs
                    + (
                        self.fuel_cost_OC if fuel is not None else 0 * CE_index_units
                    )  # 60 day fuel supply
                    # Other costs
                    + (
                        self.chemical_costs_OC
                        if chemicals is not None
                        else 0 * CE_index_units
                    )  # Initial Cost for Catalyst and Chemicals
                )
                * 1
                * pyunits.year  # variable costs for 1 year
                + (self.land_cost if land_cost is not None else 0 * CE_index_units)
                + self.total_TPC
                * self.pct_TPC  # other owners costs (other + spare parts
                # + financing + other pre-production)
            )

            self.tasc_toc_factor = Param(
                initialize=1.093,
                mutable=True,
                doc="TASC/TOC factor from Exhibit 3-7 reference 1, real for"
                + " three years 1.093, real for five years 1.154, nominal for"
                + " three years 1.242, nominal for five years 1.289",
            )

            self.total_as_spent_cost = Expression(
                expr=self.total_overnight_capital * self.tasc_toc_factor
            )

            self.fixed_charge_factor = Param(
                initialize=0.0707,
                mutable=True,
                doc="Fixed charge rate from Exhibit 3-5 based on CRF values,"
                + " real for three/five years 0.0707, nominal for three/five"
                + " years 0.0886",
            )
            self.annualized_cost = Expression(
                expr=self.fixed_charge_factor * self.total_as_spent_cost
            )

            self.additional_cost_of_electricity = Var(
                initialize=70.9e-6,
                doc="additional cost to be added to the COE calculations"
                + " in millions",
                units=CE_index_units / pyunits.MWh,
            )

            # only build COE for power plants
            if net_power is not None:
                self.cost_of_electricity = Expression(
                    expr=(
                        (
                            self.annualized_cost / pyunits.year
                            + self.total_fixed_OM_cost / pyunits.year
                            + self.total_variable_OM_cost[0] * self.capacity_factor
                        )
                        / (
                            self.capacity_factor
                            * self.net_power[0]
                            * 8760
                            * pyunits.hr
                            / pyunits.year
                        )
                        + self.additional_cost_of_electricity
                    )
                )
            if tonne_CO2_capture is not None:
                if not hasattr(self, "tonne_CO2_capture"):
                    self.tonne_CO2_capture = Param(
                        initialize=tonne_CO2_capture, mutable=True, units=pyunits.ton
                    )
                self.cost_of_capture = Expression(
                    expr=(
                        (
                            self.annualized_cost / pyunits.year
                            + self.total_fixed_OM_cost / pyunits.year
                            + self.total_variable_OM_cost[0] * self.capacity_factor
                        )
                        / self.tonne_CO2_capture
                    )
                )

                if transport_cost is not None:
                    if isinstance(transport_cost, (Expression, ScalarExpression)) or (
                        isinstance(transport_cost, (Param, Var))
                        and transport_cost.get_units is None
                    ):
                        self.transport_cost = (
                            transport_cost
                            * CE_index_units
                            / pyunits.ton
                            * self.tonne_CO2_capture
                        )
                    else:
                        self.transport_cost = transport_cost * self.tonne_CO2_capture

            else:  # except the case where transport_cost is passed but tonne_CO2_capture is not passed
                if transport_cost is not None:
                    raise Exception(
                        "If a transport_cost is not None, "
                        "tonne_CO2_capture cannot be None."
                    )

    @staticmethod
    def initialize_build(*args, **kwargs):
        """
        Here we can add initialization steps for the things we built in
        build_process_costs.

        Note that the aggregate costs will be initialized by the framework.
        """
        # TODO: For now,  no additional process level costs to initialize

    def report(self, export=False):
        var_dict = {}

        if hasattr(self, "total_TPC"):
            var_dict["Total TPC [$MM]"] = value(self.total_TPC)

        if hasattr(self, "total_overnight_capital"):
            var_dict["Total Overnight Cost [$MM]"] = value(self.total_overnight_capital)

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
            var_dict["Total variable O&M cost full capacity [$MM/year]"] = value(
                self.total_variable_OM_cost[0]
            )
            var_dict["Total variable O&M cost operating capacity [$MM/year]"] = value(
                self.total_variable_OM_cost[0] * self.capacity_factor
            )

        if (
            hasattr(self, "annualized_cost")
            and hasattr(self, "total_fixed_OM_cost")
            and hasattr(self, "total_variable_OM_cost")
        ):
            var_dict["Total Annualized Cost [$MM/year]"] = (
                value(self.annualized_cost)
                + value(self.total_fixed_OM_cost)
                + self.capacity_factor * value(self.total_variable_OM_cost[0])
            )

        if hasattr(self, "cost_of_electricity"):
            var_dict["Cost of Electricity [$/MWh]"] = value(
                self.cost_of_electricity * 1e6
            )

        if hasattr(self, "cost_of_capture"):
            var_dict["Cost of Capture [$/tonne]"] = value(self.cost_of_capture * 1e6)

        if hasattr(self, "transport_cost"):
            var_dict["Total Transport Cost [$MM]"] = value(self.transport_cost)

        # REE-specific results
        if self.library == "REE":
            if hasattr(self, "total_plant_cost"):
                var_dict["Total Plant Cost [$MM]"] = value(self.total_plant_cost)

            if hasattr(self, "bare_erected_cost"):
                var_dict["Total Bare Erected Cost [$MM]"] = value(
                    self.bare_erected_cost
                )

            if hasattr(self, "total_installation_cost"):
                var_dict["Total Installation Cost [$MM]"] = value(
                    self.total_installation_cost
                )

            if hasattr(self, "ancillary_costs"):
                var_dict["Total Ancillary Installation Cost [$MM]"] = value(
                    self.ancillary_costs
                )

            if hasattr(self, "piping_materials_and_labor_costs"):
                var_dict[
                    "Total Ancillary Piping, Materials and Labor Installation Cost [$MM]"
                ] = value(self.piping_materials_and_labor_costs)

            if hasattr(self, "electrical_materials_and_labor_costs"):
                var_dict[
                    "Total Ancillary Electrical, Materials and Labor Installation Cost [$MM]"
                ] = value(self.electrical_materials_and_labor_costs)

            if hasattr(self, "instrumentation_costs"):
                var_dict[
                    "Total Ancillary Instrumentation Installation Cost [$MM]"
                ] = value(self.instrumentation_costs)

            if hasattr(self, "plant_services_costs"):
                var_dict[
                    "Total Ancillary Plant Services Installation Cost [$MM]"
                ] = value(self.plant_services_costs)

            if hasattr(self, "buildings_costs"):
                var_dict["Total Buildings Installation Cost [$MM]"] = value(
                    self.buildings_costs
                )

            if hasattr(self, "process_buildings_costs"):
                var_dict["Total Process Buildings Installation Cost [$MM]"] = value(
                    self.process_buildings_costs
                )

            if hasattr(self, "auxiliary_buildings_costs"):
                var_dict["Total Auxiliary Buildings Installation Cost [$MM]"] = value(
                    self.auxiliary_buildings_costs
                )

            if hasattr(self, "site_improvements_costs"):
                var_dict[
                    "Total Site Improvements Buildings Installation Cost [$MM]"
                ] = value(self.site_improvements_costs)

            if hasattr(self, "epcm_costs"):
                var_dict["Total EPCM Installation Cost [$MM]"] = value(self.epcm_costs)

            if hasattr(self, "equipment_installation_costs"):
                var_dict[
                    "Total Equipment Installation EPCM Installation Cost [$MM]"
                ] = value(self.equipment_installation_costs)

            if hasattr(self, "field_expenses_costs"):
                var_dict["Total Field Expenses EPCM Cost [$MM]"] = value(
                    self.field_expenses_costs
                )

            if hasattr(self, "project_management_and_construction_costs"):
                var_dict[
                    "Total Project Management and Construction EPCM Installation Cost [$MM]"
                ] = value(self.project_management_and_construction_costs)

            if hasattr(self, "process_contingency_costs"):
                var_dict["Total Process Contingency Installation Cost [$MM]"] = value(
                    self.process_contingency_costs
                )

            if hasattr(self, "contingency_costs"):
                var_dict["Total Contingency Installation Cost [$MM]"] = value(
                    self.contingency_costs
                )

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
        tech,
        ccs="B",
        n_equip=1,
        scale_down_parallel_equip=False,
        CE_index_year="2018",
        additional_costing_params=None,
        use_additional_costing_params=False,
    ):
        """
        Power Plant Costing Method
        This method relies on the capital cost scaling methodologies developed
        by NETL. Report #DOE/NETL-341/013113
        Multiple vendors quotes have been used to determine the cost of several
        plant equipment (i.e. boiler, pumps, heat exchangers, etc.), other
        cost incurred during the plant operation (i.e. solids handling, etc.)

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

        The scaled cost is computed using ref values for different
        technologies.
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

        This method computes the capital cost of units and main components of
        the power plant, and requires a few arguments to build a constraint as
        part of your main model.

        Args:
            blk: A unit-level costing block where costing variables and
                constraints can be added to
            cost_accounts: A list of accounts to be included in the total cost
            scaled_param: the process parameter for the system(s) being costed;
                this is the total flow for all parallel trains of the system(s)
            tech: integer 1-7 representing the above categories
            ccs: 'A' or 'B' representing no CCS or CCS
            n_equip: Integer number of parallel equipment trains for unit
                operations; for example, enter '5' if a feed will be split
                among 5 identical units and then re-mixed
            scale_down_parallel_equip: Boolean flag whether to scale down
                parallel equipment trains, e.g. two trains scaled down will
                each be half the size/capacity of a single train, and two
                trains not scaled down will each be the same size as a single
                train (twice the capacity).
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars
            additional_costing_params: user-defined dictionary to append to
                existing cost accounts dictionary
            use_additional_costing_params: Boolean flag to use additional
                costing parameters when account names conflict with existing
                accounts data

        The appropriate scaling parameters for various cost accounts can be
        found in the QGESS on capital cost scaling (Report
        #DOE/NETL-341/013113).
        The correct units for the reference parameters are found in the BBR4
        COE spreadsheet.

        """
        # ------------------------ Power Plant Cost ------------------------

        try:
            CE_index_units = getattr(
                pyunits, "MUSD_" + CE_index_year
            )  # millions of USD, for base year
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

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
        if isinstance(cost_accounts, str):
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
        reference_cost_units = {}
        reference_costs_init = {}
        reference_params = {}
        cost_scaling_fractions = {}
        engineering_fees = {}
        process_contingencies = {}
        project_contingencies = {}

        # load the cost dictionaries and add custom accounts

        # load BB costing dictionary
        BB_costing_params = load_BB_costing_dictionary()

        # load generic ccs costing dictionary
        generic_ccs_costing_params = load_generic_ccs_costing_dictionary()

        # for compatibility with potential custom accounts, the loop handles
        # new technologies, new CCS types for existing technologies, and new
        # accounts for existing technology-CCS type pairs
        # Users should not be adding new entries for existing accounts

        if additional_costing_params is not None and additional_costing_params != {}:
            accounts_to_merge = [generic_ccs_costing_params, additional_costing_params]
        else:
            accounts_to_merge = [generic_ccs_costing_params]

        costing_params = BB_costing_params  # initialize with baseline accounts
        for (
            new_costing_params
        ) in accounts_to_merge:  # merge new dictionaries sequentially
            # adding ccs and any provided custom params to the ccs dictionary
            # need to "freeze" dict so it is hashable for merging keys
            frozen_dict = {**costing_params}
            for techkey, techval in new_costing_params.items():
                if (
                    techkey in frozen_dict.keys()
                ):  # if techkey already exists, append any new ccs types
                    for ccskey, ccsval in new_costing_params[techkey].items():
                        if (
                            ccskey in frozen_dict[techkey].keys()
                        ):  # if ccskey already exists, append any new accounts
                            for accountkey, accountval in new_costing_params[techkey][
                                ccskey
                            ].items():
                                if (
                                    accountkey in frozen_dict[techkey][ccskey].keys()
                                ) and not use_additional_costing_params:
                                    if accountkey not in cost_accounts:
                                        pass  # not the current account, don't fail here
                                    else:  # this is not allowed
                                        raise ValueError(
                                            "Data already exists for Account {} "
                                            "using technology {} with CCS {}. "
                                            "Please confirm that the custom "
                                            "account dictionary is correct, or "
                                            "add the new parameters as a new "
                                            "account. To use the custom account "
                                            "dictionary for all conflicts, please "
                                            "pass the argument use_additional_costing_params "
                                            "as True.".format(
                                                accountkey, str(techkey), ccskey
                                            )
                                        )
                                else:  # conflict is the account passed, and overwrite it
                                    frozen_dict[techkey][ccskey][
                                        accountkey
                                    ] = accountval
                        else:  # it's a new type, append the entry
                            frozen_dict[techkey][ccskey] = ccsval
                else:
                    frozen_dict[techkey] = techval
            costing_params = {k: frozen_dict[k] for k in sorted(frozen_dict)}

        for account in cost_accounts:
            try:  # look for data in json file info
                process_params[account] = costing_params[str(tech)][ccs][account][
                    "Process Parameter"
                ]
                reference_units[account] = costing_params[str(tech)][ccs][
                    cost_accounts[0]
                ]["Units"]
                account_names[account] = costing_params[str(tech)][ccs][account][
                    "Account Name"
                ]
                exponents[account] = float(
                    costing_params[str(tech)][ccs][account]["Exponent"]
                )
                reference_costs[account] = costing_params[str(tech)][ccs][account][
                    "BEC"
                ]
                reference_cost_units[account] = costing_params[str(tech)][ccs][account][
                    "BEC_units"
                ]
                reference_costs_init[account] = (
                    costing_params[str(tech)][ccs][account]["BEC"] * 1e-3
                )

                if isinstance(process_params[account], list):
                    for i, processparam in enumerate(process_params[account]):
                        reference_params[account, processparam] = costing_params[
                            str(tech)
                        ][ccs][account]["RP Value"][i]
                        cost_scaling_fractions[account, processparam] = costing_params[
                            str(tech)
                        ][ccs][account]["Cost scaling fraction"][i]

                elif isinstance(process_params[account], str):
                    reference_params[account] = costing_params[str(tech)][ccs][account][
                        "RP Value"
                    ]

                engineering_fees[account] = costing_params[str(tech)][ccs][account][
                    "Eng Fee"
                ]
                process_contingencies[account] = costing_params[str(tech)][ccs][
                    account
                ]["Process Contingency"]
                project_contingencies[account] = costing_params[str(tech)][ccs][
                    account
                ]["Project Contingency"]
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

        # check that the user passed the correct units type and try to convert

        for account in cost_accounts:
            ref_units = reference_units[account]
            if "/" in ref_units:
                ref_units = ref_units.split("/")
                if "**" in ref_units[0]:
                    ref_units[0] = ref_units[0].split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0][0]) ** int(
                            ref_units[0][1]
                        ) / getattr(pyunits, ref_units[1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (
                                cost_accounts[0],
                                ref_units[0][0]
                                + "**"
                                + ref_units[0][1]
                                + "/"
                                + ref_units[1],
                            )
                        )
                elif "**" in ref_units[1]:
                    ref_units[1] = ref_units[1].split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1][0]
                        ) ** int(ref_units[1][1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (
                                cost_accounts[0],
                                ref_units[0]
                                + "/"
                                + ref_units[1][0]
                                + "**"
                                + ref_units[1][1],
                            )
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1]
                        )
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )

            else:
                if "**" in ref_units:
                    ref_units = ref_units.split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) ** int(ref_units[1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units)
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )

            if isinstance(scaled_param, list):
                for sp in scaled_param:
                    if sp.get_units() is None:
                        raise ValueError(
                            "Account %s uses units of %s. "
                            "Units of %s were passed. "
                            "Scaled_param must have units."
                            % (cost_accounts[0], ref_units, sp.get_units())
                        )
                    else:
                        try:
                            pyunits.convert(sp, ref_units)
                        except InconsistentUnitsError:
                            raise Exception(
                                "Account %s uses units of %s. "
                                "Units of %s were passed. "
                                "Cannot convert unit containers."
                                % (
                                    cost_accounts[0],
                                    ref_units,
                                    sp.get_units(),
                                )
                            )
            else:
                try:
                    if pyunits.get_units(scaled_param) is None:
                        raise UnitsError(
                            "Account %s uses units of %s. "
                            "Units of %s were passed. "
                            "Scaled_param must have units."
                            % (
                                cost_accounts[0],
                                ref_units,
                                pyunits.get_units(scaled_param),
                            )
                        )
                    else:
                        try:
                            pyunits.convert(scaled_param, ref_units)
                        except InconsistentUnitsError:
                            raise UnitsError(
                                "Account %s uses units of %s. "
                                "Units of %s were passed. "
                                "Cannot convert unit containers."
                                % (
                                    cost_accounts[0],
                                    ref_units,
                                    pyunits.get_units(scaled_param),
                                )
                            )
                except InconsistentUnitsError:
                    raise UnitsError(
                        f"The expression {scaled_param.name} has inconsistent units."
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
            # units not defined here, since every account could have different
            # currency units
        )

        if isinstance(process_params[cost_accounts[0]], list):
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
        elif isinstance(process_params[cost_accounts[0]], str):
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
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        blk.total_plant_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, 1e4),
            doc="total plant cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        # rule for scaling BEC
        # reference cost is in 2018 dollars, 671.1 is CE index for 2018
        def bare_erected_cost_rule(costing, i):
            ref_units = reference_units[i]
            if "/" in ref_units:
                ref_units = ref_units.split("/")
                ref_units = getattr(pyunits, ref_units[0]) / getattr(
                    pyunits, ref_units[1]
                )
            else:
                ref_units = getattr(pyunits, ref_units)

            ref_cost_units = reference_cost_units[i]
            ref_cost_units = ref_cost_units.split("$")
            if ref_cost_units[0] == "":  # no prefix
                ref_cost_units = getattr(pyunits, "USD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "K":  # thousands of $
                ref_cost_units = getattr(pyunits, "kUSD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "M":  # millions of $
                ref_cost_units = getattr(pyunits, "MUSD_" + ref_cost_units[1])

            # determine reference parameter scaler based on train scaling
            if scale_down_parallel_equip:
                scaler = n_equip
            else:
                scaler = 1

            if isinstance(process_params[i], list):
                if len(process_params[i]) > 1:
                    return costing.bare_erected_cost[i] == (
                        n_equip
                        * pyunits.convert(
                            costing.ref_cost[i] * ref_cost_units, CE_index_units
                        )
                        * sum(
                            costing.cost_scaling_fracs[i, p]
                            * (
                                pyunits.convert(scaled_param[j], ref_units)
                                / (scaler * costing.ref_param[i, p] * ref_units)
                            )
                            ** costing.exp[i]
                            for j, p in enumerate(process_params[i])
                        )
                    )
            elif isinstance(process_params[i], str):
                return costing.bare_erected_cost[i] == (
                    n_equip
                    * pyunits.convert(
                        costing.ref_cost[i] * ref_cost_units, CE_index_units
                    )
                    * (
                        pyunits.convert(scaled_param, ref_units)
                        / (scaler * costing.ref_param[i] * ref_units)
                    )
                    ** costing.exp[i]
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
        self,
        equipment,
        scaled_param,
        temp_C=None,
        n_equip=1,
        CE_index_year="2018",
        custom_accounts=None,
    ):
        """
        Args:
            self: pyomo Block where constraints will be made
            equipment: the name of the SCO2 equipment to cost
            scaled_param: the scaling parameter (in appropriate units) for the
                elected equipment
            temp: the maximum temperature of the equipment. Not all types of
                equipment use a temperature correction factor, so it is
                optional
            n_equip: the number of pieces of equipment to cost
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars
            custom_accounts: user-defined dictionary to append to
                existing cost accounts dictionary

        Cost is in M$
        """
        # check to see if a costing block already exists
        if (
            self.parent_block().name
            in self.config.flowsheet_costing_block._registered_unit_costing  # pylint: disable=protected-access
        ):
            raise AttributeError(
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(self.name)
            )

        # define costing library
        self.library = "sCO2"
        self.equipment = equipment

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # load sCO2 costing dictionary
        sCO2_costing_params = load_sCO2_costing_dictionary()

        param_dict = sCO2_costing_params[equipment]

        if custom_accounts is not None:
            # adding custom params to the dictionary
            for key, val in custom_accounts.items():
                param_dict[key] = val

        # define parameters
        self.ref_cost = Param(
            mutable=True,
            initialize=param_dict["a"],
            doc="Reference cost",
            units=pyunits.USD_2017,
        )

        self.exp = Param(
            mutable=True,
            initialize=param_dict["b"],
            doc="Scaling exponent",
            units=pyunits.dimensionless,
        )

        self.c = Param(
            mutable=True,
            initialize=param_dict["c"],
            doc="coefficient for temperature correction",
            units=pyunits.C**-1,
        )

        self.d = Param(
            mutable=True,
            initialize=param_dict["d"],
            doc="coefficient for temperature correction",
            units=pyunits.C**-2,
        )

        self.material_cost = Param(
            mutable=True,
            doc="material installation cost factor",
            initialize=param_dict["Material Cost"],
            units=pyunits.dimensionless,
        )

        self.labor_cost = Param(
            mutable=True,
            initialize=param_dict["Labor Cost"],
            doc="labor installation cost factor",
            units=pyunits.dimensionless,
        )

        # estimates for the percentages of TPC will be added later

        self.eng_fee = Param(
            mutable=True,
            initialize=0,
            doc="engineering fee percentage",
            units=pyunits.dimensionless,
        )

        self.process_conting = Param(
            mutable=True,
            initialize=0,
            doc="process contingency percentage",
            units=pyunits.dimensionless,
        )

        self.project_conting = Param(
            mutable=True,
            initialize=0,
            doc="project contingency percentage",
            units=pyunits.dimensionless,
        )

        # define variables
        # n_equip is left as a fixed variable to support MINLP optimization
        self.n_equip = Var(
            initialize=n_equip,
            doc="number of pieces of equipment",
            units=pyunits.dimensionless,
        )

        self.n_equip.fix(n_equip)

        self.scaled_param = Var(
            initialize=scaled_param,
            bounds=(0, 1e12),
            doc="scaled parameter",
            units=scaled_param.get_units(),
        )

        self.temp_factor = Var(
            initialize=1,
            bounds=(0.9, 100),
            doc="temperature correction factor",
            units=pyunits.dimensionless,
        )

        self.equipment_cost = Var(
            initialize=self.ref_cost,
            bounds=(0, 1e4),
            doc="equipment cost of sCO2 unit in $MM",
            units=CE_index_units,
        )

        self.bare_erected_cost = Var(
            initialize=self.ref_cost,
            bounds=(0, 1e4),
            doc="bare erected cost of sCO2 unit" "in $MM",
            units=CE_index_units,
        )

        self.total_plant_cost = Var(
            initialize=self.ref_cost,
            bounds=(0, 1e4),
            doc="total plant cost of sCO2 unit" "in $MM",
            units=CE_index_units,
        )

        # divides the scaled parameter by the number of pieces of equipment
        def scaled_param_rule(costing):
            return costing.scaled_param * costing.n_equip == scaled_param

        self.scaled_param_eq = Constraint(rule=scaled_param_rule)

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
                self.temperature = Var(
                    initialize=500,
                    bounds=(0, 1e6),
                    doc="dummy var for temperature",
                    units=pyunits.C,
                )

                self.temp_eq = Constraint(expr=(self.temperature == temp_C))

                def temp_correction_rule(costing):  # rule for temp correction
                    return (
                        Expr_if(
                            costing.temperature < 550 * pyunits.C,
                            1e-6 * costing.temperature / pyunits.C + 1,
                            1
                            + costing.c * (costing.temperature - 550 * pyunits.C)
                            + costing.d * (costing.temperature - 550 * pyunits.C) ** 2,
                        )
                        == costing.temp_factor
                    )

                self.temp_correction_eq = Constraint(rule=temp_correction_rule)
        else:
            self.temp_factor.fix(1)

        # rule for equipment cost
        def equipment_cost_rule(costing):
            return (
                costing.equipment_cost
                == costing.n_equip
                * pyunits.convert(costing.ref_cost, CE_index_units)
                * (
                    (costing.scaled_param / costing.scaled_param.get_units())
                    ** costing.exp
                )
                * costing.temp_factor
            )

        self.equipment_cost_eq = Constraint(rule=equipment_cost_rule)

        # rule for bare erected cost
        def bare_erected_cost_rule(costing):
            return costing.bare_erected_cost == costing.equipment_cost * (
                1 + costing.material_cost + costing.labor_cost
            )

        self.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

        # rule for calculating total plant cost
        def total_plant_cost_rule(costing):
            return costing.total_plant_cost == costing.bare_erected_cost * (
                1 + costing.eng_fee + costing.process_conting + costing.project_conting
            )

        self.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # add variable and constraint scaling
        if equipment in ["Recuperator", "Direct air cooler"]:
            iscale.set_scaling_factor(self.scaled_param, 1e-5)
        else:
            iscale.set_scaling_factor(self.scaled_param, 1)

        iscale.set_scaling_factor(self.equipment_cost, 1e3)
        iscale.set_scaling_factor(self.bare_erected_cost, 1e3)
        iscale.set_scaling_factor(self.total_plant_cost, 1e3)
        iscale.constraint_scaling_transform(
            self.equipment_cost_eq, 1e-6, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.bare_erected_cost_eq, 1e3, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.bare_erected_cost_eq, 1e3, overwrite=False
        )

    # -----------------------------------------------------------------------------
    # Air Separation Unit Costing Library
    # -----------------------------------------------------------------------------
    def get_ASU_cost(self, scaled_param, CE_index_year="2018"):
        # scaled parameter is O2 flowrate in TPD
        # only one set of parameters used, ref is hard coded for TPD and 2017

        params = {
            "Reference Cost": 3.26e6,
            "Reference Parameter": 13078,
            "Exponent": 0.7,
            "Eng Fee": 0.097,
            "Process": 0,
            "Project": 0.110,
        }

        # check to see if a costing block already exists
        if (
            self.parent_block().name
            in self.config.flowsheet_costing_block._registered_unit_costing  # pylint: disable=protected-access
        ):
            raise AttributeError(
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(self.name)
            )

        # define costing library
        self.library = "ASU"

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # define parameters
        self.ref_cost = Param(
            initialize=params["Reference Cost"],
            mutable=True,
            doc="ASU reference cost",
            units=pyunits.USD_2008_Nov,
        )

        self.ref_param = Param(
            initialize=params["Reference Parameter"],
            mutable=True,
            doc="ASU reference parameter value",
            units=pyunits.ton / pyunits.d,
        )

        self.exp = Param(
            initialize=params["Exponent"],
            mutable=True,
            doc="ASU scaling exponent",
            units=pyunits.dimensionless,
        )

        self.eng_fee = Param(
            mutable=True,
            initialize=params["Eng Fee"],
            doc="engineering fee percentage",
            units=pyunits.dimensionless,
        )

        self.process_conting = Param(
            mutable=True,
            initialize=params["Process"],
            doc="process contingency percentage",
            units=pyunits.dimensionless,
        )

        self.project_conting = Param(
            mutable=True,
            initialize=params["Project"],
            doc="project contingency percentage",
            units=pyunits.dimensionless,
        )

        # define variables
        self.bare_erected_cost = Var(
            initialize=params["Reference Cost"],
            bounds=(0, 1e4),
            doc="scaled bare erected cost in $MM",
            units=CE_index_units,
        )

        self.total_plant_cost = Var(
            initialize=params["Reference Cost"],
            bounds=(0, 1e4),
            doc="total plant cost in $MM",
            units=CE_index_units,
        )

        # rule for scaling BEC

        def bare_erected_cost_rule(costing):
            return (
                costing.bare_erected_cost
                == pyunits.convert(
                    costing.ref_cost,
                    CE_index_units,
                )
                * (scaled_param / costing.ref_param) ** costing.exp
            )

        self.bare_erected_cost_eq = Constraint(rule=bare_erected_cost_rule)

        # rule for calculating TPC
        def total_plant_cost_rule(costing):
            return costing.total_plant_cost == costing.bare_erected_cost * (
                1 + costing.eng_fee + costing.process_conting + costing.project_conting
            )

        self.total_plant_cost_eq = Constraint(rule=total_plant_cost_rule)

        # add variable and constraint scaling
        iscale.set_scaling_factor(self.bare_erected_cost, 1)
        iscale.set_scaling_factor(self.total_plant_cost, 1)

        iscale.constraint_scaling_transform(
            self.bare_erected_cost_eq, 1e-3, overwrite=False
        )
        iscale.constraint_scaling_transform(
            self.total_plant_cost_eq, 1, overwrite=False
        )

    # -----------------------------------------------------------------------------
    # REE Recovery Costing Library
    # -----------------------------------------------------------------------------
    def get_REE_costing(
        blk,
        cost_accounts,
        scaled_param,
        tech,
        Lang_factor=None,
        n_equip=1,
        scale_down_parallel_equip=False,
        CE_index_year="2018",
        additional_costing_params=None,
        use_additional_costing_params=False,
    ):
        """
        The scaled cost is computed using reference values for different
        technologies as listed below:
            1. University of Kentucky Fire Clay Seam (Hazard No. 4) Rejects
        Args:
            blk: A unit-level costing block where costing variables and
                constraints can be added to
            cost_accounts: A list of accounts to be included in the total cost
            scaled_param: the process parameter for the system(s) being costed;
                this is the total flow for all parallel trains of the system(s)
            tech: integer representing the above categories
            Lang_factor: optional single Lang factor value used to calculate
                commercial-scale installation costs. If None, accounts must
                include necessary component factors.
            n_equip: Integer number of parallel equipment trains for unit
                operations; for example, enter '5' if a feed will be split
                among 5 identical units and then re-mixed
            scale_down_parallel_equip: Boolean flag whether to scale down
                parallel equipment trains, e.g. two trains scaled down will
                each be half the size/capacity of a single train, and two
                trains not scaled down will each be the same size as a single
                train (twice the capacity).
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars
            additional_costing_params: user-defined dictionary to append to
                existing cost accounts dictionary
            use_additional_costing_params: Boolean flag to use additional
                costing parameters when account names conflict with existing
                accounts data


        Cost is in M$
        """
        # check to see if a costing block already exists
        if (
            blk.parent_block().name
            in blk.config.flowsheet_costing_block._registered_unit_costing  # pylint: disable=protected-access
        ):
            raise AttributeError(
                "{} already has an attribute costing. "
                "Check that you are not calling get_costing"
                " twice on the same model".format(blk.name)
            )

        # define costing library
        blk.library = "REE"

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # pull data for each account into dictionaries
        process_params = {}
        reference_units = {}
        account_names = {}
        exponents = {}
        reference_costs = {}
        reference_cost_units = {}
        reference_costs_init = {}
        reference_params = {}
        piping_materials_and_labor_percentage = {}
        electrical_materials_and_labor_percentage = {}
        instrumentation_percentage = {}
        plants_services_percentage = {}
        process_buildings_percentage = {}
        auxiliary_buildings_percentage = {}
        site_improvements_percentage = {}
        equipment_installation_percentage = {}
        field_expenses_percentage = {}
        project_management_and_construction_percentage = {}
        process_contingency_percentage = {}

        # load the cost dictionaries and add custom accounts

        # load ree costing dictionary
        REE_costing_params = load_REE_costing_dictionary()

        # for compatibility with potential custom accounts, the loop handles
        # new technologies, and new accounts for existing technologies
        # Users should not be adding new entries for existing accounts

        costing_params = REE_costing_params  # initialize with baseline accounts
        if additional_costing_params is not None and additional_costing_params != {}:
            for (
                new_costing_params
            ) in additional_costing_params:  # merge new dictionaries sequentially
                # adding any provided custom params to the base dictionary
                # need to "freeze" dict so it is hashable for merging keys
                frozen_dict = {**costing_params}
                for techkey, techval in new_costing_params.items():
                    if (
                        techkey in frozen_dict.keys()
                    ):  # if techkey already exists, append any new accounts
                        for accountkey, accountval in new_costing_params[
                            techkey
                        ].items():
                            if (
                                accountkey in frozen_dict[techkey].keys()
                            ) and not use_additional_costing_params:
                                if accountkey not in cost_accounts:
                                    pass  # not the current account, don't fail here
                                else:  # this is not allowed
                                    raise ValueError(
                                        "Data already exists for Account {} "
                                        "using technology {} with CCS {}. "
                                        "Please confirm that the custom "
                                        "account dictionary is correct, or "
                                        "add the new parameters as a new "
                                        "account. To use the custom account "
                                        "dictionary for all conflicts, please "
                                        "pass the argument use_additional_costing_params "
                                        "as True.".format(accountkey, str(techkey))
                                    )
                            else:  # conflict is the account passed, and overwrite it
                                frozen_dict[techkey][accountkey] = accountval
                    else:
                        frozen_dict[techkey] = techval
                costing_params = {k: frozen_dict[k] for k in sorted(frozen_dict)}

        for account in cost_accounts:
            try:  # look for data in json file info
                process_params[account] = costing_params[str(tech)][account][
                    "Process Parameter"
                ]
                reference_units[account] = costing_params[str(tech)][cost_accounts[0]][
                    "Units"
                ]
                account_names[account] = costing_params[str(tech)][account][
                    "Account Name"
                ]
                exponents[account] = float(
                    costing_params[str(tech)][account]["Exponent"]
                )
                reference_costs[account] = costing_params[str(tech)][account]["BEC"]
                reference_cost_units[account] = costing_params[str(tech)][account][
                    "BEC_units"
                ]
                reference_costs_init[account] = (
                    costing_params[str(tech)][account]["BEC"] * 1e-3
                )

                if isinstance(process_params[account], list):
                    for i, processparam in enumerate(process_params[account]):
                        reference_params[account, processparam] = costing_params[
                            str(tech)
                        ][account]["RP Value"][i]

                elif isinstance(process_params[account], str):
                    reference_params[account] = costing_params[str(tech)][account][
                        "RP Value"
                    ]

                piping_materials_and_labor_percentage[account] = costing_params[
                    str(tech)
                ][account]["Piping, Materials and Labor"]
                electrical_materials_and_labor_percentage[account] = costing_params[
                    str(tech)
                ][account]["Electrical, Materials and Labor"]
                instrumentation_percentage[account] = costing_params[str(tech)][
                    account
                ]["Instrumentation"]
                plants_services_percentage[account] = costing_params[str(tech)][
                    account
                ]["Plant Services"]
                process_buildings_percentage[account] = costing_params[str(tech)][
                    account
                ]["Process Buildings"]
                auxiliary_buildings_percentage[account] = costing_params[str(tech)][
                    account
                ]["Auxiliary Buildings"]
                site_improvements_percentage[account] = costing_params[str(tech)][
                    account
                ]["Site Improvements"]
                equipment_installation_percentage[account] = costing_params[str(tech)][
                    account
                ]["Equipment Installation"]
                field_expenses_percentage[account] = costing_params[str(tech)][account][
                    "Field Expenses"
                ]
                project_management_and_construction_percentage[
                    account
                ] = costing_params[str(tech)][account][
                    "Project Management and Construction"
                ]
                process_contingency_percentage[account] = costing_params[str(tech)][
                    account
                ]["Process Contingency"]
            except KeyError:
                print(
                    "KeyError: Account {} could not be found in the "
                    "dictionary for technology {} with CCS {}".format(
                        account, str(tech)
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

        # check that the user passed the correct units type and try to convert

        for account in cost_accounts:
            ref_units = reference_units[account]
            if "/" in ref_units:
                ref_units = ref_units.split("/")
                if "**" in ref_units[0]:
                    ref_units[0] = ref_units[0].split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0][0]) ** int(
                            ref_units[0][1]
                        ) / getattr(pyunits, ref_units[1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (
                                cost_accounts[0],
                                ref_units[0][0]
                                + "**"
                                + ref_units[0][1]
                                + "/"
                                + ref_units[1],
                            )
                        )
                elif "**" in ref_units[1]:
                    ref_units[1] = ref_units[1].split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1][0]
                        ) ** int(ref_units[1][1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (
                                cost_accounts[0],
                                ref_units[0]
                                + "/"
                                + ref_units[1][0]
                                + "**"
                                + ref_units[1][1],
                            )
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) / getattr(
                            pyunits, ref_units[1]
                        )
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )

            else:
                if "**" in ref_units:
                    ref_units = ref_units.split("**")
                    try:
                        ref_units = getattr(pyunits, ref_units[0]) ** int(ref_units[1])
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )
                else:
                    try:
                        ref_units = getattr(pyunits, ref_units)
                    except AttributeError:
                        raise AttributeError(
                            "Account %s uses references units of %s. Cannot "
                            "parse reference units as Pyomo unit containers. "
                            "Check that source uses correct syntax for Pyomo "
                            "unit containers, for example gpm should be "
                            "gal/min, tpd should be ton/d and MMBtu should be "
                            "MBtu (using Pyomo prefix)."
                            % (cost_accounts[0], ref_units[0] + "/" + ref_units[1])
                        )

            if isinstance(scaled_param, list):
                for sp in scaled_param:
                    if sp.get_units() is None:
                        raise ValueError(
                            "Account %s uses units of %s. "
                            "Units of %s were passed. "
                            "Scaled_param must have units."
                            % (cost_accounts[0], ref_units, sp.get_units())
                        )
                    else:
                        try:
                            pyunits.convert(sp, ref_units)
                        except InconsistentUnitsError:
                            raise Exception(
                                "Account %s uses units of %s. "
                                "Units of %s were passed. "
                                "Cannot convert unit containers."
                                % (
                                    cost_accounts[0],
                                    ref_units,
                                    sp.get_units(),
                                )
                            )
            else:
                try:
                    if pyunits.get_units(scaled_param) is None:
                        raise UnitsError(
                            "Account %s uses units of %s. "
                            "Units of %s were passed. "
                            "Scaled_param must have units."
                            % (
                                cost_accounts[0],
                                ref_units,
                                pyunits.get_units(scaled_param),
                            )
                        )
                    else:
                        try:
                            pyunits.convert(scaled_param, ref_units)
                        except InconsistentUnitsError:
                            raise UnitsError(
                                "Account %s uses units of %s. "
                                "Units of %s were passed. "
                                "Cannot convert unit containers."
                                % (
                                    cost_accounts[0],
                                    ref_units,
                                    pyunits.get_units(scaled_param),
                                )
                            )
                except InconsistentUnitsError:
                    raise UnitsError(
                        f"The expression {scaled_param.name} has inconsistent units."
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
            # units not defined here, since every account could have different
            # currency units
        )

        if isinstance(process_params[cost_accounts[0]], list):
            if len(process_params[cost_accounts[0]]) > 1:
                blk.ref_param = Param(
                    cost_accounts,
                    process_params[cost_accounts[0]],
                    mutable=True,
                    initialize=reference_params,
                    doc="reference parameter for account",
                )
        elif isinstance(process_params[cost_accounts[0]], str):
            blk.ref_param = Param(
                cost_accounts,
                mutable=True,
                initialize=reference_params,
                doc="reference parameter for account",
            )

        blk.piping_materials_and_labor_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=piping_materials_and_labor_percentage,
            doc="piping, materials and labor",
        )

        blk.electrical_materials_and_labor_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=electrical_materials_and_labor_percentage,
            doc="electrical, materials and labor",
        )

        blk.instrumentation_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=instrumentation_percentage,
            doc="instrumentation",
        )

        blk.plant_services_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=plants_services_percentage,
            doc="plant services",
        )

        blk.process_buildings_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=process_buildings_percentage,
            doc="process buildings",
        )

        blk.auxiliary_buildings_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=auxiliary_buildings_percentage,
            doc="auxiliary buildings",
        )

        blk.site_improvements_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=site_improvements_percentage,
            doc="site improvements",
        )

        blk.equipment_installation_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=equipment_installation_percentage,
            doc="equipment installation",
        )

        blk.field_expenses_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=field_expenses_percentage,
            doc="field expenses",
        )

        blk.project_management_and_construction_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=project_management_and_construction_percentage,
            doc="project management and construction",
        )

        blk.process_contingency_percentage = Param(
            cost_accounts,
            mutable=True,
            initialize=process_contingency_percentage,
            doc="process contingency",
        )

        # define variables
        blk.bare_erected_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, 1e4),
            doc="scaled bare erected cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        if Lang_factor is None:
            # ancillary costs
            blk.ancillary_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.piping_materials_and_labor_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="piping, materials and labor ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.electrical_materials_and_labor_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="electrical, materials and labor ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.instrumentation_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.plant_services_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="ancillary cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # buildings costs
            blk.buildings_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.process_buildings_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="process buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.auxiliary_buildings_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="auxiliary buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.site_improvements_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="site improvements buildings cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # engineering, procurement and construction management costs
            blk.epcm_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.equipment_installation_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="equipment installation epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.field_expenses_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="field expenses epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.project_management_and_construction_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="project management and construction epcm cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            # contingency costs - generic support more contingency cost types in the future
            blk.contingency_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="contingency cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )

            blk.process_contingency_costs = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, 1e4),
                doc="contingency cost in $MM",
                units=getattr(pyunits, "MUSD_" + CE_index_year),
            )
        else:
            blk.Lang_factor = Param(
                initialize=Lang_factor,
                mutable=True,
                doc="Lang factor",
                units=pyunits.dimensionless,
            )

        # total cost variables
        blk.total_installation_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, 1e4),
            doc="total installation cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        blk.total_plant_cost = Var(
            cost_accounts,
            initialize=reference_costs_init,
            bounds=(0, 1e4),
            doc="total plant cost in $MM",
            units=getattr(pyunits, "MUSD_" + CE_index_year),
        )

        # rule for scaling BEC
        # reference cost is in 2018 dollars, 671.1 is CE index for 2018
        def bare_erected_cost_rule(costing, i):
            ref_units = reference_units[i]
            if "/" in ref_units:
                ref_units = ref_units.split("/")
                ref_units = getattr(pyunits, ref_units[0]) / getattr(
                    pyunits, ref_units[1]
                )
            else:
                ref_units = getattr(pyunits, ref_units)

            ref_cost_units = reference_cost_units[i]
            ref_cost_units = ref_cost_units.split("$")
            if ref_cost_units[0] == "":  # no prefix
                ref_cost_units = getattr(pyunits, "USD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "K":  # thousands of $
                ref_cost_units = getattr(pyunits, "kUSD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "M":  # millions of $
                ref_cost_units = getattr(pyunits, "MUSD_" + ref_cost_units[1])

            # determine reference parameter scaler based on train scaling
            if scale_down_parallel_equip:
                scaler = n_equip
            else:
                scaler = 1

            if isinstance(process_params[i], list):
                if len(process_params[i]) > 1:
                    return costing.bare_erected_cost[i] == (
                        n_equip
                        * pyunits.convert(
                            costing.ref_cost[i] * ref_cost_units, CE_index_units
                        )
                        * sum(
                            (
                                pyunits.convert(scaled_param[j], ref_units)
                                / (scaler * costing.ref_param[i, p] * ref_units)
                            )
                            ** costing.exp[i]
                            for j, p in enumerate(process_params[i])
                        )
                    )
            elif isinstance(process_params[i], str):
                return costing.bare_erected_cost[i] == (
                    n_equip
                    * pyunits.convert(
                        costing.ref_cost[i] * ref_cost_units, CE_index_units
                    )
                    * (
                        pyunits.convert(scaled_param, ref_units)
                        / (scaler * costing.ref_param[i] * ref_units)
                    )
                    ** costing.exp[i]
                )

        blk.bare_erected_cost_eq = Constraint(
            cost_accounts, rule=bare_erected_cost_rule
        )

        if Lang_factor is None:
            # rules for calculating Ancillary costs
            def piping_materials_and_labor_cost_rule(blk, i):
                return blk.piping_materials_and_labor_costs[i] == (
                    blk.bare_erected_cost[i]
                    * blk.piping_materials_and_labor_percentage[i]
                )

            blk.piping_materials_and_labor_cost_eq = Constraint(
                cost_accounts, rule=piping_materials_and_labor_cost_rule
            )

            def electrical_materials_and_labor_cost_rule(blk, i):
                return blk.electrical_materials_and_labor_costs[i] == (
                    blk.bare_erected_cost[i]
                    * blk.electrical_materials_and_labor_percentage[i]
                )

            blk.electrical_materials_and_labor_cost_eq = Constraint(
                cost_accounts, rule=electrical_materials_and_labor_cost_rule
            )

            def instrumentation_cost_rule(blk, i):
                return blk.instrumentation_costs[i] == (
                    blk.bare_erected_cost[i] * blk.instrumentation_percentage[i]
                )

            blk.instrumentation_cost_eq = Constraint(
                cost_accounts, rule=instrumentation_cost_rule
            )

            def plant_services_cost_rule(blk, i):
                return blk.plant_services_costs[i] == (
                    blk.bare_erected_cost[i] * blk.plant_services_percentage[i]
                )

            blk.plant_services_cost_eq = Constraint(
                cost_accounts, rule=plant_services_cost_rule
            )

            def ancillary_cost_rule(blk, i):
                return blk.ancillary_costs[i] == (
                    blk.piping_materials_and_labor_costs[i]
                    + blk.electrical_materials_and_labor_costs[i]
                    + blk.instrumentation_costs[i]
                    + blk.plant_services_costs[i]
                )

            blk.ancillary_cost_eq = Constraint(cost_accounts, rule=ancillary_cost_rule)

            # rules for calculating Buildings costs

            def process_buildings_cost_rule(blk, i):
                return blk.process_buildings_costs[i] == (
                    blk.bare_erected_cost[i] * blk.process_buildings_percentage[i]
                )

            blk.process_buildings_cost_eq = Constraint(
                cost_accounts, rule=process_buildings_cost_rule
            )

            def auxiliary_buildings_cost_rule(blk, i):
                return blk.auxiliary_buildings_costs[i] == (
                    blk.bare_erected_cost[i] * blk.auxiliary_buildings_percentage[i]
                )

            blk.auxiliary_buildings_cost_eq = Constraint(
                cost_accounts, rule=auxiliary_buildings_cost_rule
            )

            def site_improvements_cost_rule(blk, i):
                return blk.site_improvements_costs[i] == (
                    blk.bare_erected_cost[i] * blk.site_improvements_percentage[i]
                )

            blk.site_improvements_cost_eq = Constraint(
                cost_accounts, rule=site_improvements_cost_rule
            )

            def buildings_cost_rule(blk, i):
                return blk.buildings_costs[i] == (
                    blk.process_buildings_costs[i]
                    + blk.auxiliary_buildings_costs[i]
                    + blk.site_improvements_costs[i]
                )

            blk.buildings_cost_eq = Constraint(cost_accounts, rule=buildings_cost_rule)

            # rules for calculating Engineering, Procurement and Construction Management costs

            def equipment_installation_cost_rule(blk, i):
                return blk.equipment_installation_costs[i] == (
                    blk.bare_erected_cost[i] * blk.equipment_installation_percentage[i]
                )

            blk.equipment_installation_cost_eq = Constraint(
                cost_accounts, rule=equipment_installation_cost_rule
            )

            def field_expenses_cost_rule(blk, i):
                return blk.field_expenses_costs[i] == (
                    blk.bare_erected_cost[i] * blk.field_expenses_percentage[i]
                )

            blk.field_expenses_cost_eq = Constraint(
                cost_accounts, rule=field_expenses_cost_rule
            )

            def project_management_and_construction_cost_rule(blk, i):
                return blk.project_management_and_construction_costs[i] == (
                    blk.bare_erected_cost[i]
                    * blk.project_management_and_construction_percentage[i]
                )

            blk.project_management_and_construction_cost_eq = Constraint(
                cost_accounts, rule=project_management_and_construction_cost_rule
            )

            def epcm_cost_rule(blk, i):
                return blk.epcm_costs[i] == (
                    blk.equipment_installation_costs[i]
                    + blk.field_expenses_costs[i]
                    + blk.project_management_and_construction_costs[i]
                )

            blk.epcm_cost_eq = Constraint(cost_accounts, rule=epcm_cost_rule)

            # rules for calculating Contingency costs
            def process_contingency_cost_rule(blk, i):
                return blk.contingency_costs[i] == (
                    blk.bare_erected_cost[i] * blk.process_contingency_percentage[i]
                )

            blk.process_contingency_cost_eq = Constraint(
                cost_accounts, rule=process_contingency_cost_rule
            )

            def contingency_cost_rule(blk, i):
                return blk.contingency_costs[i] == (blk.process_contingency_costs[i])

            blk.contingency_cost_eq = Constraint(
                cost_accounts, rule=contingency_cost_rule
            )

            def total_installation_cost_rule(blk, i):
                return blk.total_installation_cost[i] == (
                    blk.ancillary_costs[i]
                    + blk.buildings_costs[i]
                    + blk.epcm_costs[i]
                    + blk.contingency_costs[i]
                )

            blk.total_installation_cost_eq = Constraint(
                cost_accounts, rule=total_installation_cost_rule
            )
        else:

            def total_installation_cost_rule(blk, i):
                return blk.total_installation_cost[i] == blk.bare_erected_cost[i] * (
                    blk.Lang_factor - 1
                )

            blk.total_installation_cost_eq = Constraint(
                cost_accounts, rule=total_installation_cost_rule
            )

        # rule for calculating TPC
        def total_plant_cost_rule(blk, i):
            return blk.total_plant_cost[i] == (
                blk.bare_erected_cost[i] + blk.total_installation_cost[i]
            )

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

    # -----------------------------------------------------------------------------
    # Operation & Maintenance Costing Library
    # -----------------------------------------------------------------------------

    def get_fixed_OM_costs(
        b,
        net_power=None,
        nameplate_capacity=650,
        capacity_factor=0.85,
        labor_rate=38.50,
        labor_burden=30,
        operators_per_shift=6,
        tech=1,
        fixed_TPC=None,
        CE_index_year="2018",
    ):
        """
        Creates constraints for the following fixed O&M costs in $MM/yr:
            1. Annual operating labor
            2. Maintenance labor
            3. Admin and support labor
            4. Property taxes and insurance
            5. Other fixed costs
            6. Total fixed O&M cost
            7. Maintenance materials (actually a variable cost,
                but scales off TPC)

        These costs apply to the project as a whole and are scaled based on the
        total TPC.

        Args:
            b: costing block to add fixed cost variables and constraints to
            net_power: actual plant output in MW, only required if calculating
                variable costs
            nameplate_capacity: rated plant output in MW
            capacity_factor: multiplicative factor for normal operating
                capacity
            labor_rate: hourly rate of plant operators in project dollar year
            labor_burden: a percentage multiplier used to estimate non-salary
                labor expenses
            operators_per_shift: average number of operators per shift
            tech: int 1-7 representing the categories in get_PP_costing, used
                to determine maintenance costs
            fixed_TPC: The TPC in $MM that will be used to determine fixed O&M,
                costs. If the value is None, the function will try to use the
                TPC calculated from the individual units.
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars

        Returns:
            None
        """

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # create and fix total_TPC if it does not exist yet
        if not hasattr(b, "total_TPC"):
            b.total_TPC = Var(
                initialize=100,
                bounds=(0, 1e4),
                doc="total TPC in $MM",
                units=CE_index_units,
            )
            if fixed_TPC is None:
                b.total_TPC.fix(100 * CE_index_units)
                _log.warning(
                    "b.costing.total_TPC does not exist and a value "
                    "for fixed_TPC was not specified, total_TPC will be "
                    "fixed to 100 MM$"
                )
            else:
                b.total_TPC.fix(value(fixed_TPC) * CE_index_units)
        else:
            if fixed_TPC is not None:
                _log.warning(
                    "b.costing.total_TPC already exists, the value "
                    "passed for fixed_TPC will be ignored."
                )

        # make params
        b.labor_rate = Param(
            initialize=labor_rate,
            mutable=True,
            units=getattr(pyunits, "USD_" + CE_index_year) / pyunits.hr,
        )
        b.labor_burden = Param(
            initialize=labor_burden, mutable=True, units=pyunits.dimensionless
        )
        b.operators_per_shift = Param(
            initialize=operators_per_shift, mutable=True, units=pyunits.dimensionless
        )
        b.nameplate_capacity = Param(
            initialize=nameplate_capacity, mutable=True, units=pyunits.MW
        )

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
            initialize=1,
            bounds=(0, 1e4),
            doc="annual labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.maintenance_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="maintenance labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.admin_and_support_labor_cost = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="admin and support labor cost in $MM/yr",
            units=CE_index_units,
        )
        b.property_taxes_and_insurance = Var(
            initialize=1,
            bounds=(0, 1e4),
            doc="property taxes and insurance cost in $MM/yr",
            units=CE_index_units,
        )
        b.total_fixed_OM_cost = Var(
            initialize=4,
            bounds=(0, 1e4),
            doc="total fixed O&M costs in $MM/yr",
            units=CE_index_units,
        )

        # variable for user to assign other fixed costs to,
        # fixed to 0 by default
        b.other_fixed_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="other fixed costs in $MM/yr",
            units=CE_index_units,
        )
        b.other_fixed_costs.fix(0)

        # maintenance material cost is technically a variable cost, but it
        # makes more sense to include with the fixed costs becuase it uses TPC
        b.maintenance_material_cost = Var(
            initialize=2e-7,
            bounds=(0, 1e4),
            doc="cost of maintenance materials in $MM/year",
            units=CE_index_units / pyunits.year,
        )

        # create constraints
        TPC = b.total_TPC  # quick reference to total_TPC

        # calculated from labor rate, labor burden, and operators per shift
        @b.Constraint()
        def annual_labor_cost_rule(c):
            return c.annual_operating_labor_cost == pyunits.convert(
                (
                    c.operators_per_shift
                    * c.labor_rate
                    * (1 + c.labor_burden / 100)
                    * 8760
                    * pyunits.hr
                ),
                CE_index_units,
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

        # technology specific percentage of TPC
        @b.Constraint()
        def maintenance_material_cost_rule(c):
            if net_power is not None:
                return c.maintenance_material_cost == (
                    TPC
                    * c.maintenance_material_TPC_split
                    * c.maintenance_material_percent
                    / (capacity_factor * pyunits.year)
                    * net_power[0]
                    / c.nameplate_capacity
                )
            else:
                return c.maintenance_material_cost == (
                    TPC
                    * c.maintenance_material_TPC_split
                    * c.maintenance_material_percent
                    / (capacity_factor * pyunits.year)
                    * (capacity_factor)
                )

    def get_variable_OM_costs(
        b,
        resources,
        rates,
        prices=None,
        CE_index_year="2018",
        capacity_factor=0.85,
    ):
        """
        This function is used to calculate the variable cost of producing
        electricity in $/MWh. The function is structured to allow the user to
        generate all fuel, consumable, and waste disposal costs at once.
        A total variable cost is created for each point in fs.time.

        Args:
            b: costing block to add fixed cost variables and constraints to
            resources: list of strings for the resources to be costed
            rates: list of pyomo vars for resource consumption rates
            prices: dict of resource prices to be added to the premade
                dictionary
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars
            capacity_factor: multiplicative factor for normal operating
                capacity

        Returns:
            None.

        """
        if prices is None:
            prices = {}

        if not hasattr(b.parent_block(), "time"):  # flowsheet is not dynamic
            b.parent_block().time = [0]
        if not hasattr(b.parent_block(), "time_units"):  # no time units set
            b.parent_block().time_units = pyunits.s

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # assert arguments are correct types
        if not isinstance(resources, list):
            raise TypeError("resources argument must be a list")
        if not isinstance(rates, list):
            raise TypeError("rates argument must be a list")
        if not isinstance(prices, dict):
            raise TypeError("prices argument must be a dictionary")

        # assert lists are the same length
        if len(resources) != len(rates):
            raise Exception("resources and rates must be lists of the same" " length")

        # dictionary of default prices
        # the currency units are millions of USD, so all prices need a 1e-6
        # multiplier to get USD
        default_prices = {
            "natural_gas": 4.42 * 1e-6 * CE_index_units / pyunits.MBtu,  # $/MMbtu
            "coal": 51.96 * 1e-6 * CE_index_units / pyunits.ton,
            "water": 1.90e-3 * 1e-6 * CE_index_units / pyunits.gallon,
            "water_treatment_chemicals": 550 * 1e-6 * CE_index_units / pyunits.ton,
            "ammonia": 300 * 1e-6 * CE_index_units / pyunits.ton,
            "SCR_catalyst": 150 * 1e-6 * CE_index_units / pyunits.ft**3,
            "triethylene_glycol": 6.80 * 1e-6 * CE_index_units / pyunits.gallon,
            "SCR_catalyst_waste": 2.50 * 1e-6 * CE_index_units / pyunits.ft**3,
            "triethylene_glycol_waste": 0.35 * 1e-6 * CE_index_units / pyunits.gallon,
            "amine_purification_unit waste": 38 * 1e-6 * CE_index_units / pyunits.ton,
            "thermal_reclaimer_unit_waste": 38 * 1e-6 * CE_index_units / pyunits.ton,
        }

        # add entries from prices to default_prices
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

        # make vars
        b.variable_operating_costs = Var(
            b.parent_block().time,
            resources,
            initialize=2e-7,
            doc="variable operating costs in $MM/year",
            units=CE_index_units / pyunits.year,
        )

        b.other_variable_costs = Var(
            b.parent_block().time,
            initialize=0,
            doc="a variable to include non-standard O&M costs in $MM/year",
            units=CE_index_units / pyunits.year,
        )
        # assume the user is not using this
        b.other_variable_costs.fix(0)

        b.total_variable_OM_cost = Var(
            b.parent_block().time,
            initialize=4e-6,
            doc="total variable operating and maintenance costs in $MM/year",
            units=CE_index_units / pyunits.year,
        )

        @b.Constraint(b.parent_block().time, resources)
        def variable_cost_rule_power(c, t, r):
            return c.variable_operating_costs[t, r] == (
                pyunits.convert(
                    resource_prices[r] * resource_rates[r][t],
                    CE_index_units / pyunits.year,
                )
            )

        if hasattr(b, "maintenance_material_cost"):  # from get_fixed_OM_costs

            @b.Constraint(b.parent_block().time)
            def total_variable_cost_rule_power(c, t):
                return (
                    c.total_variable_OM_cost[t]
                    == sum(c.variable_operating_costs[t, r] for r in resources)
                    + c.maintenance_material_cost
                    + c.other_variable_costs[t]
                )

        else:  # not built yet, don't include that variable

            @b.Constraint(b.parent_block().time)
            def total_variable_cost_rule_power(c, t):
                return (
                    c.total_variable_OM_cost[t]
                    == sum(c.variable_operating_costs[t, r] for r in resources)
                    + c.other_variable_costs[t]
                )

    def initialize_fixed_OM_costs(b):
        # b is the flowsheet-level costing block
        if hasattr(b, "total_fixed_OM_cost"):
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

    def initialize_variable_OM_costs(b):
        # b is the flowsheet-level costing block
        # initialization for power generation costs
        if hasattr(b, "variable_operating_costs"):
            for i in b.variable_operating_costs.keys():
                if hasattr(b, "variable_cost_rule_power"):
                    calculate_variable_from_constraint(
                        b.variable_operating_costs[i],
                        b.variable_cost_rule_power[i],
                    )

            for i in b.total_variable_OM_cost.keys():
                calculate_variable_from_constraint(
                    b.total_variable_OM_cost[i],
                    b.total_variable_cost_rule_power[i],
                )

    def get_REE_plant_costs(
        b,
        blocks_to_cost=None,
        # add more arguments for detailed O&M calculations later
        CE_index_year="2018",
    ):
        """
        Creates constraints for the following plant-level costs in $MM/yr:
            1. Total ancillary
            2. Piping materials and labor ancillary
            2. Electrical materials and labor ancillary
            3. Instrumentation ancillary
            4. Plant services ancillary
            5. Total buildings
            6. Process buildings
            7. Auxiliary buildings
            8. Site improvements buildings
            9. Total engineering procurement and construction management (EPCM)
            10. Equipment installation EPCM
            11. Field expenses EPCM
            12. Project management and construction EPCM
            13. Total contingency
            14. Process contingency
            (space for more costs, including contingency costs, in the future)

        These costs apply to the project as a whole and are scaled based on the
        total TPC.

        Args:
            b: flowsheet-level costing block to add installation cost to
            blocks_to_cost: let the user specify which blocks should contribute
                to the plant-wide total plant cost and total installation cost.
                Each block in the list must be a UnitModelCostingBlock(). If
                this argument is not set, the method assumes that "b" is a
                flowsheet-level costing block and defines the list based on the
                registered costing blocks.
            CE_index_year: year for cost basis, e.g. "2018" to use 2018 dollars

        Returns:
            None
        """

        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        # define costing library
        b.library = "REE"

        # make vars
        b.total_plant_cost = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Plant Cost [$MM]",
            units=CE_index_units,
        )
        b.bare_erected_cost = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Bare Erected Cost [$MM]",
            units=CE_index_units,
        )
        b.total_installation_cost = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.ancillary_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Ancillary Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.piping_materials_and_labor_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Ancillary Piping, Materials and Labor Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.electrical_materials_and_labor_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Ancillary Electrical, Materials and Labor Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.instrumentation_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Ancillary Instrumentation Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.plant_services_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Ancillary Plant Services Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.buildings_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Buildings Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.process_buildings_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Process Buildings Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.auxiliary_buildings_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Auxiliary Buildings Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.site_improvements_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Site Improvements Buildings Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.epcm_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total EPCM Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.equipment_installation_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Equipment Installation EPCM Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.field_expenses_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Field Expenses EPCM Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.project_management_and_construction_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Project Management and Construction EPCM Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.process_contingency_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Process Contingency Installation Cost [$MM]",
            units=CE_index_units,
        )
        b.contingency_costs = Var(
            initialize=0,
            bounds=(0, 1e4),
            doc="Total Contingency Installation Cost [$MM]",
            units=CE_index_units,
        )

        # define list of blocks to cost
        if blocks_to_cost is None:  # need to search for blocks to cost
            blocks_to_cost = []
            for o in b.parent_block().component_objects(descend_into=True):
                # look for costing blocks
                if o.name in [block.name for block in b._registered_unit_costing]:
                    blocks_to_cost.append(o)

        # create empty dictionary to concisely save and call block variables
        b.installation_cost_list = [
            "total_plant_cost",
            "bare_erected_cost",
            "total_installation_cost",
            "ancillary_costs",
            "piping_materials_and_labor_costs",
            "electrical_materials_and_labor_costs",
            "instrumentation_costs",
            "plant_services_costs",
            "buildings_costs",
            "process_buildings_costs",
            "auxiliary_buildings_costs",
            "site_improvements_costs",
            "epcm_costs",
            "equipment_installation_costs",
            "field_expenses_costs",
            "project_management_and_construction_costs",
            "process_contingency_costs",
            "contingency_costs",
        ]
        plant_cost_dict = dict.fromkeys(b.installation_cost_list, [])
        plant_cost_dict = dict((key, []) for key in b.installation_cost_list)

        # retrieve and save block variables
        for var in b.installation_cost_list:  # build one var row at a time
            for block in blocks_to_cost:  # get var from each block
                if hasattr(block, var):
                    for key in getattr(block, var).keys():  # get all var keys
                        plant_cost_dict[var].append(getattr(block, var)[key])

        # create constraints
        for var in b.installation_cost_list:
            if len(plant_cost_dict[var]) > 0:  # picked up vars, sum them
                total_var = getattr(b, var)  # aggregate variable
                sum_var = sum(plant_cost_dict[var])  # sum of components
                constraint = Constraint(
                    expr=(
                        total_var == pyunits.convert(sum_var, to_units=CE_index_units)
                    )
                )
                setattr(b, var + "_eq", constraint)  # generate the constraint
            elif len(plant_cost_dict[var]) == 0:  # no vars exists
                # could occur if Lang factor set and individual
                # installation costs do not exist - delete aggregate vars
                b.del_component(var)

    def initialize_REE_plant_costs(b):
        # b is the flowsheet-level costing block
        # initialization for REE plant-level costs
        for var in b.installation_cost_list:
            if hasattr(b, var):
                calculate_variable_from_constraint(
                    getattr(b, var), getattr(b, var + "_eq")
                )

    # -----------------------------------------------------------------------------
    # Costing Library Utility Functions
    # -----------------------------------------------------------------------------

    def costing_initialization(b):
        # b is the flowsheet-level costing block
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in b._registered_unit_costing and hasattr(o, "library"):
                if o.library == "sCO2":
                    if o.equipment in [
                        "Axial turbine",
                        "Radial turbine",
                        "Coal-fired heater",
                        "Natural gas-fired heater",
                        "Recouperator",
                    ]:
                        calculate_variable_from_constraint(o.temperature, o.temp_eq)
                        calculate_variable_from_constraint(
                            o.temp_factor, o.temp_correction_eq
                        )

                    calculate_variable_from_constraint(
                        o.scaled_param, o.scaled_param_eq
                    )
                    calculate_variable_from_constraint(
                        o.equipment_cost, o.equipment_cost_eq
                    )
                    calculate_variable_from_constraint(
                        o.bare_erected_cost, o.bare_erected_cost_eq
                    )
                    calculate_variable_from_constraint(
                        o.total_plant_cost, o.total_plant_cost_eq
                    )

                elif o.library in ["PP", "ASU", "REE"]:
                    for key in o.costing.bare_erected_cost.keys():
                        calculate_variable_from_constraint(
                            o.bare_erected_cost[key],
                            o.bare_erected_cost_eq[key],
                        )
                        calculate_variable_from_constraint(
                            o.total_plant_cost[key],
                            o.total_plant_cost_eq[key],
                        )

    def display_total_plant_costs(b):
        print("-----Total Plant Costs-----")
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "total_plant_cost"):
                print(
                    "%s: $%.2f Million"
                    % (
                        value(o.name),
                        value(
                            sum(o.total_plant_cost[key] for key in o.total_plant_cost)
                        ),
                    )
                )

    def display_bare_erected_costs(b):
        print("-----Bare Erected Costs-----")
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "bare_erected_costs"):
                print(
                    "%s: $%.2f Million"
                    % (
                        value(o.name),
                        value(
                            sum(
                                o.bare_erected_costs[key]
                                for key in o.bare_erected_costs
                            )
                        ),
                    )
                )

    def display_equipment_costs(b):
        print("-----Equipment Costs-----")
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "equipment_costs"):
                print(
                    "%s: $%.2f Million"
                    % (
                        value(o.name),
                        value(sum(o.equipment_costs[key] for key in o.equipment_costs)),
                    )
                )

    def get_total_TPC(b, CE_index_year):
        # This method accepts a flowsheet-level costing block

        try:
            CE_index_units = getattr(
                pyunits, "MUSD_" + CE_index_year
            )  # millions of USD, for base year
        except AttributeError:
            raise AttributeError(
                "CE_index_year %s is not a valid currency base option. "
                "Valid CE index options include CE500, CE394 and years from "
                "1990 to 2020." % (CE_index_year)
            )

        TPC_list = []

        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name for block in b._registered_unit_costing
            ] and hasattr(o, "total_plant_cost"):
                for key in o.total_plant_cost.keys():
                    TPC_list.append(o.total_plant_cost[key])

        b.total_TPC = Var(
            initialize=100,
            bounds=(0, 1e4),
            doc="total TPC in $MM",
            # assume that total_plant_cost is in millions of
            # USD_year, where year is the CE_index_year users set
            units=CE_index_units,
        )

        @b.Constraint()
        def total_TPC_eq(c):
            return c.total_TPC == sum(TPC_list)

    def display_flowsheet_cost(b):
        # This method accepts a flowsheet-level costing block
        print("\n")
        print("Total flowsheet cost: $%.3f Million" % value(b.total_TPC))

    def check_sCO2_costing_bounds(b):
        # load sCO2 costing dictionary
        sCO2_costing_params = load_sCO2_costing_dictionary()
        # This method accepts a flowsheet-level costing block
        for o in b.parent_block().component_objects(descend_into=False):
            # look for costing blocks
            if o.name in [block.name for block in b._registered_unit_costing]:
                if o.library == "sCO2":
                    equipment = o.equipment
                    lower_bound = sCO2_costing_params[equipment]["Lower Bound"]
                    upper_bound = sCO2_costing_params[equipment]["Upper Bound"]
                    if value(o.scaled_param) < lower_bound:
                        print(
                            """%s: The scaled parameter (%f) is below the lower
                            bound (%f)."""
                            % (value(o.name), value(o.scaled_param), lower_bound)
                        )
                    elif value(o.scaled_param) > upper_bound:
                        print(
                            """%s: The scaled parameter (%f) is above the upper
                            bound (%f)."""
                            % (value(o.name), value(o.scaled_param), upper_bound)
                        )
                    else:
                        print(
                            """%s: The scaled parameter is within the
                            bounds."""
                            % value(o.name)
                        )
