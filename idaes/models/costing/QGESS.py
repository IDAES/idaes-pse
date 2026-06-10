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
Costing package based on methods from:

    Quality Guidelines for Energy System Studies: Capital Cost Scaling Methodology
    Revision 4 Report) https://doi.org/10.2172/1573493

    Quality Guidelines for Energy Systems Studies: Cost Estimation Methodology for
    NETL Assessments of Power Plant Performance https://doi.org/10.2172/1567736

This costing package includes common methods for equipment cost scaling using
economies of scale, cost model initialization and scaling, fixed and variable
operating costs and estimated revenue over the plant lifetime, tax estimation,
net present value calculation, and economy of numbers.
"""

# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

# all authors
__author__ = "Costing Team (B. Paul, L. Deng, A. Fritz, A. Ojo, A. Dasgupta, A. Noring, A. Deshpande, D. Caballero, and M. Zamarripa)"


__version__ = "1.0.0"

import textwrap
from sys import stdout

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core import FlowsheetCostingBlockData, declare_process_block_class
from idaes.core.util.exceptions import BurntToast, ConfigurationError
from idaes.core.util.math import smooth_max
from idaes.core.util.tables import stream_table_dataframe_to_string
from idaes.models_extra.power_generation.costing.generic_ccs_capcost_custom_dict import (
    load_generic_ccs_costing_dictionary,
)
from idaes.models_extra.power_generation.costing.power_plant_costing_dictionaries import (
    define_preloaded_accounts,
    load_BB_costing_dictionary,
    load_default_resource_prices,
    load_fixed_OM_data,
    register_power_plant_currency_units,
)
from pandas import DataFrame
from pyomo.common.config import ConfigValue, ListOf
from pyomo.core.base.check import BuildCheck
from pyomo.core.base.units_container import InconsistentUnitsError, UnitsError
from pyomo.environ import Expression, Param, Var, log10
from pyomo.environ import units as pyunits
from pyomo.environ import value
from pyomo.util.calc_var_value import calculate_variable_from_constraint

_log = idaeslog.getLogger(__name__)

# EPS for calculations where needed
EPS = 1e-4


@declare_process_block_class("QGESSCosting")
class QGESSCostingData(FlowsheetCostingBlockData):
    """
    Costing class for building QGESS methods.
    """

    # Register currency and conversion rates based on CE Index

    if (
        not hasattr(pyunits, "USD_2008_Nov")
        and not hasattr(pyunits, "USD_2019_Sep")
        and not hasattr(pyunits, "USD_2018_Dec")
    ):
        register_power_plant_currency_units()

    # set CONFIG

    CONFIG = FlowsheetCostingBlockData.CONFIG()

    CONFIG.declare(
        "Lang_factor",
        ConfigValue(
            default=None,
            domain=float,
            description="Installation cost factor. If not set, uses installation factors based on selected tech.",
        ),
    )

    CONFIG.declare(
        "has_fixed_OM",
        ConfigValue(
            default=False,
            domain=bool,
            description="Boolean to enable fixed operating & maintenance cost calculations.",
        ),
    )

    CONFIG.declare(
        "has_variable_OM",
        ConfigValue(
            default=False,
            domain=bool,
            description="Boolean to enable variable operating & maintenance cost calculations.",
        ),
    )

    CONFIG.declare(
        "has_taxes_and_credits",
        ConfigValue(
            default=False,
            domain=bool,
            description="Boolean to enable tax and production credit calculations.",
        ),
    )

    CONFIG.declare(
        "has_production_credit_phaseout",
        ConfigValue(
            default=False,
            domain=bool,
            description="Boolean for whether to consider phaseout of production incentive tax credit.",
        ),
    )

    CONFIG.declare(
        "phaseout_fractions",
        ConfigValue(
            default=None,
            domain=dict,
            description="dictionary of phaseout fractiosn indexed by years during which production incentive tax credit is phased out. "
            "MUST be in chronological (ascending) order. The list MUST have the same length as "
            "phaseout_fractions. Defaults to [2031, 2032, 2033] based on the One Big Beautiful "
            "Bill Amendment (OBBBA) https://www.congress.gov/bill/119th-congress/house-bill/1/text "
            "(Volume 139 STAT. 274).",
        ),
    )

    CONFIG.declare(
        "has_net_present_value",
        ConfigValue(
            default=False,
            domain=bool,
            description="Boolean to enable net present value calculations.",
        ),
    )

    CONFIG.declare(
        "has_capital_expenditure_period",
        ConfigValue(
            default=False,
            domain=bool,
            description="Boolean for whether a capital expenditure period occurs.",
        ),
    )

    CONFIG.declare(
        "capital_expenditure_percentages",
        ConfigValue(
            default=None,
            domain=ListOf(float),
            description="A list of values that sum to 100 representing how "
            "capital costs are spread over a capital expenditure period; for "
            "example, an input of [10, 60, 30] is parsed as a 3-year period "
            "where capital costs are spread as 10% in year 1, 60% in year 2, "
            "and 30% in year 3. The capital period precedes the operating "
            "period, for example an input of [100] means that 100% of capital "
            "expenses occur in the year preceding the operating period (t=-1). Set "
            "to None to indicate no expenditure period, which means that all "
            "capital expenses occur at the start of the plant lifetime (t=0).",
        ),
    )

    CONFIG.declare(
        "has_economy_of_numbers",
        ConfigValue(
            default=False,
            domain=bool,
            description="Boolean to enable economy of numbers calculations.",
        ),
    )

    CONFIG.declare(
        "CE_index_year",
        ConfigValue(
            default="2021",
            domain=str,
            description="Basis year for costing. Must be a supported value from 1990 to 2023 or a user-defined value. For details, see the IDAES 'costing_base.py' module",
        ),
    )

    CONFIG.declare(
        "tech",
        ConfigValue(
            default=None,
            domain=int,
            description="Integer corresponding to supported technology libraries, where 1-7 are various power plant types, 8-9 are specific case studies, and 10 is UKy REE.",
        ),
    )

    # Categories:
    #     1. Supercritical PC, air-fired, with and without CO2 capture,
    #     Illinois No. 6 coal
    #     2. Subcritical PC, air-fired, with and without CO2 capture,
    #     Illinois No. 6 coal
    #     3. Two-stage, slurry-feed, oxygen-blown gasifier with and without
    #     CO2 capture, Illinois No. 6 coal
    #     4. Single-stage, slurry-feed, oxygen-blown gasifier with and without
    #     CO2 capture, Illinois No. 6 coal
    #     5. Single-stage, dry-feed, oxygen-blown, up-flow gasifier with
    #     and without CO2 capture, Illinois No. 6 coal
    #     6. Natural gas, air-fired, with and without CO2 capture
    #     7. Advanced Ultrasupercritical PC
    #     8. Polymer Layers accounts
    #     9. Sensors & Controls accounts
    #     10. University of Kentucky Fire Clay Seam (Hazard No. 4) Rejects

    # TODO add location factor to calculations
    # CONFIG.declare(
    #     "location",
    #     ConfigValue(
    #         default="United States; Washington, DC",
    #         domain=str,
    #         description="Basis location for costing. Must be a supported value passed as 'country; city';"
    #         "see the IDAES 'location_factors.json' dictionary. For entries with a country and no city, must "
    #         "be of the form 'country; None'.",
    #     ),
    # )

    def build_global_params(self):
        """
        This is where we can declare any global parameters such as default year (2021 per the CONFIG),
        currency (currently only USD is supported), Lang factor (TIC = Lang * BEC, TPC = BEC + TIC),
        location factor (Washington, D.C. = 1), and so on.
        """

        # check that the user selected a technology
        if self.config.tech is None:
            raise ValueError(
                "Must set a technology type. Valid options include: \n"
                "1. Supercritical PC, air-fired, with and without CO2 capture, Illinois No. 6 coal \n"
                "2. Subcritical PC, air-fired, with and without CO2 capture, Illinois No. 6 coal \n"
                "3. Two-stage, slurry-feed, oxygen-blown gasifier with and without CO2 capture, Illinois No. 6 coal \n"
                "4. Single-stage, slurry-feed, oxygen-blown gasifier with and without CO2 capture, Illinois No. 6 coal \n"
                "5. Single-stage, dry-feed, oxygen-blown, up-flow gasifier with and without CO2 capture, Illinois No. 6 coal \n"
                "6. Natural gas, air-fired, with and without CO2 capture \n"
                "7. Advanced Ultrasupercritical PC \n"
                "8. Polymer Layers accounts \n"
                "9. Sensors & Controls accounts \n"
                "10. University of Kentucky Fire Clay Seam (Hazard No. 4) Rejects \n"
            )

        # Set the base year for all costs

        # check that the currency units are built-in or user-defined
        try:
            self.base_currency = getattr(pyunits, "USD_" + self.config.CE_index_year)
        except AttributeError:
            raise AttributeError(
                f"CE_index_year {self.config.CE_index_year} is not a valid currency base option. "
                f"Valid CE index options include CE500, CE394, years from 1990 to 2023, or user-defined "
                f"units such as 2019_Sep and UKy_2019."
            )

        if self.config.has_production_credit_phaseout:
            # validate: must start or end with a 4-digit year
            if self.config.CE_index_year[:4].isdigit():
                current_year = int(self.config.CE_index_year[:4])
            elif self.config.CE_index_year[-4:].isdigit():
                current_year = int(self.config.CE_index_year[-4:])
            else:
                raise ValueError(
                    "CE_index_year must contain a valid 4‑digit year at the start or end "
                    "(e.g. 2019, 2019_Sep, UKy_2019)"
                )

            # define this for later
            self.current_year = Param(
                initialize=int(current_year),
                mutable=True,
                doc="Current year for calculating production incentive charge phaseout",
                units=pyunits.dimensionless,
            )

        self.CE_index_units = getattr(pyunits, "MUSD_" + self.config.CE_index_year)

        # Set a base period for all operating costs
        self.base_period = pyunits.year

        # Set capacity factor for overnight and levelized costs
        if self.config.tech == 10:  # UKy REE, assume operating 336 days per year
            self.capacity_factor = Param(
                initialize=0.92,
                mutable=True,
                doc="capacity factor of the plant",
                units=pyunits.dimensionless,
            )
        else:  # power plant accounts default to 85%
            self.capacity_factor = Param(
                initialize=0.85,
                mutable=True,
                doc="capacity factor of the plant",
                units=pyunits.dimensionless,
            )

        # if Lang_factor not None, add as param
        # this will be used to estimate installation costs of cost accounts
        # that do not already calculate construction fees and other costs, e.g.
        # many power plant accounts have their own construction cost factors
        if self.config.Lang_factor is not None:

            self.Lang_factor = Param(
                initialize=self.config.Lang_factor,
                mutable=True,
                doc="Lang factor",
                units=pyunits.dimensionless,
            )

        # else, if tech == 10, use % factors, Lang = sum of % factors
        elif self.config.tech == 10:

            installation_components = {
                "piping_materials_and_labor_percentage": 20,
                "electrical_materials_and_labor_percentage": 20,
                "instrumentation_percentage": 8,
                "plants_services_percentage": 10,
                "process_buildings_percentage": 40,
                "auxiliary_buildings_percentage": 15,
                "site_improvements_percentage": 10,
                "equipment_installation_percentage": 17,
                "field_expenses_percentage": 12,
                "project_management_and_construction_percentage": 30,
                "process_contingency_percentage": 15,
            }

            self.installation_components = Param(
                installation_components,
                mutable=True,
                initialize=installation_components,
                doc="Percentages of bare erected cost used to estimate installation costs by plant component",
                units=pyunits.percent,
            )

            self.Lang_factor = Expression(
                expr=pyunits.convert(
                    sum(
                        self.installation_components[k] for k in installation_components
                    ),
                    to_units=pyunits.dimensionless,
                )
            )

        else:  # assume there is no Lang factor and TPC = BEC

            self.Lang_factor = Param(
                initialize=1,
                mutable=True,
                doc="Lang factor",
                units=pyunits.dimensionless,
            )

        # fixed O&M params
        if self.config.has_fixed_OM:

            labor_types, labor_rates, maintenance_percentages = load_fixed_OM_data()

            if self.config.tech == 10:
                operators_per_shift = [2, 5, 2, 3, 1, 2, 0]
                labor_burden = 25
            else:
                operators_per_shift = [0, 0, 0, 0, 0, 0, 6]
                labor_burden = 30

                # power plant-specific params
                self.nameplate_capacity = Param(
                    initialize=650,
                    mutable=True,
                    units=pyunits.MW,
                    doc="baseline power generation",
                )

                self.maintenance_labor_TPC_split = Param(
                    initialize=maintenance_percentages[self.config.tech][0],
                    mutable=True,
                    doc="Percent of maintenance cost allocated for labor",
                )
                self.maintenance_labor_percent = Param(
                    initialize=maintenance_percentages[self.config.tech][1],
                    mutable=True,
                    doc="Percent of TPC used to estimate maintenance labor",
                )
                self.maintenance_material_TPC_split = Param(
                    initialize=(1 - maintenance_percentages[self.config.tech][0]),
                    mutable=True,
                    doc="Percent of maintenance cost allocated for materials",
                )
                self.maintenance_material_percent = Param(
                    initialize=maintenance_percentages[self.config.tech][2],
                    mutable=True,
                    doc="Percent of TPC used to estimate maintenance materials",
                )

            # make common labor params
            self.labor_rates = Param(
                labor_types,
                initialize=dict(zip(labor_types, labor_rates)),
                mutable=True,
                units=self.base_currency / pyunits.hr,
                doc="Pay rates for each labor type",
            )
            self.labor_burden = Param(
                initialize=labor_burden,
                mutable=True,
                units=pyunits.percent,
                doc="Fringe labor benefits percentage",
            )
            self.operators_per_shift = Param(
                labor_types,
                initialize=dict(zip(labor_types, operators_per_shift)),
                mutable=True,
                units=pyunits.dimensionless,
                doc="Operator staffing per shift",
            )
            self.mixed_product_sale_price_realization_factor = Param(
                initialize=0.65,
                mutable=True,
                units=pyunits.dimensionless,
                doc="Value estimation for saleable mixed basket products relative to pure products",
            )

            # tax params - this should only be calculated if has_fixed_OM is True so that sales revenue is calculated
            # this could be True even is has_variable_OM is False, as labor is required for production but not necessarily chemicals, electricity, etc.
            if self.config.has_taxes_and_credits:

                self.income_tax_percentage = Param(
                    initialize=26,
                    mutable=True,
                    doc="Combined federal and state income tax percentage"
                    "usually between 26 - 40%",
                    units=pyunits.percent,
                )
                self.mineral_depletion_percentage = Param(
                    initialize=14,
                    mutable=True,
                    doc="tax deduction percentage for mineral depletion."
                    "default value of 14% is used, as reported in the UKy report",
                    units=pyunits.percent,
                )
                self.production_incentive_percentage = Param(
                    initialize=10,
                    mutable=True,
                    doc="tax deduction percentage for producing critical minerals"
                    "default value of 10% of total production cost",
                    units=pyunits.percent,
                )
                self.royalty_charge_percentage_of_revenue = Param(
                    initialize=6.5,
                    mutable=True,
                    doc="Percentage of revenue charged as royalties",
                    units=pyunits.percent,
                )

                # phaseout is only relevant if tax calculations are enabled
                if self.config.has_production_credit_phaseout:

                    if not isinstance(
                        self.config.phaseout_fractions, dict
                    ):  # includes not set as default is None
                        raise ValueError(
                            "Must set phaseout_fractions as dict of fractions indexed by integer years."
                        )

                    if not list(self.config.phaseout_fractions.keys()) == sorted(
                        self.config.phaseout_fractions.keys()
                    ):
                        raise ValueError(
                            "Years for phaseout_fractions must be in ascending order."
                        )

                    if not list(self.config.phaseout_fractions.values()) == sorted(
                        self.config.phaseout_fractions.values(), reverse=True
                    ):
                        raise ValueError(
                            "Fractions for phaseout_fractions must be in descending order."
                        )

                    self.phaseout_fractions = Param(
                        self.config.phaseout_fractions.keys(),
                        initialize=self.config.phaseout_fractions,
                        mutable=True,
                        doc="dict of fractions of the production incentive tax credit that are applied during the "
                        "phaseout years. Each value MUST be between 0 and 1, and the list MUST have the same length "
                        "as phaseout_years. Defaults to [75, 50, 25] based on the One Big Beautiful Bill Amendment "
                        "(OBBBA) https://www.congress.gov/bill/119th-congress/house-bill/1/text (Volume 139 STAT. 274), "
                        "which specifies a 25% reduction in the first year of phaseout, 50% reduction in the second year, "
                        "and 75% reduction in the third year before the tax credit is fully phased out.",
                    )

                    self.phaseout_factor = Param(
                        initialize=1,
                        mutable=True,
                        doc="Phaseout fraction applied to project year",
                        units=pyunits.dimensionless,
                    )

                    if value(self.current_year) < min(
                        [int(y) for y in self.phaseout_fractions.keys()]
                    ):
                        self.phaseout_factor.set_value(1)
                    elif value(self.current_year) <= max(
                        [int(y) for y in self.phaseout_fractions.keys()]
                    ):
                        self.phaseout_factor.set_value(
                            value(
                                self.phaseout_fractions[str(value(self.current_year))]
                            )
                            / 100
                        )
                    else:
                        self.phaseout_factor.set_value(0)

        else:
            if self.config.has_taxes_and_credits:
                # this is not allowed, there's no income or production to tax
                raise ConfigurationError(
                    "Cannot set has_taxes_and_credits to True if has_fixed_OM is False. Calculations of income tax, royalties, "
                    "and mineral depletion require sales revenue, and production incentive is only applicable to plants "
                    "producing saleable product streams. Sale revenue and associated administrative costs are included in "
                    "the fixed O&M calculations."
                )

        # net present value params
        if self.config.has_net_present_value:

            self.discount_percentage = Param(
                initialize=10,
                units=pyunits.percent,
                doc="Rate of return used to discount future cash flows "
                "back to their present value. The value should be a percentage, "
                "for example 10 for a 10% discount. The NETL QGESS recommends "
                "setting the discount rate as the calculated after-tax weighted "
                "average cost of capital (ATWACC).",
            )
            self.plant_lifetime = Param(
                initialize=20,
                units=pyunits.year,
                doc="Length of operating period in years.",
            )

            if self.config.has_capital_expenditure_period:

                if not isinstance(self.config.capital_expenditure_percentages, list):
                    raise ValueError(
                        "Must set capital_expenditure_percentages as list of integer values on [0, 100]."
                    )

                if not sum(self.config.capital_expenditure_percentages) == 100:
                    raise ValueError(
                        f"Argument capital_expenditure_percentages has a sum of "
                        f"{sum(self.config.capital_expenditure_percentages)}. List must sum to 100 percent."
                    )

                self.capital_expenditure_percentages = Param(
                    range(len(self.config.capital_expenditure_percentages)),
                    initialize=self.config.capital_expenditure_percentages,
                    mutable=True,
                    doc="A list of values that sum to 100 representing how "
                    "capital costs are spread over a capital expenditure period; for "
                    "example, an input of [10, 60, 30] is parsed as a 3-year period "
                    "where capital costs are spread as 10% in year 1, 60% in year 2, "
                    "and 30% in year 3. The capital period precedes the operating "
                    "period, for example an input of [100] means that 100% of capital "
                    "expenses occur in the year preceding the operating period (t=-1). Set "
                    "to None to indicate no expenditure period, which means that all "
                    "capital expenses occur at the start of the plant lifetime (t=0).",
                )

            self.capital_escalation_percentage = Param(
                initialize=3.6,
                units=pyunits.percent,
                doc="Rate at which capital costs escalate during the "
                "capital expenditure period. The value should be a percentage, "
                "for example 10 for a 10% escalation rate. Set to 0 to indicate "
                "there is no cost escalation in the expenditure period.",
            )

            self.capital_loan_interest_percentage = Param(
                initialize=6,
                units=pyunits.percent,
                doc="Interest rate for capital equipment loan repayment."
                "The value should be a percentage, for example 10 for a 10% "
                "interest rate.",
            )

            self.capital_loan_repayment_period = Param(
                initialize=10,
                units=pyunits.year,
                doc="Length of loan repayment period in years.",
            )

            self.debt_percentage_of_capex = Param(
                initialize=50,
                units=pyunits.percent,
                doc="Percentage of CAPEX financed by debt; the value should be "
                "set as a percentage, for example 10 for a debt corresponding to 10% of "
                "the CAPEX. Set to zero to indicate no loans are taken out on capital.",
            )

            self.operating_inflation_percentage = Param(
                initialize=3,
                units=pyunits.percent,
                doc="Inflation rate for operating costs during the "
                "operating period. The value should be a percentage, for example "
                "10 for a 10% inflation rate. Set to 0 to indicate no inflation.",
            )

            self.revenue_inflation_percentage = Param(
                initialize=3,
                units=pyunits.percent,
                doc="Inflation rate for revenue during the operating "
                "period. The value should be a percentage, for example 10 for a "
                "10% inflation rate. Set to 0 to indicate no inflation.",
            )

        # economy of numbers params
        if self.config.has_economy_of_numbers:

            self.cum_num_units = Param(
                initialize=5,
                mutable=True,
                units=pyunits.dimensionless,
                doc="Cumulative number of units produced",
            )
            self.learning_rate = Param(
                initialize=0.04,
                mutable=True,
                units=pyunits.dimensionless,
                doc="The learning factor reflects the level of maturity of the unit/technology",
            )
            self.learning_rate_exponent = Expression(
                expr=(-log10(1 - self.learning_rate) / log10(2))
            )

        # annualization factors that depend on which tech is being used
        # TODO is there another way to do this besides setting strings manually?
        # pylint: disable=pointless-string-statement

        if self.config.tech == 10:

            tasc_toc_factor = 1.144
            tasc_toc_doc = "TASC/TOC factor calculated from UKy report using 3 year "
            "expenditure period with 10/60/30 % expenditure at 3.6% "
            "escalation at 2.94% debt interest rate with 7.84% return on "
            "equity, 26% combined federal/state tax, and 50/50 % debt and "
            "equity financed."

            fixed_charge_factor = 0.1002
            fixed_charge_doc = (
                "Fixed charge rate calculated from UKy report using a 26% "
            )
            "effective tax rate, a tax depreciation fraction of 2.231 over "
            "21 years of depreciation, a nominal capital recovery factor of "
            "0.0856, an after-tax weighted average cost of capital of 5.77%, "
            "= Present value of tax depreciation expense of 0.237"

        elif self.config.tech == 9:

            # 15% TPC for other owner's costs, 2.7% TPC for financing, 0.5% TPC for spare parts, 2% for preproduction
            self.pct_TPC = Param(
                initialize=(15 + 2.7 + 0.5 + 2) / 100,
                doc="Fixed percentage for other owners cost",
            )

            tasc_toc_factor = 1.070
            tasc_toc_doc = (
                "TASC/TOC factor from E. Lewis, S. McNaul, et al., Comparison"
            )
            " of Commercial, State-of-the-Art, Fossil-based Hydrogen Production"
            " Technologies, National Energy Technology Laboratory, Pittsburgh,"
            " September 30, 2022. Exhibit 3-32"

            fixed_charge_factor = 0.0690
            fixed_charge_doc = (
                "Fixed charge rate from E. Lewis, S. McNaul, et al., Comparison"
            )
            " of Commercial, State-of-the-Art, Fossil-based Hydrogen Production"
            " Technologies, National Energy Technology Laboratory, Pittsburgh,"
            " September 30, 2022. based on CRF values, real for three/five years"
            " X, nominal for three/five years Y"

        elif self.config.tech == 8:

            # 15% TPC for other owner's costs, 2.7% TPC for financing, 0.5% TPC for spare parts, 2% for preproduction
            self.pct_TPC = Param(
                initialize=(15 + 2.7 + 0.5 + 2) / 100,
                doc="Fixed percentage for other owners cost",
            )

            tasc_toc_factor = 1.047
            tasc_toc_doc = "TASC/TOC factor from S. McNaul, Screening Techno-economic"
            " Analysis of NETL Reactive Capture Technology, National Energy"
            " Technology Laboratory, Pittsburgh, September 30, 2022. Exhibit"
            " 4-8"

            fixed_charge_factor = 0.0664
            fixed_charge_doc = (
                "Fixed charge rate from S. McNaul, Screening Techno-economic"
            )
            " Analysis of NETL Reactive Capture Technology, National Energy"
            " Technology Laboratory, Pittsburgh, September 30, 2022. Exhibit"
            " 4-8"

        else:

            # 20.2% for total misc cost
            self.pct_TPC = Param(
                initialize=20.2 / 100, doc="Fixed percentage for other owners cost"
            )

            tasc_toc_factor = 1.093
            tasc_toc_doc = "TASC/TOC factor from Exhibit 3-7 reference 1, real for"
            " three years 1.093, real for five years 1.154, nominal for"
            " three years 1.242, nominal for five years 1.289"

            fixed_charge_factor = 0.0707
            fixed_charge_doc = "Fixed charge rate from Exhibit 3-5 based on CRF values,"
            " real for three/five years 0.0707, nominal for three/five years 0.0886"

        self.tasc_toc_factor = Param(
            initialize=tasc_toc_factor,
            mutable=True,
            doc=tasc_toc_doc,
        )

        self.fixed_charge_factor = Param(
            initialize=fixed_charge_factor,
            mutable=True,
            doc=fixed_charge_doc,
        )

        # report built params with default values, and built Expressions
        print("Costing class built with the following global parameters:")
        print("---------------------------------------------------------")
        param_dict = {}
        param_dict["base_currency"] = self.base_currency
        param_dict["CE_index_units"] = self.CE_index_units
        param_dict["base_period"] = self.base_period
        for p in self.component_data_objects(Param, descend_into=True):
            param_dict[p.name] = p
            # param_dict[p.name, ", ", pyunits.get_units(p)] = value(p)

        self.param_dir = {}
        self.param_dir["Value"] = {}
        self.param_dir["Units"] = {}
        self.param_dir["pos"] = {}

        count = 1
        for k, v in param_dict.items():
            self.param_dir["Value"][k] = value(v)
            self.param_dir["Units"][k] = pyunits.get_units(v)
            self.param_dir["pos"][k] = count
            count += 1

        df = DataFrame.from_dict(self.param_dir, orient="columns")
        del df["pos"]

        print("\n" + "=" * 84)
        print(f"{self.local_name}")
        print("-" * 84)
        stdout.write(textwrap.indent(stream_table_dataframe_to_string(df), " " * 4))
        print("\n" + "=" * 84 + "\n")

        if isinstance(self.Lang_factor, Expression):
            print("Lang factor built as Expression: \n", self.Lang_factor.expr)
        if self.config.has_economy_of_numbers:
            print(
                "EON learning rate exponent built as Expression: \n",
                self.learning_rate_exponent.expr,
            )

    def build_process_costs(
        self,
        # optional arguments that directly fix cost variables and bypass calculations
        total_purchase_cost=None,
        annual_fixed_operating_cost=None,
        annual_revenue=None,
        debt_expression=None,  # for NPV method
        # optional arguments for additional calculations and reporting
        feedstock_rate=None,
        production_rate=None,  # this could be power, total REE, water, CO2
        # required arguments for fixed_OM calculations
        pure_product_output_rates=None,
        mixed_product_output_rates=None,
        sale_prices=None,
        # required arguments for variable_OM calculations
        resources=None,  # for annual OPEX costs
        resource_prices=None,
        # optional arguments related to overnight costs
        land_cost=None,  # for startup and/or annual leasing costs
        fuel=None,  # extra inventory required as part of overnight costs
        feedstock=None,  # extra inventory required as part of overnight costs
        waste=None,  # extra storage required as part of overnight costs
        additional_waste_cost=None,
        chemicals=None,  # extra inventory required as part of overnight costs
        additional_chemicals_cost=None,
        chemicals_inventory=None,
        transport_per_unit_feedstock_cost=None,
        transport_per_unit_production_cost=None,
    ):
        """
        This method builds process-wide costing, including fixed and variable
        operating & maintenance costs, costs of production, cost of
        electricity and cost of production where production can be power or
        output of a produced or recovered component.

        Args:
            total_purchase_cost: user-defined value (MUSD in base year), Var,
                Param, or Expression for the total equipment purchase cost (not
                including installation or other plant costs). To use as the
                total plant cost, including installation, set the Lang_factor
                to 1.
            annual_fixed_operating_cost: user-defined value (MUSD/year in base year),
                Var, Param, or Expression for the total annual fixed operating cost,
                including all fixed, maintenance, and material costs.
            annual_revenue: user-defined value (MUSD/year in base year),
                Var, Param, or Expression for the total annual revenue from
                saleable products.
            debt_expression: user-defined value (MUSD in base year), Var,
                Param, or Expression for the total debt assumed to finance
                capital costs. This expression will override any value set for
                debt_percentage_of_capex if both are passed.
            feedstock_rate: Var, Param, or Expression for rate of feedstock
                input to the process. If passed, cost of per unit feedstock
                will be calculated and reported.
            production_rate: Var, Param, or Expression for the production rate of power
                or saleable product generated by the plant. If passed, cost of production
                will be calculated and reported.
            pure_product_output_rates: production rates of each pure product.
                Input should have the form {"name1": rate1, ...} where the name
                keys are strings and the rates are Vars, Params, or Expressions
                with Pyomo units of flow.
            mixed_product_output_rates: production rates of each product in the
                mixed basket. Input should have the form {"name1": rate1, ...}
                where the name keys are strings and the rates are Vars, Params,
                or Expressions with Pyomo units of flow.
            sale_prices: dictionary of products and prices to use in lieu of
                built-in sale prices. Input should have the form {"product1":
                price1, "product2": price2, ...} where the product keys match
                some or all keys in the `pure_product_output_rates` and/or
                `mixed_product_output_rates` arguments and the prices are Vars,
                Params, or Expressions with Pyomo units of currency per unit.
                Users should consult the built-in price dictionary in
                costing_dictionaries.py for available prices; note that the
                built-prices do not require input to this argument to be used
                as long as an entry in the `resources` argument matches a
                built-in price name.
            resources: dictionary of resources and rates to calculate variable
                operating costs for. Input should have the form {"resource1":
                rate1, "resource2": rate2, ...} where the resource keys are
                strings and the rates are Vars, Params, or Expressions with
                Pyomo units of flow.
            resource_prices: dictionary of resources and prices to use in lieu
                of built-in prices. Input should have the form {"resource1":
                price1, "resource2": price2, ...} where the resource keys match
                some or all keys in the `resources` argument and the prices
                are Vars, Params, or Expressions with Pyomo units of currency
                per unit. Users should consult the built-in price dictionary
                in costing_dictionaries.py for available prices; note that
                the built-prices do not require input to this argument to be
                used as long as an entry in the `resources` argument matches
                a built-in price name.
            land_cost: Expression, Var or Param to calculate land costs. If this
                argument has units of cost per time, an annual land cost will be
                included in the total variable OM cost, and plant overhead for
                tech=10 (REE) plants. If this argument has units of cost, a
                one-time land cost will be included in the overnight cost. If this
                argument has any other units or no units, an error will be returned.
            fuel: list of strings to define which resources require fuel storage and
                handling.
            feedstock: list of strings to define which resources require feedstock storage and
                handling.
            waste: list of strings to define which resources require waste
                storage and handling.
            additional_waste_cost: Expression, Var or Param to calculate additional waste
                costs, e.g. initial fills, additional resources, and so. If this argument
                has units of cost per time, an annual additional waste cost will be included
                in the total variable OM cost, and plant overhead for tech=10 (REE) plants.
                If this argument has units of cost, a one-time additional waste cost will be
                included in the overnight cost. If this argument has any other units or no
                units, an error will be returned.
            chemicals: list of strings to define which resources are chemicals
                that require special storage and handling.
            additional_chemicals_expression: Expression, Var or Param to calculate
                additional chemical costs, e.g. initial fills, additional resources,
                and so. If this argument has units of cost per time, an annual additional
                chemicals cost will be included in the total variable OM cost, and plant
                overhead for tech=10 (REE) plants. If this argument has units of cost, a
                one-time additional chemicals cost will be included in the overnight cost.
                If this argument has any other units or no units, an error will be returned.
            chemicals_inventory: list of strings to define which resources are
                chemicals that require inventory stock.
            transport_per_unit_feedstock_cost: Expression, Var or Param to calculate transport costs of feedstock
                per unit of feedstock.
            transport_per_unit_ production_cost: Expression, Var or Param to calculate transport costs of product
                per unit of product.

        """

        if hasattr(self, "built"):  # costing already exists
            raise AttributeError(
                f"Costing for the block {self} already exists. Please ensure that "
                f"the costing build method is not called twice on the same "
                f"model. Create a new flowsheet costing block, or if needed "
                f"use delattr() to remove the preexisting costing."
            )
        else:
            self.built = True

        # calculate BEC, TPC by summing over all the blocks
        self.get_total_BEC_and_TPC(total_purchase_cost)

        # create chemicals list to reference in variable operating cost method
        if chemicals is None:
            self.chemicals_list = []
        else:
            self.chemicals_list = chemicals

        # create waste list to reference in variable operating cost method
        if waste is None:
            self.waste_list = []
        else:
            self.waste_list = waste

        # define feed input, if passed
        if feedstock_rate is not None:
            if (
                pyunits.get_units(feedstock_rate) == pyunits.dimensionless
            ):  # require units
                raise UnitsError(
                    "The argument feedstock_rate was passed as a dimensionless "
                    "quantity with no units. Please ensure that the feed "
                    "rate is passed in units of consumption / time."
                )

        # define land cost, if passed
        if land_cost is not None:
            expr_units = pyunits.get_units(land_cost)

            if expr_units is None:
                raise ValueError("Expression land_cost has no units defined.")

            try:
                # flag to include land cost in overnight cost
                self.land_cost = Expression(
                    expr=pyunits.convert(land_cost, to_units=self.CE_index_units)
                )
                self.land_cost_reoccurrence = "one_time"

            except UnitsError:
                try:
                    # flag to include land cost in total variable OM cost
                    self.land_cost = Expression(
                        expr=pyunits.convert(
                            land_cost, to_units=self.CE_index_units / pyunits.year
                        )
                    )
                    self.land_cost_reoccurrence = "annual"

                except UnitsError:
                    raise ValueError(
                        f"Expression land_cost units {expr_units} are not compatible with "
                        f"{self.CE_index_units} or {self.CE_index_units/pyunits.year}. The "
                        f"expression must be compatible with cost units to be included in "
                        f"the overnight cost, or cost per time units to be included in the "
                        f"total variable operating cost."
                    )
        else:
            self.land_cost_reoccurrence = "one_time"
            self.land_cost = Expression(expr=0 * self.CE_index_units)

        # define additional chemicals cost, if passed

        if additional_chemicals_cost is not None:

            expr_units = pyunits.get_units(additional_chemicals_cost)

            if expr_units is None:
                raise ValueError(
                    "Expression additional_chemicals_cost has no units defined."
                )

            try:
                # flag to include land cost in overnight cost
                self.additional_chemicals_cost = Expression(
                    expr=pyunits.convert(
                        additional_chemicals_cost, to_units=self.CE_index_units
                    )
                )
                self.additional_chemicals_cost_reoccurrence = "one_time"

            except UnitsError:
                try:
                    # flag to include land cost in total variable OM cost
                    self.additional_chemicals_cost = Expression(
                        expr=pyunits.convert(
                            additional_chemicals_cost,
                            to_units=self.CE_index_units / pyunits.year,
                        )
                    )
                    self.additional_chemicals_cost_reoccurrence = "annual"

                except UnitsError:
                    raise ValueError(
                        f"Expression additional_chemicals_cost units {expr_units} are not compatible with "
                        f"{self.CE_index_units} or {self.CE_index_units/pyunits.year}. The "
                        f"expression must be compatible with cost units to be included in "
                        f"the overnight cost, or cost per time units to be included in the "
                        f"total variable operating cost."
                    )
        else:
            self.additional_chemicals_cost_reoccurrence = "one_time"
            self.additional_chemicals_cost = Expression(expr=0 * self.CE_index_units)

        # define waste cost, if passed

        if additional_waste_cost is not None:

            expr_units = pyunits.get_units(additional_waste_cost)

            if expr_units is None:
                raise ValueError(
                    "Expression additional_waste_cost has no units defined."
                )

            try:
                # flag to include land cost in overnight cost
                self.additional_waste_cost = Expression(
                    expr=pyunits.convert(
                        additional_waste_cost, to_units=self.CE_index_units
                    )
                )
                self.additional_waste_cost_reoccurrence = "one_time"

            except UnitsError:
                try:
                    # flag to include land cost in total variable OM cost
                    self.additional_waste_cost = Expression(
                        expr=pyunits.convert(
                            additional_waste_cost,
                            to_units=self.CE_index_units / pyunits.year,
                        )
                    )
                    self.additional_waste_cost_reoccurrence = "annual"

                except UnitsError:
                    raise ValueError(
                        f"Expression additional_waste_cost units {expr_units} are not compatible with "
                        f"{self.CE_index_units} or {self.CE_index_units/pyunits.year}. The "
                        f"expression must be compatible with cost units to be included in "
                        f"the overnight cost, or cost per time units to be included in the "
                        f"total variable operating cost."
                    )
        else:
            self.additional_waste_cost_reoccurrence = "one_time"
            self.additional_waste_cost = Expression(expr=0 * self.CE_index_units)

        # call fixed OM method
        if self.config.has_fixed_OM:
            self.get_fixed_OM_costs(
                annual_fixed_operating_cost=annual_fixed_operating_cost,
                annual_revenue=annual_revenue,
                pure_product_output_rates=pure_product_output_rates,
                mixed_product_output_rates=mixed_product_output_rates,
                sale_prices=sale_prices,
                production_rate=production_rate,
            )

        # call variable OM method
        if self.config.has_variable_OM:
            self.get_variable_OM_costs(
                resources=resources,
                resource_prices=resource_prices,
            )

        if self.config.tech == 10:
            # plant overhead variable OPEX accounts for all preproduction
            #      including annual land, chemicals, and waste costs
            # one-time land cost should be included in overnight cost
            # one-time chemicals or waste, e.g. initial fills, go here too
            self.total_overnight_capital = Expression(
                expr=(
                    self.total_TPC
                    + (
                        self.land_cost
                        if self.land_cost_reoccurrence == "one_time"
                        else 0 * self.CE_index_units
                    )
                    + (
                        self.additional_chemicals_cost
                        if self.additional_chemicals_cost_reoccurrence == "one_time"
                        else 0 * self.CE_index_units
                    )
                    + (
                        self.additional_waste_cost
                        if self.additional_waste_cost_reoccurrence == "one_time"
                        else 0 * self.CE_index_units
                    )
                )
            )
        else:
            # power plants break overnight cost down by plant components

            # some components relate to fixed operating costs
            if self.config.has_fixed_OM and annual_fixed_operating_cost is None:
                self.six_month_labor = Expression(
                    expr=(
                        self.annual_labor_cost
                        + self.maintenance_labor_cost
                        + self.admin_and_support_labor_cost
                    )
                    * pyunits.year
                    / 2
                )

            # some components relate to variable operating costs
            if self.config.has_variable_OM:

                non_fuel_feedstock_waste_resources = {}
                for i in resources:
                    non_fuel_feedstock_waste_resources[i] = resources[i]
                if fuel is not None:
                    if (
                        self.config.tech == 8 or self.config.tech == 9
                    ):  # at 100% capacity, 25% of a month
                        self.fuel_cost_OC = Expression(
                            expr=sum(self.variable_operating_costs[0, i] for i in fuel)
                            / 12
                            * 0.25
                            / self.capacity_factor,
                            doc="Owner's costs - 0.25 months of fuel costs",
                        )
                    else:  # at operating capacity, 2.25 months
                        self.fuel_cost_OC = Expression(
                            expr=sum(self.variable_operating_costs[0, i] for i in fuel)
                            / 12
                            * 2.25,
                            doc="Owner's costs - 2.25 months of fuel costs",
                        )
                    for i in fuel:
                        if i in non_fuel_feedstock_waste_resources:
                            del non_fuel_feedstock_waste_resources[
                                i
                            ]  # remove fuel from the list

                if feedstock is not None:
                    if (
                        self.config.tech == 8 or self.config.tech == 9
                    ):  # at 100% capacity, 25% of a month
                        self.feedstock_cost_OC = Expression(
                            expr=sum(
                                self.variable_operating_costs[0, i] for i in feedstock
                            )
                            / 12
                            * 0.25
                            / self.capacity_factor,
                            doc="Owner's costs - 0.25 months of feedstock costs",
                        )
                    else:  # at operating capacity, 2.25 months
                        self.feedstock_cost_OC = Expression(
                            expr=sum(
                                self.variable_operating_costs[0, i] for i in feedstock
                            )
                            / 12
                            * 0.25
                            * 0,  # exclude feedstock cost
                            doc="Owner's costs - power plants don't include feedstock costs",
                        )
                    for i in feedstock:
                        if i in non_fuel_feedstock_waste_resources:
                            del non_fuel_feedstock_waste_resources[
                                i
                            ]  # remove feedstock from the list

                if waste is not None:
                    if (
                        self.config.tech == 8 or self.config.tech == 9
                    ):  # at 100% capacity, 1 month
                        self.waste_cost_OC = Expression(
                            expr=(
                                sum(self.variable_operating_costs[0, i] for i in waste)
                                / 12
                                / self.capacity_factor
                            )
                        )
                    else:  # at operating capacity, 1 month
                        self.waste_cost_OC = Expression(
                            expr=(
                                sum(self.variable_operating_costs[0, i] for i in waste)
                                / 12
                            )
                        )
                    for i in waste:
                        if i in non_fuel_feedstock_waste_resources:
                            del non_fuel_feedstock_waste_resources[
                                i
                            ]  # remove waste from the list

                if self.config.tech == 8 or self.config.tech == 9:  # at 100% capacity
                    self.non_fuel_feedstock_waste_OC = Expression(
                        expr=(
                            sum(
                                self.variable_operating_costs[0, i]
                                for i in non_fuel_feedstock_waste_resources
                            )
                            / 12
                            / self.capacity_factor
                        )
                    )
                else:  # at operating capacity
                    self.non_fuel_feedstock_waste_OC = Expression(
                        expr=(
                            sum(
                                self.variable_operating_costs[0, i]
                                for i in non_fuel_feedstock_waste_resources
                            )
                            / 12
                        )
                    )

                if chemicals is not None:
                    self.chemicals_cost_OC = Expression(
                        expr=(
                            sum(self.variable_operating_costs[0, i] for i in chemicals)
                            / 2  # six months of chemicals
                        )
                    )

                if chemicals_inventory is not None:
                    if (
                        self.config.tech == 8 or self.config.tech == 9
                    ):  # at 100% capacity, spare parts already included in other owner's costs
                        self.chemicals_inventory_cost_OC = Expression(
                            expr=(
                                (
                                    sum(
                                        self.variable_operating_costs[0, i]
                                        for i in chemicals_inventory
                                    )
                                    / 6
                                    / self.capacity_factor
                                )
                                * pyunits.year  # two months of chemicals inventory
                            )
                        )
                    else:  # at operating capacity, need to include spare parts here
                        self.chemicals_inventory_cost_OC = Expression(
                            expr=(
                                (
                                    sum(
                                        self.variable_operating_costs[0, i]
                                        for i in chemicals_inventory
                                    )
                                    / 6
                                )
                                * pyunits.year  # two months of chemicals inventory
                                + 0.005 * self.total_TPC  # inventory spare parts
                            )
                        )

            # define overnight cost
            # if fixed OM is not defined or user defines as a lumped value, no way of knowing the labor cost so skip it
            include_labor_components = (
                self.config.has_fixed_OM and annual_fixed_operating_cost is None
            )

            self.total_overnight_capital = Expression(
                expr=self.total_TPC
                # include labor costs if present
                + (
                    self.six_month_labor
                    if include_labor_components
                    else 0 * self.CE_index_units
                )
                + (  # include material and resource costs if present
                    (
                        (
                            self.maintenance_material_cost / 12 / self.capacity_factor
                            if self.config.tech == 8 or self.config.tech == 9
                            else self.maintenance_material_cost / 12
                        )
                        if include_labor_components
                        else 0 * self.CE_index_units / pyunits.year
                    )  # 1 month materials
                    + (
                        self.non_fuel_feedstock_waste_OC
                        if self.config.has_variable_OM
                        else 0 * self.CE_index_units / pyunits.year
                    )  # 1 month nonfuel consumables
                    + (
                        self.waste_cost_OC
                        if waste is not None
                        else 0 * self.CE_index_units / pyunits.year
                    )  # 1 month waste
                    # inventory capital costs
                    + (
                        self.fuel_cost_OC
                        if fuel is not None
                        else 0 * self.CE_index_units / pyunits.year
                    )  # initial fuel supply
                    + (
                        self.feedstock_cost_OC
                        if feedstock is not None
                        else 0 * self.CE_index_units / pyunits.year
                    )  # initial feedstock supply
                    # Other costs
                    + (
                        self.chemicals_cost_OC
                        if chemicals is not None
                        else 0 * self.CE_index_units / pyunits.year
                    )  # Initial Cost for Catalyst and Chemicals
                    + (
                        self.chemicals_inventory_cost_OC / pyunits.year
                        if chemicals_inventory is not None
                        else 0 * self.CE_index_units / pyunits.year
                    )  # Initial Cost for Catalyst and Chemicals Inventory
                )
                * 1
                * pyunits.year  # overhead costs for 1 year
                + (
                    self.land_cost
                    if self.land_cost_reoccurrence == "one_time"
                    else 0 * self.CE_index_units
                )
                + (
                    self.additional_chemicals_cost
                    if self.additional_chemicals_cost_reoccurrence == "one_time"
                    else 0 * self.CE_index_units
                )
                + (
                    self.additional_waste_cost
                    if self.additional_waste_cost_reoccurrence == "one_time"
                    else 0 * self.CE_index_units
                )
                + self.total_TPC
                * self.pct_TPC  # other owners costs (other + spare parts + financing + other pre-production)
            )

        # build expressions - factors are defined as global params based on tech
        self.total_as_spent_cost = Expression(
            expr=self.total_overnight_capital * self.tasc_toc_factor
        )

        self.annualized_cost = Expression(
            expr=self.fixed_charge_factor * self.total_as_spent_cost
        )

        # calculate taxes now that all plant costs have been defined
        if self.config.has_taxes_and_credits:
            self.calculate_taxes()

        # call NPV method now that all CAPEX, OPEX, revenue, and taxes are defined
        # if any are user-supplied values, those are set in their respective methods prior to this point
        # if any required objects are not built, the method called below will raise an appropriate exception
        if self.config.has_net_present_value:
            self.calculate_net_present_value(debt_expression)

        # transport cost of production and feedstock

        if (
            transport_per_unit_production_cost is not None
            and production_rate is not None
        ):
            try:
                self.transport_production_cost = Expression(
                    expr=pyunits.convert(
                        transport_per_unit_production_cost * production_rate,
                        to_units=self.CE_index_units / pyunits.year,
                    )
                )
            except UnitsError as e:
                raise InconsistentUnitsError(
                    "production_rate",
                    "transport_per_unit_production_cost",
                    f"Expression transport_production_cost failed to build with error: {e} ",
                )
        else:
            self.transport_production_cost = Expression(
                expr=0 * self.CE_index_units / pyunits.year
            )

        if (
            transport_per_unit_feedstock_cost is not None
            and production_rate is not None
        ):
            try:
                self.transport_feedstock_cost = Expression(
                    expr=pyunits.convert(
                        transport_per_unit_feedstock_cost * feedstock_rate,
                        to_units=self.CE_index_units / pyunits.year,
                    )
                )
            except UnitsError as e:
                raise InconsistentUnitsError(
                    "feedstock_rate",
                    "transport_per_unit_feedstock_cost",
                    f"Expression transport_feedstock_cost failed to build with error: {e} ",
                )
        else:
            self.transport_feedstock_cost = Expression(
                expr=0 * self.CE_index_units / pyunits.year
            )

        # levelized cost per unit production
        # store valid units to check production, feedstock rates against
        ureg = pyunits.get_units(
            pyunits.kg
        )._pint_unit._REGISTRY  # pylint: disable=protected-access
        valid_units = {}
        for name, unit in ureg._units.items():  # pylint: disable=protected-access
            try:
                dim = (1 * ureg.parse_units(name)).dimensionality
                valid_units.setdefault(dim, []).append(name)
            except Exception:
                pass

        if production_rate is not None:
            # TODO look into protected access issue with pint
            # check that production rate units are valid_units/time
            dim = pyunits.get_units(
                production_rate
            )._pint_unit.dimensionality  # pylint: disable=protected-access
            dim_numerator = pyunits.get_units(
                production_rate * pyunits.s
            )._pint_unit.dimensionality  # pylint: disable=protected-access

            allowed_production_units = (
                True  # check if production units are compatible with LCOP
            )

            if (
                "[time]" not in dim
                or dim["[time]"] >= 0
                or dim_numerator not in valid_units
            ):  # units not in form valid_units/time
                allowed_production_units = False

            if (
                dim
                == pyunits.get_units(
                    pyunits.MW
                )._pint_unit.dimensionality  # pylint: disable=protected-access
            ):  # objects has UOM production/time, but production is not a valid Pyomo UOM container
                allowed_production_units = True

            if not allowed_production_units:
                raise UnitsError(
                    f"The argument production_rate was passed with the units {pyunits.get_units(production_rate)} "
                    f"with dimensions {dim}. Please ensure that the production rate has units of valid_units/time."
                )

            self.additional_cost_of_production = Var(
                initialize=0,
                doc="Additional cost to be added to the LCOP calculations; use 70.9 for NGCC baseline "
                + "report and 43.3 for NGCC without carbon capture, both for production rate in MWh/[time]",
                units=self.base_currency
                / pyunits.get_units(production_rate * pyunits.year),
            )
            self.additional_cost_of_production.fix()

            self.cost_of_production = Expression(
                expr=(
                    pyunits.convert(
                        (
                            self.annualized_cost / pyunits.year
                            + (
                                self.total_fixed_OM_cost
                                if self.config.has_fixed_OM
                                else 0 * self.CE_index_units / pyunits.year
                            )
                            + (
                                self.total_variable_OM_cost[0] * self.capacity_factor
                                if self.config.has_variable_OM
                                else 0 * self.CE_index_units / pyunits.year
                            )
                            + (
                                self.net_tax_owed
                                if self.config.has_taxes_and_credits
                                else 0 * self.CE_index_units / pyunits.year
                            )
                            + self.transport_production_cost
                            + self.transport_feedstock_cost
                        )
                        * pyunits.year
                        / (production_rate * pyunits.year * self.capacity_factor),
                        to_units=self.base_currency
                        / pyunits.get_units(production_rate * pyunits.year),
                    )
                    + self.additional_cost_of_production
                )
            )

        # levelized cost per unit feedstock

        if feedstock_rate is not None:
            # TODO look into protected access issue with pint
            # check that feedstock rate units are valid_units/time
            dim = pyunits.get_units(
                feedstock_rate
            )._pint_unit.dimensionality  # pylint: disable=protected-access
            dim_numerator = pyunits.get_units(
                feedstock_rate * pyunits.s
            )._pint_unit.dimensionality  # pylint: disable=protected-access

            allowed_feedstock_units = (
                True  # check if feedstock units are compatible with LCOP
            )

            if (
                "[time]" not in dim
                or dim["[time]"] >= 0
                or dim_numerator not in valid_units
            ):  # units not in form valid_units/time
                allowed_feedstock_units = False

            if (
                dim
                == pyunits.get_units(
                    pyunits.MW
                )._pint_unit.dimensionality  # pylint: disable=protected-access
            ):  # objects has UOM feedstock/time, but feedstock is not a valid Pyomo UOM container
                allowed_feedstock_units = True

            if not allowed_feedstock_units:
                raise UnitsError(
                    f"The argument feedstock_rate was passed with the units {pyunits.get_units(feedstock_rate)} "
                    f"with dimensions {dim}. Please ensure that the feedstock rate has units of valid_units/time."
                )

            self.additional_cost_of_feedstock = Var(
                initialize=0,
                doc="Additional cost to be added to the LCOF calculations",
                units=self.base_currency
                / pyunits.get_units(feedstock_rate * pyunits.year),
            )
            self.additional_cost_of_feedstock.fix()

            self.cost_of_feedstock = Expression(
                expr=(
                    pyunits.convert(
                        (
                            self.annualized_cost / pyunits.year
                            + (
                                self.total_fixed_OM_cost
                                if self.config.has_fixed_OM
                                else 0 * self.CE_index_units / pyunits.year
                            )
                            + (
                                self.total_variable_OM_cost[0] * self.capacity_factor
                                if self.config.has_variable_OM
                                else 0 * self.CE_index_units / pyunits.year
                            )
                            + (
                                self.net_tax_owed
                                if self.config.has_taxes_and_credits
                                else 0 * self.CE_index_units / pyunits.year
                            )
                            + self.transport_production_cost
                            + self.transport_feedstock_cost
                        )
                        * pyunits.year
                        / (feedstock_rate * pyunits.year * self.capacity_factor),
                        to_units=self.base_currency
                        / pyunits.get_units(feedstock_rate * pyunits.year),
                    )
                    + self.additional_cost_of_feedstock
                )
            )

    def get_equipment_costing(
        blk,
        cost_accounts,
        scaled_param,
        tech,  # a flag to search through the JSON cost account dictionary
        ccs,  # flag whether technology does not include CCS ("A") or does ("B")
        n_equip=1,
        scale_down_parallel_equip=False,
        CE_index_year="2021",
        additional_costing_params=None,  # this is for user-defined or external cost accounts
        use_additional_costing_params=False,  # guard against user accidentally overwriting built-in accounts
        multiply_project_conting=True,
    ):
        """
        Equipment Costing Method
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
        8. Polymer Layers accounts
        9. Sensors & Controls accounts
        10. University of Kentucky Fire Clay Seam (Hazard No. 4) Rejects

        This method computes the capital cost of units and main components of
        the power plant, and requires a few arguments to build a constraint as
        part of your main model.

        Args:
            blk: A unit-level costing block where costing variables and
                constraints can be added to
            cost_accounts: A list of accounts to be included in the total cost
            scaled_param: the process parameter for the system(s) being costed;
                this is the total flow for all parallel trains of the system(s)
            tech: integer representing the above categories
            ccs: 'A' or 'B' representing no CCS or CCS
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
                train (twice the capacity). If duplicating a fixed operation
                into multiple parallel trains, use the default value ('False').
            CE_index_year: year for cost basis, e.g. "2021" to use 2021 dollars
            additional_costing_params: user-defined dictionary to append to
                existing cost accounts dictionary, e.g. UKy REE from PrOMMiS
            use_additional_costing_params: Boolean flag to use additional costing
                parameters when account names conflict with existing accounts data
            multiply_project_conting: whether project contingency should be multiplied
                by the entire total plant cost, or added with the other contigencies, e.g.
                "True" would yield
                TPC = BEC * (1 + eng_fee + process_conting) * (1 + project_conting)
                whereas "False" would yield
                TPC = BEC * (1 + eng_fee + process_conting + project_conting)

        The appropriate scaling parameters for various cost accounts can be
        found in the QGESS on capital cost scaling (Report
        #DOE/NETL-341/013113).
        The correct units for the reference parameters are found in the BBR4
        COE spreadsheet.

        """

        # check to see if a costing block already exists
        if (
            blk.parent_block().name
            in blk.config.flowsheet_costing_block._registered_unit_costing  # pylint: disable=protected-access
        ):
            raise AttributeError(
                f"{blk.name} already has an attribute costing. "
                f"Check that you are not calling get_costing "
                f"twice on the same model and if needed "
                f"use delattr() to remove the preexisting costing."
            )

        # validate currency units
        try:
            CE_index_units = getattr(pyunits, "MUSD_" + CE_index_year)
        except AttributeError:
            raise AttributeError(
                f"CE_index_year {CE_index_year} is not a valid currency base option. "
                f"Valid CE index options include CE500, CE394 and years from "
                f"1990 to 2023."
            )

        # preloaded account handling
        if isinstance(cost_accounts, str):
            if tech in [1, 2, 3, 4, 5, 6, 7]:
                (
                    PC_preloaded_accounts,
                    IGCC_preloaded_accounts,
                    NGCC_preloaded_accounts,
                    AUSC_preloaded_accounts,
                ) = define_preloaded_accounts()
                if tech in [1, 2]:
                    cost_accounts = PC_preloaded_accounts[cost_accounts]
                elif tech in [3, 4, 5]:
                    cost_accounts = IGCC_preloaded_accounts[cost_accounts]
                elif tech == 6:
                    cost_accounts = NGCC_preloaded_accounts[cost_accounts]
                elif tech == 7:
                    cost_accounts = AUSC_preloaded_accounts[cost_accounts]
            elif tech in [8, 9, 10]:
                cost_accounts = (
                    dict()
                )  # empty dictionary to append these additional costing accounts to
            else:
                raise AttributeError(
                    f"{blk.name} technology not supported with tech = {tech}, please use an integer value between 1 and 10."
                )

        # Load data into cost dictionary

        # pull data for each account into dictionaries

        costing_params = (
            {}
        )  # dictionary to append all pre-built and user-defined cost accounts to
        accounts_to_merge = (
            []
        )  # place to store cost dictionaries that will be checked sequentially

        process_params = {}
        reference_units = {}
        account_names = {}
        exponents = {}
        reference_costs = {}
        reference_cost_units = {}
        reference_costs_init = {}
        reference_params = {}

        if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            cost_scaling_fractions = {}
            engineering_fees = {}
            process_contingencies = {}
            project_contingencies = {}

            # load BB costing dictionary
            BB_costing_params = load_BB_costing_dictionary()
            accounts_to_merge.append(BB_costing_params)

            # load generic ccs costing dictionary
            generic_ccs_costing_params = load_generic_ccs_costing_dictionary()
            accounts_to_merge.append(generic_ccs_costing_params)

            # append them
            # costing_params.append(BB_costing_params)

        # add additional costing params - enforce that this is a list of dictionaries for backwards compatibility
        if additional_costing_params is not None:

            if isinstance(additional_costing_params, dict):
                raise TypeError(
                    "additional_costing_params must be a list of dicts, not a single dict, "
                    "e.g. [{'1': data}, {'2': data},] or [{'1': data},] and not {'1': data}."
                )

            if not (
                isinstance(additional_costing_params, list)
                and all(isinstance(item, dict) for item in additional_costing_params)
            ):
                raise TypeError(
                    "additional_costing_params must be a list of dicts, not a single dict, "
                    "e.g. [{'1': data}, {'2': data},] or [{'1': data},] and not {'1': data}."
                )

            # add all new costing parameter dictionaries to the existing list
            accounts_to_merge.extend(additional_costing_params)

        # for compatibility with potential custom accounts, the loop handles
        # new technologies, new CCS types for existing technologies, and new
        # accounts for existing technology-CCS type pairs
        # Users should not be adding new entries for existing accounts

        for (
            new_costing_params
        ) in accounts_to_merge:  # merge new dictionaries sequentially
            # adding ccs and any provided custom params to the ccs dictionary
            # need to "freeze" dict so it is hashable for merging keys
            frozen_dict = {**costing_params}
            for techkey, techval in new_costing_params.items():
                if (
                    techkey in frozen_dict
                ):  # if techkey already exists, append any new ccs types
                    for ccskey, ccsval in new_costing_params[techkey].items():
                        if (
                            ccskey in frozen_dict[techkey]
                        ):  # if ccskey already exists, append any new accounts
                            for accountkey, accountval in new_costing_params[techkey][
                                ccskey
                            ].items():
                                if (
                                    accountkey in frozen_dict[techkey][ccskey]
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
                        if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
                            cost_scaling_fractions[account, processparam] = (
                                costing_params[str(tech)][ccs][account][
                                    "Cost scaling fraction"
                                ][i]
                            )

                elif isinstance(process_params[account], str):
                    reference_params[account] = costing_params[str(tech)][ccs][account][
                        "RP Value"
                    ]

                if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
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
                            raise ValueError(
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
                except UnitsError:
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

                if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
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

        if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
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
            bounds=(0, None),
            doc="scaled bare erected cost in millions of USD",
            units=CE_index_units,
        )

        if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            blk.total_plant_cost = Var(
                cost_accounts,
                initialize=reference_costs_init,
                bounds=(0, None),
                doc="total plant cost in millions of USD",
                units=CE_index_units,
            )

        # rule for scaling BEC
        @blk.Constraint(cost_accounts)
        def bare_erected_cost_eq(costing, i):
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
            elif ref_cost_units[0] == "K":  # thousands of USD
                ref_cost_units = getattr(pyunits, "kUSD_" + ref_cost_units[1])
            elif ref_cost_units[0] == "M":  # millions of USD
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
                                costing.cost_scaling_fracs[i, p]
                                if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]
                                else 1
                            )
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

        if tech in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
            # rule for calculating TPC
            if multiply_project_conting is True:

                @blk.Constraint(cost_accounts)
                def total_plant_cost_eq(blk, i):
                    return blk.total_plant_cost[i] == blk.bare_erected_cost[i] * (
                        1 + blk.eng_fee[i] + blk.process_conting[i]
                    ) * (1 + blk.project_conting[i])

            else:

                @blk.Constraint(cost_accounts)
                def total_plant_cost_eq(blk, i):
                    return blk.total_plant_cost[i] == blk.bare_erected_cost[i] * (
                        1
                        + blk.eng_fee[i]
                        + blk.process_conting[i]
                        + blk.project_conting[i]
                    )

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

    def get_total_BEC_and_TPC(b, total_purchase_cost):
        # This method accepts a flowsheet-level costing block
        b.BEC_list = []
        b.custom_fixed_costs_list = []
        b.custom_variable_costs_list = []
        # blocks that already include TPC factors shouldn't be multiplied by Lang factor, but should be summed with plant TPC
        b.BEC_blocks_with_TPC_list = []  # BEC to exclude when applying Lang factor
        b.TPC_blocks_with_TPC_list = []  # TPC of those blocks to add in at the end

        if total_purchase_cost is not None:
            b.total_BEC = Var(
                initialize=pyunits.convert(
                    total_purchase_cost, to_units=b.CE_index_units
                ),
                units=b.CE_index_units,
            )
            b.total_BEC.fix()
            b.BEC_list.append(b.total_BEC)
        else:  # calculate total BEC by looping over equipment costing blocks

            for o in b.parent_block().component_objects(descend_into=True):
                # look for costing blocks
                if o.name in [
                    block.name for block in b._registered_unit_costing
                ]:  # pylint: disable=protected-access
                    if hasattr(
                        o, "total_plant_cost"
                    ):  # block that should not be multiplied by Lang factor
                        for key in o.bare_erected_cost:
                            b.BEC_blocks_with_TPC_list.append(
                                pyunits.convert(
                                    o.bare_erected_cost[key], to_units=b.CE_index_units
                                )
                            )

                            # the keys should match, so just append the corresponding TPC for later
                            b.TPC_blocks_with_TPC_list.append(
                                pyunits.convert(
                                    o.total_plant_cost[key], to_units=b.CE_index_units
                                )
                            )

                    if hasattr(o, "bare_erected_cost"):  # added from cost accounts
                        for key in o.bare_erected_cost:
                            b.BEC_list.append(
                                pyunits.convert(
                                    o.bare_erected_cost[key], to_units=b.CE_index_units
                                )
                            )
                    elif hasattr(
                        o, "capital_cost"
                    ):  # added from custom model or WaterTAP model
                        b.BEC_list.append(
                            pyunits.convert(o.capital_cost, to_units=b.CE_index_units)
                        )
                    if hasattr(o, "fixed_operating_cost"):
                        b.custom_fixed_costs_list.append(
                            pyunits.convert(
                                o.fixed_operating_cost,
                                to_units=b.CE_index_units / pyunits.year,
                            )
                        )
                    if hasattr(o, "variable_operating_cost"):
                        b.custom_variable_costs_list.append(
                            pyunits.convert(
                                o.variable_operating_cost,
                                to_units=b.CE_index_units / pyunits.year,
                            )
                        )

        b.total_BEC = Var(
            initialize=100,
            bounds=(0, None),
            doc="Plant bare erected cost",
            # assume that BEC is in millions of USD_year, where year is the CE_index_year users set
            units=b.CE_index_units,
        )

        b.total_TPC = Var(
            initialize=100,
            bounds=(0, None),
            doc="Total plant cost",
            # assume that TPC is in millions of USD_year, where year is the CE_index_year users set
            units=b.CE_index_units,
        )

        # add other plant costs to catch non-equipment capital costs, e.g. reagent fills, packing to include in CAPEX
        b.other_plant_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Additional plant costs",
            units=b.CE_index_units,
        )
        b.other_plant_costs.fix(1e-12)

        @b.Constraint()
        def total_BEC_eq(c):
            if len(c.BEC_list) == 0:
                c.BEC_list.append(1e-12 * b.CE_index_units)
                return c.total_BEC == 1e-12 * b.CE_index_units
            else:
                return c.total_BEC == sum(b.BEC_list)

        # calculate economy of numbers if enabled
        # Economy of Numbers (EoN) estimates the future profitability of novel/First-of-A-Kind (FOAK)
        # equipment. This is because the cost of manufacturing a piece of equipment tends to decline
        # as the cumulative production quantity rises, resulting from a consistent improvement in
        # technical know-how.

        # Y = A(X^-b)

        # b = - log(1-R)/log(2)
        # where Y is the cost of the Nth-of-A-Kind (NOAK) of the equipment, A is the cost of the FOAK,
        # X is the cumulative number of units, b is the learning rate exponent, and R is the learning
        # rate constant.

        # The equations above are derived from Faber G, Ruttinger A, Strunge T, Langhorst T, Zimmermann A,
        # van der Hulst M, Bensebaa F, Moni S and Tao L (2022) Adapting Technology Learning Curves for
        # Prospective Techno-Economic and Life Cycle Assessments of Emerging Carbon Capture and Utilization
        # Pathways. Front. Clim. 4:820261. doi: 10.3389/fclim.2022.820261
        # Learning rates vary between 0.01 - 0.1, depending on the level of maturity (i.e., experimental,
        # growing, proven, etc.) as described in Rubin, E. S., Mantripragada, H., and Zhai, H., "An"
        # "Assessment of the NETL Cost Estimation Methodology". Department of Engineering and Public Policy,
        # Carnegie Mellon University, Pittsburgh, PA (2016). p. 31, Fig. 6-4.

        if b.config.has_economy_of_numbers:
            b.NOAK_factor = Var(
                initialize=1,
                bounds=(0, None),
                doc="Factor to calculate the Nth-of-A-Kind cost of the unit",
                units=pyunits.dimensionless,
            )

            @b.Constraint()
            def NOAK_eq(c):
                return c.NOAK_factor == pyunits.convert(
                    (c.cum_num_units) ** -(c.learning_rate_exponent)
                )

        @b.Constraint()
        def total_TPC_eq(c):
            # TPC = BEC that needs Lang factor applied to get TPC + TPC that is already calculated
            # apply location factor and economy of numbers here too
            return c.total_TPC == (
                (
                    (sum(c.BEC_list) - sum(c.BEC_blocks_with_TPC_list)) * c.Lang_factor
                    + sum(c.TPC_blocks_with_TPC_list)
                )
                *
                # apply economy of numbers if enabled
                (c.NOAK_factor if c.config.has_economy_of_numbers else 1)
                # TODO * c.location_factor
            )

    def get_fixed_OM_costs(
        b,
        annual_fixed_operating_cost=None,
        annual_revenue=None,
        pure_product_output_rates=None,
        mixed_product_output_rates=None,
        sale_prices=None,
        additional_sales_price_dictionaries=None,
        production_rate=None,
    ):
        """
        Args:
            pure_product_output_rates: production rates of each pure product.
                Input should have the form {"name1": rate1, ...} where the name
                keys are strings and the rates are Vars, Params, or Expressions
                with Pyomo units of flow.
            mixed_product_output_rates: production rates of each product in the
                mixed basket. Input should have the form {"name1": rate1, ...}
                where the name keys are strings and the rates are Vars, Params,
                or Expressions with Pyomo units of flow.
            sale_prices: dictionary of products and prices to use in lieu of
                built-in sale prices. Input should have the form {"product1":
                price1, "product2": price2, ...} where the product keys match
                some or all keys in the `pure_product_output_rates` and/or
                `mixed_product_output_rates` arguments and the prices are Vars,
                Params, or Expressions with Pyomo units of currency per unit.
                Users should consult the built-in price dictionary in
                costing_dictionaries.py for available prices; note that the
                built-prices do not require input to this argument to be used
                as long as an entry in the `resources` argument matches a
                built-in price name.
            additional_sales_price_dictionaries: external dictionaries to append
            production_rate: Var, Param, or Expression for the production rate of power
                or saleable product generated by the plant. If passed, cost of production
                will be calculated and reported.
        """

        if (
            annual_revenue is None
        ):  # then we will calculate the annual revenue from product sales
            # let product rates be optional
            if pure_product_output_rates is None:
                pure_product_output_rates = {}
            if mixed_product_output_rates is None:
                mixed_product_output_rates = {}

            # check that required product arguments were passed
            if not isinstance(pure_product_output_rates, dict):
                raise TypeError("product_output_rates argument must be a dict")
            if not isinstance(mixed_product_output_rates, dict):
                raise TypeError("product_output_rates argument must be a dict")

            # dictionary of default sale prices
            # the currency units are of USD
            # Purity, purchase quantity, purchasing time, and location, all affect the cost.
            default_sale_prices = {}

            if additional_sales_price_dictionaries is not None:
                for dictionary in additional_sales_price_dictionaries:
                    for key in dictionary:
                        default_sale_prices[key] = dictionary[key]

            if sale_prices is None:
                sale_prices = {}

            # add entries from sale_prices to default_sale_prices
            if not isinstance(sale_prices, dict):
                raise TypeError(
                    "Dictionary of custom sale_prices must be a dict object."
                )
            else:
                for key in sale_prices:
                    default_sale_prices[key] = sale_prices[key]

            # raise error if the user included a product not in default_sale_prices
            if not set(pure_product_output_rates).issubset(default_sale_prices):
                raise AttributeError(
                    f"A pure product was included that does not contain a "
                    f"sale price. Sale prices exist for the following products: "
                    f"{list(default_sale_prices)}"
                )
            elif not set(mixed_product_output_rates).issubset(default_sale_prices):
                raise AttributeError(
                    f"A mixed product was included that does not contain a "
                    f"sale price. Sale prices exist for the following products: "
                    f"{list(default_sale_prices)}"
                )

        # sales revenue
        b.total_sales_revenue = Var(
            initialize=4,
            bounds=(0, None),
            doc="Total sales revenue",
            units=b.CE_index_units / pyunits.year,
        )

        if annual_revenue is not None:  # use user-supplied value

            @b.Constraint()
            def total_sales_revenue_eq(c):
                return c.total_sales_revenue == pyunits.convert(
                    annual_revenue,
                    c.CE_index_units / pyunits.year,
                )

        else:

            @b.Constraint()
            def total_sales_revenue_eq(c):
                return c.total_sales_revenue == pyunits.convert(
                    (
                        (
                            sum(
                                pyunits.convert(
                                    pure_product_output_rates[p]
                                    * default_sale_prices[p],
                                    to_units=b.CE_index_units / pyunits.h,
                                )
                                for p in pure_product_output_rates
                            )
                            if len(pure_product_output_rates) > 0
                            else 1e-12 * b.CE_index_units / pyunits.h
                        )
                        + (
                            c.mixed_product_sale_price_realization_factor
                            * sum(
                                pyunits.convert(
                                    mixed_product_output_rates[p]
                                    * default_sale_prices[p],
                                    to_units=b.CE_index_units / pyunits.h,
                                )
                                for p in mixed_product_output_rates
                            )
                            if len(mixed_product_output_rates) > 0
                            else 1e-12 * c.CE_index_units / pyunits.h
                        )
                    )
                    * c.capacity_factor,
                    c.CE_index_units / pyunits.year,
                )

        if (
            annual_fixed_operating_cost is None
        ):  # we will calculate the fixed operating cost

            # labor costs, common to all tech types
            b.annual_operating_labor_cost = Var(
                initialize=1,
                bounds=(0, None),
                doc="Annual operating labor cost",
                units=b.CE_index_units / pyunits.year,
            )

            b.annual_technical_labor_cost = Var(
                initialize=1,
                bounds=(0, None),
                doc="Annual technical labor cost",
                units=b.CE_index_units / pyunits.year,
            )

            b.annual_labor_cost = Var(
                initialize=1,
                bounds=(0, None),
                doc="Annual labor cost",
                units=b.CE_index_units / pyunits.year,
            )

            @b.Constraint()
            def annual_operating_labor_cost_eq(c):
                return (
                    c.annual_operating_labor_cost
                    == pyunits.convert(
                        (
                            sum(
                                c.operators_per_shift[i] * c.labor_rates[i]
                                for i in [
                                    "skilled",
                                    "unskilled",
                                    "supervisor",
                                    "maintenance",
                                    "operator",
                                ]
                            )
                            * (
                                1
                                + pyunits.convert(
                                    c.labor_burden, to_units=pyunits.dimensionless
                                )
                            )
                            * c.capacity_factor
                        ),
                        c.CE_index_units / pyunits.year,
                    )
                    + 1e-12 * b.CE_index_units / pyunits.year
                )  # so annual_operating_labor_cost has a value if no operators are selected

            @b.Constraint()
            def annual_technical_labor_cost_eq(c):
                return (
                    c.annual_technical_labor_cost
                    == pyunits.convert(
                        (
                            sum(
                                c.operators_per_shift[i] * c.labor_rates[i]
                                for i in ["technician", "engineer"]
                            )
                            * (
                                1
                                + pyunits.convert(
                                    c.labor_burden, to_units=pyunits.dimensionless
                                )
                            )
                            * c.capacity_factor
                        ),
                        c.CE_index_units / pyunits.year,
                    )
                    + 1e-12 * b.CE_index_units / pyunits.year
                )  # so annual_operating_labor_cost has a value if no operators are selected

            @b.Constraint()
            def annual_labor_cost_eq(c):
                return c.annual_labor_cost == pyunits.convert(
                    (c.annual_operating_labor_cost + c.annual_technical_labor_cost),
                    b.CE_index_units / pyunits.year,
                )

            if b.config.tech == 10:  # PrOMMiS REE UKy
                # UKy-specific fixed costs
                b.maintenance_and_material_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="Maintenance and material cost",
                    units=b.CE_index_units / pyunits.year,
                )

                b.quality_assurance_and_control_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="Quality assurance and control cost",
                    units=b.CE_index_units / pyunits.year,
                )

                b.sales_patenting_and_research_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="Sales, patenting and research cost",
                    units=b.CE_index_units / pyunits.year,
                )

                b.admin_and_support_labor_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="Admin and support labor cost",
                    units=b.CE_index_units / pyunits.year,
                )

                b.property_taxes_and_insurance_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="Property taxes and insurance cost",
                    units=b.CE_index_units / pyunits.year,
                )

                # TODO apply location factor to some of these components
                # maintenance cost is 2% of TPC
                @b.Constraint()
                def maintenance_and_material_cost_eq(c):
                    return (
                        c.maintenance_and_material_cost
                        == 0.02 * c.total_TPC / pyunits.year
                    )

                # quality assurance cost is 10% of operating labor
                @b.Constraint()
                def quality_assurance_and_control_cost_eq(c):
                    return (
                        c.quality_assurance_and_control_cost
                        == 0.10
                        * pyunits.convert(
                            (c.annual_operating_labor_cost),
                            b.CE_index_units / pyunits.year,
                        )
                    )

                # sales cost is 0.5% of total revenue
                @b.Constraint()
                def sales_patenting_and_research_cost_eq(c):
                    return (
                        c.sales_patenting_and_research_cost
                        == 0.005
                        * pyunits.convert(
                            (c.total_sales_revenue),
                            c.CE_index_units / pyunits.year,
                        )
                    )

                # admin cost is 20% of direct labor
                @b.Constraint()
                def admin_and_support_labor_cost_eq(c):
                    return c.admin_and_support_labor_cost == 0.20 * pyunits.convert(
                        (c.annual_operating_labor_cost),
                        c.CE_index_units / pyunits.year,
                    )

                # taxes are 1% of TPC
                @b.Constraint()
                def property_taxes_and_insurance_cost_eq(c):
                    return (
                        c.property_taxes_and_insurance_cost
                        == 0.01 * c.total_TPC / pyunits.year
                    )

            else:
                # power plant-specific fixed costs
                b.maintenance_labor_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="maintenance labor cost in millions of USD/yr",
                    units=b.CE_index_units / pyunits.year,
                )
                # maintenance material cost is technically a variable cost, but it
                # makes more sense to include with the fixed costs because it uses TPC
                b.maintenance_material_cost = Var(
                    initialize=2e-7,
                    bounds=(0, None),
                    doc="cost of maintenance materials in $MM/year",
                    units=b.CE_index_units / pyunits.year,
                )
                b.admin_and_support_labor_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="admin and support labor cost in millions of USD/yr",
                    units=b.CE_index_units / pyunits.year,
                )
                b.property_taxes_and_insurance_cost = Var(
                    initialize=1,
                    bounds=(0, None),
                    doc="property taxes and insurance cost in millions of USD/yr",
                    units=b.CE_index_units / pyunits.year,
                )

                # technology specific percentage of TPC
                @b.Constraint()
                def maintenance_labor_cost_eq(c):
                    return c.maintenance_labor_cost == (
                        c.total_TPC
                        * c.maintenance_labor_TPC_split
                        * pyunits.convert(
                            c.maintenance_labor_percent, to_units=pyunits.dimensionless
                        )
                        / pyunits.year
                    )

                # 25% of the sum of annual operating labor and maintenance labor
                @b.Constraint()
                def admin_and_support_labor_cost_eq(c):
                    return c.admin_and_support_labor_cost == (
                        0.25
                        * (c.annual_operating_labor_cost + c.maintenance_labor_cost)
                    )

                # 2% of TPC
                @b.Constraint()
                def property_taxes_and_insurance_cost_eq(c):
                    return (
                        c.property_taxes_and_insurance_cost
                        == 0.02 * c.total_TPC / pyunits.year
                    )

                # technology specific percentage of TPC
                @b.Constraint()
                def maintenance_material_cost_eq(c):
                    # check if an annual production rate was defined in units of power
                    try:
                        assert production_rate is not None
                        power_ratio = pyunits.convert(
                            production_rate / c.nameplate_capacity,
                            to_units=pyunits.dimensionless,
                        )

                        # adjust by the power ratio
                        return c.maintenance_material_cost == (
                            c.total_TPC
                            * c.maintenance_material_TPC_split
                            * pyunits.convert(
                                c.maintenance_material_percent,
                                to_units=pyunits.dimensionless,
                            )
                            / (c.capacity_factor * pyunits.year)
                            * power_ratio
                        )

                    except (AssertionError, InconsistentUnitsError):
                        # no power rate defined, don't adjust
                        return c.maintenance_material_cost == (
                            c.total_TPC
                            * c.maintenance_material_TPC_split
                            * pyunits.convert(
                                c.maintenance_material_percent,
                                to_units=pyunits.dimensionless,
                            )
                            / pyunits.year
                        )

        # support custom models including WaterTAP models, applies to all tech types

        # variable for user to assign other fixed costs to,
        # fixed to 0 by default
        b.other_fixed_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Other fixed costs",
            units=b.CE_index_units / pyunits.year,
        )
        b.other_fixed_costs.fix(1e-12)

        # variable for user to assign custom fixed costs to,
        # constraint sets to sum of list, which is 0 for empty list
        b.custom_fixed_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Custom fixed costs",
            units=b.CE_index_units / pyunits.year,
        )

        # sum of fixed operating costs of custom units
        @b.Constraint()
        def sum_custom_fixed_costs(c):
            if len(c.custom_fixed_costs_list) == 0:
                return c.custom_fixed_costs == 1e-12 * b.CE_index_units / pyunits.year
            else:
                return c.custom_fixed_costs == sum(b.custom_fixed_costs_list)

        # total fixed OM cost

        b.total_fixed_OM_cost = Var(
            initialize=4,
            bounds=(0, None),
            doc="Total fixed operating & maintenance costs",
            units=b.CE_index_units / pyunits.year,
        )

        if (
            annual_fixed_operating_cost is not None
        ):  # use user-supplied value for estimated OPEX components

            @b.Constraint()
            def total_fixed_OM_cost_eq(c):
                return c.total_fixed_OM_cost == (
                    pyunits.convert(
                        annual_fixed_operating_cost, b.CE_index_units / pyunits.year
                    )
                    + c.other_fixed_costs
                    + c.custom_fixed_costs
                )

        else:

            @b.Constraint()
            def total_fixed_OM_cost_eq(c):
                return c.total_fixed_OM_cost == (
                    c.annual_labor_cost
                    + (
                        c.maintenance_material_cost
                        if hasattr(c, "maintenance_material_cost")
                        else 0 * c.CE_index_units / pyunits.year
                    )
                    + (
                        c.quality_assurance_and_control_cost
                        if hasattr(c, "quality_assurance_and_control_cost")
                        else 0 * c.CE_index_units / pyunits.year
                    )
                    + c.admin_and_support_labor_cost
                    + (
                        c.sales_patenting_and_research_cost
                        if hasattr(c, "sales_patenting_and_research_cost")
                        else 0 * c.CE_index_units / pyunits.year
                    )
                    + c.property_taxes_and_insurance_cost
                    + c.other_fixed_costs
                    + c.custom_fixed_costs
                )

    def get_variable_OM_costs(
        b,
        resources=None,
        resource_prices=None,
        additional_resource_price_dictionaries=None,
    ):

        if resource_prices is None:
            resource_prices = {}

        if not hasattr(b.parent_block(), "time"):  # flowsheet is not dynamic
            b.parent_block().time = [0]
        if not hasattr(b.parent_block(), "time_units"):  # no time units set
            b.parent_block().time_units = pyunits.s

        # assert arguments are correct types
        if not isinstance(resources, dict):
            raise TypeError("resources argument must be a dictionary")
        if not isinstance(resource_prices, dict):
            raise TypeError("prices argument must be a dictionary")

        # dictionary of default prices
        default_resource_prices = load_default_resource_prices()

        # add entries from prices to default_prices
        for key in resource_prices:
            default_resource_prices[key] = resource_prices[key]

        if additional_resource_price_dictionaries is not None:
            for dictionary in additional_resource_price_dictionaries:
                for key in dictionary:
                    default_resource_prices[key] = dictionary[key]

        # raise error if the user included a resource not in default_prices
        if not set(resources).issubset(default_resource_prices):
            raise AttributeError(
                f"A resource was included that does not contain a "
                f"price. Prices exist for the following resources: "
                f"{list(default_resource_prices)}"
            )

        # make vars
        b.variable_operating_costs = Var(
            b.parent_block().time,
            resources,
            initialize=2e-7,
            doc="Variable operating costs",
            units=b.CE_index_units / pyunits.year,
        )

        b.other_variable_costs = Var(
            b.parent_block().time,
            initialize=0,
            bounds=(0, None),
            doc="A variable to include non-standard operating & maintenance costs",
            units=b.CE_index_units / pyunits.year,
        )

        # assume the user is not using this
        b.other_variable_costs.fix(1e-12)

        # variable for user to assign custom variable costs to,
        # constraint sets to sum of list, which is 0 for empty list
        b.custom_variable_costs = Var(
            initialize=0,
            bounds=(0, None),
            doc="Custom variable costs",
            units=b.CE_index_units / pyunits.year,
        )

        b.total_variable_OM_cost = Var(
            b.parent_block().time,
            initialize=4e-6,
            doc="Total variable operating and maintenance costs",
            units=b.CE_index_units / pyunits.year,
        )

        @b.Constraint(b.parent_block().time, resources.keys())
        def variable_cost_eq(c, t, r):
            return c.variable_operating_costs[t, r] == (
                pyunits.convert(
                    default_resource_prices[r] * resources[r][t] * c.capacity_factor,
                    to_units=b.CE_index_units / pyunits.year,
                )
            )

        # sum of variable operating costs of custom units
        @b.Constraint()
        def sum_custom_variable_costs(c):
            if len(c.custom_variable_costs_list) == 0:
                return (
                    c.custom_variable_costs == 1e-12 * b.CE_index_units / pyunits.year
                )
            else:
                return c.custom_variable_costs == sum(b.custom_variable_costs_list)

        if b.config.tech == 10:  # REE UKy has plant overhead
            if hasattr(b, "total_fixed_OM_cost"):

                # define overhead cost
                # plant overhead, 20% of direct costs - fixed OM, power, water, lease/land, chemicals, waste
                b.plant_overhead_cost = Var(
                    b.parent_block().time,
                    initialize=0,
                    doc="Plant overhead costs",
                    units=b.CE_index_units / pyunits.year,
                )

                @b.Constraint(b.parent_block().time)
                def plant_overhead_cost_eq(c, t):
                    return c.plant_overhead_cost[t] == 0.20 * (
                        c.total_fixed_OM_cost
                        + (
                            c.variable_operating_costs[0, "power"]
                            if (0, "power")
                            in b.variable_operating_costs.id_index_map().values()
                            else 0 * c.CE_index_units / pyunits.year
                        )
                        + sum(
                            c.variable_operating_costs[0, chemical]
                            for chemical in c.chemicals_list
                        )
                        + sum(
                            c.variable_operating_costs[0, waste]
                            for waste in c.waste_list
                        )
                        + (
                            c.land_cost
                            if c.land_cost_reoccurrence == "annual"
                            else 0 * c.CE_index_units / pyunits.year
                        )
                        + (
                            c.additional_chemicals_cost
                            if c.additional_chemicals_cost_reoccurrence == "annual"
                            else 0 * c.CE_index_units / pyunits.year
                        )
                        + (
                            c.additional_waste_cost
                            if c.additional_waste_cost_reoccurrence == "annual"
                            else 0 * c.CE_index_units / pyunits.year
                        )
                        + c.custom_variable_costs
                    )

        @b.Constraint(b.parent_block().time)
        def total_variable_cost_eq(c, t):
            return (
                c.total_variable_OM_cost[t]
                == sum(c.variable_operating_costs[t, r] for r in resources)
                + c.other_variable_costs[t]
                + (
                    c.plant_overhead_cost[t]
                    if hasattr(c, "plant_overhead_cost")
                    else 0 * c.CE_index_units / pyunits.year
                )
                + (
                    c.land_cost
                    if c.land_cost_reoccurrence == "annual"
                    else 0 * c.CE_index_units / pyunits.year
                )
                + (
                    c.additional_chemicals_cost
                    if c.additional_chemicals_cost_reoccurrence == "annual"
                    else 0 * c.CE_index_units / pyunits.year
                )
                + (
                    c.additional_waste_cost
                    if c.additional_waste_cost_reoccurrence == "annual"
                    else 0 * c.CE_index_units / pyunits.year
                )
                # for power plants, include maintenance material costs here if defined
                + (
                    c.maintenance_material_cost
                    if hasattr(c, "maintenance_material_cost")
                    else 0 * c.CE_index_units / pyunits.year
                )
            )

    def calculate_taxes(b):
        """
        Method to calculate income and production tax components.
        """

        # TODO add option and calculation of depreciation
        # build variables and expressions

        b.net_tax_owed = Var(
            initialize=0.40 * b.total_sales_revenue,
            doc="Net tax owed in millions USD",
            units=b.CE_index_units / pyunits.year,
        )
        b.income_tax = Var(
            initialize=0.26 * b.total_sales_revenue,
            doc="Income tax in millions USD",
            units=b.CE_index_units / pyunits.year,
        )
        b.additional_tax_credit = Var(
            initialize=0,
            doc="Additional tax credit",
            units=b.CE_index_units / pyunits.year,
        )
        b.additional_tax_credit.fix(1e-12)

        b.additional_tax_owed = Var(
            initialize=0,
            doc="Additional tax owed",
            units=b.CE_index_units / pyunits.year,
        )
        b.additional_tax_owed.fix(1e-12)

        b.min_net_tax_owed = Var(
            initialize=0,
            doc="Minimum net tax owed in millions USD",
            units=b.CE_index_units / pyunits.year,
        )
        b.min_net_tax_owed.fix(1e-12)

        b.royalty_charge = Expression(
            expr=pyunits.convert(
                b.royalty_charge_percentage_of_revenue,
                to_units=pyunits.dimensionless,
            )
            * b.total_sales_revenue
        )
        b.mineral_depletion_charge = Expression(
            expr=pyunits.convert(
                b.mineral_depletion_percentage,
                to_units=pyunits.dimensionless,
            )
            * (b.total_sales_revenue - b.royalty_charge)
        )

        if b.config.has_production_credit_phaseout:
            b.production_incentive_charge = Expression(
                expr=pyunits.convert(
                    b.production_incentive_percentage,
                    to_units=pyunits.dimensionless,
                )
                * b.phaseout_factor
                * (
                    (
                        b.total_variable_OM_cost[0]
                        if b.config.has_variable_OM
                        else 0 * b.CE_index_units / pyunits.year
                    )
                    + b.total_fixed_OM_cost
                )
            )
        else:
            b.production_incentive_charge = Expression(
                expr=pyunits.convert(
                    b.production_incentive_percentage,
                    to_units=pyunits.dimensionless,
                )
                * (
                    (
                        b.total_variable_OM_cost[0]
                        if b.config.has_variable_OM
                        else 0 * b.CE_index_units / pyunits.year
                    )
                    + b.total_fixed_OM_cost
                )
            )

        # build constraints

        @b.Constraint()
        def income_tax_eq(c):
            return c.income_tax == pyunits.convert(
                b.income_tax_percentage, to_units=pyunits.dimensionless
            ) * (
                b.total_sales_revenue
                - (
                    (
                        b.total_variable_OM_cost[0]
                        if b.config.has_variable_OM
                        else 0 * b.CE_index_units / pyunits.year
                    )
                    + b.total_fixed_OM_cost
                    + b.annualized_cost / pyunits.year
                )
            )

        @b.Constraint()
        def net_tax_owed_eq(c):
            tax_units = (
                c.CE_index_units / pyunits.year
            )  # for unit consistency since EPS has no units
            return (
                c.net_tax_owed
                == smooth_max(
                    b.min_net_tax_owed / tax_units,
                    (
                        (b.income_tax + b.royalty_charge + b.additional_tax_owed)
                        - (
                            b.mineral_depletion_charge
                            + b.production_incentive_charge
                            + b.additional_tax_credit
                        )
                    )
                    / tax_units,
                    eps=EPS,
                )
                * tax_units
            )

    def initialize(b):

        # loop through costing block looking for equipment costs to initialize
        # b is the flowsheet-level costing block
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in b._registered_unit_costing:  # pylint: disable=protected-access
                if hasattr(o, "library"):
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

                if hasattr(o.costing, "bare_erected_cost"):
                    for key in o.costing.bare_erected_cost:
                        calculate_variable_from_constraint(
                            o.bare_erected_cost[key],
                            o.bare_erected_cost_eq[key],
                        )
                if hasattr(o.costing, "total_plant_cost"):
                    for key in o.costing.total_plant_cost:
                        calculate_variable_from_constraint(
                            o.total_plant_cost[key],
                            o.total_plant_cost_eq[key],
                        )
                if hasattr(o.costing, "capital_cost"):
                    calculate_variable_from_constraint(
                        o.capital_cost,
                        o.capital_cost_constraint,
                    )
                if hasattr(o.costing, "fixed_operating_cost"):
                    calculate_variable_from_constraint(
                        o.fixed_operating_cost,
                        o.fixed_operating_cost_constraint,
                    )
                if hasattr(o.costing, "variable_operating_cost"):
                    calculate_variable_from_constraint(
                        o.variable_operating_cost,
                        o.variable_operating_cost_constraint,
                    )

                else:  # shouldn't happen if user is using one of the pre-built costing packages
                    raise BurntToast(
                        "This should not occur. If this message displays, please contact the IDAES developers."
                    )

            # if fixed OM costs exist, initialize them here
            # variables and constraints for REE are different
            if hasattr(b, "total_fixed_OM_cost"):

                calculate_variable_from_constraint(
                    b.admin_and_support_labor_cost, b.admin_and_support_labor_cost_eq
                )

                calculate_variable_from_constraint(
                    b.property_taxes_and_insurance_cost,
                    b.property_taxes_and_insurance_cost_eq,
                )

                if hasattr(b, "annual_operating_labor_cost_eq"):
                    calculate_variable_from_constraint(
                        b.annual_operating_labor_cost, b.annual_operating_labor_cost_eq
                    )

                if hasattr(b, "annual_technical_labor_cost_eq"):
                    calculate_variable_from_constraint(
                        b.annual_technical_labor_cost, b.annual_technical_labor_cost_eq
                    )

                calculate_variable_from_constraint(
                    b.annual_labor_cost, b.annual_labor_cost_eq
                )

                calculate_variable_from_constraint(
                    b.total_sales_revenue, b.total_sales_revenue_eq
                )

                calculate_variable_from_constraint(
                    b.custom_fixed_costs, b.sum_custom_fixed_costs
                )

                if b.config.tech == 10:

                    calculate_variable_from_constraint(
                        b.maintenance_and_material_cost,
                        b.maintenance_and_material_cost_eq,
                    )

                    calculate_variable_from_constraint(
                        b.quality_assurance_and_control_cost,
                        b.quality_assurance_and_control_cost_eq,
                    )

                    calculate_variable_from_constraint(
                        b.sales_patenting_and_research_cost,
                        b.sales_patenting_and_research_cost_eq,
                    )

                else:

                    calculate_variable_from_constraint(
                        b.maintenance_labor_cost, b.maintenance_labor_cost_eq
                    )

                if hasattr(b, "maintenance_material_cost"):

                    calculate_variable_from_constraint(
                        b.maintenance_material_cost, b.maintenance_material_cost_eq
                    )

                calculate_variable_from_constraint(
                    b.total_fixed_OM_cost, b.total_fixed_OM_cost_eq
                )

            # if variable OM costs exist, initialize them here
            if hasattr(b, "variable_operating_costs"):
                for i in b.variable_operating_costs:
                    if hasattr(b, "variable_cost_eq"):
                        calculate_variable_from_constraint(
                            b.variable_operating_costs[i],
                            b.variable_cost_eq[i],
                        )

                for i in b.total_variable_OM_cost:
                    calculate_variable_from_constraint(
                        b.total_variable_OM_cost[i],
                        b.total_variable_cost_eq[i],
                    )

                if hasattr(b, "plant_overhead_cost"):

                    for i in b.plant_overhead_cost:
                        calculate_variable_from_constraint(
                            b.plant_overhead_cost[i],
                            b.plant_overhead_cost_eq[i],
                        )

    def calculate_net_present_value(b, debt_expression=None):
        """
        Equations for cash flow expressions derived from the textbook
        Engineering Economy: Applying Theory to Practice, 3rd Ed. by Ted. G. Eschenbach.

        The net present value (NPV) is a representative measure of the "current
        day" value of a chemical plant over the total lifetime, including all
        cash flows.

        This method supports capital expenditure, loan repayment, inflation,
        and taxes if previously calculated. The NPV formulation assumes that
        negative cash flows consist of capital and operating costs scaled to a
        constant present value. The general NPV formula with 100% of capital
        expenditure upfront at the start of the operating period and no capital
        or operating growth rate is

        NPV = [(REVENUE - OPEX - TAXES) * P/A(r, N)] - CAPEX

        where P/A(r, N) is the series present worth factor; this factor scales
        a future cost to its present value from a known discount rate r and
        project lifetime N based on annuity growth over time. This factor is
        calculated as

        P/A(r, N) = [ 1 - (1+r)**(-N) ] / r

        where r is expressed as a decimal and N is expressed in years. In the
        NPV expression above, REVENUE is the constant annual revenue, OPEX is
        the constant annual operating cost, CAPEX is the total capital cost,
        and TAXES are the total net tax liability including income tax,
        royalties, depletion rights, and incentives.

        Operating costs and revenues are adjusted based on predicted annuity
        growth to obtain the present value. These expressions are implemented
        if there is no capital expenditure period or additional growth rate.

        Often, uniform series of cash flows do not start in the first period
        (t=1), but in some later period, T, after a known delay. In this case
        the delay is accounted for using

        PV_year1_cashflow = Cash_Flow_Value * P/A(r, N)

        PV_yearT_cashflow = Cash_Flow_Value * [P/A(r, N) - P/A(r, T-1)]

        ----------------------------------------------------------------------

        The general NPV formulation allows capital costs, operating costs,
        and revenues to escalate over time via geometric gradient growth, e.g.
        a constant proportional growth rate expressed as an escalation or
        inflation percentage. The general formulation is given by

        NPV = PV_Revenue - PV_Operating_Cost - PV_Taxes
              - PV_Capital_Cost - PV_Loan_Interest

        For costs escalating at a constant rate for the project lifetime, the
        series present worth factor is modified to account for escalation,
        yielding a modified formula

        P/A(r, g, N) = ( 1 - [ (1+g)**(N) ] * [(1+r)**(-N)] ) / (r - g)

        where r is the discount rate expressed as a decimal, N is the project
        lifetime, and g is the escalation rate (e.g. inflation) expressed as a
        decimal.

        The general formulation considers a capital escalation period followed
        by the operating period, and capital costs may be distributed across
        the capital escalation period rather than fully paid for upfront. For
        example, if the capital expenditures are distributed across a 3-year
        capital escalation period, the PV from the capital costs are given as

        PV_Capital_Cost = Y1_% * CAPEX * [P/A(r, gCap, 1) - P/A(r, gCap, 0)]
                          + Y2_% * CAPEX * [P/A(r, gCap, 2) - P/A(r, gCap, 1)]
                          + Y3_% * CAPEX * [P/A(r, gCap, 3) - P/A(r, gCap, 2)]

        where Y1_%, Y2_%, and Y3_% are the percentages of capital expenditure
        in each year expressed as decimals, CAPEX is the total capital cost from
        equipment purchasing, gCap is the capital escalation growth rate
        expressed as a decimal. The capital costs spent in each year are handled
        separately to properly account for the value growth over time. Loan
        repayment and interest owed are calculated as

        PV_Loan_Interest_Owed = Debt * [P/A(r, 0, Nloan) / P/A(iLoan_%, 0, Nloan) - 1]

        where Debt is the loan principal (typically a percentage of the CAPEX), Nloan
        is the loan repayment period, and iLoan_% is the capital equipment loan interest
        rate expressed as a decimal.

        Revenue, operating costs and taxes based on revenue escalate with
        standard inflation. Notably, these cash flows occur after any capital
        expenditure period, meaning that the annuity growth must be offset by
        the length of the capital expenditure period. This yields the expressions

        PV_Revenue = REVENUE * [ P/A(r, gRev, NOp+NCap) - P/A(r, gRev, NCap) ]

        PV_Operating_Cost = OPEX * [ P/A(r, gOp, NOp+NCap) - P/A(r, gOp, NCap) ]

        PV_Taxes = TAXES * [ P/A(r, 0, NOp+NCap) - P/A(r, gRev, NCap) ]

        where REVENUE is the annual revenue, OPEX is the annual operating cost,
        TAXES is the annual net tax liability, gRev is the inflation or growth
        rate of revenue year-on-year expressed as a decimal, gOp is the inflation
        or growth rate of operating costs year-on-year expressed as a decimal,
        NOp is the length of the operating period or plant lifetime, and NCap is
        the length of the capital expenditure period. The expressions above take
        the annuity growth during the entire analysis period (NOp+NCap) and subtract
        the capital expenditure period (NCap) as there is no operation or production
        during that time. Note that the taxes do not include inflation growth, as future
        tax rates are difficult to predict and do not rise proportionally with inflation.

        Args:
            b: costing block to retrieve total plant cost (capital), total plant
                operating cost, and total revenue from, and add net present
                value (NPV) calculations to
        """

        # if OPEX methods were not called, assume the user purposefully did not define fixed OPEX or revenue
        # and they are zero
        if not b.config.has_fixed_OM:

            b.total_fixed_OM_cost = Var(
                initialize=0,
                bounds=(0, None),
                doc="Total fixed operating & maintenance costs",
                units=b.CE_index_units / pyunits.year,
            )
            b.total_fixed_OM_cost.fix(0)

            b.total_sales_revenue = Var(
                initialize=0,
                bounds=(0, None),
                doc="Total sales revenue",
                units=b.CE_index_units / pyunits.year,
            )
            b.total_sales_revenue.fix(0)

        if not b.config.has_variable_OM:

            b.total_variable_OM_cost = Var(
                b.parent_block().time,
                initialize=0,
                doc="Total variable operating and maintenance costs",
                units=b.CE_index_units / pyunits.year,
            )
            b.total_variable_OM_cost.fix(0)

        # input verification

        try:
            b.capex = Expression(
                expr=(
                    b.total_TPC
                    + b.other_plant_costs
                    + (
                        b.land_cost
                        if b.land_cost_reoccurrence == "one_time"
                        else 0 * b.CE_index_units
                    )
                )
            )
            b.opex = Expression(
                expr=(
                    (
                        (
                            b.total_fixed_OM_cost
                            if b.config.has_fixed_OM
                            else 0 * b.CE_index_units / pyunits.year
                        )
                        + (
                            b.total_variable_OM_cost[0]
                            if b.config.has_variable_OM
                            else 0 * b.CE_index_units / pyunits.year
                        )
                        + (
                            b.land_cost
                            if b.land_cost_reoccurrence == "annual"
                            else 0 * b.CE_index_units / pyunits.year
                        )
                    )
                )
            )

        except AttributeError:
            raise AttributeError(
                "Expected FlowsheetCostingBlockData object "
                "with attributes total_BEC,"
                "total_fixed_OM_cost, total_variable_OM_cost, "
                "other_plant_costs, land_cost, and total_sales_revenue. "
                "Please confirm that b is a FlowsheetCostingBlockData object "
                "and that all expected attributes exist."
            )

        # TODO add sum to 100 check on capital expenditure percentages in case user changes from default values
        # if b.config.has_capital_expenditure_period:

        # @b.BuildCheck()
        # def capital_expenditure_percentages_sum_to_100(b, tol=EPS):
        #     total = sum(b.capital_expenditure_percentages[i] for i in b.capital_expenditure_percentages)
        #     return abs(total - 100) < tol

        # b.capital_expenditure_percentages_sum_to_100 = capital_expenditure_percentages_sum_to_100.__get__(b)

        # check debt expression, if defined
        if debt_expression is not None:

            expr_units = pyunits.get_units(debt_expression)

            if expr_units is None:
                raise ValueError("Expression debt_expression has no units defined.")

            try:
                pyunits.convert(expr_units, to_units=b.CE_index_units)
            except UnitsError:
                raise ValueError(
                    f"Expression debt_expression with units {expr_units} are not compatible with "
                    f"{b.CE_index_units}. The  expression must be compatible with cost units to be "
                    f"included in the NPV calculation."
                )

        # build variables

        b.pv_capital_cost = Var(
            initialize=-b.capex,
            bounds=(None, 0),
            doc="Present value of total lifetime capital costs; negative cash flow",
            units=b.CE_index_units,
        )

        b.loan_debt = Var(
            initialize=b.capex,
            bounds=(0, 1e4),
            doc="total debt from loans in millions of USD",
            units=b.CE_index_units,
        )

        b.pv_loan_interest = Var(
            initialize=-b.capex,
            bounds=(-1e4, 1e4),
            doc="present value of total lifetime loan interest in millions of USD; normally a negative cash flow, but can be positive depending on the discount and interest rates",
            units=b.CE_index_units,
        )

        b.pv_operating_cost = Var(
            initialize=-b.opex * b.plant_lifetime,
            bounds=(None, 0),
            doc="Present value of total lifetime operating costs; negative cash flow",
            units=b.CE_index_units,
        )

        b.pv_revenue = Var(
            initialize=b.total_sales_revenue * b.plant_lifetime,
            bounds=(0, None),
            doc="Present value of total lifetime sales revenue; positive cash flow",
            units=b.CE_index_units,
        )

        # TODO add option and calculation of depreciation
        if b.config.has_taxes_and_credits:
            b.pv_taxes = Var(
                initialize=-b.net_tax_owed * b.plant_lifetime,
                bounds=(None, None),
                doc="Present value of total lifetime tax owed; typically negative cash flow, positive in the case of an effective negative tax rate",
                units=b.CE_index_units,
            )

            if b.config.has_production_credit_phaseout:

                b.pv_production_incentive = Var(
                    initialize=b.production_incentive_charge * b.plant_lifetime,
                    bounds=(0, None),
                    doc="Present value of total lifetime production incentive",
                    units=b.CE_index_units,
                )

                current_year = int(value(b.current_year))
                b.plant_end_year = Param(
                    initialize=current_year + int(value(b.plant_lifetime)),
                    mutable=True,
                    doc="Year that plant ceases operation; used to determine phaseout of production incentive",
                )

                full_credit_years = int(
                    value(
                        smooth_max(
                            0,
                            min([int(y) for y in b.phaseout_fractions.keys()])
                            - current_year,
                            eps=EPS,
                        )
                    )
                )

                b.production_incentive_charge_percent_list = []

                for i in range(full_credit_years):
                    b.production_incentive_charge_percent_list.append(
                        value(
                            pyunits.convert(
                                b.production_incentive_percentage,
                                to_units=pyunits.dimensionless,
                            )
                        )
                    )

                for year, fraction in b.phaseout_fractions.items():
                    if value(b.plant_end_year) >= int(year):
                        b.production_incentive_charge_percent_list.append(
                            fraction
                            * value(
                                pyunits.convert(
                                    b.production_incentive_percentage,
                                    to_units=pyunits.dimensionless,
                                )
                            )
                        )

                if value(b.plant_end_year) > max(
                    [int(y) for y in b.phaseout_fractions.keys()]
                ):
                    zero_credit_years = value(b.plant_end_year) - max(
                        [int(y) for y in b.phaseout_fractions.keys()]
                    )

                    for i in range(value(zero_credit_years)):
                        b.production_incentive_charge_percent_list.append(0.0)

        b.npv = Var(
            initialize=(-b.capex + (b.total_sales_revenue - b.opex) * b.plant_lifetime),
            bounds=(None, None),
            doc="Present value of plant over entire capital and operation lifetime",
            units=b.CE_index_units,
        )

        # define series present worth factor as a method so it can be called

        def series_present_worth_factor(r, g, N):
            """
            Returns expression for series present worth factor where r is the discount rate
            expressed as a decimal, N is the project lifetime, and g is the escalation rate
            (e.g. inflation) expressed as a decimal.
            """
            if (
                abs(value(r - g)) < EPS / 1e2
            ):  # r and g are basically the same value, use simpler expression to avoid numerical issues
                return N / (1 + r)
            else:
                return (1 - ((1 + g) ** (N)) * ((1 + r) ** (-N))) / (r - g)

        # build constraints

        if b.config.has_capital_expenditure_period:

            @b.Constraint()
            def pv_capital_cost_constraint(c):
                # percentage of CAPEX is basis for each capital expenditure year
                # since the expenditure series restarts in each year, we need to split
                # the terms for each year out and subtract off the delayed years
                # PV_Capital_Cost = - (
                # %year1 * CAPEX * (P/A_year1 - P/A_year0)     change from year 1 only
                # + %year2 * CAPEX * (P/A_year2 - P/A_year1)   change from year 2 only
                # + %year3 * CAPEX * (P/A_year3 - P/A_year2)   change from year 2 only
                # + ...)
                # P/A_year0 = 0, which places each CAPEX expenditure at the end of each period

                return c.pv_capital_cost == -pyunits.convert(
                    sum(
                        pyunits.convert(
                            c.capital_expenditure_percentages[idx] * pyunits.percent,
                            to_units=pyunits.dimensionless,
                        )
                        * c.capex
                        * (  # P/A_year(i) - P/A_year(i-1))
                            series_present_worth_factor(
                                pyunits.convert(
                                    c.discount_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                pyunits.convert(
                                    c.capital_escalation_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                idx + 1,
                            )
                            - series_present_worth_factor(
                                pyunits.convert(
                                    c.discount_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                pyunits.convert(
                                    c.capital_escalation_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                idx,
                            )
                        )
                        for idx in range(len(c.capital_expenditure_percentages))
                    ),
                    to_units=c.CE_index_units,
                )

        else:

            @b.Constraint()
            def pv_capital_cost_constraint(c):
                # no expenditure period, so cash flow occurs at t=0 (project year)
                # PV_Capital_Cost = - CAPEX

                return c.pv_capital_cost == -pyunits.convert(
                    c.capex,
                    to_units=c.CE_index_units,
                )

        if debt_expression is None:

            @b.Constraint()
            def loan_debt_constraint(c):
                # Debt  = %debt_charge_of_CAPEX * CAPEX

                return c.loan_debt == pyunits.convert(
                    pyunits.convert(
                        c.debt_percentage_of_capex, to_units=pyunits.dimensionless
                    )
                    * c.capex,
                    to_units=c.CE_index_units,
                )

        else:

            @b.Constraint()
            def loan_debt_constraint(c):

                return c.loan_debt == pyunits.convert(
                    debt_expression,
                    to_units=c.CE_index_units,
                )

        @b.Constraint()
        def pv_loan_interest_constraint(c):
            # PV_Loan_Interest_Owed = Debt * [P/A(r, 0, loan_length) / P/A(%interest, 0, loan_length) - 1]
            # when r > %interest, this is a negative value; the loan amount borrowed devalues faster than
            # the loan amount repaid
            # when r < %interest, this is a positive value; the loan amount borrowed devalues slower than
            # the loan amount repaid

            return c.pv_loan_interest == pyunits.convert(
                c.loan_debt
                * (
                    series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        0,  # loan value does not grow over time
                        c.capital_loan_repayment_period / pyunits.year,
                    )
                    / series_present_worth_factor(
                        pyunits.convert(
                            c.capital_loan_interest_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        0,  # loan payments do not grow over time
                        c.capital_loan_repayment_period / pyunits.year,
                    )
                    - 1
                ),
                to_units=c.CE_index_units,
            )

        @b.Constraint()
        def pv_operating_cost_constraint(c):
            # OPEX starts after the capital expenditure period, so we need to account for a delay
            # PV_Operating_Cost = - OPEX * [ P/A(r, g, OPEX_end_year) - P/A(r, g, CAPEX_end_year) ]

            return c.pv_operating_cost == -pyunits.convert(
                c.opex
                * pyunits.year
                * (
                    series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.operating_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        c.plant_lifetime / pyunits.year
                        + len(c.capital_expenditure_percentages),
                    )
                    - series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.operating_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        len(c.capital_expenditure_percentages),
                    )
                ),
                to_units=c.CE_index_units,
            )

        @b.Constraint()
        def pv_revenue_constraint(c):
            # Revenue starts after the capital expenditure period, so we need to account for a delay
            # PV_Revenue = REVENUE * [ P/A(r, g, Revenue_end_year) - P/A(r, g, CAPEX_end_year) ]

            return c.pv_revenue == pyunits.convert(
                c.total_sales_revenue
                * pyunits.year
                * (
                    series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.revenue_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        c.plant_lifetime / pyunits.year
                        + len(c.capital_expenditure_percentages),
                    )
                    - series_present_worth_factor(
                        pyunits.convert(
                            c.discount_percentage, to_units=pyunits.dimensionless
                        ),
                        pyunits.convert(
                            c.revenue_inflation_percentage,
                            to_units=pyunits.dimensionless,
                        ),
                        len(c.capital_expenditure_percentages),
                    )
                ),
                to_units=c.CE_index_units,
            )

        if b.config.has_taxes_and_credits:

            @b.Constraint()
            def pv_taxes_constraint(c):
                # Taxes start after the capital expenditure period, so we need to account for a delay
                # PV_taxes = net_tax_owed * [ P/A(r, 0, Operating_end_year) - P/A(r, 0, CAPEX_end_year) ]

                return c.pv_taxes == -pyunits.convert(
                    (
                        (c.net_tax_owed + c.production_incentive_charge)
                        if b.config.has_production_credit_phaseout
                        else c.net_tax_owed
                    )
                    * pyunits.year
                    * (
                        series_present_worth_factor(
                            pyunits.convert(
                                c.discount_percentage, to_units=pyunits.dimensionless
                            ),
                            pyunits.convert(
                                0,
                                to_units=pyunits.dimensionless,
                            ),
                            c.plant_lifetime / pyunits.year
                            + len(c.capital_expenditure_percentages),
                        )
                        - series_present_worth_factor(
                            pyunits.convert(
                                c.discount_percentage, to_units=pyunits.dimensionless
                            ),
                            pyunits.convert(
                                0,
                                to_units=pyunits.dimensionless,
                            ),
                            len(c.capital_expenditure_percentages),
                        )
                    ),
                    to_units=c.CE_index_units,
                )

        if b.config.has_production_credit_phaseout:

            # The present value of the production incentive is:
            #
            #   PV_prod_incentive = sum over k=0,1,...,M-1 of [ p_k * C_opex * (SPWF(r, g, k+1) - SPWF(r, g, k)) ]

            # where:
            #   p_k    = production_incentive_charge_percent_list[k]  (incentive fraction in year k)
            #   C_opex = annual operating cost (base year dollars)
            #   M      = len(production_incentive_charge_percent_list)
            #   SPWF(r, g, N) = series present worth factor with discount rate r, inflation rate g, and N years

            @b.Constraint()
            def pv_production_incentive_constraint(c):

                return c.pv_production_incentive == pyunits.convert(
                    sum(
                        b.production_incentive_charge_percent_list[idx]
                        * c.opex
                        * pyunits.year
                        * (
                            series_present_worth_factor(
                                pyunits.convert(
                                    c.discount_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                pyunits.convert(
                                    c.operating_inflation_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                idx + 1,
                            )
                            - series_present_worth_factor(
                                pyunits.convert(
                                    c.discount_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                pyunits.convert(
                                    c.operating_inflation_percentage,
                                    to_units=pyunits.dimensionless,
                                ),
                                idx,
                            )
                        )
                        for idx in range(
                            len(b.production_incentive_charge_percent_list)
                        )
                    ),
                    to_units=c.CE_index_units,
                )

        @b.Constraint()
        def npv_constraint(c):

            return c.npv == pyunits.convert(
                c.pv_revenue
                + c.pv_capital_cost
                + c.pv_loan_interest
                + c.pv_operating_cost
                + (
                    c.pv_taxes
                    if b.config.has_taxes_and_credits
                    else 0 * c.CE_index_units
                )
                + (
                    c.pv_production_incentive
                    if b.config.has_production_credit_phaseout
                    else 0 * c.CE_index_units
                ),
                to_units=c.CE_index_units,
            )

    def display_total_plant_costs(b):
        print("-----Total Plant Costs (MUSD)-----")
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name
                for block in b._registered_unit_costing  # pylint: disable=protected-access
            ] and hasattr(o, "total_plant_cost"):
                print(
                    "%s: %.2f"
                    % (
                        value(o.name),
                        value(
                            sum(o.total_plant_cost[key] for key in o.total_plant_cost)
                        ),
                    )
                )

    def display_bare_erected_costs(b):
        print("-----Bare Erected Costs (MUSD)-----")
        for o in b.parent_block().component_objects(descend_into=True):
            # look for costing blocks
            if o.name in [
                block.name
                for block in b._registered_unit_costing  # pylint: disable=protected-access
            ] and hasattr(o, "bare_erected_cost"):
                print(
                    "%s: %.5f"
                    % (
                        value(o.name),
                        value(
                            sum(o.bare_erected_cost[key] for key in o.bare_erected_cost)
                        ),
                    )
                )
