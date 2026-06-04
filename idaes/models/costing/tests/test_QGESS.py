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

__author__ = "Costing Team (B. Paul, L. Deng, A. Fritz, A. Ojo, A. Dasgupta, A. Noring, A. Deshpande, D. Caballero, and M. Zamarripa)"
__version__ = "1.0.0"

import pytest

import pyomo.environ as pyo
from pyomo.environ import units as pyunits

from idaes.core import FlowsheetBlock
from idaes.core.util.exceptions import ConfigurationError

from idaes.models.costing.QGESS import (
    QGESSCosting,
)

from idaes.models_extra.power_generation.costing.power_plant_costing_dictionaries import (
    load_fixed_OM_data,
    )

class TestQGESSConfigParameters(object):
    @pytest.mark.component
    def test_no_config_set(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError,
            match="Must set a technology type. Valid options include: \n"
            "1. Supercritical PC, air-fired, with and without CO2 capture, Illinois No. 6 coal \n"
            "2. Subcritical PC, air-fired, with and without CO2 capture, Illinois No. 6 coal \n"
            "3. Two-stage, slurry-feed, oxygen-blown gasifier with and without CO2 capture, Illinois No. 6 coal \n"
            "4. Single-stage, slurry-feed, oxygen-blown gasifier with and without CO2 capture, Illinois No. 6 coal \n"
            "5. Single-stage, dry-feed, oxygen-blown, up-flow gasifier with and without CO2 capture, Illinois No. 6 coal \n"
            "6. Natural gas, air-fired, with and without CO2 capture \n"
            "7. Advanced Ultrasupercritical PC \n"
            "8. Polymer Layers accounts \n"
            "9. Sensors & Controls accounts \n"
            "10. University of Kentucky Fire Clay Seam \\(Hazard No. 4\\) Rejects \n",
        ):
            m.fs.costing = QGESSCosting()
    
    
    @pytest.mark.component
    def test_use_defaults(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = QGESSCosting(tech=1)
    
        expected = {
            "Value": {
                "base_currency": 1,
                "CE_index_units": 1,
                "base_period": 1,
                "fs.costing.capacity_factor": 0.85,
                "fs.costing.Lang_factor": 1,
                "fs.costing.pct_TPC": 0.2020,
                "fs.costing.tasc_toc_factor": 1.093,
                "fs.costing.fixed_charge_factor": 0.0707,
            },
            "Units": {
                "base_currency": "USD_2021",
                "CE_index_units": "MUSD_2021",
                "base_period": "a",
                "fs.costing.capacity_factor": "dimensionless",
                "fs.costing.Lang_factor": "dimensionless",
                "fs.costing.pct_TPC": "dimensionless",
                "fs.costing.tasc_toc_factor": "dimensionless",
                "fs.costing.fixed_charge_factor": "dimensionless",
            },
        }
    
        assert m.fs.costing.param_dir["Value"].keys() == expected["Value"].keys()
    
        for k in m.fs.costing.param_dir["Value"]:
            assert m.fs.costing.param_dir["Value"][k] == pytest.approx(
                expected["Value"][k], rel=1e-4
            ), f"Value mismatch for key {k}"
    
        for k in m.fs.costing.param_dir["Units"]:
            assert (
                str(m.fs.costing.param_dir["Units"][k]) == expected["Units"][k]
            ), f"Units mismatch for key {k}"
    
    
    @pytest.mark.component
    def test_invalid_CE_index_year(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            AttributeError,
            match="CE_index_year notayear is not a valid currency base option. "
            "Valid CE index options include CE500, CE394, years from 1990 to 2023, or user-defined "
            "units such as 2019_Sep and UKy_2019.",
        ):
            m.fs.costing = QGESSCosting(tech=1, CE_index_year="notayear")
    
    
    @pytest.mark.component
    def test_invalid_year_for_phaseout(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError,
            match="CE_index_year must contain a valid 4\\‑digit year at the start or end "
            "\\(e.g. 2019, 2019_Sep, UKy_2019\\)",
        ):
            m.fs.costing = QGESSCosting(
                tech=1, CE_index_year="CE500", has_production_credit_phaseout=True
            )
    
    
    @pytest.mark.parametrize(
        "CE_index_year",
        [
            "2021",
            "2023",
            "2019_Sep",
            "UKy_2019",
            "2008_Nov",
            "2018_Dec",
        ],
    )
    @pytest.mark.component
    def test_valid_years_for_phaseout(self, CE_index_year):
    
        # mock UKy_2019 for the purpose of this test
        pyunits.load_definitions_from_strings(
            [
                "USD_UKy_2019 = 500/609.495 * USD_CE500",
            ]
        )
    
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = QGESSCosting(
            tech=1, CE_index_year=CE_index_year, has_production_credit_phaseout=True
        )
    
        if CE_index_year == "UKy_2019":
            assert m.fs.costing.current_year.value == int(CE_index_year[-4:])
        else:
            assert m.fs.costing.current_year.value == int(CE_index_year[:4])
    
    
    @pytest.mark.component
    def test_taxes_no_fixed_OM(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ConfigurationError,
            match="Cannot set has_taxes_and_credits to True if has_fixed_OM is False. Calculations of income tax, royalties, "
            "and mineral depletion require sales revenue, and production incentive is only applicable to plants "
            "producing saleable product streams. Sale revenue and associated administrative costs are included in "
            "the fixed O&M calculations.",
        ):
            m.fs.costing = QGESSCosting(
                tech=1, CE_index_year="CE500", has_taxes_and_credits=True
            )
    
    
    @pytest.mark.component
    def test_no_phaseout_years_set(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError,
            match="Must set phaseout_fractions as dict of fractions indexed by integer years.",
        ):
            m.fs.costing = QGESSCosting(
                tech=1,
                has_production_credit_phaseout=True,
                has_taxes_and_credits=True,
                has_fixed_OM=True,
            )
    
    
    @pytest.mark.component
    def test_phaseout_years_not_ascending(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError, match="Years for phaseout_fractions must be in ascending order."
        ):
            m.fs.costing = QGESSCosting(
                tech=1,
                has_fixed_OM=True,
                has_taxes_and_credits=True,
                has_production_credit_phaseout=True,
                phaseout_fractions=dict(
                    zip(
                        ["2049", "2029", "2030", "2031", "2032", "2033"],
                        [90, 75, 60, 45, 30, 15],
                    )
                ),
            )
    
    
    @pytest.mark.component
    def test_phaseout_fractions_not_descending(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError,
            match="Fractions for phaseout_fractions must be in descending order.",
        ):
            m.fs.costing = QGESSCosting(
                tech=1,
                has_fixed_OM=True,
                has_taxes_and_credits=True,
                has_production_credit_phaseout=True,
                phaseout_fractions=dict(
                    zip(
                        ["2028", "2029", "2030", "2031", "2032", "2033"],
                        [70, 75, 60, 45, 30, 15],
                    )
                ),
            )
    
    
    @pytest.mark.parametrize(
        "CE_index_year",
        [
            "2021",
            "2023",
            "2027",
            "2028",
            "2029",
            "2030",
            "2031",
            "2032",
            "2033",
            "2034",
            "2035",
        ],
    )
    @pytest.mark.component
    def test_phaseout_factors(self, CE_index_year):
    
        # mock future years for the purpose of this test
        for year in range(27, 36):
            pyunits.load_definitions_from_strings(
                [
                    "USD_20" + str(year) + " = 500/500 * USD_CE500",
                ]
            )
    
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = QGESSCosting(
            tech=1,
            CE_index_year=CE_index_year,
            has_fixed_OM=True,
            has_taxes_and_credits=True,
            has_production_credit_phaseout=True,
            phaseout_fractions=dict(
                zip(
                    ["2028", "2029", "2030", "2031", "2032", "2033"],
                    [90, 75, 60, 45, 30, 15],
                )
            ),
        )
    
        expected = dict(
            zip(
                [
                    "2021",
                    "2023",
                    "2027",
                    "2028",
                    "2029",
                    "2030",
                    "2031",
                    "2032",
                    "2033",
                    "2034",
                    "2035",
                ],
                [
                    1.00,
                    1.00,
                    1.00,
                    0.90,
                    0.75,
                    0.60,
                    0.45,
                    0.30,
                    0.15,
                    0.00,
                    0.00,
                ],
            )
        )
    
        assert (
            pyo.value(m.fs.costing.phaseout_factor) == expected[CE_index_year]
        ), f"Expected {expected[CE_index_year]}, got {pyo.value(m.fs.costing.phaseout_factor)}"
    
    
    @pytest.mark.component
    def test_capital_expenditure_percentages_not_set(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError,
            match="Must set capital_expenditure_percentages as list of integer values on \\[0, 100\\].",
        ):
            m.fs.costing = QGESSCosting(
                tech=1,
                has_net_present_value=True,
                has_capital_expenditure_period=True,
            )
    
    
    @pytest.mark.component
    def test_capital_expenditure_percentages_not_sum_to_100(self):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        with pytest.raises(
            ValueError,
            match="Argument capital_expenditure_percentages has a sum of 105.0. List must sum to 100 percent.",
        ):
            m.fs.costing = QGESSCosting(
                tech=1,
                has_net_present_value=True,
                has_capital_expenditure_period=True,
                capital_expenditure_percentages=[10, 60, 35],
            )
    
    
    @pytest.mark.parametrize(
        "tech",
        [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
        ],
    )
    @pytest.mark.component
    def test_all_config_set(self, tech):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = QGESSCosting(
            Lang_factor=2,
            has_fixed_OM=True,
            has_variable_OM=True,
            has_taxes_and_credits=True,
            has_production_credit_phaseout=True,
            phaseout_fractions=dict(zip([2031, 2032, 2033], [75, 50, 25])),
            has_net_present_value=True,
            has_capital_expenditure_period=True,
            capital_expenditure_percentages=[10, 60, 30],
            has_economy_of_numbers=True,
            CE_index_year="2023",
            tech=tech,
        )

        # tech-dependent expected values to check against
        _, _, maintenance_percentages = load_fixed_OM_data()

        if tech == 10:
            operators_per_shift = [2, 5, 2, 3, 1, 2, 0]
            labor_burden = 25
            capacity_factor = 0.92
        else:
            operators_per_shift = [0, 0, 0, 0, 0, 0, 6]
            labor_burden = 30
            capacity_factor = 0.85

        if tech == 10:
            tasc_toc_factor = 1.144
            fixed_charge_factor = 0.1002
        elif tech == 9:
            pct_TPC = (15 + 2.7 + 0.5 + 2) / 100
            tasc_toc_factor = 1.070
            fixed_charge_factor = 0.0690
        elif tech == 8:
            pct_TPC = (15 + 2.7 + 0.5 + 2) / 100
            tasc_toc_factor = 1.047
            fixed_charge_factor = 0.0664
        else:
            pct_TPC = 20.2 / 100
            tasc_toc_factor = 1.093
            fixed_charge_factor = 0.0707

        if tech == 10:
            expected = {
                "Value": {
                    "base_currency": 1,
                    "CE_index_units": 1,
                    "base_period": 1,
                    "fs.costing.current_year": 2023,
                    "fs.costing.capacity_factor": capacity_factor,
                    "fs.costing.Lang_factor": 2,
                    "fs.costing.labor_rates[skilled]": 27.9,
                    "fs.costing.labor_rates[unskilled]": 23.26,
                    "fs.costing.labor_rates[supervisor]": 30.290,
                    "fs.costing.labor_rates[maintenance]": 24.06,
                    "fs.costing.labor_rates[technician]": 23.43,
                    "fs.costing.labor_rates[engineer]": 46.82,
                    "fs.costing.labor_rates[operator]": 38.5,
                    "fs.costing.labor_burden": labor_burden,
                    "fs.costing.operators_per_shift[skilled]": operators_per_shift[0],
                    "fs.costing.operators_per_shift[unskilled]": operators_per_shift[1],
                    "fs.costing.operators_per_shift[supervisor]": operators_per_shift[2],
                    "fs.costing.operators_per_shift[maintenance]": operators_per_shift[3],
                    "fs.costing.operators_per_shift[technician]": operators_per_shift[4],
                    "fs.costing.operators_per_shift[engineer]": operators_per_shift[5],
                    "fs.costing.operators_per_shift[operator]": operators_per_shift[6],
                    "fs.costing.mixed_product_sale_price_realization_factor": 0.65,
                    "fs.costing.income_tax_percentage": 26,
                    "fs.costing.mineral_depletion_percentage": 14,
                    "fs.costing.production_incentive_percentage": 10,
                    "fs.costing.royalty_charge_percentage_of_revenue": 6.5,
                    "fs.costing.phaseout_fractions[2031]": 75,
                    "fs.costing.phaseout_fractions[2032]": 50,
                    "fs.costing.phaseout_fractions[2033]": 25,
                    "fs.costing.phaseout_factor": 1,
                    "fs.costing.discount_percentage": 10,
                    "fs.costing.plant_lifetime": 20,
                    "fs.costing.capital_expenditure_percentages[0]": 10,
                    "fs.costing.capital_expenditure_percentages[1]": 60,
                    "fs.costing.capital_expenditure_percentages[2]": 30,
                    "fs.costing.capital_escalation_percentage": 3.6,
                    "fs.costing.capital_loan_interest_percentage": 6,
                    "fs.costing.capital_loan_repayment_period": 10,
                    "fs.costing.debt_percentage_of_capex": 50,
                    "fs.costing.operating_inflation_percentage": 3,
                    "fs.costing.revenue_inflation_percentage": 3,
                    "fs.costing.cum_num_units": 5,
                    "fs.costing.learning_rate": 0.04,
                    "fs.costing.tasc_toc_factor": tasc_toc_factor,
                    "fs.costing.fixed_charge_factor": fixed_charge_factor,
                },
                "Units": {
                    "base_currency": "USD_2023",
                    "CE_index_units": "MUSD_2023",
                    "base_period": "a",
                    "fs.costing.current_year": "dimensionless",
                    "fs.costing.capacity_factor": "dimensionless",
                    "fs.costing.Lang_factor": "dimensionless",
                    "fs.costing.labor_rates[skilled]": "USD_2023/h",
                    "fs.costing.labor_rates[unskilled]": "USD_2023/h",
                    "fs.costing.labor_rates[supervisor]": "USD_2023/h",
                    "fs.costing.labor_rates[maintenance]": "USD_2023/h",
                    "fs.costing.labor_rates[technician]": "USD_2023/h",
                    "fs.costing.labor_rates[engineer]": "USD_2023/h",
                    "fs.costing.labor_rates[operator]": "USD_2023/h",
                    "fs.costing.labor_burden": "%",
                    "fs.costing.operators_per_shift[skilled]": "dimensionless",
                    "fs.costing.operators_per_shift[unskilled]": "dimensionless",
                    "fs.costing.operators_per_shift[supervisor]": "dimensionless",
                    "fs.costing.operators_per_shift[maintenance]": "dimensionless",
                    "fs.costing.operators_per_shift[technician]": "dimensionless",
                    "fs.costing.operators_per_shift[engineer]": "dimensionless",
                    "fs.costing.operators_per_shift[operator]": "dimensionless",
                    "fs.costing.mixed_product_sale_price_realization_factor": "dimensionless",
                    "fs.costing.income_tax_percentage": "%",
                    "fs.costing.mineral_depletion_percentage": "%",
                    "fs.costing.production_incentive_percentage": "%",
                    "fs.costing.royalty_charge_percentage_of_revenue": "%",
                    "fs.costing.phaseout_fractions[2031]": "dimensionless",
                    "fs.costing.phaseout_fractions[2032]": "dimensionless",
                    "fs.costing.phaseout_fractions[2033]": "dimensionless",
                    "fs.costing.phaseout_factor": "dimensionless",
                    "fs.costing.discount_percentage": "%",
                    "fs.costing.plant_lifetime": "a",
                    "fs.costing.capital_expenditure_percentages[0]": "dimensionless",
                    "fs.costing.capital_expenditure_percentages[1]": "dimensionless",
                    "fs.costing.capital_expenditure_percentages[2]": "dimensionless",
                    "fs.costing.capital_escalation_percentage": "%",
                    "fs.costing.capital_loan_interest_percentage": "%",
                    "fs.costing.capital_loan_repayment_period": "a",
                    "fs.costing.debt_percentage_of_capex": "%",
                    "fs.costing.operating_inflation_percentage": "%",
                    "fs.costing.revenue_inflation_percentage": "%",
                    "fs.costing.cum_num_units": "dimensionless",
                    "fs.costing.learning_rate": "dimensionless",
                    "fs.costing.tasc_toc_factor": "dimensionless",
                    "fs.costing.fixed_charge_factor": "dimensionless",
                },
            }
        else:
            expected = {
                "Value": {
                    "base_currency": 1,
                    "CE_index_units": 1,
                    "base_period": 1,
                    "fs.costing.current_year": 2023,
                    "fs.costing.capacity_factor": capacity_factor,
                    "fs.costing.Lang_factor": 2,
                    "fs.costing.nameplate_capacity": 650,
                    "fs.costing.maintenance_labor_TPC_split": maintenance_percentages[tech][0],
                    "fs.costing.maintenance_labor_percent":maintenance_percentages[tech][1],
                    "fs.costing.maintenance_material_TPC_split": 1 - maintenance_percentages[tech][0],
                    "fs.costing.maintenance_material_percent": maintenance_percentages[tech][2],
                    "fs.costing.labor_rates[skilled]": 27.9,
                    "fs.costing.labor_rates[unskilled]": 23.26,
                    "fs.costing.labor_rates[supervisor]": 30.290,
                    "fs.costing.labor_rates[maintenance]": 24.06,
                    "fs.costing.labor_rates[technician]": 23.43,
                    "fs.costing.labor_rates[engineer]": 46.82,
                    "fs.costing.labor_rates[operator]": 38.5,
                    "fs.costing.labor_burden": labor_burden,
                    "fs.costing.operators_per_shift[skilled]": operators_per_shift[0],
                    "fs.costing.operators_per_shift[unskilled]": operators_per_shift[1],
                    "fs.costing.operators_per_shift[supervisor]": operators_per_shift[2],
                    "fs.costing.operators_per_shift[maintenance]": operators_per_shift[3],
                    "fs.costing.operators_per_shift[technician]": operators_per_shift[4],
                    "fs.costing.operators_per_shift[engineer]": operators_per_shift[5],
                    "fs.costing.operators_per_shift[operator]": operators_per_shift[6],
                    "fs.costing.mixed_product_sale_price_realization_factor": 0.65,
                    "fs.costing.income_tax_percentage": 26,
                    "fs.costing.mineral_depletion_percentage": 14,
                    "fs.costing.production_incentive_percentage": 10,
                    "fs.costing.royalty_charge_percentage_of_revenue": 6.5,
                    "fs.costing.phaseout_fractions[2031]": 75,
                    "fs.costing.phaseout_fractions[2032]": 50,
                    "fs.costing.phaseout_fractions[2033]": 25,
                    "fs.costing.phaseout_factor": 1,
                    "fs.costing.discount_percentage": 10,
                    "fs.costing.plant_lifetime": 20,
                    "fs.costing.capital_expenditure_percentages[0]": 10,
                    "fs.costing.capital_expenditure_percentages[1]": 60,
                    "fs.costing.capital_expenditure_percentages[2]": 30,
                    "fs.costing.capital_escalation_percentage": 3.6,
                    "fs.costing.capital_loan_interest_percentage": 6,
                    "fs.costing.capital_loan_repayment_period": 10,
                    "fs.costing.debt_percentage_of_capex": 50,
                    "fs.costing.operating_inflation_percentage": 3,
                    "fs.costing.revenue_inflation_percentage": 3,
                    "fs.costing.cum_num_units": 5,
                    "fs.costing.learning_rate": 0.04,
                    "fs.costing.pct_TPC": pct_TPC,
                    "fs.costing.tasc_toc_factor": tasc_toc_factor,
                    "fs.costing.fixed_charge_factor": fixed_charge_factor,
                },
                "Units": {
                    "base_currency": "USD_2023",
                    "CE_index_units": "MUSD_2023",
                    "base_period": "a",
                    "fs.costing.current_year": "dimensionless",
                    "fs.costing.capacity_factor": "dimensionless",
                    "fs.costing.Lang_factor": "dimensionless",
                    "fs.costing.nameplate_capacity": "MW",
                    "fs.costing.maintenance_labor_TPC_split": "dimensionless",
                    "fs.costing.maintenance_labor_percent": "dimensionless",
                    "fs.costing.maintenance_material_TPC_split": "dimensionless",
                    "fs.costing.maintenance_material_percent": "dimensionless",
                    "fs.costing.labor_rates[skilled]": "USD_2023/h",
                    "fs.costing.labor_rates[unskilled]": "USD_2023/h",
                    "fs.costing.labor_rates[supervisor]": "USD_2023/h",
                    "fs.costing.labor_rates[maintenance]": "USD_2023/h",
                    "fs.costing.labor_rates[technician]": "USD_2023/h",
                    "fs.costing.labor_rates[engineer]": "USD_2023/h",
                    "fs.costing.labor_rates[operator]": "USD_2023/h",
                    "fs.costing.labor_burden": "%",
                    "fs.costing.operators_per_shift[skilled]": "dimensionless",
                    "fs.costing.operators_per_shift[unskilled]": "dimensionless",
                    "fs.costing.operators_per_shift[supervisor]": "dimensionless",
                    "fs.costing.operators_per_shift[maintenance]": "dimensionless",
                    "fs.costing.operators_per_shift[technician]": "dimensionless",
                    "fs.costing.operators_per_shift[engineer]": "dimensionless",
                    "fs.costing.operators_per_shift[operator]": "dimensionless",
                    "fs.costing.mixed_product_sale_price_realization_factor": "dimensionless",
                    "fs.costing.income_tax_percentage": "%",
                    "fs.costing.mineral_depletion_percentage": "%",
                    "fs.costing.production_incentive_percentage": "%",
                    "fs.costing.royalty_charge_percentage_of_revenue": "%",
                    "fs.costing.phaseout_fractions[2031]": "dimensionless",
                    "fs.costing.phaseout_fractions[2032]": "dimensionless",
                    "fs.costing.phaseout_fractions[2033]": "dimensionless",
                    "fs.costing.phaseout_factor": "dimensionless",
                    "fs.costing.discount_percentage": "%",
                    "fs.costing.plant_lifetime": "a",
                    "fs.costing.capital_expenditure_percentages[0]": "dimensionless",
                    "fs.costing.capital_expenditure_percentages[1]": "dimensionless",
                    "fs.costing.capital_expenditure_percentages[2]": "dimensionless",
                    "fs.costing.capital_escalation_percentage": "%",
                    "fs.costing.capital_loan_interest_percentage": "%",
                    "fs.costing.capital_loan_repayment_period": "a",
                    "fs.costing.debt_percentage_of_capex": "%",
                    "fs.costing.operating_inflation_percentage": "%",
                    "fs.costing.revenue_inflation_percentage": "%",
                    "fs.costing.cum_num_units": "dimensionless",
                    "fs.costing.learning_rate": "dimensionless",
                    "fs.costing.pct_TPC": "dimensionless",
                    "fs.costing.tasc_toc_factor": "dimensionless",
                    "fs.costing.fixed_charge_factor": "dimensionless",
                },
            }
    
        assert set(m.fs.costing.param_dir["Value"].keys()) == set(
            expected["Value"].keys()
        ), f"Expected {expected['Value'].keys()}, got {m.fs.costing.param_dir['Value'].keys()}"
    
        for k in m.fs.costing.param_dir["Value"]:
            if expected["Value"][k] == 0:
                assert m.fs.costing.param_dir["Value"][k] == pytest.approx(
                    expected["Value"][k], abs=1e-4
                ), f"Value mismatch for key {k}"
            else:
                assert m.fs.costing.param_dir["Value"][k] == pytest.approx(
                    expected["Value"][k], rel=1e-4
                ), f"Value mismatch for key {k}"
    
        for k in m.fs.costing.param_dir["Units"]:
            assert (
                str(m.fs.costing.param_dir["Units"][k]) == expected["Units"][k]
            ), f"Units mismatch for key {k}"

class TestQGESSBuildProcessCosts(object):

    @pytest.mark.parametrize(
        "tech",
        [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
        ],
    )
    @pytest.mark.component
    def test_default_arguments(self, tech):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = QGESSCosting(
            tech=tech,
        )

        m.fs.costing.build_process_costs()

    @pytest.mark.parametrize(
        "tech",
        [
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
        ],
    )
    @pytest.mark.component
    def test_all_arguments(self, tech):
        # Create a Concrete Model as the top level object
        m = pyo.ConcreteModel()
    
        # Add a flowsheet object to the model
        m.fs = FlowsheetBlock(dynamic=False)
        m.fs.costing = QGESSCosting(
            Lang_factor=2,
            has_fixed_OM=True,
            has_variable_OM=True,
            has_taxes_and_credits=True,
            has_production_credit_phaseout=True,
            phaseout_fractions=dict(zip([2031, 2032, 2033], [75, 50, 25])),
            has_net_present_value=True,
            has_capital_expenditure_period=True,
            capital_expenditure_percentages=[10, 60, 30],
            has_economy_of_numbers=True,
            CE_index_year="2023",
            tech=tech,
        )

        m.fs.coal = pyo.Var(m.fs.time, initialize=1, units=pyunits.tonne/pyunits.h)
        m.fs.coal.fix()
        
        m.fs.natural_gas = pyo.Var(m.fs.time, initialize=1, units=pyunits.MBtu/pyunits.s)
        m.fs.natural_gas.fix()
        
        m.fs.water = pyo.Var(m.fs.time, initialize=1, units=pyunits.gallon/pyunits.s)
        m.fs.water.fix()
        
        m.fs.chemicals = pyo.Var(m.fs.time, initialize=1, units=pyunits.kg/pyunits.s)
        m.fs.chemicals.fix()
        
        m.fs.nonharzardous_waste_disposal = pyo.Var(m.fs.time, initialize=1, units=pyunits.kg/pyunits.s)
        m.fs.nonharzardous_waste_disposal.fix()
        
        m.fs.pure_product = pyo.Var(m.fs.time, initialize=1, units=pyunits.kg/pyunits.s)
        m.fs.pure_product.fix()
        
        m.fs.mixed_product = pyo.Var(m.fs.time, initialize=0.1, units=pyunits.kg/pyunits.s)
        m.fs.mixed_product.fix()
        

        m.fs.costing.build_process_costs(
            # optional arguments that directly fix cost variables and bypass calculations
            total_purchase_cost=100 * pyunits.MUSD_2021,
            annual_fixed_operating_cost=10 * pyunits.MUSD_2021/pyunits.year,
            annual_revenue=15 * pyunits.MUSD_2021/pyunits.year,
            debt_expression=50 * pyunits.MUSD_2021,  # for NPV method
            # optional arguments for additional calculations and reporting
            feedstock_rate=m.fs.coal,
            production_rate=m.fs.pure_product[0] + m.fs.mixed_product[0],  # this could be power, total REE, water, CO2
            # required arguments for fixed_OM calculations
            pure_product_output_rates={"main_product": m.fs.pure_product[0]},
            mixed_product_output_rates={"main_product": 0.5 * m.fs.mixed_product[0], "byproduct": 0.5 * m.fs.mixed_product[0],},
            sale_prices={"main_product": 1 * pyunits.USD_2021/pyunits.kg, "byproduct": 0.25 * pyunits.USD_2021/pyunits.kg,},
            # required arguments for variable_OM calculations
            resources={"coal": m.fs.coal, "natural_gas": m.fs.natural_gas, "water": m.fs.water, "chemicals":  m.fs.chemicals, "nonharzardous_waste_disposal": m.fs.nonharzardous_waste_disposal},  # for annual OPEX costs
            resource_prices={"chemicals": 1 * pyunits.USD_2021/pyunits.kg},
            # optional arguments related to overnight costs
            land_cost=0.30 * pyunits.USD_2021/pyunits.kg * 1 * pyunits.kg/pyunits.s,  # for startup and/or annual leasing costs
            # TODO allow more than one fuel or feedstock resource
            fuel=["natural_gas",],  # extra inventory required as part of overnight costs
            feedstock=["coal",],  # extra inventory required as part of overnight costs
            waste=["nonharzardous_waste_disposal",],  # extra storage required as part of overnight costs
            additional_waste_cost=1 * pyunits.USD_2021/pyunits.d,
            chemicals=["chemicals",],  # extra inventory required as part of overnight costs
            additional_chemicals_cost=1 * pyunits.USD_2021/pyunits.d,
            chemicals_inventory=["chemicals",],
            transport_cost=10 * pyunits.USD_2021/pyunits.kg * 0.2 * pyunits.kg/pyunits.s,
            )

    
