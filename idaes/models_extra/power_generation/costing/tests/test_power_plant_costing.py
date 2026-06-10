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

__author__ = (
    "Costing Team (A. Noring, A. Deshpande, B. Paul, D. Caballero, and M. Zamarripa)"
)
__version__ = "1.0.0"

import pyomo.environ as pyo
import pytest
from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_diagnostics import DiagnosticsToolbox
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.models_extra.power_generation.costing.power_plant_costing import (
    PowerPlantCosting,
    PowerPlantCostingData,
)
from pyomo.environ import units as pyunits


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
        10,  # this is not a power plant tech, but test for coverage and to ensure it works as expected
    ],
)
@pytest.mark.component
def test_PP_costing_with_all_methods(tech):
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = PowerPlantCosting(
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
        CE_index_year="2021",
        tech=tech,
    )

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

    ###########################################################################
    #  Create costing constraints                                             #
    ###########################################################################

    # coal flow rate
    # accounts 1.x and 2.x are coal handling, preparation and feed
    # accounts 4.x are for boiler BOP and foundations
    coal_accounts = ["1.1", "1.2", "1.3", "1.4", "2.1", "2.2", "4.11", "4.15", "4.16"]
    m.fs.boiler = UnitModelBlock()
    m.fs.boiler.coal_mass_flow = pyo.Var(
        initialize=7238.95, units=pyunits.ton / pyunits.day
    )
    m.fs.boiler.coal_mass_flow.fix()
    m.fs.boiler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=PowerPlantCostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": coal_accounts,
            "scaled_param": m.fs.boiler.coal_mass_flow,
            "tech": 2,
            "ccs": "A",
        },
    )

    # total fuel feed
    # accounts 3.x are for start up systems and miscellaneous plant equipment
    # accounts 7.x are for ductwork and stack foundations
    fuel_accounts = ["3.6", "3.9", "7.3", "7.5"]
    m.fs.fuel_feed = UnitModelBlock()
    m.fs.fuel_feed.total_fuel_feed = pyo.Var(
        initialize=603246, units=pyunits.lb / pyunits.hr
    )
    m.fs.fuel_feed.total_fuel_feed.fix()
    m.fs.fuel_feed.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=PowerPlantCostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": fuel_accounts,
            "scaled_param": m.fs.fuel_feed.total_fuel_feed,
            "tech": 2,
            "ccs": "A",
        },
    )

    # HP BFW flow rate
    # accounts 3.x are for feedwater systems
    # account 4.9 is for the boiler
    # account 8.4 is steam piping
    BFW_accounts = ["3.1", "3.3", "3.5", "4.9", "8.4"]
    m.fs.bfp = UnitModelBlock()
    m.fs.bfp.BFW_mass_flow = pyo.Var(initialize=5316158, units=pyunits.lb / pyunits.hr)
    m.fs.bfp.BFW_mass_flow.fix()
    m.fs.bfp.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=PowerPlantCostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": BFW_accounts,
            "scaled_param": m.fs.bfp.BFW_mass_flow,
            "tech": 2,
            "ccs": "A",
        },
    )

    # Steam turbine power
    # accounts 8.x are for the steam turbine and its foundations
    power_accounts = ["8.1"]
    m.fs.turb = UnitModelBlock()
    m.fs.turb.power = pyo.Var(initialize=769600, units=pyunits.kW)
    m.fs.turb.power.fix()
    m.fs.turb.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=PowerPlantCostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": power_accounts,
            "scaled_param": m.fs.turb.power,
            "tech": 2,
            "ccs": "A",
        },
    )

    # Condernser duty
    cond_accounts = ["8.3"]
    m.fs.condenser = UnitModelBlock()
    m.fs.condenser.duty_MMBtu = pyo.Var(
        initialize=2016, units=pyunits.MBtu / pyunits.hr
    )
    m.fs.condenser.duty_MMBtu.fix()
    m.fs.condenser.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=PowerPlantCostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": cond_accounts,
            "scaled_param": m.fs.condenser.duty_MMBtu,
            "tech": 2,
            "ccs": "A",
        },
    )

    # Circulating water flow rate
    # accounts 9.x are for circulating water systems
    # account 14.5 is for the pumphouse
    circ_accounts = ["9.2", "9.3", "9.4", "9.6", "9.7", "14.5"]
    m.fs.circulating_water = UnitModelBlock()
    m.fs.circulating_water.vol_flow = pyo.Var(
        initialize=463371, units=pyunits.gal / pyunits.min
    )
    m.fs.circulating_water.vol_flow.fix()
    m.fs.circulating_water.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=PowerPlantCostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": circ_accounts,
            "scaled_param": m.fs.circulating_water.vol_flow,
            "tech": 2,
            "ccs": "A",
        },
    )

    # Ash flow rate
    # accounts are for ash storage and handling
    ash_accounts = ["10.6", "10.7", "10.9"]
    m.fs.ash_handling = UnitModelBlock()
    m.fs.ash_handling.ash_mass_flow = pyo.Var(
        initialize=66903, units=pyunits.lb / pyunits.hr
    )
    m.fs.ash_handling.ash_mass_flow.fix()
    m.fs.ash_handling.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=PowerPlantCostingData.get_equipment_costing,
        costing_method_arguments={
            "cost_accounts": ash_accounts,
            "scaled_param": m.fs.ash_handling.ash_mass_flow,
            "tech": 2,
            "ccs": "A",
        },
    )

    # Fixed and Variable Costs:
    # build variable costs components
    m.fs.net_power = pyo.Var(m.fs.time, initialize=650, units=pyunits.MW)
    m.fs.net_power.fix()

    m.fs.NG_rate = pyo.Var(m.fs.time, initialize=1.2, units=pyunits.MBtu / pyunits.s)
    m.fs.NG_rate.fix()

    m.fs.solvent_rate = pyo.Var(
        m.fs.time, initialize=40, units=pyunits.ton / pyunits.day
    )
    m.fs.solvent_rate.fix()

    resources = ["natural_gas", "solvent"]
    rates = [m.fs.NG_rate, m.fs.solvent_rate]
    prices = {"solvent": 500 * pyunits.USD_2018 / pyunits.ton}

    m.fs.costing.build_process_costs(
        production_rate=m.fs.net_power[0],
        resources=dict(zip(resources, rates)),
        resource_prices=prices,
        fuel=[
            "natural_gas",
        ],
        feedstock=[
            "natural_gas",
        ],
        feedstock_rate=m.fs.NG_rate[0],
        chemicals=[
            "solvent",
        ],
        chemicals_inventory=[
            "solvent",
        ],
    )

    # add initialize
    PowerPlantCostingData.initialize(m.fs.costing)

    # check diagnostics
    dt = DiagnosticsToolbox(m)
    dt.report_structural_issues()
    dt.assert_no_structural_warnings()

    solver = get_solver()
    results = solver.solve(m, tee=True)
    pyo.assert_optimal_termination(results)
    dt.report_numerical_issues()
    dt.display_variables_at_or_outside_bounds()
    dt.assert_no_numerical_warnings()
