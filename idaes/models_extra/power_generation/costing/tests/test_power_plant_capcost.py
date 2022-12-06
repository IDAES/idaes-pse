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

__author__ = "Costing Team (A. Noring, B. Paul, D. Caballero, and M. Zamarripa)"
__version__ = "1.0.0"

import pytest

import pyomo.environ as pyo
from pyomo.core.base.constraint import ScalarConstraint, IndexedConstraint
from pyomo.core.base.expression import ScalarExpression
from pyomo.core.base.var import IndexedVar
from pyomo.environ import units as pyunits
from pyomo.core.base.units_container import UnitsError
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)


@pytest.mark.component
def test_PP_costing_exp():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

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
    m.fs.boiler.coal_mass_flow = pyo.Expression(
        expr=7238.95 * pyunits.ton / pyunits.day
    )
    m.fs.boiler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": coal_accounts,
            "scaled_param": m.fs.boiler.coal_mass_flow,
            "tech": 2,
            "ccs": "A",
        },
    )


@pytest.mark.component
def test_PP_costing_err1():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

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
    m.fs.boiler.coal_mass_flow = pyo.Expression(
        expr=7238.95 * pyunits.ton / pyunits.day + 1 * pyunits.gal
    )
    with pytest.raises(UnitsError):
        m.fs.boiler.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": coal_accounts,
                "scaled_param": m.fs.boiler.coal_mass_flow,
                "tech": 2,
                "ccs": "A",
            },
        )


@pytest.mark.component
def test_PP_costing_err2():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

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
    m.fs.boiler.coal_mass_flow = pyo.Expression(expr=7238.95 * pyunits.ton)
    with pytest.raises(UnitsError):
        m.fs.boiler.costing = UnitModelCostingBlock(
            flowsheet_costing_block=m.fs.costing,
            costing_method=QGESSCostingData.get_PP_costing,
            costing_method_arguments={
                "cost_accounts": coal_accounts,
                "scaled_param": m.fs.boiler.coal_mass_flow,
                "tech": 2,
                "ccs": "A",
            },
        )


@pytest.mark.component
def test_PP_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

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
        costing_method=QGESSCostingData.get_PP_costing,
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
        costing_method=QGESSCostingData.get_PP_costing,
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
        costing_method=QGESSCostingData.get_PP_costing,
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
        costing_method=QGESSCostingData.get_PP_costing,
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
        costing_method=QGESSCostingData.get_PP_costing,
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
        costing_method=QGESSCostingData.get_PP_costing,
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
        costing_method=QGESSCostingData.get_PP_costing,
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
        net_power=m.fs.net_power,
        fixed_OM=True,
        variable_OM=False,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
    )

    # add initialize
    QGESSCostingData.costing_initialization(m.fs.costing)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal

    #  all numbers come from the NETL excel file:
    # "201.001.001_BBR4 COE Spreadsheet_Rev0U_20190919_njk.xlsm"
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.1"]), abs=1e-1)
        == 2306 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.2"]), abs=1e-1)
        == 6385 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.3"]), abs=1e-1)
        == 59527 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.4"]), abs=1e-1)
        == 8086 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["2.1"]), abs=1e-1)
        == 4073 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["2.2"]), abs=1e-1)
        == 13976 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["4.11"]), abs=1e-1)
        == 3751 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["4.15"]), abs=1e-1)
        == 197 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["4.16"]), abs=1e-1)
        == 1014 / 1e3
    )

    assert (
        pytest.approx(
            pyo.value(m.fs.fuel_feed.costing.total_plant_cost["3.6"]), abs=1e-1
        )
        == 4864 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.fuel_feed.costing.total_plant_cost["3.9"]), abs=1e-1
        )
        == 522 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.fuel_feed.costing.total_plant_cost["7.3"]), abs=1e-1
        )
        == 1710 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.fuel_feed.costing.total_plant_cost["7.5"]), abs=1e-1
        )
        == 647 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.bfp.costing.total_plant_cost["3.1"]), abs=1e-1)
        == 19233 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.bfp.costing.total_plant_cost["3.3"]), abs=1e-1)
        == 6897 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.bfp.costing.total_plant_cost["3.5"]), abs=1e-1)
        == 2366 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.bfp.costing.total_plant_cost["4.9"]), abs=1e-1)
        == 572550 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.bfp.costing.total_plant_cost["8.4"]), abs=1e-1)
        == 81916 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.turb.costing.total_plant_cost["8.1"]), abs=1e-1)
        == 110166 / 1e3
    )

    assert (
        pytest.approx(
            pyo.value(m.fs.condenser.costing.total_plant_cost["8.3"]), abs=1e-1
        )
        == 21223 / 1e3
    )

    assert (
        pytest.approx(
            pyo.value(m.fs.circulating_water.costing.total_plant_cost["9.2"]), abs=1e-1
        )
        == 4133 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.circulating_water.costing.total_plant_cost["9.3"]), abs=1e-1
        )
        == 25518 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.circulating_water.costing.total_plant_cost["9.4"]), abs=1e-1
        )
        == 19859 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.circulating_water.costing.total_plant_cost["9.6"]), abs=1e-1
        )
        == 2870 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.circulating_water.costing.total_plant_cost["9.7"]), abs=1e-1
        )
        == 2690 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.circulating_water.costing.total_plant_cost["14.5"]), abs=1e-1
        )
        == 464 / 1e3
    )

    assert (
        pytest.approx(
            pyo.value(m.fs.ash_handling.costing.total_plant_cost["10.6"]), abs=1e-1
        )
        == 6429 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.ash_handling.costing.total_plant_cost["10.7"]), abs=1e-1
        )
        == 10725 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.ash_handling.costing.total_plant_cost["10.9"]), abs=1e-1
        )
        == 2564 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.costing.total_TPC), abs=1e-1) == 996662 / 1e3
    )  # 993753 / 1e3

    return m


@pytest.mark.component
def test_build_process_costs_emptymodel():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

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

    assert not hasattr(m.fs.costing, "total_TPC")
    assert not hasattr(m.fs.costing, "total_TPC_eq")

    m.fs.costing.build_process_costs(
        net_power=m.fs.net_power,
        fixed_OM=False,
        variable_OM=False,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
    )

    assert hasattr(m.fs.costing, "total_TPC")
    assert type(m.fs.costing.total_TPC) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "total_TPC_eq")
    print(type(m.fs.costing.total_TPC_eq))
    print(pyo.Constraint)
    assert type(m.fs.costing.total_TPC_eq) is ScalarConstraint


@pytest.mark.component
def test_build_process_costs_emptymodel_nonearguments():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # Fixed and Variable Costs:
    # build variable costs components

    assert not hasattr(m.fs.costing, "total_TPC")
    assert not hasattr(m.fs.costing, "total_TPC_eq")

    m.fs.costing.build_process_costs(
        net_power=None,
        fixed_OM=False,
        variable_OM=False,
        resources=None,
        rates=None,
        prices=None,
        fuel=None,
    )

    assert hasattr(m.fs.costing, "total_TPC")
    assert type(m.fs.costing.total_TPC) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "total_TPC_eq")
    print(type(m.fs.costing.total_TPC_eq))
    print(pyo.Constraint)
    assert type(m.fs.costing.total_TPC_eq) is ScalarConstraint


@pytest.mark.component
def test_build_process_costs_fixedonly():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

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
        net_power=m.fs.net_power,
        fixed_OM=True,
        variable_OM=False,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
    )

    assert hasattr(m.fs.costing, "annual_operating_labor_cost")
    assert type(m.fs.costing.annual_operating_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "maintenance_labor_cost")
    assert type(m.fs.costing.maintenance_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "admin_and_support_labor_cost")
    assert type(m.fs.costing.admin_and_support_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "property_taxes_and_insurance")
    assert type(m.fs.costing.property_taxes_and_insurance) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert type(m.fs.costing.total_fixed_OM_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "other_fixed_costs")
    assert type(m.fs.costing.other_fixed_costs) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "maintenance_material_cost")
    assert type(m.fs.costing.maintenance_material_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "annual_labor_cost_rule")
    assert type(m.fs.costing.annual_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "maintenance_labor_cost_rule")
    assert type(m.fs.costing.maintenance_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "admin_and_support_labor_cost_rule")
    assert type(m.fs.costing.admin_and_support_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "taxes_and_insurance_cost_rule")
    assert type(m.fs.costing.taxes_and_insurance_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "total_fixed_OM_cost_rule")
    assert type(m.fs.costing.total_fixed_OM_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "maintenance_material_cost_rule")
    assert type(m.fs.costing.maintenance_material_cost_rule) is ScalarConstraint


@pytest.mark.component
def test_build_process_costs_fixedonly_nonearguments():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

    # Fixed and Variable Costs:
    # build variable costs components

    m.fs.costing.build_process_costs(
        net_power=None,
        fixed_OM=True,
        variable_OM=False,
        resources=None,
        rates=None,
        prices=None,
        fuel=None,
    )

    assert hasattr(m.fs.costing, "annual_operating_labor_cost")
    assert type(m.fs.costing.annual_operating_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "maintenance_labor_cost")
    assert type(m.fs.costing.maintenance_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "admin_and_support_labor_cost")
    assert type(m.fs.costing.admin_and_support_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "property_taxes_and_insurance")
    assert type(m.fs.costing.property_taxes_and_insurance) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert type(m.fs.costing.total_fixed_OM_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "other_fixed_costs")
    assert type(m.fs.costing.other_fixed_costs) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "maintenance_material_cost")
    assert type(m.fs.costing.maintenance_material_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "annual_labor_cost_rule")
    assert type(m.fs.costing.annual_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "maintenance_labor_cost_rule")
    assert type(m.fs.costing.maintenance_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "admin_and_support_labor_cost_rule")
    assert type(m.fs.costing.admin_and_support_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "taxes_and_insurance_cost_rule")
    assert type(m.fs.costing.taxes_and_insurance_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "total_fixed_OM_cost_rule")
    assert type(m.fs.costing.total_fixed_OM_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "maintenance_material_cost_rule")
    assert type(m.fs.costing.maintenance_material_cost_rule) is ScalarConstraint


@pytest.mark.component
def test_build_process_costs_variableonly():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

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
        net_power=m.fs.net_power,
        fixed_OM=False,
        variable_OM=True,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
    )

    assert hasattr(m.fs.costing, "variable_operating_costs")
    assert type(m.fs.costing.variable_operating_costs) is IndexedVar
    assert hasattr(m.fs.costing, "other_variable_costs")
    assert type(m.fs.costing.other_variable_costs) is IndexedVar
    assert hasattr(m.fs.costing, "total_variable_OM_cost")
    assert type(m.fs.costing.total_variable_OM_cost) is IndexedVar
    assert hasattr(m.fs.costing, "variable_cost_rule_power")
    assert type(m.fs.costing.variable_cost_rule_power) is IndexedConstraint
    assert hasattr(m.fs.costing, "total_variable_cost_rule_power")
    assert type(m.fs.costing.total_variable_cost_rule_power) is IndexedConstraint


@pytest.mark.component
def test_build_process_costs_variableonly_nonearguments():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

    # Fixed and Variable Costs:
    # build variable costs components

    # methods don't require net_power, will fail on next error
    with pytest.raises(TypeError, match="resources argument must be a list"):
        m.fs.costing.build_process_costs(
            net_power=None,
            fixed_OM=False,
            variable_OM=True,
            resources=None,
            rates=None,
            prices=None,
            fuel=None,
        )

    # checking that same error occurs whether net_power is passed or not
    m.fs.net_power = pyo.Var(m.fs.time, initialize=650, units=pyunits.MW)
    m.fs.net_power.fix()
    with pytest.raises(TypeError, match="resources argument must be a list"):
        m.fs.costing.build_process_costs(
            net_power=m.fs.net_power,
            fixed_OM=False,
            variable_OM=True,
            resources=None,
            rates=None,
            prices=None,
            fuel=None,
        )

    resources = list()
    with pytest.raises(TypeError, match="rates argument must be a list"):
        m.fs.costing.build_process_costs(
            net_power=m.fs.net_power,
            fixed_OM=False,
            variable_OM=True,
            resources=resources,
            rates=None,
            prices=None,
            fuel=None,
        )

    rates = list()
    with pytest.raises(TypeError, match="prices argument must be a dictionary"):
        m.fs.costing.build_process_costs(
            net_power=m.fs.net_power,
            fixed_OM=False,
            variable_OM=True,
            resources=resources,
            rates=rates,
            prices=None,
            fuel=None,
        )

    prices = dict()
    # fuel being None should not raise an Exception
    m.fs.costing.build_process_costs(
        net_power=m.fs.net_power,
        fixed_OM=False,
        variable_OM=True,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
    )


@pytest.mark.component
def test_build_process_costs_allOM():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

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

    m.fs.costing.sorbent = pyo.Var(
        m.fs.time, initialize=1.0, units=pyunits.ft**3 / pyunits.d
    )  # ft**3/day (EPAT: 6971.75)

    resources = ["natural_gas", "solvent", "waste_sorbent", "sorbent"]
    rates = [
        m.fs.NG_rate,
        m.fs.solvent_rate,
        m.fs.costing.sorbent,
        m.fs.costing.sorbent,
    ]
    prices = {
        "solvent": 500 * pyunits.USD_2018 / pyunits.ton,
        "sorbent": 4 * pyunits.USD_2018 / pyunits.ft**3,
        "waste_sorbent": 0.86 * pyunits.USD_2018 / pyunits.ft**3,
    }

    m.fs.costing.land_cost = pyo.Expression(
        expr=156000 * (30 / 120) ** (0.78)
    )  # 30 is a fixed value, 30 must be replaced by a model variable
    m.fs.tonne_CO2_capture = pyo.Var(initialize=1, units=pyunits.ton)

    m.fs.costing.build_process_costs(
        net_power=m.fs.net_power,
        fixed_OM=True,
        variable_OM=True,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
        waste=["waste_sorbent"],
        chemicals=["sorbent"],
        land_cost=m.fs.costing.land_cost,
        tonne_CO2_capture=m.fs.tonne_CO2_capture,
    )

    assert hasattr(m.fs.costing, "annual_operating_labor_cost")
    assert type(m.fs.costing.annual_operating_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "maintenance_labor_cost")
    assert type(m.fs.costing.maintenance_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "admin_and_support_labor_cost")
    assert type(m.fs.costing.admin_and_support_labor_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "property_taxes_and_insurance")
    assert type(m.fs.costing.property_taxes_and_insurance) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert type(m.fs.costing.total_fixed_OM_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "other_fixed_costs")
    assert type(m.fs.costing.other_fixed_costs) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "maintenance_material_cost")
    assert type(m.fs.costing.maintenance_material_cost) is pyo.ScalarVar
    assert hasattr(m.fs.costing, "annual_labor_cost_rule")
    assert type(m.fs.costing.annual_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "maintenance_labor_cost_rule")
    assert type(m.fs.costing.maintenance_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "admin_and_support_labor_cost_rule")
    assert type(m.fs.costing.admin_and_support_labor_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "taxes_and_insurance_cost_rule")
    assert type(m.fs.costing.taxes_and_insurance_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "total_fixed_OM_cost_rule")
    assert type(m.fs.costing.total_fixed_OM_cost_rule) is ScalarConstraint
    assert hasattr(m.fs.costing, "maintenance_material_cost_rule")
    assert type(m.fs.costing.maintenance_material_cost_rule) is ScalarConstraint

    assert hasattr(m.fs.costing, "variable_operating_costs")
    assert type(m.fs.costing.variable_operating_costs) is IndexedVar
    assert hasattr(m.fs.costing, "other_variable_costs")
    assert type(m.fs.costing.other_variable_costs) is IndexedVar
    assert hasattr(m.fs.costing, "total_variable_OM_cost")
    assert type(m.fs.costing.total_variable_OM_cost) is IndexedVar
    assert hasattr(m.fs.costing, "variable_cost_rule_power")
    assert type(m.fs.costing.variable_cost_rule_power) is IndexedConstraint
    assert hasattr(m.fs.costing, "total_variable_cost_rule_power")
    assert type(m.fs.costing.total_variable_cost_rule_power) is IndexedConstraint

    assert hasattr(m.fs.costing, "waste_costs_OC")
    assert type(m.fs.costing.waste_costs_OC) is ScalarExpression
    assert hasattr(m.fs.costing, "chemical_costs_OC")
    assert type(m.fs.costing.waste_costs_OC) is ScalarExpression
    assert hasattr(m.fs.costing, "cost_of_capture")
    assert type(m.fs.costing.cost_of_capture) is ScalarExpression


@pytest.mark.component
def test_build_process_costs_allOM_nonearguments():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

    # Fixed and Variable Costs:
    # build variable costs components

    # methods don't require net_power, will fail on next error
    with pytest.raises(TypeError, match="resources argument must be a list"):
        m.fs.costing.build_process_costs(
            net_power=None,
            fixed_OM=False,
            variable_OM=True,
            resources=None,
            rates=None,
            prices=None,
            fuel=None,
        )

    # checking that same error occurs whether net_power is passed or not
    m.fs.net_power = pyo.Var(m.fs.time, initialize=650, units=pyunits.MW)
    m.fs.net_power.fix()
    with pytest.raises(TypeError, match="resources argument must be a list"):
        m.fs.costing.build_process_costs(
            net_power=m.fs.net_power,
            fixed_OM=False,
            variable_OM=True,
            resources=None,
            rates=None,
            prices=None,
            fuel=None,
        )

    resources = list()
    with pytest.raises(TypeError, match="rates argument must be a list"):
        m.fs.costing.build_process_costs(
            net_power=m.fs.net_power,
            fixed_OM=False,
            variable_OM=True,
            resources=resources,
            rates=None,
            prices=None,
            fuel=None,
        )

    rates = list()
    with pytest.raises(TypeError, match="prices argument must be a dictionary"):
        m.fs.costing.build_process_costs(
            net_power=m.fs.net_power,
            fixed_OM=False,
            variable_OM=True,
            resources=resources,
            rates=rates,
            prices=None,
            fuel=None,
        )

    prices = dict()
    # fuel being None should not raise an Exception
    m.fs.costing.build_process_costs(
        net_power=m.fs.net_power,
        fixed_OM=False,
        variable_OM=True,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel=None,
    )


@pytest.mark.component
def test_build_process_costs_invalid_currency_units():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    with pytest.raises(
        AttributeError,
        match="CE_index_year notavalidvalue is not a "
        "valid currency base option. Valid CE index options "
        "include CE500, CE394 and years from 1990 to 2020.",
    ):
        m.fs.costing.build_process_costs(
            net_power=None,
            fixed_OM=False,
            variable_OM=False,
            resources=None,
            rates=None,
            prices=None,
            fuel=None,
            CE_index_year="notavalidvalue",
        )


@pytest.mark.component
def test_power_plant_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    ###########################################################################
    #  Create costing constraints                                             #
    ###########################################################################

    # subcritical PC
    coal_accounts = ["1.1", "1.2", "1.3"]
    m.fs.subcritical_PC = UnitModelBlock()
    m.fs.subcritical_PC.coal_feed_rate = pyo.Var(
        initialize=7613.37, units=pyunits.ton / pyunits.day
    )
    m.fs.subcritical_PC.coal_feed_rate.fix()
    m.fs.subcritical_PC.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": coal_accounts,
            "scaled_param": m.fs.subcritical_PC.coal_feed_rate,
            "tech": 1,
            "ccs": "A",
        },
    )

    # two-stage, slurry-feed IGCC
    feedwater_accounts = ["3.1", "3.3", "3.5"]
    m.fs.IGCC_1 = UnitModelBlock()
    m.fs.IGCC_1.feedwater_flow_rate = pyo.Var(
        initialize=1576062.15, units=pyunits.lb / pyunits.hr
    )
    m.fs.IGCC_1.feedwater_flow_rate.fix()
    m.fs.IGCC_1.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": feedwater_accounts,
            "scaled_param": m.fs.IGCC_1.feedwater_flow_rate,
            "tech": 3,
            "ccs": "A",
        },
    )

    # single-stage, slurry-feed, IGCC
    syngas_accounts = ["6.1", "6.2", "6.3"]
    m.fs.IGCC_2 = UnitModelBlock()
    m.fs.IGCC_2.syngas_flow_rate = pyo.Var(
        initialize=182335.921, units=pyunits.lb / pyunits.hr
    )
    m.fs.IGCC_2.syngas_flow_rate.fix()
    m.fs.IGCC_2.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": syngas_accounts,
            "scaled_param": m.fs.IGCC_2.syngas_flow_rate,
            "tech": 4,
            "ccs": "A",
        },
    )

    # single-stage, dry-feed, IGCC
    HRSG_accounts = ["7.1", "7.2"]
    m.fs.IGCC_3 = UnitModelBlock()
    m.fs.IGCC_3.HRSG_duty = pyo.Var(initialize=1777.86, units=pyunits.MBtu / pyunits.hr)
    m.fs.IGCC_3.HRSG_duty.fix()
    m.fs.IGCC_3.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": HRSG_accounts,
            "scaled_param": m.fs.IGCC_3.HRSG_duty,
            "tech": 5,
            "ccs": "A",
        },
    )

    # NGCC
    steam_turbine_accounts = ["8.1", "8.2", "8.5"]
    m.fs.NGCC = UnitModelBlock()
    m.fs.NGCC.turbine_power = pyo.Var(initialize=212500, units=pyunits.kW)
    m.fs.NGCC.turbine_power.fix()
    m.fs.NGCC.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": steam_turbine_accounts,
            "scaled_param": m.fs.NGCC.turbine_power,
            "tech": 6,
            "ccs": "A",
        },
    )

    # AUSC PC
    AUSC_accounts = ["4.9", "8.4"]
    m.fs.AUSC = UnitModelBlock()
    m.fs.AUSC.feedwater_flow = pyo.Var(
        initialize=3298815.58, units=pyunits.lb / pyunits.hr
    )
    m.fs.AUSC.feedwater_flow.fix()
    m.fs.AUSC.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": AUSC_accounts,
            "scaled_param": m.fs.AUSC.feedwater_flow,
            "tech": 7,
            "ccs": "B",
        },
    )

    # custom carbon capture
    CCS_accounts = ["5.1.a.epri"]
    m.fs.CCS = UnitModelBlock()
    m.fs.CCS.CO2_flow = pyo.Var(initialize=493587.88, units=pyunits.lb / pyunits.hr)
    m.fs.CCS.CO2_flow.fix()
    m.fs.CCS.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": CCS_accounts,
            "scaled_param": m.fs.CCS.CO2_flow,
            "tech": 6,
            "ccs": "B",
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
        net_power=m.fs.net_power,
        fixed_OM=True,
        variable_OM=False,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
    )

    # add initialize
    QGESSCostingData.costing_initialization(m.fs.costing)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal
    #  all numbers come from the NETL excel file
    # "201.001.001_BBR4 COE Spreadsheet_Rev0U_20190919_njk.xlsm"
    assert (
        pytest.approx(
            pyo.value(m.fs.subcritical_PC.costing.total_plant_cost["1.1"]), abs=1e-1
        )
        == 2379 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.subcritical_PC.costing.total_plant_cost["1.2"]), abs=1e-1
        )
        == 6588 / 1e3
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.subcritical_PC.costing.total_plant_cost["1.3"]), abs=1e-1
        )
        == 61409 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.IGCC_1.costing.total_plant_cost["3.1"]), abs=1e-1)
        == 10807 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.IGCC_1.costing.total_plant_cost["3.3"]), abs=1e-1)
        == 2564 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.IGCC_1.costing.total_plant_cost["3.5"]), abs=1e-1)
        == 923 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.IGCC_2.costing.total_plant_cost["6.1"]), abs=1e-1)
        == 110873 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.IGCC_2.costing.total_plant_cost["6.2"]), abs=1e-1)
        == 3207 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.IGCC_2.costing.total_plant_cost["6.3"]), abs=1e-1)
        == 3770 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.IGCC_3.costing.total_plant_cost["7.1"]), abs=1e-1)
        == 53530 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.IGCC_3.costing.total_plant_cost["7.2"]), abs=1e-1)
        == 19113 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.NGCC.costing.total_plant_cost["8.1"]), abs=1e-1)
        == 49468 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.NGCC.costing.total_plant_cost["8.2"]), abs=1e-1)
        == 565 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.NGCC.costing.total_plant_cost["8.5"]), abs=1e-1)
        == 4094 / 1e3
    )

    assert (
        pytest.approx(pyo.value(m.fs.AUSC.costing.bare_erected_cost["4.9"]), abs=1e-1)
        == 295509 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.AUSC.costing.bare_erected_cost["8.4"]), abs=1e-1)
        == 57265 / 1e3
    )

    # test various report utilities - "smoke tests" only
    m.fs.costing.report()
    QGESSCostingData.display_total_plant_costs(m.fs.costing)
    QGESSCostingData.display_bare_erected_costs(m.fs.costing)
    QGESSCostingData.display_equipment_costs(m.fs.costing)
    QGESSCostingData.display_flowsheet_cost(m.fs.costing)

    return m


@pytest.mark.component
def test_sCO2_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # ######################################################
    # Primary Heater
    m.fs.boiler = UnitModelBlock()
    m.fs.boiler.heat_duty = pyo.Var(initialize=1461.5, units=pyunits.MW)
    m.fs.boiler.heat_duty.fix()
    m.fs.boiler.temp = pyo.Var(initialize=620, units=pyunits.C)  # C
    m.fs.boiler.temp.fix()
    m.fs.boiler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Coal-fired heater",
            "scaled_param": m.fs.boiler.heat_duty,
            "temp_C": m.fs.boiler.temp,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # CO2 Turbine
    m.fs.turbine = UnitModelBlock()
    m.fs.turbine.work_isentropic = pyo.Var(initialize=1006.2, units=pyunits.MW)
    m.fs.turbine.work_isentropic.fix()
    m.fs.turbine.temp = pyo.Var(initialize=620, units=pyunits.C)
    m.fs.turbine.temp.fix()
    m.fs.turbine.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Axial turbine",
            "scaled_param": m.fs.turbine.work_isentropic,
            "temp_C": m.fs.turbine.temp,
            "n_equip": 1,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # Generator
    m.fs.generator = UnitModelBlock()
    m.fs.generator.work_isentropic = pyo.Var(initialize=1006.2, units=pyunits.MW)
    m.fs.generator.work_isentropic.fix()
    m.fs.generator.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Generator",
            "scaled_param": m.fs.generator.work_isentropic,
            "n_equip": 1,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # High Temperature Recuperator
    m.fs.HTR = UnitModelBlock()
    m.fs.HTR.heat_duty = pyo.Var(initialize=1461e6, units=pyunits.W)
    m.fs.HTR.heat_duty.fix()
    m.fs.HTR.LMTD = pyo.Var(initialize=21.45, units=pyunits.C)
    m.fs.HTR.LMTD.fix()
    m.fs.HTR.temp = pyo.Var(initialize=453, units=pyunits.C)
    m.fs.HTR.temp.fix()

    m.fs.HTR.UA = pyo.Var(initialize=1e8, units=pyunits.W / pyunits.C)

    # gives units of W/C = W/K
    @m.fs.Constraint()
    def HTR_UA_rule(b):
        return b.HTR.UA * b.HTR.LMTD == b.HTR.heat_duty

    m.fs.HTR.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Recuperator",
            "scaled_param": m.fs.HTR.UA,
            "temp_C": m.fs.HTR.temp,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # Low Temperature Recuperator
    m.fs.LTR = UnitModelBlock()
    m.fs.LTR.heat_duty = pyo.Var(initialize=911.7e6, units=pyunits.W)
    m.fs.LTR.heat_duty.fix()
    m.fs.LTR.LMTD = pyo.Var(initialize=5.21, units=pyunits.C)
    m.fs.LTR.LMTD.fix()
    m.fs.LTR.temp = pyo.Var(initialize=216, units=pyunits.C)
    m.fs.LTR.temp.fix()
    m.fs.LTR.UA = pyo.Var(initialize=1e8, units=pyunits.W / pyunits.C)

    @m.fs.Constraint()
    def LTR_UA_rule(b):
        return b.LTR.UA * b.LTR.LMTD == b.LTR.heat_duty

    m.fs.LTR.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Recuperator",
            "scaled_param": m.fs.LTR.UA,
            "temp_C": m.fs.LTR.temp,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # CO2 Cooler, costed using the recouperator not dry cooler
    m.fs.co2_cooler = UnitModelBlock()
    m.fs.co2_cooler.heat_duty = pyo.Var(initialize=739.421217e6, units=pyunits.W)
    m.fs.co2_cooler.heat_duty.fix()
    m.fs.co2_cooler.temp = pyo.Var(initialize=81, units=pyunits.C)
    m.fs.co2_cooler.temp.fix()

    # Estimating LMTD
    # Cost from report: $27,780 thousand
    # Back-calculated UA: 41819213 W/K
    # Heat duty from report: 2523 MMBTu/hr --> 739421217 W
    # Estimated LMTD: 17.68 K
    m.fs.co2_cooler.LMTD = pyo.Var(initialize=5, units=pyunits.C)
    m.fs.co2_cooler.UA = pyo.Var(initialize=1e5, units=pyunits.W / pyunits.C)
    m.fs.co2_cooler.LMTD.fix(17.68 * pyunits.C)

    @m.fs.Constraint()
    def co2_cooler_UA_rule(b):
        return b.co2_cooler.UA * b.co2_cooler.LMTD == b.co2_cooler.heat_duty

    m.fs.co2_cooler.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Recuperator",
            "scaled_param": m.fs.co2_cooler.UA,
            "temp_C": m.fs.co2_cooler.temp,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # Main Compressor - 5.99 m^3/s in Baseline620
    m.fs.main_compressor = UnitModelBlock()
    m.fs.main_compressor.flow_vol = pyo.Var(
        initialize=5.99, units=pyunits.m**3 / pyunits.s
    )
    m.fs.main_compressor.flow_vol.fix()
    m.fs.main_compressor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Barrel type compressor",
            "scaled_param": m.fs.main_compressor.flow_vol,
            "n_equip": 5.0,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # Main Compressor Motor
    m.fs.main_compressor_motor = UnitModelBlock()
    m.fs.main_compressor_motor.work_isentropic = pyo.Var(
        initialize=159.7, units=pyunits.MW
    )
    m.fs.main_compressor_motor.work_isentropic.fix()
    m.fs.main_compressor_motor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Open drip-proof motor",
            "scaled_param": m.fs.main_compressor_motor.work_isentropic,
            "n_equip": 5.0,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # Recompressor - 6.89 m^3/s in Baseline620
    m.fs.bypass_compressor = UnitModelBlock()
    m.fs.bypass_compressor.flow_vol = pyo.Var(
        initialize=6.89, units=pyunits.m**3 / pyunits.s
    )
    m.fs.bypass_compressor.flow_vol.fix()
    m.fs.bypass_compressor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Barrel type compressor",
            "scaled_param": m.fs.bypass_compressor.flow_vol,
            "n_equip": 4.0,
            "CE_index_year": "2017",
        },
    )

    # ######################################################
    # Recompressor Motor
    m.fs.bypass_compressor_motor = UnitModelBlock()
    m.fs.bypass_compressor_motor.work_isentropic = pyo.Var(
        initialize=124.3, units=pyunits.MW
    )
    m.fs.bypass_compressor_motor.work_isentropic.fix()

    m.fs.bypass_compressor_motor.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_sCO2_unit_cost,
        costing_method_arguments={
            "equipment": "Open drip-proof motor",
            "scaled_param": m.fs.bypass_compressor_motor.work_isentropic,
            "n_equip": 4.0,
            "CE_index_year": "2017",
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
        net_power=m.fs.net_power,
        fixed_OM=True,
        variable_OM=False,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
        CE_index_year="2017",
    )

    # add initialize
    QGESSCostingData.costing_initialization(m.fs.costing)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal

    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.equipment_cost), abs=1e-1)
        == 216.291
    )
    assert (
        pytest.approx(pyo.value(m.fs.turbine.costing.equipment_cost), abs=1e-1)
        == 13.180
    )
    assert (
        pytest.approx(pyo.value(m.fs.generator.costing.equipment_cost), abs=1e-1)
        == 4.758
    )
    assert pytest.approx(pyo.value(m.fs.HTR.costing.equipment_cost), abs=1e-1) == 40.156
    assert pytest.approx(pyo.value(m.fs.LTR.costing.equipment_cost), abs=1e-1) == 81.807
    assert (
        pytest.approx(pyo.value(m.fs.co2_cooler.costing.equipment_cost), abs=1e-1)
        == 27.784
    )
    assert (
        pytest.approx(pyo.value(m.fs.main_compressor.costing.equipment_cost), abs=1e-1)
        == 31.732
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.bypass_compressor.costing.equipment_cost), abs=1e-1
        )
        == 26.434
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.main_compressor_motor.costing.equipment_cost)
            + pyo.value(m.fs.bypass_compressor_motor.costing.equipment_cost),
            abs=1e-1,
        )
        == 29.133
    )
    assert (
        pytest.approx(
            pyo.value(m.fs.bypass_compressor.costing.equipment_cost), abs=1e-1
        )
        == 26.434
    )

    # test bound check utility
    QGESSCostingData.check_sCO2_costing_bounds(m.fs.costing)

    return m


@pytest.mark.component
def test_ASU_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    m.fs.ASU = UnitModelBlock()
    m.fs.ASU.O2_flow = pyo.Var(initialize=13078, units=pyunits.ton / pyunits.d)
    m.fs.ASU.O2_flow.fix()

    m.fs.ASU.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_ASU_cost,
        costing_method_arguments={
            "scaled_param": m.fs.ASU.O2_flow,
            "CE_index_year": "2017",
        },
    )

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal

    assert pytest.approx(pyo.value(m.fs.ASU.costing.bare_erected_cost), abs=1) == 3.4725

    return m


@pytest.mark.component
def test_OM_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=True, time_units=pyunits.s)
    m.fs.costing = QGESSCosting()

    # build fixed costs
    nameplate_capacity = 650  # MW
    labor_rate = 40
    labor_burden = 30
    operators_per_shift = 10
    tech = 1
    fixed_TPC = 800  # MM$

    QGESSCostingData.get_fixed_OM_costs(
        m.fs.costing,
        nameplate_capacity=nameplate_capacity,
        labor_rate=labor_rate,
        labor_burden=labor_burden,
        operators_per_shift=operators_per_shift,
        tech=tech,
        fixed_TPC=fixed_TPC,
    )

    # build variable costs
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

    QGESSCostingData.get_variable_OM_costs(
        m.fs.costing,
        resources=resources,
        rates=rates,
        prices=prices,  # pass a flowsheet object
    )

    QGESSCostingData.initialize_fixed_OM_costs(m.fs.costing)
    QGESSCostingData.initialize_variable_OM_costs(m.fs.costing)

    assert hasattr(m.fs, "costing")
    assert hasattr(m.fs.costing, "total_fixed_OM_cost")
    assert hasattr(m.fs.costing, "total_variable_OM_cost")

    assert degrees_of_freedom(m) == 0

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)

    assert results.solver.termination_condition == pyo.TerminationCondition.optimal

    assert pytest.approx(28.094, abs=0.1) == (
        pyo.value(m.fs.costing.total_fixed_OM_cost)
    )

    assert pytest.approx(182.367, abs=0.1) == (
        pyo.value(m.fs.costing.total_variable_OM_cost[0])
    )
