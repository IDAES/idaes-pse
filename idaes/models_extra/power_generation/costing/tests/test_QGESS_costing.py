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
Tests for costing package based on methods from:

    Process and Product Design Principles: Synthesis, Analysis, and
    Evaluation
    Seider, Seader, Lewin, Windagdo, 3rd Ed. John Wiley and Sons
    Chapter 22. Cost Accounting and Capital Cost Estimation
    22.2 Cost Indexes and Capital Investment
"""
import pytest

from pyomo.environ import (
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Param,
    units as pyunits,
    value,
    Var,
    Expression,
)
from pyomo.util.check_units import assert_units_consistent
from pyomo.common.config import ConfigValue

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom


from idaes.models.properties import iapws95


from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)
import pyomo.environ as pyo

# Some more information about this module
__author__ = "Costing Team (A. Noring and M. Zamarripa)"


solver = get_solver()


@pytest.fixture
def model():
    m = ConcreteModel()

    m.fs = FlowsheetBlock()

    m.fs.costing = QGESSCosting()

    # Add a placeholder to represent a unit model
    m.fs.boiler = UnitModelBlock()
    coal_accounts = ["1.1", "1.2", "1.3", "1.4", "2.1",
                     "2.2", "4.11", "4.15", "4.16"]
    # m.fs.boiler = pyo.Block()
    m.fs.boiler.coal_mass_flow = pyo.Var(initialize=7238.95)  # tpd
    m.fs.boiler.coal_mass_flow.fix()
    # get_PP_costing(m.fs.boiler, coal_accounts,
    #                m.fs.boiler.coal_mass_flow, "tpd", 2)

    m.fs.boiler.costing = UnitModelCostingBlock(
        default={
            "flowsheet_costing_block": m.fs.costing,
            "costing_method": QGESSCostingData.get_PP_costing,
            "costing_method_arguments": {
                "cost_accounts": coal_accounts,
                "scaled_param": m.fs.boiler.coal_mass_flow,
                "units": "tpd",
                "tech": 2,
                "ccs": "A",
            },
        }
    )

    # Fixed and Variable Costs:
    # build variable costs components
    m.fs.net_power = pyo.Var(m.fs.time, initialize=650, units=pyunits.MW)
    m.fs.net_power.fix()

    m.fs.NG_rate = pyo.Var(m.fs.time, initialize=1.2,
                           units=pyunits.MBtu / pyunits.s)
    m.fs.NG_rate.fix()

    m.fs.solvent_rate = pyo.Var(m.fs.time, initialize=40,
                                units=pyunits.ton / pyunits.day)
    m.fs.solvent_rate.fix()

    resources = ["natural_gas", "solvent"]
    rates = [m.fs.NG_rate, m.fs.solvent_rate]
    prices = {"solvent": 500 * pyunits.USD_2018 / pyunits.ton}

    m.fs.costing.land_cost = Expression(
        expr=156000 * (30 / 120) ** (0.78)
    )  # 30 is a fixed value, 30 must be replaced by a model variable
    m.fs.costing.display()
    m.fs.costing.build_process_costs(
        net_power=m.fs.net_power,
        fixed_OM=True,
        variable_OM=False,
        resources=resources,
        rates=rates,
        prices=prices,
        fuel="natural_gas",
    )

    return m


@pytest.mark.unit
def test_solve_model(model):
    m = model  # shorter alias
    solver = get_solver()
    results = solver.solve(m, tee=True)

    assert results.solver.termination_condition == \
        pyo.TerminationCondition.optimal

    # #  all numbers come from the NETL excel file:
    # # "201.001.001_BBR4 COE Spreadsheet_Rev0U_20190919_njk.xlsm"
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.1"]),
                      abs=1e-1)
        == 2306 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.2"]),
                      abs=1e-1)
        == 6385 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.3"]),
                      abs=1e-1)
        == 59527 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["1.4"]),
                      abs=1e-1)
        == 8086 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["2.1"]),
                      abs=1e-1)
        == 4073 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["2.2"]),
                      abs=1e-1)
        == 13976 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["4.11"]),
                      abs=1e-1)
        == 3751 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["4.15"]),
                      abs=1e-1)
        == 197 / 1e3
    )
    assert (
        pytest.approx(pyo.value(m.fs.boiler.costing.total_plant_cost["4.16"]),
                      abs=1e-1)
        == 1014 / 1e3
    )
