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

__author__ = (
    "Costing Team (B. Paul and M. Zamarripa)"
)
__version__ = "1.0.0"

import pytest

import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)


@pytest.mark.component
def test_REE_costing():
    # Create a Concrete Model as the top level object
    m = pyo.ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.costing = QGESSCosting()

    # check that the model solved properly and has 0 degrees of freedom
    assert degrees_of_freedom(m) == 0

    # create new account to test, this is a new entry for a type 8 REE process
    additional_costing_params = {
        "8": {
            "A": {
                "1.1": {
                    "Account Name": "Crushing & Screening - Front End Loader",
                    "BEC": 737000/2.97/5,
                    "BEC_units": "$2019",
                    "Eng Fee": 1.82,
                    "Exponent": 0.55,
                    "Process Contingency": 0.15,
                    "Process Parameter": "Total load (TPH)",
                    "Project Contingency": 0,
                    "RP Value": 500/5,
                    "Units": "ton/hr"
                    }
                }
            }
        }

    ###########################################################################
    #  Create costing constraints                                             #
    ###########################################################################

    # ree load rate
    # accounts 1.x are crushing and screening
    ree_accounts = ["1.1"]
    m.fs.front_end_loader = UnitModelBlock()
    m.fs.front_end_loader.load = pyo.Var(
        initialize=500,
        units=pyunits.ton / pyunits.hr
    )
    m.fs.front_end_loader.load.fix()
    m.fs.front_end_loader.costing = UnitModelCostingBlock(
        flowsheet_costing_block=m.fs.costing,
        costing_method=QGESSCostingData.get_PP_costing,
        costing_method_arguments={
            "cost_accounts": ree_accounts,
            "scaled_param": m.fs.front_end_loader.load,
            "tech": 8,
            "ccs": "A",
            "additional_costing_params": additional_costing_params,
            "use_additional_costing_params": True,
            "number_parallel_trains": 5
        },
    )

    # add initialize
    QGESSCostingData.costing_initialization(m.fs.costing)

    # try solving
    solver = get_solver()
    results = solver.solve(m, tee=True)

    # check unit consistency
    assert_units_consistent(m)
    print(pyo.value(m.fs.front_end_loader.costing.total_plant_cost["1.1"]))
    # check that existing account was replaced, and the new data was used
    assert (
        pytest.approx(pyo.value(
            m.fs.front_end_loader.costing.total_plant_cost["1.1"]
            ), rel=1e-3)
        == 1.773144  # 731.662 / 1e3
    )
