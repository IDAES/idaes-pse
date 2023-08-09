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

import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.util.check_units import assert_units_consistent
from pyomo.environ import check_optimal_termination

from idaes.core import FlowsheetBlock, UnitModelBlock, UnitModelCostingBlock
from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import degrees_of_freedom

from idaes.models_extra.power_generation.costing.power_plant_capcost import (
    QGESSCosting,
    QGESSCostingData,
)

# Create a Concrete Model as the top level object
m = pyo.ConcreteModel()

# Add a flowsheet object to the model
m.fs = FlowsheetBlock(dynamic=False)
m.fs.costing = QGESSCosting()
CE_index_year = "2016"

# check that the model solved properly and has 0 degrees of freedom
assert degrees_of_freedom(m) == 0

###########################################################################
#  Create costing constraints                                             #
###########################################################################

# accounts 1.x are crushing and screening
# 1.1 is Front End Loader (2 cuyd)
# this is a constant-cost unit, where n_equip is the scaling parameter
front_end_loader_2yd3_accounts = ["1.1"]

m.fs.front_end_loader_2yd3 = UnitModelBlock()
m.fs.front_end_loader_2yd3.n_equip = pyo.Var(
    initialize=5,
    units=pyunits.dimensionless
)
m.fs.front_end_loader_2yd3.n_equip.fix()

m.fs.front_end_loader_2yd3.costing = UnitModelCostingBlock(
    flowsheet_costing_block=m.fs.costing,
    costing_method=QGESSCostingData.get_REE_costing,
    costing_method_arguments={
        "cost_accounts": front_end_loader_2yd3_accounts,
        "scaled_param": m.fs.front_end_loader_2yd3.n_equip, # n_equip = 5
        "tech": 1,
        # n_equip is the scaling parameter for constant-cost units,
        # so here "n_equip" should be entered as 1 set of loaders
        "n_equip": 1,
        "scale_down_parallel_equip": False,
        "CE_index_year": CE_index_year,
    },
)

# 1.2 is Jaw Crusher
jaw_crusher_accounts = ["1.2"]
m.fs.jaw_crusher = UnitModelBlock()
m.fs.jaw_crusher.power = pyo.Var(
    initialize=589,
    units=pyunits.hp
)
m.fs.jaw_crusher.power.fix()
m.fs.jaw_crusher.costing = UnitModelCostingBlock(
    flowsheet_costing_block=m.fs.costing,
    costing_method=QGESSCostingData.get_REE_costing,
    costing_method_arguments={
        "cost_accounts": jaw_crusher_accounts,
        "scaled_param": m.fs.jaw_crusher.power,
        "tech": 1,
        "n_equip": 1,
        "scale_down_parallel_equip": False,
        "CE_index_year": CE_index_year,
    },
)

# add plant-level cost constraints
QGESSCostingData.get_REE_plant_costs(m.fs.costing,
                                     CE_index_year=CE_index_year)

# add initialize
QGESSCostingData.costing_initialization(m.fs.costing)

# add initialize for plant-level costs
QGESSCostingData.initialize_REE_plant_costs(m.fs.costing)

# try solving
solver = get_solver()
results = solver.solve(m, tee=True)
assert check_optimal_termination(results)

# check unit consistency
assert_units_consistent(m)

# report results
QGESSCostingData.report(m.fs.costing)
