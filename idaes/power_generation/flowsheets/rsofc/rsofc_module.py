###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University,
# West Virginia University Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
###############################################################################

__author__ = "Chinedu Okoli"


# Import Python libraries
import sys
import os

import numpy as np
import matplotlib.pyplot as plt


# Import Pyomo modules
import pyomo.environ as pyo


# Import IDAES modules
import idaes.core.plugins
from idaes.core.solvers import use_idaes_solver_configuration_defaults

from idaes.power_generation.flowsheets.rsofc import (
    rsofc_costing as rsofc_cost,
    rsofc_sofc_mode_flowsheet_v3 as rsofc_sofc,
    rsofc_soec_mode_flowsheet_v6 as rsofc_soec
    )


def add_sofc_mode_flowsheet(m, name="sofc_mode_flowsheet"):
    # if m is None:
    #     m = pyo.ConcreteModel(name)
    # if not hasattr(m, "sofc_fs"):
    #     m.sofc_fs = FlowsheetBlock(default={"dynamic": False})

    rsofc_sofc.get_model(m)


def add_soec_mode_flowsheet(m, name="soec_mode_flowsheet"):
    # if m is None:
    #     m = pyo.ConcreteModel(name)
    # if not hasattr(m, "soec_fs"):
    #     m.soec_fs = FlowsheetBlock(default={"dynamic": False})

    rsofc_soec.get_model(m)


def get_solver():
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    idaes.cfg.ipopt["options"]["tol"] = 1e-7
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-9
    idaes.cfg.ipopt["options"]["linear_solver"] = "ma27"
    idaes.cfg.ipopt["options"]["max_iter"] = 400
    idaes.cfg.ipopt["options"]["ma27_pivtol"] = 1e-3
    # idaes.cfg.ipopt["options"]["ma57_pivtol"] = 1e-1
    return pyo.SolverFactory("ipopt")


def get_model(use_DNN=True):
    # Create model and add flowsheets
    m = pyo.ConcreteModel()
    add_sofc_mode_flowsheet(m)
    add_soec_mode_flowsheet(m)
    return m


def cost_rsofc(m):
    # Capital and fixed costs
    # Add both sofc and soec mode flowsheets before computing cap&fixed costs
    rsofc_cost.get_rsofc_capital_costing(m)
    design_rsofc_netpower = 650  # MW
    rsofc_cost.get_rsofc_fixed_OM_costing(m, design_rsofc_netpower)

    # # Solve for capital and fixed costs
    # solver = get_solver()
    # solver.solve(m, tee=True)  # cap&fixed costs calculated in fs flowsheet

    rsofc_cost.lock_capital_cost(m)

    # Variable O&M costs
    rsofc_cost.get_rsofc_soec_variable_OM_costing(m.soec_fs)
    rsofc_cost.get_rsofc_sofc_variable_OM_costing(m.sofc_fs)

    # # Display results
    # m.fs.costing.display()
    # m.sofc_fs.costing.display()
    # m.soec_fs.H2_costing.display()

    return m


# %%       # ------------------------------------------------------------------
if __name__ == "__main__":
    m = get_model()
    cost_rsofc(m)
