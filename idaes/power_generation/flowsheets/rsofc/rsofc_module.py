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

"""
Example simulation and costing module for the reversible sofc

"""

__author__ = "Chinedu Okoli"


# Import Python libraries
# import sys
# import os
# import numpy as np
# import matplotlib.pyplot as plt


# Import Pyomo modules
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES modules
import idaes.core.plugins
from idaes.core.solvers import use_idaes_solver_configuration_defaults

from idaes.power_generation.flowsheets.rsofc import (
    rsofc_costing as rsofc_cost,
    rsofc_sofc_flowsheet as rsofc_sofc,
    rsofc_soec_flowsheet as rsofc_soec
    )


def add_sofc_mode_flowsheet(m, name="sofc_mode_flowsheet"):
    rsofc_sofc.get_model(m)


def add_soec_mode_flowsheet(m, name="soec_mode_flowsheet"):
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
    rsofc_sofc.get_model(m)
    rsofc_soec.get_model(m)
    # add_sofc_mode_flowsheet(m)
    # add_soec_mode_flowsheet(m)
    return m


def cost_rsofc(m):
    # Capital and fixed costs
    # Add both sofc and soec mode flowsheets before computing cap&fixed costs
    rsofc_cost.get_rsofc_sofc_capital_cost(m.sofc_fs)
    rsofc_cost.get_rsofc_soec_capital_cost(m.soec_fs)

    design_rsofc_netpower = 650  # MW
    fixed_TPC = 693.019  # TPC $MM for sofc
    rsofc_cost.get_rsofc_sofc_fixed_OM_costing(m.sofc_fs,
                                               design_rsofc_netpower,
                                               fixed_TPC)
    rsofc_cost.get_rsofc_soec_fixed_OM_costing(m.soec_fs,
                                               design_rsofc_netpower)

    # # Solve for capital and fixed costs
    # solver = get_solver()
    # solver.solve(m.sofc_fs, tee=True)
    # solver.solve(m.soec_fs, tee=True)

    rsofc_cost.lock_rsofc_sofc_capital_cost(m.sofc_fs)
    rsofc_cost.lock_rsofc_soec_capital_cost(m.soec_fs)

    # Variable O&M costs
    rsofc_cost.get_rsofc_soec_variable_OM_costing(m.soec_fs)
    rsofc_cost.get_rsofc_sofc_variable_OM_costing(m.sofc_fs)

    m.rsofc_total_TPC = pyo.Var(
        initialize=1e-6,
        doc="rsofc total plant cost in $MM",
    )

    @m.Constraint()
    def rsofc_total_TPC_eq(c):
        return c.rsofc_total_TPC == (
            c.sofc_fs.costing.total_plant_cost + c.soec_fs.costing.total_TPC)

    calculate_variable_from_constraint(
        m.rsofc_total_TPC,
        m.rsofc_total_TPC_eq
    )

    m.rsofc_total_fixed_OM_cost = pyo.Var(
        initialize=1e-6,
        doc="rsofc total plant cost in $MM/yr",
    )

    @m.Constraint()
    def rsofc_total_fixed_OM_cost_eq(c):
        return c.rsofc_total_fixed_OM_cost == (
            c.sofc_fs.costing.total_fixed_OM_cost +
            c.soec_fs.costing.total_fixed_OM_cost)

    calculate_variable_from_constraint(
        m.rsofc_total_fixed_OM_cost,
        m.rsofc_total_fixed_OM_cost_eq
    )

    # Display results
    # m.fs.costing.display()
    m.sofc_fs.costing.display()
    m.soec_fs.H2_costing.display()

    return m


# %%       # ------------------------------------------------------------------
if __name__ == "__main__":
    m = get_model()
    cost_rsofc(m)
