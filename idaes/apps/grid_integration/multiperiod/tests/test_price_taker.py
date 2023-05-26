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


from pyomo.environ import (
    Var,
    Binary,
    NonNegativeReals,
    Constraint,
    Expression,
)

from idaes.apps.grid_integration import DesignModel, OperationModel


def dfc_design(m, params, capacity_range=(650, 900)):
    _dfc_capacity = params["dfc_capacity"]
    _ng_flow = params["ng_flow"]
    _capex = params["capex"]
    _fom_factor = params["fom_factor"]

    m.capacity = Var(
        within=NonNegativeReals,
        initialize=_dfc_capacity,
        bounds=(0, capacity_range[1]),
        doc="Capacity of the power plant [in MW]",
    )

    # Define a variable that informs whether the plant needs to built or not. 
    m.build_unit = Var(
        within=Binary,
        doc="1: Plant is built, 0: Plant is not built",
    )

    # Bound the capcity of the plant in the specified range
    m.capacity_lb_con = Constraint(expr=m.capacity >= m.build_unit * capacity_range[0])
    m.capacity_ub_con = Constraint(expr=m.capacity <= m.build_unit * capacity_range[1])

    # Compute the natural gas flowrate required at maximum capacity. 
    m.ng_flow = Expression(
        expr=m.capacity * (_ng_flow / _dfc_capacity),
        doc="Computes the natural flowrate required [in kg/s] at full load",
    )

    # Define an expression for CAPEX and FOM
    m.capex = Expression(
        expr=_capex * (m.capacity / _dfc_capacity),
        doc="CAPEX of the power cycle [in 1000$]",
    )
    m.fom = Expression(
        expr=_fom_factor * m.capex,
        doc="Fixed O&M cost [in 1000$/year]",
    )