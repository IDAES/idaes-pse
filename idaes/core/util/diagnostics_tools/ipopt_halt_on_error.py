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
This module contains a helper function for running IPOPT in a halt-on-error mode.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"


from pyomo.environ import (
    SolverFactory,
)


def ipopt_solve_halt_on_error(model, options=None):
    """
    Run IPOPT to solve model with debugging output enabled.

    This function calls IPOPT to solve the model provided with settings
    to halt on AMPL evaluation errors and report these with symbolic names.

    Args:
        model: Pyomo model to be solved.
        options: solver options to be passed to IPOPT

    Returns:
        Pyomo solver results dict

    """
    if options is None:
        options = {}

    solver = SolverFactory("ipopt")
    solver.options = options
    solver.options["halt_on_ampl_error"] = "yes"

    return solver.solve(
        model, tee=True, symbolic_solver_labels=True, export_defined_variables=False
    )
