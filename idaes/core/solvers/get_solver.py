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
This module contains the IDAES get_solver method.
"""

import idaes.logger as idaeslog
import idaes.core.solvers

_log = idaeslog.getLogger(__name__)


# Author: Andrew Lee
def get_solver(solver=None, options=None):
    """
    General method for getting a solver object which defaults to the standard
    IDAES solver (defined in the IDAES configuration).

    Args:
        solver: string name for desired solver. Default=None, use default solver
        options: dict of solver options to use, overwrites any settings
                 provided by IDAES configuration. Default = None, use default
                 solver options.

    Returns:
        A Pyomo solver object
    """
    if solver is None:
        solver = "default"
    solver_obj = idaes.core.solvers.SolverWrapper(solver, register=False)()

    if options is not None:
        solver_obj.options.update(options)

    return solver_obj
