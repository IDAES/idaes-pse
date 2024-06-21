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

"""
This module contains the IDAES get_solver method.
"""
from pyomo.contrib.solver.base import LegacySolverWrapper

import idaes.logger as idaeslog
import idaes.core.solvers

_log = idaeslog.getLogger(__name__)


# Author: Andrew Lee
def get_solver(
    solver=None,
    solver_options: dict = None,
    writer_config: dict = None,
    options: dict = None,
):
    """
    General method for getting a solver object which defaults to the standard
    IDAES solver (defined in the IDAES configuration).

    Args:
        solver: string name for desired solver. Default=None, use default solver
        solver_options: dict of solver options to use, overwrites any settings
                 provided by IDAES configuration. Default = None, use default
                 solver options.
        writer_config: dict of configuration options for solver writer, overwrites
                 ny settings provided by IDAES configuration. Default = None, use
                 default solver options.
        options: DEPRECATED. Alias of solver_options.

    Returns:
        A Pyomo solver object
    """
    if solver_options is not None:
        if options is not None:
            raise ValueError(
                "Cannot provide both the 'options' and 'solver_options' argument. "
                "'options' has been deprecated in favor of 'solver_options'."
            )
        options = solver_options

    if solver is None:
        solver = "default"
    solver_obj = idaes.core.solvers.SolverWrapper(solver, register=False)()

    if isinstance(solver_obj, LegacySolverWrapper):
        # New solver interface.
        # LegacySolverWrapper is a wrapper for the new solver interface that makes it
        # backward compatible.
        if options is not None:
            for k, v in options.items():
                solver_obj.options[k] = v
        if writer_config is not None:
            for k, v in writer_config.items():
                solver_obj.config.writer_config[k] = v
    else:
        # Old solver interface
        if options is not None:
            solver_obj.options.update(options)
        if writer_config is not None:
            _log.info(
                "Older Pyomo solver interface does not support writer_config argument: ignoring."
            )

    return solver_obj
