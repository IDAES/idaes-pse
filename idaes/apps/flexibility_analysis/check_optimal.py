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
This module provides a function like Pyomo's assert_optimal_termination but that 
works for both APPSI solver interfaces and non-appsi solver interfaces.
"""
import pyomo.environ as pe
from pyomo.contrib import appsi


def assert_optimal_termination(results):
    """
    Raise an exception if the termination condition was not optimal.

    Parameters
    ----------
    results: pyomo results object from calling solve()
    """
    if hasattr(results, "termination_condition"):
        assert results.termination_condition == appsi.base.TerminationCondition.optimal
    else:
        pe.assert_optimal_termination(results)
