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
Initializer class for implementing initialization from a data source
"""
from pyomo.environ import check_optimal_termination, Constraint, Var
from pyomo.common.config import ConfigValue
from pyomo.contrib.incidence_analysis.util import solve_strongly_connected_components

from idaes.core.initialization.initializer_base import InitializerBase
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

__author__ = "Andrew Lee"


class FromDataInitializer(InitializerBase):
    def initialization_routine(self, model):
        # No action required
        pass
