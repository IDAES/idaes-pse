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

from pyomo.solvers.plugins.solvers.IPOPT import IPOPT
from pyomo.common import Executable
from pyomo.opt.base.solvers import SolverFactory
from pyomo.opt.solver import SystemCallSolver

import logging

logger = logging.getLogger("pyomo.solvers")


@SolverFactory.register("ipopt_l1", doc="The Ipopt NLP solver")
class IPOPT_L1(IPOPT):
    def _default_executable(self):
        executable = Executable("ipopt_l1")
        if not executable:
            logger.warning(
                "Could not locate the 'ipopt_l1' executable, "
                "which is required for solver %s" % self.name
            )
            self.enable = False
            return None
        return executable.path()
