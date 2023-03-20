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
Class for users who deliberately do not want to assign a method of calculation for a property in a particular phase.
Upon constructing a parameter block, a warning is logged that no method is assigned for a property-phase combination.
When the Expression is constructed, an Expression.Skip is returned for this phase.
"""

from pyomo.environ import Expression
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class NoMethod(object):
    """Class to pass when user wants to deliberately not assign a property to a given phase"""

    class visc_d_phase(object):
        """No method for dynamic viscosity"""

        @staticmethod
        def build_parameters(pobj):
            """Logs that property is not being constructed. We only need one log for the property package"""
            _log.warning(
                f"Skipping construction of dynamic viscosity for phase {pobj.local_name}"
            )

        @staticmethod
        def return_expression(b, p):
            """All we need is the expression skip"""
            return Expression.Skip

    class therm_cond_phase(object):
        """No method for thermal conductivity"""

        @staticmethod
        def build_parameters(pobj):
            """Logs that property is not being constructed. We only need one log for the property package"""
            _log.warning(
                f"Skipping construction of thermal conductivity for phase {pobj.local_name}"
            )

        @staticmethod
        def return_expression(b, p):
            """All we need is the expression skip"""
            return Expression.Skip
