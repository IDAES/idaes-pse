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
Deprecation path for renamed module.
"""
from pyomo.common.deprecation import deprecation_warning

deprecation_warning("The generic_models.properties.iapws95 has been "
                    "moved to idaes.models.properties.iapws95",
                    version="2.0.0.alpha0")

from idaes.models.properties.iapws95 import *
