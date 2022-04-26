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
Deprecation path for relocated unit_model module.
"""
from pyomo.common.deprecation import relocated_module_attribute


relocated_module_attribute(
    "UnitModelBlock",
    "idaes.core.base.unit_model.UnitModelBlock",
    version="2.0.0.alpha0",
)

relocated_module_attribute(
    "UnitModelBlockData",
    "idaes.core.base.unit_model.UnitModelBlockData",
    version="2.0.0.alpha0",
)
