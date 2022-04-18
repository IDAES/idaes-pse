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
Deprecation path for relocated phases module.
"""
from pyomo.common.deprecation import relocated_module_attribute


relocated_module_attribute(
    "PhaseType", "idaes.core.base.phases.PhaseType", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "Phase", "idaes.core.base.phases.Phase", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "PhaseData", "idaes.core.base.phases.PhaseData", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "LiquidPhase", "idaes.core.base.phases.LiquidPhase", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "SolidPhase", "idaes.core.base.phases.SolidPhase", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "VaporPhase", "idaes.core.base.phases.VaporPhase", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "AqueousPhase", "idaes.core.base.phases.AqueousPhase", version="2.0.0.alpha0"
)
