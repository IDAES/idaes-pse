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
Deprecation path for relocated Components module.
"""
from pyomo.common.deprecation import relocated_module_attribute


relocated_module_attribute(
    "Component", "idaes.core.base.components.Component", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "ComponentData", "idaes.core.base.components.ComponentData", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "Solute", "idaes.core.base.components.Solute", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "Solvent", "idaes.core.base.components.Solvent", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "Ion", "idaes.core.base.components.Ion", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "Anion", "idaes.core.base.components.Anion", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "Cation", "idaes.core.base.components.Cation", version="2.0.0.alpha0"
)

relocated_module_attribute(
    "Apparent", "idaes.core.base.components.Apparent", version="2.0.0.alpha0"
)
