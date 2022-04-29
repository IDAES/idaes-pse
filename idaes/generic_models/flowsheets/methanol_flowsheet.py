##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Deprecation path for renamed module.
"""
from pyomo.common.deprecation import deprecation_warning

deprecation_warning("The generic_models.flowsheets.methanol_flowsheet module "
                    "has been moved to idaes.models.flowsheets.methanol_flowsheet",
                    version="2.0.0.alpha0")

from idaes.models.flowsheets.methanol_flowsheet import *
