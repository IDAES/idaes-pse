#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
IDAES Data Management Framework (DMF)

The DMF lets you save, search, and retrieve provenance related
to your models.
"""
__author__ = "Dan Gunter"

from pyomo.common.deprecation import deprecation_warning


deprecation_warning(
    msg=(
        "idaes.core.dmf is no longer supported and has been moved into a dedicated repository for archival."
        " Visit https://github.com/IDAES/dmf for more information."
    ),
    version="2.6.0dev0",
    remove_in="2.6.0rc0",
)
