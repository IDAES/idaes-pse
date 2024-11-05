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
API for accessing core IDAES datasets

Usage, e.g., for the Pitzer(1984) data::

    from idaes.core.datasets import Pitzer
    pitzer = Pitzer()
    gibbs_data = pitzer.get_table("Standard G").data

"""
__authors__ = ["Dan Gunter (LBNL)"]
__author__ = __authors__[0]

from pyomo.common.deprecation import deprecation_warning


deprecation_warning(
    msg=(
        "idaes.core.datasets is no longer available due to the DMF being unsupported."
        " Visit https://github.com/IDAES/dmf for more information."
    ),
    version="2.6.0dev0",
    remove_in="2.6.0rc0",
)
