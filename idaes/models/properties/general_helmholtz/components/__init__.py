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
"""This module provides a list of suppored component and Pyomo expressions for
some properties not implimented as external functions.
"""

__author__ = "John Eslick"

import idaes.models.properties.general_helmholtz.components.h2o as h2o
import idaes.models.properties.general_helmholtz.components.co2 as co2
import idaes.models.properties.general_helmholtz.components.r134a as r134a
import idaes.models.properties.general_helmholtz.components.r1234ze as r1234ze

components = {
    "h2o": h2o,
    "co2": co2,
    "r134a": r134a,
    "r1234ze": r1234ze,
}
