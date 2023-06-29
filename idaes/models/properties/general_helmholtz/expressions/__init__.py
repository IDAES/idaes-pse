#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

__author__ = "John Eslick"

from .phi_ideal_type01 import phi_ideal_expressions_type01
from .phi_ideal_type02 import phi_ideal_expressions_type02
from .phi_ideal_type03 import phi_ideal_expressions_type03
from .phi_residual_type01 import phi_residual_expressions_type01
from .phi_residual_type02 import phi_residual_expressions_type02
from .phi_residual_type03 import phi_residual_expressions_type03
from .phi_residual_type04 import phi_residual_expressions_type04
from .surface_tension_type01 import surface_tension_type01
from .sat_delta_approx import sat_delta_type01, sat_delta_type02

phi_residual_types = {
    0: None,  # custom
    1: phi_residual_expressions_type01,
    2: phi_residual_expressions_type02,
    3: phi_residual_expressions_type03,
    4: phi_residual_expressions_type04,
}

phi_ideal_types = {
    0: None,  # custom
    1: phi_ideal_expressions_type01,
    2: phi_ideal_expressions_type02,
    3: phi_ideal_expressions_type03,
}

delta_sat_types = {
    0: None,  # custom
    1: sat_delta_type01,
    2: sat_delta_type02,
}

surface_tension_types = {
    0: None,
    1: surface_tension_type01,
}
