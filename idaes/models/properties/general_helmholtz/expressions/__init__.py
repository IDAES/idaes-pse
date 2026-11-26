#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

__author__ = "John Eslick"

from .surface_tension_type01 import surface_tension_type01
from .sat_delta_approx import sat_delta_type01, sat_delta_type02, sat_delta_type03
from .phi_ideal_expressions import (
    phi_ideal_expressions_lead,
    phi_ideal_expressions_logtau,
    phi_ideal_expressions_planck_einstein1,
    phi_ideal_expressions_planck_einstein2,
    phi_ideal_expressions_planck_einstein3,
    phi_ideal_expressions_cp_constant,
    phi_ideal_expressions_power,
    phi_ideal_expressions_enth_entr_offset,
    phi_ideal_expressions_GERG_Cosh,
    phi_ideal_expressions_GERG_Sinh,
)
from .phi_residual_expressions import (
    phi_residual_expressions_exponential_reduced_density,
    phi_residual_expressions_exponential_delta_tau,
    phi_residual_expressions_gaussian,
    phi_residual_expressions_power,
    phi_residual_expressions_double_exponential,
    phi_residual_expressions_gaob,
    phi_residual_expressions_gaussian_GERG2008,
)
from .phi_ideal_type01 import phi_ideal_expressions_type01
from .phi_ideal_type02 import phi_ideal_expressions_type02
from .phi_ideal_type03 import phi_ideal_expressions_type03
from .phi_ideal_type04 import phi_ideal_expressions_type04
from .phi_residual_type01 import phi_residual_expressions_type01
from .phi_residual_type02 import phi_residual_expressions_type02
from .phi_residual_type03 import phi_residual_expressions_type03
from .phi_residual_type04 import phi_residual_expressions_type04
from .phi_residual_type05 import phi_residual_expressions_type05

phi_residual_types = {
    0: None,  # custom
    1: phi_residual_expressions_type01,
    2: phi_residual_expressions_type02,
    3: phi_residual_expressions_type03,
    4: phi_residual_expressions_type04,
    5: phi_residual_expressions_type05,
}

phi_ideal_types = {
    0: None,  # custom
    1: phi_ideal_expressions_type01,
    2: phi_ideal_expressions_type02,
    3: phi_ideal_expressions_type03,
    4: phi_ideal_expressions_type04,
}

phi_ideal_modular_parts = {
    0: None,  # custom
    1: phi_ideal_expressions_lead,
    2: phi_ideal_expressions_logtau,
    3: phi_ideal_expressions_planck_einstein1,
    4: phi_ideal_expressions_planck_einstein2,
    5: phi_ideal_expressions_planck_einstein3,
    6: phi_ideal_expressions_cp_constant,
    7: phi_ideal_expressions_power,
    8: phi_ideal_expressions_enth_entr_offset,
    9: phi_ideal_expressions_GERG_Cosh,
    10: phi_ideal_expressions_GERG_Sinh,
}

phi_residual_modular_parts = {
    0: None,  # custom
    1: phi_residual_expressions_power,
    2: phi_residual_expressions_gaussian,
    3: phi_residual_expressions_gaussian_GERG2008,
    4: phi_residual_expressions_gaob,
    5: phi_residual_expressions_exponential_delta_tau,
    6: phi_residual_expressions_double_exponential,
    7: phi_residual_expressions_exponential_reduced_density,
}

delta_sat_types = {
    0: None,  # custom
    1: sat_delta_type01,
    2: sat_delta_type02,
    3: sat_delta_type03,
}

surface_tension_types = {
    0: None,
    1: surface_tension_type01,
}
