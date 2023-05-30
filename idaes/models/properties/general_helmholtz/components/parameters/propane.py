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
"""Generate parameter and expression files for propane
"""

__author__ = "John Eslick"

import os
import math
import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)
from idaes.core.util.math import smooth_max
from idaes.core.util.constants import Constants


def thermal_conductivity_rule(m):
    """
    Thermal conductivity expression rule.  See:

    Marsh, K., Perkins, R., and Ramires, M.L.V., "Measurement and Correlation
        of the Thermal Conductivity of Propane from 86 to 600 K at Pressures
        to 70 MPa, J. Chem. Eng. Data, 47(4):932-940, 2002.
    """
    a = [
        -0.00124778,
        0.00816371,
        0.0199374,
    ]
    b = {
        (1, 1): -0.03695,
        (1, 2): 0.0482798,
        (2, 1): 0.148658,
        (2, 2): -0.135636,
        (3, 1): -0.119986,
        (3, 2): 0.117588,
        (4, 1): 0.0412431,
        (4, 2): -0.0436911,
        (5, 1): -0.00486905,
        (5, 2): 0.00616079,
    }

    T = m.T_star / m.tau
    MW = m.MW / 1000.0  # kg/mol
    rho_star = m.rho_star / MW
    rho = m.delta * rho_star  # mol/m^3
    tau = m.tau
    delta = m.delta
    rho_crit = 5000  # mol/m^3
    # lamb_red = 1  # reducing tc W/m/K
    Tref = 554.73  # reference T K
    k = Constants.boltzmann_constant
    Pc = 4.2512  # MPa
    big_gam = 0.0496
    R0 = 1.03
    gamma = 1.239
    qd = 1 / 7.16635e-10
    xi0 = 1.94e-10  # m
    nu = 0.63
    Tred = 1.0 / m.tau
    l0 = sum(ai * Tred**i for i, ai in enumerate(a))
    lr = sum(m.delta**i * (b[i, 1] + b[i, 2] * Tred) for i in range(1, 6))

    m.cp = pyo.ExternalFunction(library="", function="cp")
    m.cv = pyo.ExternalFunction(library="", function="cv")
    m.mu = pyo.ExternalFunction(library="", function="mu")
    m.itc = pyo.ExternalFunction(library="", function="itc")

    drho_dp = m.itc("propane", delta, m.T_star / T) * rho_star * m.delta
    drho_dp_ref = m.itc("propane", delta, m.T_star / Tref) * rho_star * m.delta

    deltchi = smooth_max(
        Pc * rho / rho_crit**2 * (drho_dp - drho_dp_ref * Tref / T), 0, eps=1e-8
    )
    xi = xi0 * (deltchi / big_gam) ** (nu / gamma)

    cp = m.cp("propane", delta, tau)
    cv = m.cv("propane", delta, tau)
    mu = m.mu("propane", delta, tau) / 1e6
    y = qd * xi
    kappa_inv = cv / cp
    pi = math.pi
    Omega = 2.0 / pi * ((1.0 - kappa_inv) * pyo.atan(y) + kappa_inv * y)
    Omega0 = (
        2.0
        / pi
        * (1.0 - pyo.exp(-1.0 / (1.0 / y + 1.0 / 3.0 * (y * rho_crit / rho) ** 2)))
    )
    lambda_c = (
        1000 * MW * cp * rho * R0 * k * T / 6.0 / math.pi / mu / xi * (Omega - Omega0)
    )
    return 1000 * (l0 + lr + lambda_c)


def viscosity_rule(m):
    """
    Viscosity expression rule.  See:

    Vogel E, C Küchenmeister, E Bich, A Laesecke, Reference Correlation of the
        Viscosity of Propane. Journal of Physical and Chemical Reference Data
        27, 947–970 (1998)
    """
    a = {
        0: 0.25104574,
        1: -0.47271238,
        3: 0.060836515,
    }
    b = {
        0: -19.572881,
        1: 219.73999,
        2: -1015.3226,
        3: 2471.01251,
        4: -3375.1717,
        5: 2491.6597,
        6: -787.26086,
        7: 14.085455,
        8: -0.34664158,
    }
    t = {
        0: 0.0,
        1: -0.25,
        2: -0.50,
        3: -0.75,
        4: -1.00,
        5: -1.25,
        6: -1.50,
        7: -2.50,
        8: -5.50,
    }
    e = {
        (2, 0): 35.9873030195,
        (2, 1): -180.512188564,
        (2, 2): 87.7124888223,
        (3, 0): -105.773052525,
        (3, 1): 205.319740877,
        (3, 2): -129.210932610,
        (4, 0): 58.9491587759,
        (4, 1): -129.740033100,
        (4, 2): 76.6280419971,
        (5, 0): -9.59407868475,
        (5, 1): 21.0726986598,
        (5, 2): -14.3971968187,
    }
    f = {1: 1616.88405374}
    g = {1: 2.50053938863, 2: 0.860516059264}

    M = m.MW  # kmol/kg
    T = m.T_star / m.tau  # K
    tau = T / m.Tc
    rho = m.delta * m.rho_star / M  # mol/l
    sigma = 0.49748  # nm
    eok = 263.88  # K
    NA = 6.0221408e23  # number/mol
    Ts = T / eok  # dimensionless
    vs = pyo.exp(sum(ai * pyo.log(Ts) ** i for i, ai in a.items()))
    etas = 0.021357 * pyo.sqrt(M * T) / (sigma**2 * vs)

    Bs = sum(bi * Ts ** t[i] for i, bi in b.items())
    B = NA * sigma**3 * Bs * 1e-24  # l/mol
    delta0 = g[1] * (1 + g[2] * pyo.sqrt(tau))
    delta = rho / 5.0
    eta = sum(e[i, j] * delta**i / tau**j for i, j in e) + f[1] * (
        delta / (delta0 - delta) - delta / delta0
    )
    return etas * (1 + B * rho) + eta


def main(dry_run=False):
    """Generate parameter and expression files.

    Args:
        dry_run (bool): If dry run don't generate files

    Returns:
        None
    """
    main_param_file = os.path.join(this_file_dir(), "propane.json")
    we = WriteParameters(parameters=main_param_file)
    we.add(
        {
            "viscosity": viscosity_rule,
            "thermal_conductivity": thermal_conductivity_rule,
        }
    )
    we.write(dry_run=dry_run)
    return we


if __name__ == "__main__":
    main()
