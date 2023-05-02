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
"""Generate parameter and expression files for r134a
"""

__author__ = "John Eslick"

import os
import math
from pyomo.common.fileutils import this_file_dir
import pyomo.environ as pyo
from idaes.core.util.math import smooth_max
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def thermal_conductivity_rule(m):
    """Thermal Conductivity Rule

    Perkins, R.A.; Laesecke, A.; Howley, J.; Ramires, M.L.V.; Gurova, A.N.; Cusco, L.,
        Experimental thermal conductivity values for the IUPAC round-robin sample of
        1,1,1,2-tetrafluoroethane (R134a), NIST Interagency/Internal Report (NISTIR)
        - 6605, 2000, https://doi.org/10.6028/NIST.IR.6605.
    """

    a = [
        -1.05248e-2,
        8.00982e-5,
    ]
    b = {
        1: 1.836526,
        2: 5.126143,
        3: -1.436883,
        4: 6.261441e-1,
    }

    T = m.T_star / m.tau
    MW = m.MW / 1000.0  # kg/mol
    rho_star = m.rho_star / MW
    rho = m.delta * rho_star  # mol/m^3
    tau = m.tau
    delta = m.delta
    rho_crit = 5049.886  # mol/m^3
    lamb_red = 2.055e-3  # reducing tc W/m/K
    Tref = 561.411  # reference T K
    k = 1.380649e-23
    Pc = 4.05928  # MPa
    big_gam = 0.0496
    R0 = 1.03
    gamma = 1.239
    qd = 1892020000.0
    xi0 = 1.94e-10  # m
    nu = 0.63

    m.cp = pyo.ExternalFunction(library="", function="cp")
    m.cv = pyo.ExternalFunction(library="", function="cv")
    m.mu = pyo.ExternalFunction(library="", function="mu")
    m.itc = pyo.ExternalFunction(library="", function="itc")

    drho_dp = m.itc("r134a", delta, m.T_star / T) * rho_star * m.delta
    drho_dp_ref = m.itc("r134a", delta, m.T_star / Tref) * rho_star * m.delta

    deltchi = smooth_max(
        Pc * rho / rho_crit**2 * (drho_dp - drho_dp_ref * Tref / T), 0, eps=1e-8
    )
    xi = xi0 * (deltchi / big_gam) ** (nu / gamma)

    cp = m.cp("r134a", delta, tau)
    cv = m.cv("r134a", delta, tau)
    mu = m.mu("r134a", delta, tau) / 1e6
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
    lambda_dg = a[0] + a[1] * T
    lambda_r = lamb_red * sum(bi * (rho / rho_crit) ** i for i, bi in b.items())
    return 1000 * (lambda_dg + lambda_r + lambda_c)


def viscosity_rule(m):
    """Viscosity Rule

    Perkins, R.A.; Laesecke, A.; Howley, J.; Ramires, M.L.V.; Gurova, A.N.; Cusco, L.,
        Experimental thermal conductivity values for the IUPAC round-robin sample of
        1,1,1,2-tetrafluoroethane (R134a), NIST Interagency/Internal Report (NISTIR)
        - 6605, 2000, https://doi.org/10.6028/NIST.IR.6605.
    """
    a = [
        0.355404,
        -0.464337,
        0.257353e-1,
    ]
    b = [
        -19.572881,
        219.73999,
        -1015.3226,
        2471.01251,
        -3375.1717,
        2491.6597,
        -787.26086,
        14.085455,
        -0.34664158,
    ]
    te = [
        0.0,
        -0.25,
        -0.50,
        -0.75,
        -1.00,
        -1.25,
        -1.50,
        -2.50,
        -5.50,
    ]
    c = {
        1: -20.6900719,
        2: 0.356029549,
        3: 2.11101816,
        4: 13.9601415,
        5: -4.5643502,
        6: -3.51593275,
        7: 214.76332,
        8: -0.890173375e-1,
        9: 0.100035295,
        10: 3.163695636,
    }
    T = m.T_star / m.tau
    rho = m.delta * m.rho_star / m.MW * 1000
    M = 102.031
    sigma = 0.46893
    eok = 299.363
    NA = 6.0221408e23
    Ts = T / eok
    vs = pyo.exp(sum(ai * pyo.log(Ts) ** i for i, ai in enumerate(a)))
    Bs = sum(bi * Ts**ti for bi, ti in zip(b, te))
    B = NA * sigma**3 * Bs / 1e9**3
    etas = 0.021357 * pyo.sqrt(M * T) / (sigma**2 * vs)
    tau = T / m.Tc
    delta0 = c[10] / (1 + c[8] * tau + c[9] * tau**2)
    delta = rho / 5017.053
    eta = (
        c[1] * delta
        + (c[2] / tau**6 + c[3] / tau**2 + c[4] / tau**0.5 + c[5] * tau**2)
        * delta**2
        + c[6] * delta**3
        + c[7] / (delta0 - delta)
        - c[7] / delta0
    )
    return etas * (1 + B * rho) + eta


def main(dry_run=False):
    """Generate parameter and expression files.

    Args:
        dry_run (bool): If dry run don't generate files

    Returns:
        None
    """
    main_param_file = os.path.join(this_file_dir(), "r134a.json")
    we = WriteParameters(parameters=main_param_file)
    we.add(
        {
            "viscosity": viscosity_rule,
            "thermal_conductivity": thermal_conductivity_rule,
        }
    )
    we.write(dry_run=dry_run)

    print("ASHRAE Offset")
    print(we.calculate_reference_offset(2.79075439914, 1.60488955608, 0, 0))
    return we


if __name__ == "__main__":
    main()
