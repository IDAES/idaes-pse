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
"""Generate parameter and expression files for NH3
"""

__author__ = "Stephen Burroughs"

import os
import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def thermal_conductivity_rule(m):
    """
    Thermal conductivity rule

    Monogenidou, S. A., Marc J. Assael, and Marcia L. Huber.
    "Reference correlation for the thermal conductivity of ammonia
    from the triple-point temperature to 680 K and pressures up to 80 MPa."
    Journal of Physical and Chemical Reference Data 47.4 (2018).
    """
    b = {
        1: 0.07152,
        2: 130.228,
        3: 9569.817,
    }
    u = {
        1: 1646,
        2: 3965,
        3: 7231,
    }
    v = {
        1: 2.224,
        2: 3.148,
        3: 0.9579,
    }
    B = {
        1: {
            1: 0.103432e0,
            2: -0.112597e0,
            3: 0.233301e0,
            4: -0.112536e0,
            5: 0.141129e-1,
        },
        2: {
            1: -0.283976e-1,
            2: 0.482520e-1,
            3: -0.644124e-1,
            4: 0.529376e-2,
            5: 0.891203e-2,
        },
    }
    T = m.T_star / m.tau
    # Ts = T / 251.196
    rho = m.rho_star * m.delta
    G = sum(bval / T**i for i, bval in b.items())
    cint_over_k = 4.0 + sum(
        v[i] * (u[i] / T) ** 2 * pyo.exp(u[i] / T) / (pyo.exp(u[i] / T) - 1) ** 2
        for i in range(1, len(u.items()) + 1)
    )
    return 0.1351767 * pyo.sqrt(T) * cint_over_k / G + sum(
        (B[1][i] + B[2][i] * T) * rho**i for i in range(1, 6)
    )


def viscosity_rule(m):
    """
    Viscosity rule

    Fenghour, A., et al. "The viscosity of ammonia."
    Journal of Physical and Chemical Reference Data 24.5 (1995): 1649-1667.
    """
    a = {
        0: 4.99318220,
        1: -0.61122364,
        2: 0,
        3: 0.18535124,
        4: -0.11160946,
    }
    c = {
        0: -0.17999496e1,
        1: 0.46692621e2,
        2: -0.53460794e3,
        3: 0.33604074e4,
        4: -0.13019164e5,
        5: 0.33414230e5,
        6: -0.58711743e5,
        7: 0.71426686e5,
        8: -0.59834012e5,
        9: 0.33652741e5,
        10: -0.12027350e5,
        11: 0.24348205e4,
        12: -0.20807957e3,
    }
    d = {
        2: {
            0: 0,
            1: 0,
            2: 2.19664285e-1,
            3: 0,
            4: -0.83651107e-1,
        },
        3: {
            0: 0.17366936e-2,
            1: -0.64250359e-2,
            2: 0,
            3: 0,
            4: 0,
        },
        4: {
            0: 0,
            1: 0,
            2: 1.67668649e-4,
            3: -1.49710093e-4,
            4: 0.77012274e-4,
        },
    }

    T = m.T_star / m.tau
    rho = m.delta * m.rho_star
    Ts = T / 386
    return (
        0.021357
        * pyo.sqrt(17.03026 * T)
        / (0.2957**2 * pyo.exp(sum(aval * pyo.log(Ts) ** i for i, aval in a.items())))
        + 0.021357
        * pyo.sqrt(17.03026 * T)
        / (0.2957**2 * pyo.exp(sum(aval * pyo.log(Ts) ** i for i, aval in a.items())))
        * sum((c[i] * pyo.sqrt(Ts) ** -i) for i in range(0, len(c.items())))
        / (0.6022137 * 0.2957**3)
        * rho
        + sum(
            sum(d[i][j] / Ts**j for j in range(0, len(d[i].items()))) * rho**i
            for i in range(2, 5)
        )
    )


def main(dry_run=False):
    """Generate parameter and expression files.

    Args:
        dry_run (bool): If dry run don't generate files

    Returns:
        None
    """
    main_param_file = os.path.join(this_file_dir(), "nh3.json")
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
