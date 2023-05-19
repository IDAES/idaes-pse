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
"""Generate parameter and expression files for water
"""

__author__ = "John Eslick"

import os
import math
import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir
from pyomo.common.fileutils import find_library
from idaes.core.util.math import smooth_max
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def thermal_conductivity_rule(m):
    """Thermal conductivity rule

    International Association for the Properties of Water and Steam (2011).
        IAPWS R15-11, "Release on the IAPWS Formulation 2011 for the
        Thermal Conductivity of Ordinary Water Substance,"
        URL: http://iapws.org/relguide/ThCond.pdf
    """

    L0 = {
        0: 2.443221e-3,
        1: 1.323095e-2,
        2: 6.770357e-3,
        3: -3.454586e-3,
        4: 4.096266e-4,
    }

    L1 = {
        (0, 0): 1.60397357,
        (1, 0): 2.33771842,
        (2, 0): 2.19650529,
        (3, 0): -1.21051378,
        (4, 0): -2.7203370,
        (0, 1): -0.646013523,
        (1, 1): -2.78843778,
        (2, 1): -4.54580785,
        (3, 1): 1.60812989,
        (4, 1): 4.57586331,
        (0, 2): 0.111443906,
        (1, 2): 1.53616167,
        (2, 2): 3.55777244,
        (3, 2): -0.621178141,
        (4, 2): -3.18369245,
        (0, 3): 0.102997357,
        (1, 3): -0.463045512,
        (2, 3): -1.40944978,
        (3, 3): 0.0716373224,
        (4, 3): 1.1168348,
        (0, 4): -0.0504123634,
        (1, 4): 0.0832827019,
        (2, 4): 0.275418278,
        (3, 4): 0.0,
        (4, 4): -0.19268305,
        (0, 5): 0.00609859258,
        (1, 5): -0.00719201245,
        (2, 5): -0.0205938816,
        (3, 5): 0.0,
        (4, 5): 0.012913842,
    }
    delta = m.delta
    tau = m.tau
    big_lam = 177.8514
    qdb = 1.0 / 0.40
    nu = 0.630
    gamma = 1.239
    xi0 = 0.13
    big_gam0 = 0.06
    Tbr = 1.5
    Tb = 1 / tau
    lib_path = find_library("general_helmholtz_external")
    m.cp = pyo.ExternalFunction(library=lib_path, function="cp")
    m.cv = pyo.ExternalFunction(library=lib_path, function="cv")
    m.mu = pyo.ExternalFunction(library=lib_path, function="mu")
    m.itc = pyo.ExternalFunction(
        library=lib_path, function="itc"
    )  # isothermal compressibility [1/MPa]
    deltchi = smooth_max(
        delta**2
        * (m.itc("h2o", delta, tau) - m.itc("h2o", delta, 1 / Tbr) * Tbr / Tb)
        * m.Pc
        / 1000,
        0,
        eps=1e-8,
    )
    xi = xi0 * (deltchi / big_gam0) ** (nu / gamma)
    y = qdb * xi
    kappa = m.cp("h2o", delta, tau) / m.cv("h2o", delta, tau)
    Z = (
        2.0
        / math.pi
        / y
        * (
            ((1 - 1.0 / kappa) * pyo.atan(y) + 1.0 / kappa * y)
            - (1 - pyo.exp(-1 / (1.0 / y + y**2 / 3 / delta**2)))
        )
    )
    cpb = m.cp("h2o", delta, tau) / m.R
    mub = m.mu("h2o", delta, tau)

    lam2 = big_lam * delta * Tb * cpb / mub * Z
    return lam2 + pyo.sqrt(1.0 / m.tau) / sum(
        L0val * m.tau**i for i, L0val in L0.items()
    ) * pyo.exp(
        m.delta
        * sum(
            (m.tau - 1) ** i * sum(L1[i, j] * (m.delta - 1) ** j for j in range(0, 6))
            for i in range(0, 5)
        )
    )


def viscosity_rule(mvisc):
    """Viscosity calculation for water

    International Association for the Properties of Water and Steam (2008).
        IAPWS R12-08, "Release on the IAPWS Formulation 2008 for the Viscosity
        of Ordinary Water Substance, URL: http://iapws.org/relguide/visc.pdf
    """

    H0 = {0: 1.67752, 1: 2.20462, 2: 0.6366564, 3: -0.241605}
    H1 = {
        (0, 0): 5.20094e-1,
        (1, 0): 8.50895e-2,
        (2, 0): -1.08374,
        (3, 0): -2.89555e-1,
        (4, 0): 0.0,
        (5, 0): 0.0,
        (0, 1): 2.22531e-1,
        (1, 1): 9.99115e-1,
        (2, 1): 1.88797,
        (3, 1): 1.26613,
        (4, 1): 0.0,
        (5, 1): 1.20573e-1,
        (0, 2): -2.81378e-1,
        (1, 2): -9.06851e-1,
        (2, 2): -7.72479e-1,
        (3, 2): -4.89837e-1,
        (4, 2): -2.57040e-1,
        (5, 2): 0.0,
        (0, 3): 1.61913e-1,
        (1, 3): 2.57399e-1,
        (2, 3): 0.0,
        (3, 3): 0.0,
        (4, 3): 0.0,
        (5, 3): 0.0,
        (0, 4): -3.25372e-2,
        (1, 4): 0.0,
        (2, 4): 0.0,
        (3, 4): 6.98452e-2,
        (4, 4): 0.0,
        (5, 4): 0.0,
        (0, 5): 0.0,
        (1, 5): 0.0,
        (2, 5): 0.0,
        (3, 5): 0.0,
        (4, 5): 8.72102e-3,
        (5, 5): 0.0,
        (0, 6): 0.0,
        (1, 6): 0.0,
        (2, 6): 0.0,
        (3, 6): -4.35673e-3,
        (4, 6): 0.0,
        (5, 6): -5.93264e-4,
    }

    return (
        1e2
        * pyo.sqrt(1.0 / mvisc.tau)
        / sum(H0val * mvisc.tau**i for i, H0val in H0.items())
        * pyo.exp(
            mvisc.delta
            * sum(
                (mvisc.tau - 1) ** i
                * sum(H1[i, j] * (mvisc.delta - 1) ** j for j in range(0, 7))
                for i in range(0, 6)
            )
        )
    )


def main(dry_run=False):
    """Generate parameter and expression files.

    Args:
        dry_run (bool): If dry run don't generate files

    Returns:
        WriteParameters
    """
    main_param_file = os.path.join(this_file_dir(), "h2o.json")
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
