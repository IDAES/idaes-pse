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
"""Generate parameter and expression files for CO2
"""

__author__ = "John Eslick"

import os
import pyomo.environ as pyo
from pyomo.common.fileutils import this_file_dir
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def thermal_conductivity_rule(m):
    """Thermal conductivity rule

    Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). "The Viscosity of Carbon
        Dioxide." J. Phys. Chem. Ref. Data, 27, 31-44.
    """
    b = {
        0: 0.4226159,
        1: 0.6280115,
        2: -0.5387661,
        3: 0.6735941,
        4: 0,
        5: 0,
        6: -0.4362677,
        7: 0.2255388,
    }
    c = {
        1: 2.387869e-2,
        2: 4.350794,
        3: -10.33404,
        4: 7.981590,
        5: -1.940558,
    }
    d = {
        1: 24.47164,
        2: 8.705605e-2,
        3: -6.547950e-5,
        4: 6.594919e-8,
    }
    T = m.T_star / m.tau
    Ts = T / 251.196
    rho = m.rho_star * m.delta
    G = sum(bval / Ts**i for i, bval in b.items())
    cint_over_k = 1.0 + pyo.exp(-183.5 / T) * sum(
        cval * (T / 100) ** (2 - i) for i, cval in c.items()
    )
    return (
        475.598 * pyo.sqrt(T) * (1 + 2.0 / 5.0 * cint_over_k) / G
        + d[1] * rho
        + d[2] * rho**2
        + d[3] * rho**3
        + d[4] * rho**4
    ) / 1e3


def viscosity_rule(m):
    """Viscosity rule

    Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). "The Viscosity of Carbon
        Dioxide." J. Phys. Chem. Ref. Data, 27, 31-44.
    """
    a = {
        0: 0.235156,
        1: -0.491266,
        2: 5.211155e-2,
        3: 5.347906e-2,
        4: -1.537102e-2,
    }
    d = {
        1: 0.4071119e-2,
        2: 0.7198037e-4,
        3: 0.2411697e-16,
        4: 0.2971072e-22,
        5: -0.1627888e-22,
    }
    T = m.T_star / m.tau
    rho = m.delta * m.rho_star
    Ts = T / 251.196
    return (
        1.00697
        * pyo.sqrt(T)
        / pyo.exp(sum(aval * pyo.log(Ts) ** i for i, aval in a.items()))
        + d[1] * rho
        + d[2] * rho**2
        + d[3] * rho**6 / Ts**3
        + d[4] * rho**8
        + d[5] * rho**8 / Ts
    )


def main(dry_run=False):
    """Generate parameter and expression files.

    Args:
        dry_run (bool): If dry run don't generate files

    Returns:
        None
    """
    main_param_file = os.path.join(this_file_dir(), "co2.json")
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
