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
"""Generate parameter and expression files for r1234ze(e)
"""

__author__ = "John Eslick"

import os
from pyomo.common.fileutils import this_file_dir
from idaes.models.properties.general_helmholtz.helmholtz_parameters import (
    WriteParameters,
)


def thermal_conductivity_rule(m):
    """Thermal Conductivity

    Richard A. Perkins and Marcia L. Huber. Measurement and Correlation of the Thermal
        Conductivity of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf) and
        trans-1,3,3,3-Tetrafluoropropene (R1234ze(E)). J. Chem. Eng. Data, 56:4868–4874,
        2011. doi:10.1021/je200811n.
    """
    a = [
        -0.0103589,
        0.0308929,
        0.000230348,
    ]
    b = {
        (1, 1): -0.0428296,
        (1, 2): 0.0434288,
        (2, 1): 0.0927099,
        (2, 2): -0.0605844,
        (3, 1): -0.0702107,
        (3, 2): 0.0440187,
        (4, 1): 0.0249708,
        (4, 2): -0.0155082,
        (5, 1): -0.00301838,
        (5, 2): 0.00210190,
    }
    Tred = 1.0 / m.tau
    l0 = sum(ai * Tred**i for i, ai in enumerate(a))
    lr = sum(m.delta**i * (b[i, 1] + b[i, 2] * Tred) for i in range(1, 6))
    lc = 0
    return 1000 * (l0 + lr + lc)


def viscosity_rule(m):
    """Viscosity rule

    Huber ML, Assael MJ. Correlations for the Viscosity of 2,3,3,3-Tetrafluoroprop-1-ene
        (R1234yf) and trans-1,3,3,3-Tetrafluoropropene (R1234ze(E)). Int J Refrig.
        2016; 71:39-45. doi:10.1016/j.ijrefrig.2016.08.007
    """

    a = [
        -963382,
        9614.09,
        -13.233,
        0.0360562,
        122059,
        -224.741,
    ]
    b = [
        -19.572881,
        219.73999,
        -1015.3226,
        2471.0125,
        -3375.1717,
        2491.6597,
        -787.26086,
        14.085455,
        -0.34664158,
    ]
    c = [
        8.61691913,
        0,
        20.83024738,
        0,
        0.54243690,
        -10.49684841,
        -1.38137689,
        1,
        0,
    ]
    T = m.T_star / m.tau
    rho = m.delta * m.rho_star / m.MW * 1000
    eok = 340
    sigma = 5.017e-1
    NA = 6.0221408e23
    Ts = T / eok
    eta0 = (a[0] + a[1] * T + a[2] * T**2 + a[3] * T**3) / (
        a[4] + a[5] * T + T**2
    )
    Bs = (
        sum(b[i] * Ts ** (-0.25 * i) for i in range(0, 7))
        + b[7] * Ts**-2.5
        + b[8] * Ts**-5.5
    )
    B = NA * Bs * sigma**3 / 1e9**3
    eta1 = eta0 * B
    Tr = 1 / m.tau
    etar = (
        m.delta ** (2.0 / 3.0)
        * Tr**0.5
        * (
            c[0]
            + c[1] * m.delta
            + c[2] * m.delta**2
            + (
                c[3] * m.delta
                + c[4] * m.delta**6
                + c[5] * m.delta * Tr**2
                + c[6] * m.delta**5 * Tr
            )
            / (c[7] * Tr + c[8] * m.delta * Tr)
        )
    )
    return eta0 + eta1 * rho + etar


def main(dry_run=False):
    """Generate parameter and expression files.

    Args:
        dry_run (bool): If dry run don't generate files

    Returns:
        None
    """
    main_param_file = os.path.join(this_file_dir(), "r1234ze.json")
    we = WriteParameters(parameters=main_param_file)

    we.add(
        {
            "viscosity": viscosity_rule,
            "thermal_conductivity": thermal_conductivity_rule,
        }
    )
    we.write(dry_run=dry_run)

    print("ASHRAE Offset")
    print(we.calculate_reference_offset(2.7584034882, 1.64063049539, 0, 0))
    return we


if __name__ == "__main__":
    main()
