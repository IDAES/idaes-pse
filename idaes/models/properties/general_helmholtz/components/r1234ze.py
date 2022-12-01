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
"""This module provides r1234ze property expressions

Richard A. Perkins and Marcia L. Huber. Measurement and Correlation of the Thermal
    Conductivity of 2,3,3,3-Tetrafluoroprop-1-ene (R1234yf) and 
    trans-1,3,3,3-Tetrafluoropropene (R1234ze(E)). J. Chem. Eng. Data, 56:4868–4874, 
    2011. doi:10.1021/je200811n.

Huber ML, Assael MJ. Correlations for the Viscosity of 2,3,3,3-Tetrafluoroprop-1-ene 
    (R1234yf) and trans-1,3,3,3-Tetrafluoropropene (R1234ze(E)). Int J Refrig. 
    2016; 71:39-45. doi:10.1016/j.ijrefrig.2016.08.007

"""

import pyomo.environ as pyo


def _thermal_conductivity(blk, delta, tau, on_blk=None):
    """Thermal condutivity (for now) ommiting critical enhancment"""
    T = blk.temperature_star / tau / pyo.units.K
    rho = delta * blk.dens_mass_star / pyo.units.kg * pyo.units.m**3
    Tred = 1.0 / tau

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
    l0 = sum(ai * Tred**i for i, ai in enumerate(a))
    lr = sum(delta**i * (b[i, 1] + b[i, 2] * Tred) for i in range(1, 6))
    lc = 0
    return (l0 + lr + lc) * pyo.units.W / pyo.units.m / pyo.units.K


def _viscosity(blk, delta, tau, on_blk=None):
    """Thermal condutivity (for now) ommiting critical enhancment"""

    T = blk.temperature_star / tau / pyo.units.K
    rho = delta * blk.dens_mol_star / pyo.units.mol * pyo.units.m**3

    eok = 340
    sigma = 5.017e-1
    M = 114.0415928
    NA = 6.0221408e23
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
    Tr = 1 / tau
    etar = (
        delta ** (2.0 / 3.0)
        * Tr**0.5
        * (
            c[0]
            + c[1] * delta
            + c[2] * delta**2
            + (
                c[3] * delta
                + c[4] * delta**6
                + c[5] * delta * Tr**2
                + c[6] * delta**5 * Tr
            )
            / (c[7] * Tr + c[8] * delta * Tr)
        )
    )

    return (eta0 + eta1 * rho + etar) / 1e6 * pyo.units.Pa * pyo.units.s
