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
"""This module provides co2 property expressions

Vesovic, V., W.A. Wakeham, G.A. Olchowy, J.V. Sengers, J.T.R. Watson, J.
   Millat, (1990). "The transport properties of carbon dioxide." J. Phys.
   Chem. Ref. Data, 19, 763-808.
Fenghour, A., W.A. Wakeham, V. Vesovic, (1998). "The Viscosity of Carbon
   Dioxide." J. Phys. Chem. Ref. Data, 27, 31-44.
"""

import pyomo.environ as pyo


def _thermal_conductivity(blk, delta, tau, on_blk=None):
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
        1: 2.447164e-5,
        2: 8.705605e-8,
        3: -6.547950e-11,
        4: 6.594919e-14,
    }
    T = blk.temperature_star / tau
    Ts = T / (251.196 * pyo.units.K)
    G = sum(b[i] / Ts**i for i in b)
    cint_over_k = 1.0 + pyo.exp(-183.5 * pyo.units.K / T) * sum(
        c[i] * (T / 100 / pyo.units.K) ** (2 - i) for i in c
    )
    return (
        (
            0.475598 * pyo.sqrt(T / pyo.units.K) * (1 + 2.0 / 5.0 * cint_over_k) / G
            + d[1] * delta
            + d[2] * delta**2
            + d[3] * delta**3
            + d[4] * delta**4
        )
        / 1e3
        * pyo.units.W
        / pyo.units.K
        / pyo.units.m
    )


def _viscosity(blk, delta, tau, on_blk=None):
    a = {
        0: 0.235156,
        1: -0.491266,
        2: 5.211155e-2,
        3: 5.347906e-2,
        4: -1.537102e-2,
    }
    d = {
        1: 0.4071119e-8 * pyo.value(blk.dens_mass_crit),
        2: 0.7198037e-10 * pyo.value(blk.dens_mass_crit) ** 2,
        3: 0.2411697e-22 * pyo.value(blk.dens_mass_crit) ** 6,
        4: 0.2971072e-28 * pyo.value(blk.dens_mass_crit) ** 8,
        5: -0.1627888e-28 * pyo.value(blk.dens_mass_crit) ** 8,
    }
    T = blk.temperature_star / tau
    Ts = T / (251.196 * pyo.units.K)
    return (
        (
            1.00697
            * pyo.sqrt(T / pyo.units.K)
            / pyo.exp(sum(a[i] * pyo.log(Ts) ** i for i in a))
            / 1e6
            + d[1] * delta
            + d[2] * delta**2
            + d[3] * delta**6 / Ts**3
            + d[4] * delta**8
            + d[5] * delta**8 / Ts
        )
        * pyo.units.Pa
        * pyo.units.s
    )
