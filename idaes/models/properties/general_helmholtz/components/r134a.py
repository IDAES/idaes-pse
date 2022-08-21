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
"""This module provides r134a property expressions

Perkins, R.A.; Laesecke, A.; Howley, J.; Ramires, M.L.V.; Gurova, A.N.; Cusco, L., 
    Experimental thermal conductivity values for the IUPAC round-robin sample of 
    1,1,1,2-tetrafluoroethane (R134a), NIST Interagency/Internal Report (NISTIR) 
    - 6605, 2000, https://doi.org/10.6028/NIST.IR.6605.

Huber, M.L.; Laesecke, A.; Perkins, R.A., Model for the Viscosity and Thermal 
    Conductivity of Refrigerants, Including a New Correlation for the Viscosity 
    of R134a, Ind. Eng. Chem. Res., 2003, 42, 13, 3163-3178, 
    https://doi.org/10.1021/ie0300880.

Thermal conductivity parameter errata correction from CoolProp parameter file.
lambda^d.g. a1: 8.00982 -> 8.00982e-5
"""

import pyomo.environ as pyo


def _thermal_conductivity(blk, delta, tau, on_blk=None):
    """Thermal condutivity (for now) ommiting critical enhancment"""
    T = blk.temperature_star / tau / pyo.units.K
    rho = delta * blk.dens_mass_star / pyo.units.kg * pyo.units.m**3
    tau = T / blk.temperature_crit * pyo.units.K
    a = [
        -1.05248e-2,
        8.00982e-05,
    ]
    b = {
        1: 0.0037740609300000003,
        2: 0.010534223865,
        3: -0.002952794565,
        4: 0.00128672592,
    }
    delta = rho / (5.049886 * 102.032)
    lambda_dg = a[0] + a[1] * T
    lambda_r = sum(bi * delta**i for i, bi in b.items())

    return (lambda_dg + lambda_r) * pyo.units.W / pyo.units.m / pyo.units.K


def _viscosity(blk, delta, tau, on_blk=None):
    T = blk.temperature_star / tau / pyo.units.K
    rho = delta * blk.dens_mol_star / pyo.units.mol * pyo.units.m**3
    M = 102.031
    sigma = 0.46893
    eok = 299.363
    NA = 6.0221408e23
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
    c = [
        0,
        -2.06900719e-05,
        3.56029549e-07,
        2.11101816e-06,
        1.39601415e-05,
        -4.5643502e-06,
        -3.51593275e-06,
        0.00021476332,
        -0.890173375e-1,
        0.100035295,
        3.163695636,
    ]
    Ts = T / eok
    vs = pyo.exp(sum(ai * pyo.log(Ts) ** i for i, ai in enumerate(a)))
    Bs = sum(bi * Ts**ti for bi, ti in zip(b, te))
    B = NA * sigma**3 * Bs / 1e9**3
    etas = 0.021357 * pyo.sqrt(M * T) / (sigma**2 * vs)
    tau = T / blk.temperature_crit * pyo.units.K
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
    return (etas * (1 + B * rho) / 1e6 + eta) * pyo.units.Pa * pyo.units.s
