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

import numpy as np

# _all__ = ["lin", "linjac", "arr", "arrjac", "refarr", "refarrjac"]
# This subroutine contains definitions for the arrhenious form accepted by RIPE
# Currently, the forms cannot be defined by the user
# The lin,arr, and refarra subroutines return squared residuals
# the linjac, arrjac, and refarrjac subroutiens return the jocbian of the
# esitmated parameter


def lin(y, x, a):
    n, s, ln = np.shape(x)
    guess = 0
    out = np.zeros(n * s)
    q = 0
    for i in range(n):
        for j in range(s):
            guess = 0.0
            for k in range(ln):
                guess += a[k] * x[i, j, k]
            out[q] = (y[i, j] - guess) ** 2
            q += 1
    return out


def linjac(y, x, a):
    n, s, ln = np.shape(x)
    out = np.zeros([n * s, ln])
    q = 0
    for i in range(n):
        for j in range(s):
            for k in range(ln):
                out[q, k] = x[i, j, k]
            q += 1
    return out


def arr(y, T, x, p, gc):
    n, s, ln = np.shape(x)
    out = np.zeros(n * s)
    q = 0
    for i in range(n):
        for j in range(s):
            guess = 0.0
            for k in range(ln):
                guess += p[k] * np.exp((-p[k + ln] / (gc * T[i]))) * x[i, j, k]
            out[q] = (y[i, j] - guess) ** 2  # / sig[i,j]
            q += 1
    return out


def arrjac(y, T, x, p, gc):
    n, s, ln = np.shape(x)
    out = np.zeros([n * s, 2 * ln])
    q = 0
    for i in range(n):
        for j in range(s):
            for k in range(ln):
                out[q, k] = np.exp((-p[k + ln] / (gc * T[i]))) * x[i, j, k]
                out[q, k + ln] = (
                    p[k]
                    * (-1 / (gc * T[i]))
                    * np.exp((-p[k + ln] / (gc * T[i])))
                    * x[i, j, k]
                )
            q += 1
    return out


def refarr(y, T, Tref, x, p, gc):
    n, s, ln = np.shape(x)
    out = np.zeros(n * s)
    q = 0
    for i in range(n):
        for j in range(s):
            guess = 0.0
            for k in range(ln):
                guess += (
                    p[k]
                    * np.exp(-1.0 * p[k + ln] / (gc) * (1 / (T[i]) - 1 / (Tref)))
                    * x[i, j, k]
                )
            out[q] = (y[i, j] - guess) ** 2  # / sig[i,j]
            q += 1
    return out


def refarrjac(y, T, Tref, x, p, gc):
    n, s, ln = np.shape(x)
    out = np.zeros([n * s, 2 * ln])
    q = 0
    for i in range(n):
        for j in range(s):
            for k in range(ln):
                out[q, k] = (
                    np.exp(((-1.0 * p[k + ln]) / (gc)) * (1 / (T[i]) - 1 / (Tref)))
                    * x[i, j, k]
                )
                out[q, k + ln] = (
                    p[k]
                    * ((-1.0 / (gc)) * (1 / (T[i]) - 1 / (Tref)))
                    * np.exp(((-1 * p[k + ln]) / (gc)) * ((1 / (T[i])) - (1 / (Tref))))
                    * x[i, j, k]
                )
            q += 1
    return out
