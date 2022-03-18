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
"""
Basis functions for generating the multiparameter equation of state
"""
import numpy as np
import sympy as sy


drd_vals = []
ar_vals = []
ai_vals = []
itt_val = 0
rtt_vals = []
dtrdt_vals = []
d2rd_vals = []
# PYLINT-TODO check if fix OK: adding d2rdt_vals variable here to address the undefined-variable error on L286
d2rdt_vals = []
d3rd_vals = []
d4rd_vals = []
d5rd_Vals = []
critT, critD, critP, acc, R, M, Rm = [0, 0, 0, 0, 0, 0, 0]
molecule = ""
coeffs = []
indexes = []


def molData(fluidData, Dmolecule, RVal):
    """Passing of the Data from the main module ::module:: MPEOSDeveloperModule

    :param fluidData: (critT, critP, critD, M, triple, acc).
    :type fluidData: array
    :param Dmolecule: Name of the molecule of interest.
    :type Dmolecule: str.
    :param RVal: Gas Constant.
    :type RVal: int.
    """

    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    molecule = Dmolecule
    (critT, critP, critD, M, triple, acc) = fluidData
    R = RVal
    Rm = RVal


def formCustomBasis(LemJac=False):
    """Basis Functions developed a bank of terms based on literature (Lemmon, Span, Wagner)"""
    global coeffs, indexes
    coeffs = []
    if LemJac:
        for i in range(1, 9):
            for j in range(1, 13):  # 12
                coeffs.append([i, j / 8.0, 0, 0])
        indexes.append(len(coeffs))
        for i in range(1, 6):
            for j in range(1, 24):  # 24
                coeffs.append([i, j / 8.0, 1, 0])
        for i in range(1, 6):
            for j in range(1, 30):  # 30
                coeffs.append([i, j / 8.0, 2, 0])
        for i in range(2, 5):
            for j in range(24, 38):  # 38
                coeffs.append([i, j / 2.0, 3, 0])
        indexes.append(len(coeffs))
        for i in range(1, 6):
            for j in range(1, 24):  # 24
                for m in range(1, 7):
                    coeffs.append([i, j / 8.0, 1, m / 2])
        for i in range(1, 6):
            for j in range(1, 30):  # 24
                for m in range(1, 7):
                    coeffs.append([i, j / 8.0, 2, m / 2])
        for i in range(2, 5):
            for j in range(24, 38):  # 38
                for m in range(1, 7):
                    coeffs.append([i, j / 2.0, 3, m / 2])
        indexes.append(len(coeffs))
    else:
        for i in range(1, 9):
            for j in range(1, 13):  # 12
                coeffs.append([i, j / 8.0, 0, 0])
        indexes.append(len(coeffs))
        for i in range(1, 6):
            for j in range(1, 24):  # 24
                coeffs.append([i, j / 8.0, 1, 0])
        for i in range(1, 6):
            for j in range(1, 30):  # 24
                coeffs.append([i, j / 8.0, 2, 0])
        for i in range(2, 5):
            for j in range(24, 38):  # 38
                coeffs.append([i, j / 2.0, 3, 0])
        indexes.append(len(coeffs))


def getTerm(Y):
    """Prints index and basis function based on Y index"""
    for y in Y:
        print(y)
        d, t, c, m = coeffs[y - 1]
        print("Delta^%f Tau^%f np.exp(-Delta^%f) np.exp(-Tau^%f)" % (d, t, c, m))


def idealBY(D, T, Y, Beta):
    """Ideal Helmholtz contribution"""

    # TOLUENE
    if molecule == "TOL":
        a = [3.5241174832, 1.1360823464]
        c = [4.0, 0, 0]
        # B = [
        #     0.96464,
        #     -2.7855,
        #     0.86712,
        #     -0.18860,
        #     0.11804,
        #     0.0025181,
        #     0.57196,
        #     -0.029287,
        #     -0.43351,
        #     -0.12540,
        #     -0.028207,
        #     0.014076,
        # ]
        v = [1.6994, 8.0577, 17.059, 8.4567, 8.6423]
        u = [190, 797, 1619, 3072, 7915]

    # Carbon Monoxide
    if molecule == "CO":
        a = [-3.3728318564, 3.3683460039]
        c = [3.5, 0.22311e-6, 1.5]
        # B = [
        #     0.90554,
        #     -2.4515,
        #     0.53149,
        #     0.024173,
        #     0.072156,
        #     0.00018818,
        #     0.19405,
        #     -0.043268,
        #     -0.12778,
        #     -0.027896,
        #     -0.03414,
        #     0.016329,
        # ]
        v = [1.0128]
        u = [3089]

    # Water
    if molecule == "H2O":
        a = [-8.32044648201, 6.6832105268]
        c = [3.00632 + 1, 0, 0]
        v = [0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
        u = [
            1.2878967 * float(critT),
            3.53734222 * float(critT),
            7.74073708 * float(critT),
            9.24437796 * float(critT),
            27.5075105 * float(critT),
        ]

    if c[2] > 0:
        idealA = (
            a[0]
            + a[1] * T
            + np.log(D)
            + (c[0] - 1) * np.log(T)
            - c[1] * (critT) ** c[2] / (c[2] * (c[2] + 1)) * T ** (-c[2])
        )
    else:
        idealA = a[0] + a[1] * T + np.log(D) + (c[0] - 1) * np.log(T)
    for uit, vi in zip(u, v):
        ui = uit / float(critT)
        idealA = idealA + vi * np.log(1 - np.exp(-ui * T))
    ai_vals = idealA
    return ai_vals


def arBY(D, T, Y, Beta):
    """Residual Helmholtz contribution"""
    global drd_vals
    global coeffs
    # drd_vals = [];
    val = 0
    # print len(Y)
    if type(Y) is not int:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                try:
                    val += y * (D**d) * (T**t)
                    pass
                except Exception:
                    print("Term", x, d, t, T, D)
            elif m == 0:
                val = np.exp(-(D**c)) * (D ** (d)) * (T**t)
            else:
                val = np.exp(-(D**c)) * (D ** (d)) * (T**t) * np.exp(-(T**m))

        return val


def drd(D, T):
    """Partial derivative with respect to density"""
    global drd_vals
    global coeffs

    drd_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (D**d) * (T**t)
        elif m == 0:
            val = np.exp(-(D**c)) * ((D ** (d)) * (T**t) * (d - c * (D**c)))
        else:
            val = (
                np.exp(-(D**c))
                * (D ** (d))
                * (T**t)
                * np.exp(-(T**m))
                * (d - c * (D**c))
            )
        drd_vals.append(val)


def d2rd(D, T):
    """Partial derivative with respect to density twice"""
    global d2rd_vals
    global coeffs
    d2rd_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (D**d) * (T**t)
        elif m == 0:
            val = np.exp(-(D**c)) * (
                (D ** (d))
                * (T**t)
                * ((d - c * (D**c)) * (d - 1 - c * (D**c)) - (c**2) * (D**c))
            )
        else:
            val = np.exp(-(D**c)) * (
                (D ** (d))
                * (T**t)
                * np.exp(-(T**m))
                * ((d - c * (D**c)) * (d - 1 - c * (D**c)) - (c**2) * (D**c))
            )
        d2rd_vals.append(val)


def d2rdt(D, T):
    """Partial derivative with respect to density twice and temperature"""
    global d2rdt_vals
    global coeffs

    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * t * (D**d) * (T**t)
        elif m == 0:
            val = np.exp(-(D**c)) * (
                (D ** (d))
                * (T**t)
                * ((d - c * (D**c)) * (d - 1 - c * (D**c)) - (c**2) * (D**c))
            )

        else:
            val = (
                np.exp(-(D**c))
                * (
                    (D ** (d))
                    * (T**t)
                    * np.exp(-(T**m))
                    * (
                        (d - c * (D**c)) * (d - 1 - c * (D**c))
                        - (c**2) * (D**c)
                    )
                )
                * (t - m * (T**m))
            )
        d2rdt_vals.append(val)


def d3rd(D, T):
    """Third partial derivative with respect to density"""
    global d3rd_vals
    global coeffs
    d3rd_vals = []
    # i = 1
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (d - 2) * (D**d) * (T**t)
        elif m == 0:
            # val = (D**(d))*(T**t)*exp(-(D**c))*(d*(d-1)*(d-2) + (D**c)*(-2*c + 6*d*c-3*(d**2)*c - 3*d*(c**2) + 3*(c**2) - (c**3)) + (D**(2*c))*(3*d*(c**2) -3*(c**2) + 3*(c**3)) - (c**3)*(D**(3*c)))
            val = (
                (D ** (d))
                * (T ** (t))
                * np.exp(-(D**c))
                * (
                    (d - c * (D**c)) * (d - 1 - c * (D**c)) * (d - 2 - c * (D**c))
                    + (3 * d - 3 + c - 3 * c * (D**c)) * (-(c**2) * (D ** (c)))
                )
            )

        else:
            val = (
                (D ** (d))
                * (T**t)
                * np.exp(-(T**m))
                * np.exp(-(D**c))
                * (
                    d * (d - 1) * (d - 2)
                    + (D**c)
                    * (
                        -2 * c
                        + 6 * d * c
                        - 3 * (d**2) * c
                        - 3 * d * (c**2)
                        + 3 * (c**2)
                        - (c**3)
                    )
                    + (D ** (2 * c)) * (3 * d * (c**2) - 3 * (c**2) + 3 * (c**3))
                    - (c**3) * (D ** (3 * c))
                )
            )
        d3rd_vals.append(val)


def d4rd(D, T):
    """Fourth partial derivative with respect to density"""
    global d4rd_vals
    global coeffs
    d4rd_vals = []
    # i = 1
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (d - 2) * (d - 3) * (D**d) * (T**t)
        elif m == 0:
            # val = (D**(d))*(T**t)*exp(-(D**c))*(d*(d-1)*(d-2) + (D**c)*(-2*c + 6*d*c-3*(d**2)*c - 3*d*(c**2) + 3*(c**2) - (c**3)) + (D**(2*c))*(3*d*(c**2) -3*(c**2) + 3*(c**3)) - (c**3)*(D**(3*c)))
            val = (
                (D ** (d))
                * (T ** (t))
                * np.exp(-(D**c))
                * (
                    (d - c * (D**c))
                    * (d - 1 - c * (D**c))
                    * (d - 2 - c * (D**c))
                    * (d - 3 - c * (D**c))
                    + (
                        6 * (c**2) * D ** (2 * c)
                        - 7 * (c**2) * (D**c)
                        + (c**2)
                        - 12 * c * d * (D**c)
                        + 4 * c * d
                        + 18 * c * (D**c)
                        - 6 * c
                        + 6 * d**2
                        - 18 * d
                        + 11
                    )
                    * (-(c**2) * (D ** (c)))
                )
            )

        else:
            val = (
                (D ** (d))
                * (T**t)
                * np.exp(-(T**m))
                * np.exp(-(D**c))
                * (
                    d * (d - 1) * (d - 2)
                    + (D**c)
                    * (
                        -2 * c
                        + 6 * d * c
                        - 3 * (d**2) * c
                        - 3 * d * (c**2)
                        + 3 * (c**2)
                        - (c**3)
                    )
                    + (D ** (2 * c)) * (3 * d * (c**2) - 3 * (c**2) + 3 * (c**3))
                    - (c**3) * (D ** (3 * c))
                )
            )
            val = 0
        d4rd_vals.append(val)


def d5rd(D, T):
    """Fifth partial derivative with respect to density"""
    global d5rd_vals
    global coeffs
    d5rd_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * (d - 1) * (d - 2) * (d - 3) * (d - 4) * (D**d) * (T**t)
        elif m == 0:
            x, y = sy.symbols("x, y")
            expr = x**d * y**t * sy.exp(-(x**c))
            f_prime = expr.diff(x, 5)
            val = f_prime.subs([(x, D), (y, T)]) * (D**5)
        else:
            val = (
                (D ** (d))
                * (T**t)
                * np.exp(-(T**m))
                * np.exp(-(D**c))
                * (
                    d * (d - 1) * (d - 2)
                    + (D**c)
                    * (
                        -2 * c
                        + 6 * d * c
                        - 3 * (d**2) * c
                        - 3 * d * (c**2)
                        + 3 * (c**2)
                        - (c**3)
                    )
                    + (D ** (2 * c)) * (3 * d * (c**2) - 3 * (c**2) + 3 * (c**3))
                    - (c**3) * (D ** (3 * c))
                )
            )
            val = 0
        d5rd_vals.append(val)


def dtrdt(D, T):
    """Second partial derivative with respect to density and temperature"""
    global dtrdt_vals
    global coeffs
    dtrdt_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = d * t * (D**d) * (T**t)
        elif m == 0:
            val = np.exp(-(D**c)) * ((D ** (d)) * (T**t) * (d - c * (D**c)))
        else:
            val = (
                np.exp(-(D**c))
                * ((D ** (d)) * (T**t) * (d - c * (D**c)))
                * np.exp(-(T**m))
                * (t - m * (T**m))
            )
        dtrdt_vals.append(val)


def rTT(D, T):
    """Second partial derivative with respect to temperature"""
    global rtt_vals
    global coeffs
    rtt_vals = []
    for d, t, c, m in coeffs:
        val = 0
        if c == 0:
            val = t * (t - 1) * (D**d) * (T ** (t))
        elif m == 0:
            val = np.exp(-(D**c)) * ((D ** (d)) * (T ** (t))) * t * (t - 1)
        else:  # NOT DONE
            val = np.exp(-(D**c)) * (D ** (d)) * (T**t) * np.exp(-(T**t))
        rtt_vals.append(val)


def d3rdRes(D, T, Y, Beta):
    """Third partial derivative with respect to density"""
    global d3rd_vals
    global coeffs
    d3rd_vals = []
    val = 0
    for x, y in zip(Y, Beta):
        x = x - 1
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val += y * d * (d - 1) * (d - 2) * (D**d) * (T**t)
        elif m == 0:  # TODO
            val += (
                y
                * np.exp(-(D**c))
                * (D ** (d))
                * (T**t)
                * (
                    d * (d - 1) * (d - 2)
                    + c
                    * (D**c)
                    * (-2 + 6 * d - 3 * (d**2) - 3 * d * c + 3 * c - (c**2))
                    + (D ** (2 * c)) * (3 * d * (c**2) - 3 * (c**2) + 3 * (c**3))
                    - (c**3) * (D ** (3 * c))
                )
            )
        else:
            val += (
                y
                * np.exp(-(D**c))
                * (D ** (d))
                * (T**t)
                * np.exp(-(T**m))
                * (
                    d * (d - 1) * (d - 2)
                    + c
                    * (D**c)
                    * (-2 + 6 * d - 3 * (d**2) - 3 * d * c + 3 * c - (c**2))
                    + (D ** (2 * c)) * (3 * d * (c**2) - 3 * (c**2) + 3 * (c**3))
                    - (c**3) * (D ** (3 * c))
                )
            )
        # print val
    return val


# PVT derivatives
def drdRes(D, T, Y, Beta):
    """
    Calculates the partial derivaties w.r.t. density

    Args:
        D: Delta
        T: Tau
        Y: index of basis Function (int or array)
        Beta: weighting (float or array)
    """
    global drd_vals
    global coeffs
    val = 0
    if type(Y) is not int:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                try:
                    val += y * d * (D**d) * (T**t)
                    pass
                except Exception:
                    print("Term", x, d, t, T, D)
            elif m == 0:
                val += (
                    y * np.exp(-(D**c)) * ((D ** (d)) * (T**t) * (d - c * (D**c)))
                )
            else:
                val += (
                    y
                    * np.exp(-(D**c))
                    * (D ** (d))
                    * (T**t)
                    * np.exp(-(T**m))
                    * (d - c * (D**c))
                )
        return val
    else:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val += y * d * (D**d) * (T**t)
        elif m == 0:
            val += y * np.exp(-(D**c)) * ((D ** (d)) * (T**t) * (d - c * (D**c)))
        else:
            val += (
                y
                * np.exp(-(D**c))
                * (D ** (d))
                * (T**t)
                * np.exp(-(T**m))
                * (d - c * (D**c))
            )
        return val


# CV derivatives
def iTT(D, T):
    """Ideal helmholtz contribution second partial derivative with respect to temperature"""
    global itt_val

    # TOLUENE
    if molecule == "TOL":
        # a = [3.5241174832, 1.1360823464]
        c = [4.0, 0, 0]
        # B = [.96464, -2.7855, 0.86712, -0.18860, 0.11804, 0.0025181, 0.57196, -0.029287, -0.43351, -0.12540, -0.028207, 0.014076];
        v = [1.6994, 8.0577, 17.059, 8.4567, 8.6423]
        u = [190, 797, 1619, 3072, 7915]

    # Carbon Monoxide
    if molecule == "CO":
        # a = [-3.3728318564, 3.3683460039]
        c = [3.5, 0.22311e-6, 1.5]
        # B = [0.90554, -2.4515, 0.53149, 0.024173, 0.072156, 0.00018818, 0.19405, -0.043268, -0.12778, -0.027896, -0.03414, 0.016329]
        v = [1.0128]
        u = [3089]

    # Carbon Dioxide
    if molecule == "CO2":
        # a = [8.37304456, -3.70454304, 2.5]
        c = [2.5 + 1, 0, 0]
        v = [1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678]
        u = [
            3.15163 * float(critT),
            6.11190 * float(critT),
            6.77708 * float(critT),
            11.32384 * float(critT),
            27.0 * float(critT),
        ]

    # Water
    if molecule == "H2O":
        # a = [-8.32044648201, 6.6832105268]
        c = [3.00632 + 1, 0, 0]
        v = [0.012436, 0.97315, 1.27950, 0.96956, 0.24873]
        u = [
            1.2878967 * float(critT),
            3.53734222 * float(critT),
            7.74073708 * float(critT),
            9.24437796 * float(critT),
            27.5075105 * float(critT),
        ]

    itt_val = 0
    if c[2] != 0:
        itt_val = -(c[0] - 1) - c[1] * (critT ** c[2]) / (c[2] * (c[2] + 1)) * (
            T ** (-c[2])
        ) * (-c[2]) * (-c[2] - 1)
        # /T**2;
    else:
        itt_val = -(c[0] - 1)
    for uit, vi in zip(u, v):
        ui = uit / float(critT)
        itt_val = itt_val - vi * ((ui) ** 2) * (T**2) * np.exp(-ui * T) * (
            (1 - np.exp(-ui * T)) ** (-2)
        )
    return itt_val


def rTTRes(D, T, Y, Beta):
    """Residual helmholtz contribution second partial derivative with respect to temperature"""
    global rtt_vals
    global coeffs
    val = 0
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            value = t * (t - 1) * (D**d) * (T ** (t))
        elif m == 0:
            value = np.exp(-(D**c)) * ((D ** (d)) * (T ** (t))) * t * (t - 1)
        else:
            value = (
                np.exp(-(D**c))
                * ((D ** (d)) * (T ** (t)))
                * ((t - m * (T**m)) * (t - 1 - m * (T**m)) - (m**2) * (T**m))
            )
        val = y * value
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                value = t * (t - 1) * (D**d) * (T ** (t))
            elif m == 0:
                value = np.exp(-(D**c)) * ((D ** (d)) * (T ** (t))) * t * (t - 1)
            else:
                value = (
                    np.exp(-(D**c))
                    * ((D ** (d)) * (T ** (t)))
                    * (
                        (t - m * (T**m)) * (t - 1 - m * (T**m))
                        - (m**2) * (T**m)
                    )
                )
            val += y * value
    return val


def d2rdRes(D, T, Y, Beta):
    """Residual helmholtz contribution second partial derivative with respect to density"""
    global d2rd_vals
    global coeffs
    val = 0
    d2rd_vals = []
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val = y * d * (d - 1) * (D**d) * (T**t)
        elif m != 0:  # TODO
            val = np.exp(-(D**c)) * (D ** (d)) * (T**t) * np.exp(-(T**t))
        else:
            val = (
                y
                * np.exp(-(D**c))
                * (
                    (D ** (d))
                    * (T**t)
                    * (
                        (d - c * (D**c)) * (d - 1 - c * (D**c))
                        - (c**2) * (D**c)
                    )
                )
            )
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                val += y * d * (d - 1) * (D**d) * (T**t)
            elif m == 0:  # TODO
                val += (
                    y
                    * np.exp(-(D**c))
                    * (
                        (D ** (d))
                        * (T**t)
                        * (
                            (d - c * (D**c)) * (d - 1 - c * (D**c))
                            - (c**2) * (D**c)
                        )
                    )
                )
            else:
                val += (
                    y
                    * np.exp(-(D**c))
                    * (
                        (D ** (d))
                        * (T**t)
                        * np.exp(-(T**m))
                        * (
                            (d - c * (D**c)) * (d - 1 - c * (D**c))
                            - (c**2) * (D**c)
                        )
                    )
                )
    return val


def dtrdtRes(D, T, Y, Beta):
    """Residual helmholtz contribution second partial derivative with respect to density and temperature"""
    global dtrdt_vals
    global coeffs
    val = 0
    dtrdt_vals = []
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val = y * d * t * (D**d) * (T**t)
        elif m == 0:  # TODO
            val = (
                y * t * np.exp(-(D**c)) * ((D ** (d)) * (T**t) * (d - c * (D**c)))
            )
        else:
            val = (
                y
                * t
                * np.exp(-(D**c))
                * ((D ** (d)) * (T**t) * np.exp(-(T**m)) * (d - c * (D**c)))
                * (t - m * (T**m))
            )
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                val += y * d * t * (D**d) * (T**t)
            elif m == 0:
                val += (
                    y
                    * t
                    * np.exp(-(D**c))
                    * ((D ** (d)) * (T**t) * (d - c * (D**c)))
                )
            else:
                val += (
                    y
                    * t
                    * np.exp(-(D**c))
                    * ((D ** (d)) * (T**t) * np.exp(-(T**m)) * (d - c * (D**c)))
                    * (t - m * (T**m))
                )
    return val


def d2rdrtRes(D, T, Y, Beta):
    """Residual helmholtz contribution third partial derivative with respect to density(2) and temperature(1)"""
    global d2rd_vals
    global coeffs
    val = 0
    d2rd_vals = []
    if type(Y) is int:
        x = Y - 1
        y = Beta
        [d, t, c, m] = coeffs[x]
        if c == 0:
            val = y * d * (d - 1) * (D**d) * t * (T**t)
        elif m != 0:
            val = (
                np.exp(-(D**c))
                * (D ** (d))
                * (T**t)
                * np.exp(-(T**t))
                * (t - m * (T**m))
            )
        else:
            val = (
                y
                * np.exp(-(D**c))
                * (
                    (D ** (d))
                    * (T**t)
                    * (
                        (d - c * (D**c)) * (d - 1 - c * (D**c))
                        - (c**2) * (D**c)
                    )
                )
                * t
            )
    else:
        for x, y in zip(Y, Beta):
            x = x - 1
            [d, t, c, m] = coeffs[x]
            if c == 0:
                val += y * d * (d - 1) * (D**d) * t * (T**t)
            elif m == 0:
                val += (
                    y
                    * np.exp(-(D**c))
                    * (
                        (D ** (d))
                        * (T**t)
                        * (
                            (d - c * (D**c)) * (d - 1 - c * (D**c))
                            - (c**2) * (D**c)
                        )
                    )
                    * t
                )
            else:
                val += (
                    y
                    * np.exp(-(D**c))
                    * (
                        (D ** (d))
                        * (T**t)
                        * np.exp(-(T**m))
                        * (
                            (d - c * (D**c)) * (d - 1 - c * (D**c))
                            - (c**2) * (D**c)
                        )
                    )
                    * (t - m * (T**m))
                )
    return val
