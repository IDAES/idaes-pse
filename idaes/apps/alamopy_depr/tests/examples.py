#!/usr/bin/python
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
import sys


def sixcamel(*x):
    x1, x2 = x
    t1 = np.multiply(
        4.0 - 2.1 * np.power(x1, 2) + np.divide(np.power(x1, 4), 3.0), np.power(x1, 2)
    )
    t2 = np.multiply(4 * np.power(x2, 2) - 4, np.power(x2, 2))
    z = t1 + np.multiply(x1, x2) + t2
    return z


def ackley(*x):
    import numpy as np

    x1, x2 = x
    a = 20
    b = 0.2
    c = 2 * 3.14159
    z = (
        -a * np.exp(-b * np.sqrt(0.5 * (x1**2 + x2**2)))
        - np.exp(0.5 * (np.cos(c * x1) + np.cos(c * x2)))
        + a
        + np.exp(1)
    )
    return z


def branin(*x):
    import numpy as np

    x1, x2 = x
    pi = 3.14159
    z = (
        (x2 - (5.1 / (4 * pi**2)) * x1**2 + (5 / pi) * x1 - 6) ** 2
        + 10 * (1 - (1 / (8 * pi)) * np.cos(x1) + 10)
        + np.random.normal(0, 0.1)
    )
    return z


if __name__ == "__main__":
    sys.stdout.write(" ALAMOpy example functions ")
    sys.stdout.write(" call functions with : ")
    sys.stdout.write(" examples.<name>")
    sys.stdout.write(" <name> = branin ")
    sys.stdout.write("          sixcamel ")
    sys.stdout.write("          ackley ")
