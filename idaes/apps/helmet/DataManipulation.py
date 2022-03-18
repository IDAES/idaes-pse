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
Calculates dimensionless data
"""

import math

molecule = ""
critT, critD, critP, acc, R, M, Rm = [0, 0, 0, 0, 0, 0, 0]


def molData(fluidData, mol, RVal):
    """
    Sets up important global values
    """
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    (critT, critP, critD, M, triple, acc) = fluidData
    molecule = mol
    R = RVal
    Rm = RVal


def PVT(x):
    """
    Calculate dimensionless compressibility

    Args:
        X: [Pressure, Density, Temperature]

    Returns:
        [Delta, Tau, Compressibility]
    """
    Pressure = float(x[0])
    Density = float(x[1])
    Temperature = float(x[2])

    Tau = float(critT) / Temperature
    Delta = Density / float(critD)

    if Density != 0:
        Z = (Pressure * 1000 / R / Temperature / Density) - 1.00
    else:
        Z = 0
    return [Delta, Tau, Z]


def P(x):
    """
    Calculate reduced density and inverse reduced temperature
    Return array of Delta, Tau, Pressure
    """
    Temperature = float(x[2])
    Density = float(x[1])
    Pressure = float(x[0])
    Tau = float(critT) / Temperature
    Delta = Density / float(critD)
    return [Delta, Tau, Pressure]


def CP(x):
    """
    Calculate dimensionless isobaric heat capacity

    Args:
        X: [Density, Temperature, Isobaric Heat Capacity]

    Returns:
        [Delta, Tau, CP]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    CP = float(x[2]) / R

    return [Delta, Tau, CP]


def CV(x):
    """
    Calculate dimensionless isochoric heat capacity

    Args:
        X: [Density, Temperature, Isochoric Heat Capacity]

    Returns:
        [Delta, Tau, CV]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    CV = float(x[2]) / R
    return [Delta, Tau, CV]


def SND(x):
    """
    Calculate dimensionless speed of sound

    Args:
        X: [Density, Temperature, Speed of Sound]

    Returns:
        [Delta, Tau, W]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    w = (float(x[2]) ** 2) / R / 1000 * M / float(x[1])
    return [Delta, Tau, w]


def CP0(x):
    """
    Calculate dimensionless ideal isobaric heat capacity
    """
    Cp0 = float(x[1]) / Rm / 1000.0
    Temp = float(x[0]) / 1000.0
    return [Temp, Cp0]


def DL(x):
    """
    Calculate Theta and Delta for saturated liquid density

    Args:
         X: [Density, Temperature]

    Returns:
        [Theta, Delta]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Delta = float(x[0]) / float(critD) - 1
    return [Theta, Delta]


def DV(x):
    """
    Calculate Theta and Delta for saturated vapor density

    Args:
         X: [Density, Temperature]

    Returns:
        [Theta, Delta]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Delta = math.log(float(x[0]) / float(critD))
    return [Theta, Delta]


def Dsat(x):
    """
    Calculate dimensionless terms
    """
    Theta = float(critT) / float(x[1])
    Delta = float(x[0]) / float(critD)
    return [Delta, Theta]


def PV(x):
    """
    Calculate Theta, Tau, and Psi for saturated liquid density

    Args:
         X: [Pressure, Temperature]

    Returns:
        [Tau, Theta, Psi]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Psi = math.log(float(x[0]) / float(critP))
    Tau = float(critT) / float(x[1])
    return [Tau, Theta, Psi]
