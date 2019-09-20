##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
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
    Inputs:
        X = [Pressure, Density, Temperature]
    OutputS:
        X = [Delta, Tau, Compressibility]
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
    Inputs:
        X = [Density, Temperature, Isobaric Heat Capacity]
    Outputs:
        X = [Delta, Tau, CP]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    CP = float(x[2]) / R

    return [Delta, Tau, CP]


def CV(x):
    """
    Calculate dimensionless isochoric heat capacity
    Inputs:
        X = [Density, Temperature, Isochoric Heat Capacity]
    Outputs:
        X = [Delta, Tau, CV]
    """
    Delta = float(x[0]) / float(critD)
    Tau = float(critT) / float(x[1])
    CV = float(x[2]) / R  
    return [Delta, Tau, CV]


def SND(x):
    """
    Calculate dimensionless speed of sound
    Inputs:
        X = [Density, Temperature, Speed of Sound]
    Outputs:
        X = [Delta, Tau, W]
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
    Inputs:
         X = [Density, Temperature]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Delta = float(x[0]) / float(critD) - 1
    return [Theta, Delta]


def DV(x):
    """
    Calculate Theta and Delta for saturated vapor density
    Inputs:
         X = [Density, Temperature]
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
    Inputs:
         X = [Pressure, Temperature]
    """
    Theta = 1 - float(x[1]) / float(critT)
    Psi = math.log(float(x[0]) / float(critP))
    Tau = float(critT) / float(x[1])
    return [Tau, Theta, Psi]
