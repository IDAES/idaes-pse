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
Modeling for saturated densities and vapor pressure
"""

from idaes.apps import alamopy_depr as alamopy
import importlib

from . import DataImport
from . import DataManipulation


molecule = ""
max_time = 1000


global DLfun, DVfun, PVfun
DLfun, DVfun, PVfun = None, None, None

DataToWrite = []


def DL():
    "ALAMO regression of Saturated Liquid Density"
    global DLfun
    DataImport.DL(molecule)
    Values = DataImport.DLValues

    xval, zval = [], []
    for x in Values:
        manVal = DataManipulation.DL(x)
        xval.append(manVal[0])
        zval.append(manVal[1])

    res = alamopy.alamo(
        xval,
        zval,
        zlabels=["%sDL" % molecule],
        almname="%sDL" % molecule,
        monomialpower=(-1, 1, 2, 3),
        expfcns=1,
        savetrace=True,
        savescratch=True,
    )

    if res is None:
        raise Exception("Model does not exist for Saturated Liquid Density")

    DLfun = importlib.import_module("%sDL" % molecule, "f")


def getDL():
    """
    Imports the regressed saturated liquid density function
    """
    global DLfun
    if DLfun is None:
        DLfun = importlib.import_module("%sDL" % molecule, "f")
    return DLfun


def DV():
    """
    ALAMO regression of saturated vapor density
    """
    global DVfun
    DataImport.DV(molecule)
    Values = DataImport.DVValues

    xval, zval = [], []
    for x in Values:
        manVal = DataManipulation.DV(x)
        xval.append(manVal[0])
        zval.append(manVal[1])

    res = alamopy.alamo(
        xval,
        zval,
        zlabels=["%sDV" % molecule],
        almname="%sDV" % molecule,
        monomialpower=(1, 2, 3, 4, 5, 6),
        savetrace=True,
    )
    # print(res['model'], res['ssr'], res['R2'])

    if res is None:
        raise Exception("Model does not exist for Saturated Vapor Density")

    DVfun = importlib.import_module("%sDV" % molecule, "f")


def getDV():
    """
    Imports the regressed saturated vapor density function
    """
    global DVfun
    if DVfun is None:
        DVfun = importlib.import_module("%sDV" % molecule, "f")
    return DVfun


def PV():
    "ALAMO regression of vapor pressure"
    global PVfun
    DataImport.PV(molecule)
    Values = DataImport.Values

    xval, zval = [], []
    for x in Values:
        manVal = DataManipulation.PV(x)
        xval.append([manVal[0], manVal[1]])
        zval.append(manVal[2])

    res = alamopy.alamo(
        xval,
        zval,
        zlabels=["%sPV" % molecule],
        almname="%sPV" % molecule,
        monomialpower=(1, 2, 3, 4, 5, 6),
    )
    # print(res['model'], res['ssr'], res['R2'])

    if res is None:
        raise Exception("Model does not exist for Vapor Pressure")

    PVfun = importlib.import_module("%sPV" % molecule, "f")


def getPV():
    """
    Imports the regressed vapor pressure function
    """
    global PVfun
    if PVfun is None:
        PVfun = importlib.import_module("%sPV" % molecule, "f")
    return PVfun
