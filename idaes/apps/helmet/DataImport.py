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
Importing thermodynamic data, specific structures for text files
"""

import numpy as np


filename, sampleRatio = "", ""
(
    Values,
    DVValues,
    DLValues,
    CPValues,
    PVTValues,
    PVTSamples,
    CP0Values,
    BValues,
    CVValues,
    CVValuesn,
    SNDValues,
) = (
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
    [],
)
InSatValues = []
PVTindexes, CVindexes, CPindexes, SNDindexes = [], [], [], []
isothermIndex = 0

critT, critD, critP, acc, R, M, Rm = [0, 0, 0, 0, 0, 0, 0]


def molData(fluidData, RVal):
    """
    Molecular data passed to the module
    """
    global critT, critP, critD, M, triple, acc, R, Rm
    (critT, critP, critD, M, triple, acc) = fluidData
    R = RVal
    Rm = RVal


def regionsOfData(molecule, DataValues, PVT=False, CV=False):
    """
    Organization of data into regions
    """

    Reg1, Reg2, Reg3, Reg4, Reg5, Reg6 = [], [], [], [], [], []

    if PVT:
        for P, D, T in DataValues:
            T = float(T)
            P = float(P)
            D = float(D)

            vals = [P, D, T]
            if T < float(critT):
                if D < float(critD):
                    Reg1.append(vals)
                else:
                    Reg2.append(vals)
            else:
                if (
                    T / float(critT) < 1.1
                    and T / float(critT) > 0.98
                    and D / float(critD) > 0.7
                    and D / float(critD) < 1.4
                ):
                    Reg3.append(vals)
                else:
                    Rho = D / float(critD)
                    if Rho < 0.6:
                        Reg4.append(vals)
                    elif Rho < 1.5:
                        Reg5.append(vals)
                    else:
                        Reg6.append(vals)
    elif CV:
        for D, T, CV in CVValues:
            T = float(T)
            CV = float(CV)
            D = float(D)

            vals = [D, T, CV]
            if D > 0:
                if T < float(critT):
                    if D < float(critD):
                        Reg1.append(vals)
                    else:
                        Reg2.append(vals)
                else:
                    if (
                        T / float(critT) < 1.1
                        and T / float(critT) > 0.98
                        and D / float(critD) > 0.7
                        and D / float(critD) < 1.4
                    ):
                        Reg3.append(vals)
                    else:
                        Rho = D / float(critD)
                        if Rho < 0.6:
                            Reg4.append(vals)
                        elif Rho < 1.5:
                            Reg5.append(vals)
                        else:
                            Reg6.append(vals)
    else:
        for D, T, Prop in DataValues:
            T = float(T)
            Prop = float(Prop)
            D = float(D)

            vals = [D, T, Prop]
            if T < float(critT):
                if D < float(critD):
                    Reg1.append(vals)
                else:
                    Reg2.append(vals)
            else:
                if (
                    T / float(critT) < 1.1
                    and T / float(critT) > 0.98
                    and D / float(critD) > 0.7
                    and D / float(critD) < 1.4
                ):
                    Reg3.append(vals)
                else:
                    Rho = D / float(critD)
                    if Rho < 0.6:
                        Reg4.append(vals)
                    elif Rho < 1.5:
                        Reg5.append(vals)
                    else:
                        Reg6.append(vals)

    return Reg1, Reg2, Reg3, Reg4, Reg5, Reg6


def sampleData(Regions, ratio):
    """
    Sampling of the data regions
    """
    global isothermIndex
    Reg1, Reg2, Reg3, Reg4, Reg5, Reg6 = Regions
    Indexes = []

    ind1, ind2, ind3, ind4, ind5, ind6 = 0, 0, 0, 0, 0, 0

    ind1_size = int(len(Reg1) / ratio)
    if len(Reg1) < ratio:
        ind1_size = len(Reg1)
    ind2_size = int(len(Reg2) / ratio)
    if len(Reg2) < ratio:
        ind2_size = len(Reg2)
    ind3_size = int(len(Reg3) / ratio)
    if len(Reg3) < ratio:
        ind3_size = len(Reg3)
    ind4_size = int(len(Reg4) / ratio)
    if len(Reg4) < ratio:
        ind4_size = len(Reg4)
    ind5_size = int(len(Reg5) / ratio)
    if len(Reg5) < ratio:
        ind5_size = len(Reg5)
    ind6_size = int(len(Reg6) / ratio)
    if len(Reg6) < ratio:
        ind6_size = len(Reg6)

    ind1 = np.random.choice(len(Reg1), size=ind1_size, replace=False)
    ind1 = np.asarray(ind1)
    Reg1s = []
    for i in ind1:
        Reg1s.append(Reg1[i])

    ind2 = np.random.choice(len(Reg2), size=ind2_size, replace=False)
    ind2 = np.asarray(ind2)
    Reg2s = []
    for i in ind2:
        Reg2s.append(Reg2[i])

    ind3 = np.random.choice(len(Reg3), size=ind3_size, replace=False)
    ind3 = np.asarray(ind3)
    Reg3s = []
    for i in ind3:
        Reg3s.append(Reg3[i])

    ind4 = np.random.choice(len(Reg4), size=ind4_size, replace=False)
    ind4 = np.asarray(ind4)
    Reg4s = []
    for i in ind4:
        Reg4s.append(Reg4[i])

    ind5 = np.random.choice(len(Reg5), size=ind5_size, replace=False)
    ind5 = np.asarray(ind5)
    Reg5s = []
    for i in ind5:
        Reg5s.append(Reg5[i])

    ind6 = np.random.choice(len(Reg6), size=ind6_size, replace=False)
    ind6 = np.asarray(ind6)
    Reg6s = []
    for i in ind6:
        Reg6s.append(Reg6[i])

    # ind = np.concatenate((ind1, ind2, ind3, ind4, ind5, ind6))

    Total = 0
    for x in [
        len(Reg1s),
        len(Reg2s),
        len(Reg3s),
        len(Reg4s),
        len(Reg5s),
        len(Reg6s),
    ]:
        Total = Total + x
        Indexes.append(Total)

    if len(Reg1s) > 0:
        Reg1s.sort(key=lambda row: row[2:], reverse=False)
        # PUT BACK  AFTER CHECK ENGLE
        try:
            isothermIndex = next(
                i for i, x in enumerate(Reg1s) if x[2] > float(critT) - 20
            )
        except Exception:
            pass

    Data = []
    [Data.append(x) for x in Reg1s[:]]
    [Data.append(x) for x in Reg2s[:]]
    [Data.append(x) for x in Reg3s[:]]
    [Data.append(x) for x in Reg4s[:]]
    [Data.append(x) for x in Reg5s[:]]
    [Data.append(x) for x in Reg6s[:]]

    return Data, Indexes


def PVT(molecule, sample=False, ratio=5):
    """
    Import pressure-volume-temperature data
    """

    global PVTValues, PVTindexes, isothermIndex
    PVTValues, PVTindexes, isothermIndex = [], [], 0

    rawfile = molecule + "PVT.txt"
    dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i > 1:
            line = line.strip()
            b = line.split()
            PVTValues.append(b[0:3])
            i += 1
        else:
            i += 1

    Reg1, Reg2, Reg3, Reg4, Reg5, Reg6 = regionsOfData(molecule, PVTValues, PVT=True)

    if sample:
        PVTValues, PVTindexes = sampleData((Reg1, Reg2, Reg3, Reg4, Reg5, Reg6), ratio)
    else:
        Reg1.sort(key=lambda row: row[2:], reverse=False)

        if len(Reg1) > 0:
            try:
                isothermIndex = next(
                    i for i, x in enumerate(Reg1) if x[2] > (float(critT) - 30)
                )
            except Exception:
                pass

        Total = 0
        for x in [len(Reg1), len(Reg2), len(Reg3), len(Reg4), len(Reg5), len(Reg6)]:
            Total = Total + x
            PVTindexes.append(Total)

        PVTValues = []
        [PVTValues.append(x) for x in Reg1[:]]
        [PVTValues.append(x) for x in Reg2[:]]
        [PVTValues.append(x) for x in Reg3[:]]
        [PVTValues.append(x) for x in Reg4[:]]
        [PVTValues.append(x) for x in Reg5[:]]
        [PVTValues.append(x) for x in Reg6[:]]


def CP(molecule, sample=False, ratio=5):
    """
    Import isobaric heat capacity data
    """
    global CPValues, CPindexes
    CPValues, CPindexes = [], []

    rawfile = molecule + "CP.txt"
    dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i != 0 and i != 1:
            line = line.strip()
            b = line.split()
            CPValues.append(b[0:3])
        else:
            i += 1

    Reg1, Reg2, Reg3, Reg4, Reg5, Reg6 = regionsOfData(molecule, CPValues)

    if sample:
        CPValues, CPindexes = sampleData((Reg1, Reg2, Reg3, Reg4, Reg5, Reg6), ratio)

    else:
        Total = 0
        for x in [len(Reg1), len(Reg2), len(Reg3), len(Reg4), len(Reg5), len(Reg6)]:
            if x == 0:
                CPindexes.append(0)
            else:
                Total = Total + x
                CPindexes.append(Total)

        CPValues = []
        [CPValues.append(x) for x in Reg1[:]]
        [CPValues.append(x) for x in Reg2[:]]
        [CPValues.append(x) for x in Reg3[:]]
        [CPValues.append(x) for x in Reg4[:]]
        [CPValues.append(x) for x in Reg5[:]]
        [CPValues.append(x) for x in Reg6[:]]


def CV(molecule, sample=False, ratio=5):
    """
    Import isochoric heat capacity
    """
    global CVValues, CVindexes
    CVValues, CVindexes = [], []

    rawfile = molecule + "CV.txt"
    dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i != 0 and i != 1:
            line = line.strip()
            b = line.split()
            CVValues.append(b[0:3])
        else:
            i += 1

    Reg1, Reg2, Reg3, Reg4, Reg5, Reg6 = regionsOfData(molecule, CVValues, CV=True)

    if sample:
        CVValues, CVindexes = sampleData((Reg1, Reg2, Reg3, Reg4, Reg5, Reg6), ratio)
    else:
        Total = 0
        for x in [len(Reg1), len(Reg2), len(Reg3), len(Reg4), len(Reg5), len(Reg6)]:
            if x == 0:
                CVindexes.append(0)
            else:
                Total = Total + x
                CVindexes.append(Total)

        CVValues = []
        [CVValues.append(x) for x in Reg1[:]]
        [CVValues.append(x) for x in Reg2[:]]
        [CVValues.append(x) for x in Reg3[:]]
        [CVValues.append(x) for x in Reg4[:]]
        [CVValues.append(x) for x in Reg5[:]]
        [CVValues.append(x) for x in Reg6[:]]


def SND(molecule, sample=False, ratio=5):
    """
    Import speed of sound data
    """

    global SNDValues, SNDindexes
    SNDValues, SNDindexes = [], []

    rawfile = molecule + "SND.txt"
    dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i != 0 and i != 1:
            line = line.strip()
            b = line.split()
            SNDValues.append(b[0:3])
            i += 1
        else:
            i += 1

    Reg1, Reg2, Reg3, Reg4, Reg5, Reg6 = regionsOfData(molecule, SNDValues)

    if sample:
        SNDValues, SNDindexes = sampleData((Reg1, Reg2, Reg3, Reg4, Reg5, Reg6), ratio)
    else:
        Total = 0
        for x in [len(Reg1), len(Reg2), len(Reg3), len(Reg4), len(Reg5), len(Reg6)]:
            if x == 0:
                SNDindexes.append(0)
            else:
                Total = Total + x
                SNDindexes.append(Total)

        SNDValues = []
        [SNDValues.append(x) for x in Reg1[:]]
        [SNDValues.append(x) for x in Reg2[:]]
        [SNDValues.append(x) for x in Reg3[:]]
        [SNDValues.append(x) for x in Reg4[:]]
        [SNDValues.append(x) for x in Reg5[:]]
        [SNDValues.append(x) for x in Reg6[:]]


def CP0(molecule):
    """
    Import ideal isobaric heat capacity
    """
    rawfile = molecule + "CP0.RAW"
    dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i > 1:
            line = line.strip()
            b = line.split()
            CP0Values.append(b[0:2])
            i += 1
        else:
            i += 1


def DL(molecule):
    """
    Import saturated liquid density
    """
    rawfile = ""
    dataFile = ""
    try:
        rawfile = molecule + "DL.RAW"
        dataFile = open(filename + "/" + rawfile)
    except Exception:
        rawfile = molecule + "DL.txt"
        dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i != 0 and i != 1:
            line = line.strip()
            b = line.split()
            DLValues.append(b[0:2])
            i += 1
        else:
            i += 1


def DV(molecule):
    """
    Import of saturated vapor density
    """
    rawfile = ""
    dataFile = ""
    try:
        rawfile = molecule + "DV.RAW"
        dataFile = open(filename + "/" + rawfile)
    except Exception:
        rawfile = molecule + "DV.txt"
        dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i != 0 and i != 1:
            line = line.strip()
            b = line.split()
            DVValues.append(b[0:2])
            i += 1
        else:
            i += 1


def PV(molecule):
    """
    Import saturated vapor pressure
    """
    rawfile = ""
    dataFile = ""
    try:
        rawfile = molecule + "PV.RAW"
        dataFile = open(filename + "/" + rawfile)
    except Exception:
        rawfile = molecule + "PV.txt"
        dataFile = open(filename + "/" + rawfile)
    i = 0
    for line in dataFile:
        if i != 0 and i != 1:
            line = line.strip()
            b = line.split()
            Values.append(b[0:2])
            i += 1
        else:
            i += 1


# def InSat(molecule):
#     # Cubic-ness
#     Tsat = float(critT) * 0.75
#     DL = SoaveDensity.Sat_Liq_Density(Tsat)
#     DV = SoaveDensity.Sat_Vap_Density(Tsat)
#     for d in np.linspace(DV, DL, 10):
#         Dc = d / float(critD)
#         InSatValues.append([Dc, Tsat])
#     return Tsat, DV, DL
