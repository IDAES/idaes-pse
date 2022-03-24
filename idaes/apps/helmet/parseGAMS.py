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
Parses and prints the solutions of the multiparameter equation of state solution
"""
# stdlib
import re

# pkg
from . import BasisFunctions

indexes = []
betas = []


def parser(filename, num=2):
    """
    Parse solution files for the muliparameter equation of state
    """
    global indexes, betas
    print("Parsing %s" % filename)
    dataFile = open(filename)

    # Y variables regex
    regexY = r"(?<=(VARIABLE y.L))(\s*)(([0-9,. \s]+)(\n))"
    regexB = r"(?<=(VARIABLE beta.L))(\s*)(([-E0-9,. \s]+)(\n))"

    dataFile = open(filename)

    Yb = re.findall(regexY, dataFile.read())

    indexes = []
    for ybi in Yb:
        indexes_sub = []
        ybi = " ".join(ybi[2].split())
        groups = ybi.split(",")
        for i in groups:
            trimmed = i.strip()
            splits = trimmed.split(" ")
            indexes_sub.append(int(splits[0]))
            if len(splits) > 3:
                indexes_sub.append(int(splits[2]))
        indexes.append(indexes_sub)

    dataFile2 = open(filename)
    Yb = re.findall(regexB, dataFile2.read())

    betas = []
    ind = 0
    for ybi in Yb:
        betas_sub = []
        ybi = " ".join(ybi[2].split())
        groups = ybi.split(",")
        for i in groups:
            i = i.strip()
            splits = i.split(" ")
            if len(splits) <= 2:
                beta = splits[1]
                index = int(splits[0])
            else:
                beta = splits[1]
                index = int(splits[0])
                if index in indexes[ind]:
                    betas_sub.append(beta)
                beta = splits[3]
                index = int(splits[2])

            if index in indexes[ind]:
                betas_sub.append(beta)
        ind = ind + 1
        betas.append(betas_sub)


def getIndexes():
    """
    Returns indexes of the basis function terms
    """
    global indexes
    return indexes


def getBetas():
    """
    Returns the weights of the basis functions
    """
    global betas
    return betas


def writeTerm(index):
    """
    Writes the basis function term with the given index
    """
    global coeffs
    coeffs = []
    eqtn = ""
    term_index = 0
    for i in range(1, 9):
        for j in range(1, 13):  # 12
            term_index += 1
            coeffs.append([i, j / 8.0, 0])
            if index == term_index:
                eqtn = eqtn + " %s^%d * %s^%.3f" % ("D", i, "T", j / 8.0)
    for i in range(1, 6):
        for j in range(1, 24):  # 24
            term_index += 1
            coeffs.append([i, j / 8.0, 1])
            if index == term_index:
                eqtn = eqtn + " %s^%d * %s^%.3f * exp(-D^1)" % ("D", i, "T", j / 8.0)
    for i in range(1, 6):
        for j in range(1, 30):  # 24
            term_index += 1
            coeffs.append([i, j / 8.0, 2])
            if index == term_index:
                eqtn = eqtn + " %s^%d * %s^%.3f * exp(-D^2)" % ("D", i, "T", j / 8.0)
    for i in range(2, 5):
        for j in range(24, 38):  # 38
            term_index += 1
            coeffs.append([i, j / 2.0, 3])
            if index == term_index:
                eqtn = eqtn + " %s^%d * %s^%.3f * exp(-D^3)" % ("D", i, "T", j / 8.0)
    return eqtn


def writeEquation(Y, Beta=None):
    """
    Write full multiparameter equation
    """
    global coeffs
    global indexes

    eqtn = "Equation for Helmholtz Energy \n ="
    coeffs = BasisFunctions.coeffs
    indexes = Y
    lastInd = indexes[-1]
    for ind, b in zip(indexes, Beta):
        d, t, c, m = coeffs[ind - 1]
        if not c == 0:
            eqtn = eqtn + "%f * %s^%d * %s^%.3f * exp(-D^%d)" % (b, "D", d, "T", t, c)
        else:
            eqtn = eqtn + "%f * %s^%d * %s^%.3f" % (b, "D", d, "T", t)
        if not ind == lastInd:
            eqtn = eqtn + " +"

    print(eqtn)
