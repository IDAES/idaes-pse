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
HELMholtz Energy Thermodynamics (HELMET)

Main capabilities of HELMET
default HELMET use
"""

__author__ = "Marissa Engle <mengle@andrew.cmu.edu>"

import os
import platform
import subprocess


from idaes.apps import alamopy_depr as alamopy
from . import AncillaryEquations  # , Certainty
from . import Plotting, DataImport, DataManipulation
from . import GAMSWrite, BasisFunctions

from matplotlib import cm


# global R
R = 8.314472  # kJ mol^-1 K^-1

# global molecule, filename, gamsname, data_name
molecule, filename, gamsname, data_name = None, None, None, None

# global sample, sample_ratio
sample = False
sample_ratio = 5

# global critT, critP, critD, M, triple, acc
critT, critP, critD, M, triple, acc = 0, 0, 0, 0, 0, 0

# global max_time, num_terms
max_time = 500
num_terms = 12

# global props
props = []

# global has_alamo
has_alamo = None

# global flag_dirty
flag_dirty = False


def initialize(**kwargs):
    """
    filename - location of data
    gamsname - name of the gams file made
    molecule - name of the molecule/compound
    data_name - name of the data
    fluid data - [critT, critP, critD, M, triple, acentric factor]
    R - gas constant value
    """

    global R, molecule, filename, gamsname, data_name
    global critT, critP, critD, M, triple, acc
    global max_time, num_terms, props
    global sample, sample_ratio
    global flag_dirty
    global has_alamo

    k_dict = {
        "R": R,
        "filename": filename,
        "gamsname": gamsname,
        "molecule": molecule,
        "fluid_data": (critT, critP, critD, M, triple, acc),
        "max_time": max_time,
        "num_terms": num_terms,
        "props": props,
        "sample": 1,
    }

    for arg in k_dict:
        if arg in kwargs:
            if arg == "R":
                R = kwargs[arg]
            elif arg == "filename":
                filename = kwargs[arg]
            elif arg == "gamsname":
                gamsname = kwargs[arg]
            elif arg == "molecule":
                molecule = kwargs[arg]
            elif arg == "fluid_data":
                (critT, critP, critD, M, triple, acc) = kwargs[arg]
            elif arg == "max_time":
                max_time = kwargs[arg]
            elif arg == "num_terms":
                num_terms = kwargs[arg]
            elif arg == "props":
                props = kwargs[arg]
            elif arg == "sample":
                sample = True
                sample_ratio = kwargs[arg]
            else:
                raise Exception("Not a keyword argument")

    if has_alamo is None:
        has_alamo = alamopy.multos.has_alamo()
        if not has_alamo:
            print("No ALAMO software found.")

    updateModelSettings()

    flag_dirty = True


def updateModelSettings():
    """
    Settings of the model based on the chemical passed to the
    different python methods
    """

    global R, molecule, filename, gamsname, data_name
    global critT, critP, critD, M, triple, acc
    global max_time, num_terms, props
    global flag_dirty

    # SoaveDensity.molData((critT, critP, critD, M, triple, acc), molecule, R)
    Plotting.molData((critT, critP, critD, M, triple, acc), molecule, R)
    Plotting.props = props
    DataImport.molData((critT, critP, critD, M, triple, acc), R)
    DataImport.filename = filename
    DataManipulation.molData((critT, critP, critD, M, triple, acc), molecule, R)
    AncillaryEquations.molecule = molecule
    AncillaryEquations.max_time = max_time
    GAMSWrite.molData(
        (critT, critP, critD, M, triple, acc), molecule, data_name, num_terms, max_time
    )
    BasisFunctions.molData((critT, critP, critD, M, triple, acc), molecule, R)

    flag_dirty = False


def prepareAncillaryEquations(plot=False, keepFiles=False):
    """
    Develops ancillary equations of state using ALAMOPY
        DL - saturated liquid density
        DV - saturated vapor density
        PV - vapor pressure

    Dependent on ALAMO
    """
    global has_alamo

    if has_alamo:
        AncillaryEquations.DL()
        AncillaryEquations.DV()
        AncillaryEquations.PV()

        if plot:
            Plotting.viewAnc()

        if not keepFiles:
            for p in ["DL", "DV", "PV"]:
                os.remove("%s%s" % (molecule, p))
                os.remove("%s%s.lst" % (molecule, p))
    else:
        if plot:
            Plotting.viewAnc()
        print("Couldn't regress ancillary equations. ALAMO executable not found")


def viewPropertyData():
    """
    Plot imported data
    """
    Plotting.viewData()


def setupRegression(numTerms=14, gams=False, pyomo=False):
    """
    setup gams regression
    """
    global props, sample, sample_ratio
    GAMSWrite.num_terms = numTerms
    GAMSWrite.props = props
    GAMSWrite.sample = sample
    GAMSWrite.sample = sample_ratio

    GAMSWrite.importData()

    GAMSWrite.runFile = "gdx"
    GAMSWrite.GenerateGDXGamsFiledtlmv()

    GAMSWrite.runFile = "main"
    GAMSWrite.GenerateGamsShell()


def runRegression(gams=False, pyomo=False):
    """
    Runs the gdx and main regression gams file
    """

    GAMSWrite.runFile = "gdx"
    command = "gams %s%s.gms" % (molecule, GAMSWrite.runFile)
    process = subprocess.check_call(command, shell=True)

    GAMSWrite.runFile = "main"
    command = "gams %s%s.gms" % (molecule, GAMSWrite.runFile)
    process = subprocess.check_call(command, shell=True)


def getFlag():
    """
    Returns flag, marks a change in the construction of the model
    """
    return flag_dirty


def viewResults(lstFile=None, plot=False, report=False, surface=cm.coolwarm):
    """
    Plot results from gams or pyomo
        lstFile - gams listing file
        surface - colormapping color eg. cm.coolwarm
    """
    Plotting.sseCombo(lstFile=lstFile, plot=plot, report=report, surface=surface)


def viewMultResults(lstFile, numTerms=0):
    """
    View mutliple results from a lst file
    """
    # PYLINT-TODO-FIX the multSSECombo function doesn't seem exist in the Plotting module
    # pylint: disable=no-member
    Plotting.multSSECombo(lstFile, numTerms)


def deletefile(*fname):
    """
    Deletes files
    """
    tos = platform.platform()
    if "Windows" in tos:
        for name in fname:
            os.system("del " + name)
    else:
        for name in fname:
            os.system("rm " + name)
