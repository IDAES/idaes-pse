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
GAMS writer for the regression
"""
import sys
from . import (
    DataImport,
    BasisFunctions,
    GAMSDataWrite,
    DataManipulation,
    # SoaveDensity
)

import numpy as np

global props, sample, sample_ratio
props = []
sample = False
sample_ratio = 5

global num_terms
molecule, data_name, runFile = "", "", ""
constraints = 0
Combination = False
max_time = 1000
num_terms = 12
Rm = 8.314472
critT, critP, critD, M, triple, acc = [0, 0, 0, 0, 0, 0]


# General
def molData(Dfluids, Dmolecule, Ddata_name, Dterms, Dmax_time):
    """Passing of the Data from the main module ::module:: MPEOSDeveloperModule

    :param Dmolecule: Name of molecule of interest
    :type Dmolecule: str
    :param Ddata_name: Name of Data files
    :type Ddata_name: str.
    :param Dmax_time: Running time limit for the GAMS file.
    :type Dmax_time: int.
    """
    global molecule, data_name, max_time, constraints
    global critT, critP, critD, M, triple, acc
    global num_terms
    (critT, critP, critD, M, triple, acc) = Dfluids

    molecule = Dmolecule
    max_time = Dmax_time
    num_terms = Dterms
    data_name = Ddata_name
    # print "Constraints", Dconstraints
    # constraints = Dconstraints;


def openFile(data_name, ending=".gms"):
    """Opens the GAMS file based on the data set

    :param data_name: name of the molecule
    :type data_name: str
    """
    global textFile
    name = molecule + data_name + ending
    print("Creating file: %s" % name)
    try:
        textFile = open(name, "w")
    except Exception:
        print("couldn't open the file Something wrong")
        sys.exit(0)


def closeFile():
    """Closes the GAMS file"""
    textFile.close()
    print("Closing File")


def getTextFile():
    "Returns name of the current textFile"
    global textFile
    return textFile


def setCombination(isCombination):
    "Sets the regression as multiple properties"
    global Combination
    Combination = isCombination


def getRunFile():
    "Returns the gams file name"
    global runFile
    return runFile


def setNumberTerms(numterms):
    "Set number of basis functions allowed in regression"
    global num_terms
    num_terms = numterms


# Data Import
def importData():
    global props
    global sample, sample_ratio

    if "PVT" in props:
        DataImport.PVT(molecule, sample=sample)
    # DataImport.samplePVT6(molecule); #TOL
    if "CV" in props:
        DataImport.CV(molecule)

    if "CP" in props:
        DataImport.CP(molecule, sample=sample)
        # DataImport.sampleCP1(molecule);
    # DataImport.sampleCP2(molecule);
    # DataImport.sampleCP3(molecule);  #TOL
    if "SND" in props:
        DataImport.SND(molecule, sample=sample)


# GDX file
def GenerateGDXGamsFiledtlmv():
    """
    Generates a gams file that creates a data .gdx file
    Creates Combination of PVT, CV, CP, and SND GAMS file. (Titles precoded in).
    """
    global Combination, data_name, props

    Combination = True
    LemJac = False
    BasisFunctions.formCustomBasis(LemJac)
    terms = len(BasisFunctions.coeffs)

    openFile(runFile)

    indexes = {}
    PVTValues, CVValues, CPValues, SNDValues = [], [], [], []
    dataPoints = 0
    if "PVT" in props:
        PVTValues = DataImport.PVTValues
        indexes["PVT"] = DataImport.PVTindexes

    if "CV" in props:
        CVValues = DataImport.CVValues
        indexes["CV"] = DataImport.CVindexes

    if "CP" in props:
        CPValues = DataImport.CPValues
        indexes["CP"] = DataImport.CPindexes

    if "SND" in props:
        SNDValues = DataImport.SNDValues
        indexes["SND"] = DataImport.SNDindexes

    props = indexes.keys()
    dataPoints = [len(PVTValues), len(CVValues), len(CPValues), len(SNDValues)]

    writeGamsHeaderdtl(dataPoints, terms, "", ["PVT", "CV", "CP", "SND"], indexes)

    textFile.write("$Offlisting\n")

    GAMSDataWrite.writeExp(textFile)

    DataToWrite = []
    if "PVT" in props:
        for x in PVTValues:
            manVal = DataManipulation.PVT(x)
            DataToWrite.append(manVal)
        GAMSDataWrite.PVTdt(textFile, DataToWrite, Combination)

    if "CV" in props:
        DataToWrite = []
        for x in CVValues:
            manVal = DataManipulation.CV(x)
            DataToWrite.append(manVal)
        GAMSDataWrite.CVdt(textFile, DataToWrite, Combination)

    if "CP" in props:
        DataToWrite = []
        for x in CPValues:
            manVal = DataManipulation.CP(x)
            DataToWrite.append(manVal)
        GAMSDataWrite.CPdt(textFile, DataToWrite, Combination)

    if "SND" in props:
        DataToWrite = []
        for x in SNDValues:
            manVal = DataManipulation.SND(x)
            DataToWrite.append(manVal)
        GAMSDataWrite.SNDdt(textFile, DataToWrite, Combination)

    textFile.write("$Onlisting\n")

    listProps = ""
    for p in props:
        listProps = "%s %s," % (listProps, p)
        for val, r in zip(indexes[p], ["G", "L", "C", "LD", "MD", "HD"]):
            if val > 0:
                listProps = "%s %s%s," % (listProps, p, r)

    if DataImport.isothermIndex != 0:
        listProps = "%s IT," % (listProps)

    textFile.write(
        "execute_unload '%sData.gdx', p, i, j, itt, delta,tau,d,t,l, %s z;\n"
        % (molecule, listProps)
    )


# GDX Header
def writeGamsHeaderdtl(num_points, terms, kset, pset, regions=None):
    # TODO: What is the kset and pset?
    """Writes the GAMS file Header including the number of terms and data points as well as different thermodynamic properties.

    :param num_points: number of data points
    :type num_points: int
    :param terms: number of basis functions
    :type terms: int
    :param kset: what is this
    :param pset: what is this

    """
    # print num_points
    global Combination
    global props
    # print "Combination", Combination
    textFile.write("$offdigit\n$offsymxref offsymlist\n")
    if len(num_points) >= 1:
        textFile.write("set i /1 * %d/" % max(num_points))
        for data, prop in zip(num_points, pset):
            if prop in props:
                textFile.write("\n %s(i) /1 * %d/" % (prop, data))
                if regions[prop] is not None:
                    prevx = 0
                    # print regions[prop]
                    for x, r in zip(regions[prop], ["G", "L", "C", "LD", "MD", "HD"]):
                        if x > 1:
                            if r == "HD":
                                textFile.write(
                                    "\n %s%s(i) /%d * %d/" % (prop, r, prevx + 1, x)
                                )
                            else:
                                textFile.write(
                                    "\n %s%s(i) /%d * %d/" % (prop, r, prevx + 1, x)
                                )
                            prevx = x
                        elif x > 0:
                            if r == "HD":
                                textFile.write("\n %s%s(i) /%d/" % (prop, r, x))
                            else:
                                textFile.write("\n %s%s(i) /%d/" % (prop, r, x))
                            prevx = x
                    if DataImport.isothermIndex != 0 and prop == "PVT":
                        textFile.write(
                            "\n %s(i) /%d * %d/" % ("IT", 1, DataImport.isothermIndex)
                        )
                    if len(DataImport.InSatValues) > 0 and prop == "PVT":  # 5/9
                        textFile.write(
                            "\n %s(i) /%d * %d/"
                            % ("IS", 1, len(DataImport.InSatValues))
                        )
        textFile.write(";\n")
    else:
        textFile.write("set i /1 * %d/" % max(num_points))
        prevx = 0
        if prop in props and regions is not None:
            for x, r in zip(regions, ["G", "L", "C", "LD", "MD", "HD"]):
                if x > 0:
                    # print x
                    # print regions
                    if r == "HD":
                        textFile.write("\n %s(i) /%d * %d/" % (r, prevx + 1, x))
                    else:
                        textFile.write("\n %s(i) /%d * %d/" % (r, prevx + 1, x))
                    prevx = x
            if DataImport.isothermIndex != 0:
                textFile.write(
                    "\n %s(i) /%d * %d/" % ("IT", 1, DataImport.isothermIndex)
                )
        textFile.write(";\n")

    textFile.write("set j /1* %d/" % terms)
    prevx = 1
    textFile.write(";\n")

    if Combination:
        propString = "Crit"
        for p in pset:
            if propString == "":
                propString = propString + "%s" % p
            else:
                propString = propString + ", %s" % p
        textFile.write("set p /%s/;\n" % propString)
        # textFile.write("parameters z(p, i), xijk(p, i,j, k),crit(p,i,j,k), isoT(p,i,j,k);\n parameters betalo(j), betaup(j);\n")
        textFile.write(
            "parameters z(p, i), d(j), l(j), t(j), delta(p,i), tau(p,i), itt(p,i);\n"
        )
    else:
        textFile.write(
            "parameters z(i),  d(j), l(j), t(j), delta(i), tau(i), itt(i);\n"
        )


# GAMS Shells
def GenerateGamsShell():
    """
    Generates the multiparameter equation of state regression through gams.
    Imports the molecule data.gdx file
    Creates Combination of PVT, CV, CP, and SND GAMS file. (Title precoded in).
    """
    global Combination
    global data_name
    global props

    Combination = True
    BasisFunctions.formCustomBasis()
    terms = len(BasisFunctions.coeffs)

    openFile(runFile)

    indexes = {}
    PVTValues, CVValues, CPValues, SNDValues = [], [], [], []

    if "PVT" in props:
        PVTValues = DataImport.PVTValues
        indexes["PVT"] = DataImport.PVTindexes

    if "CV" in props:
        CVValues = DataImport.CVValues
        indexes["CV"] = DataImport.CVindexes

    if "CP" in props:
        CPValues = DataImport.CPValues
        indexes["CP"] = DataImport.CPindexes

    if "SND" in props:
        SNDValues = DataImport.SNDValues
        indexes["SND"] = DataImport.SNDindexes

    writeGamsShellHeaderB(props, indexes)

    DataToWrite = []
    print(indexes["SND"])

    PVTVar = []
    if "PVT" in props:
        for x in PVTValues:
            manVal = DataManipulation.PVT(x)
            DataToWrite.append(manVal)
        last = 0
        PVTVar = []
        for y in indexes["PVT"]:
            if not y == 0:
                Z = [x[2] for x in DataToWrite[last:y]]
                if len(Z) > 1:
                    AvgZ = sum(Z) / len(Z)
                    Z_1 = [(AvgZ - x) ** 2 for x in Z]
                    VarZ = sum(Z_1) / len(Z)
                    last = y
                    PVTVar.append(VarZ)
                else:
                    PVTVar.append(1)
            else:
                PVTVar.append(0)

    DataToWrite = []
    CVVar = []
    if "CV" in props:
        for x in CVValues:
            manVal = DataManipulation.CV(x)
            DataToWrite.append(manVal)

        last = 0
        for y in indexes["CV"]:
            if not y == 0:
                CV = [x[2] for x in DataToWrite[last:y]]
                AvgCV = sum(CV) / len(CV)
                CV_1 = [(AvgCV - x) ** 2 for x in CV]
                VarCV = sum(CV_1) / len(CV)
                CVaad = [x for x in CV]
                VarCV = np.var(CVaad)
                last = y
                CVVar.append(VarCV)
            else:
                CVVar.append(0)

    DataToWrite = []
    CPVar = []
    if "CP" in props:
        for x in CPValues:
            manVal = DataManipulation.CP(x)
            DataToWrite.append(manVal)

        last = 0
        for y in indexes["CP"]:
            if not y == 0:
                CP = [x[2] for x in DataToWrite[last:y]]
                VarCP = 1
                if len(CP) > 1:
                    AvgCP = sum(CP) / len(CP)
                    CP_1 = [(AvgCP - x) ** 2 for x in CP]
                    VarCP = sum(CP_1) / len(CP)
                last = y
                if VarCP == 0:
                    VarCP = 1
                CPVar.append(VarCP)
            else:
                CPVar.append(0)

    DataToWrite = []
    SNDVar = []
    if "SND" in props:
        for x in SNDValues:
            manVal = DataManipulation.SND(x)
            DataToWrite.append(manVal)

        SNDVar = []
        last = 0
        for y in indexes["SND"]:
            if not y == 0:
                W = [x[2] for x in DataToWrite[last:y]]
                if len(W) > 1:
                    AvgW = sum(W) / len(W)
                    W_1 = [(AvgW - x) ** 2 for x in W]
                    VarW = sum(W_1) / len(W)
                last = y
                if VarW < 1e-5:
                    VarW = 1
                SNDVar.append(VarW)
            else:
                SNDVar.append(0)
        print(SNDVar)

    # Write Variances
    textFile.write("Parameter variance(p,i);\n")
    for var, r in zip(PVTVar, ["G", "L", "C", "LD", "MD", "HD"]):
        if var > 0:
            textFile.write("variance('PVT',PVT%s) = %f;\n" % (r, var))
    for var, r in zip(CVVar, ["G", "L", "C", "LD", "MD", "HD"]):
        if var > 0:
            textFile.write("variance('CV',CV%s) = %f;\n" % (r, var))
    for var, r in zip(CPVar, ["G", "L", "C", "LD", "MD", "HD"]):
        if var > 0:
            textFile.write("variance('CP',CP%s) = %f;\n" % (r, var))
    for var, r in zip(SNDVar, ["G", "L", "C", "LD", "MD", "HD"]):
        if var > 0:
            textFile.write("variance('SND',SND%s) = %f;\n" % (r, var))

    textFile.write("$Offlisting\n")
    writeConstants(terms)
    textFile.write("$Onlisting\n")

    writeDerivatives(props)
    writeBasisFunctions()

    writeEquationsAndVariablesB(props)

    textFile.write("\n\n******** OBJECTIVES\n\n")
    if "PVT" in props and len(props) > 1:
        writeObjectivesB(["PVT"])
    writeObjectivesB(props)
    writeBoundsB(props)

    writeGamsShellFooterB(data_name, False)

    textFile.write("\n\n******** SOLVES\n\n")
    if "PVT" in props and len(props) > 1:
        writeModelB(500, ["PVT"])
    writeCalculateIntermediates(props)
    writeModelB(3600, props)

    writeModelPostEvaluations(props)

    closeFile()


# GAMS Shell imports in GDX
def writeGamsShellHeaderB(pset, regions=None):
    """Writes the GAMS file Header including the number of terms and data points as well as different thermodynamic properties.

    :param pset: list of thermodynamic properties
    :param regions: regions of thermodynamic property data
    """
    textFile.write("$offdigit\n$offsymxref offsymlist\n")

    textFile.write("**************** IMPORT\n")
    textFile.write("Set jl /0*398/\n j(jl) /1*398/\n jld(j) /72*96/;\n\n")

    textFile.write("$Gdxin %sData\n" % molecule)

    setText = "Sets i"
    loadSetText = "$load i"
    for prop in pset:
        if regions[prop] is not None:
            setText = "%s, %s(i)" % (setText, prop)
            loadSetText = "%s, %s" % (loadSetText, prop)
            for x, r in zip(regions[prop], ["G", "L", "C", "LD", "MD", "HD"]):
                if x > 0:
                    setText = "%s, %s%s(%s)" % (setText, prop, r, prop)
                    loadSetText = "%s, %s%s" % (loadSetText, prop, r)

            if DataImport.isothermIndex != 0 and prop == "PVT":
                setText = "%s, IT(i)" % (setText)
                loadSetText = "%s, IT" % (loadSetText)
            # if len(DataImport.InSatValues)>0 and prop =="PVT": # 5/9
            #   setText = "%s, IS(i)"%(setText)
            #   loadSetText = "%s, IS"%(loadSetText);
    textFile.write("%s;\n" % setText)
    textFile.write("%s\n\n" % loadSetText)

    textFile.write("Sets p /'Crit', 'PVT', 'CV', 'CP', 'SND'/;")

    textFile.write("\n\n**************** IMPORT BASIS FUNCTION VALUES\n")
    textFile.write(
        "parameters z(p,i), delta(p,i), tau(p,i), itt(p,i), d(j), l(j), t(j);\n"
    )
    textFile.write("$loaddc z delta tau itt d l t\n")
    textFile.write("$gdxin\n\n")

    textFile.write("Set Crit(i) /1/;\ndelta('Crit','1') = 1;\ntau('Crit','1') = 1;\n\n")

    textFile.write("Parameter data(p,i);\n")
    for prop in pset:
        textFile.write("data('%s', %s(i)) = 1;\n" % (prop, prop))
    textFile.write("data('Crit',Crit(i)) =1;\n\n")


def writeConstants(terms):
    """Writes down the ranges of the fitting Beta value and sets regular to one.
    param terms: number of basis functions
    param terms: int
    """
    textFile.write("Parameter betalo(j), betaup(j);\n")
    betalo = -5
    betaup = 5
    i = "j"
    textFile.write("betalo(%s) = %d ;" % (i, betalo))
    textFile.write("betaup(%s) =  %d ;\n" % (i, betaup))


def writeGamsShellFooterB(data_name, load_in=False):
    """
    Writes the Gam Footer options, model, and display.
    """

    textFile.write("\n\n******** GENERAL OPTIONS\n")

    textFile.write("option optca = 0.0001; \n")
    textFile.write("option optcr = 0.0001 ;\n")
    textFile.write("option limrow = 0, limcol =0 ;\n")
    textFile.write("option solprint=off;\n")
    textFile.write("option reslim = %s;\n" % max_time)
    textFile.write("option sys12 = 1;\n")
    textFile.write("maxexecerror = 10000;\n")

    textFile.write("option rminlp = baron;\n")
    textFile.write("option minlp = baron;\n")
    textFile.write("option decimals = 6;\n")
    textFile.write("option savepoint = 1;\n")


# EQUATIONS
def writeDerivatives(props):
    """
    Calculate and write the important derivatives of the basis functions
    """
    textFile.write("\n\n****** Derivatives\n\n")

    if "PVT" in props:
        textFile.write("parameter PVTvec(i,j), PVTvecCV(i,j);\n")
        textFile.write("PVTvec(PVT(i),j)$(l(j) = 0) = d(j);\n")
        textFile.write(
            "PVTvec(PVT(i),j)$(l(j) <> 0) = d(j)- l(j)*(power(delta('PVT',i),l(j)));\n\n"
        )

        textFile.write("PVTvecCV(PVT(i),j) = t(j)*(t(j)-1);\n")

    if "CV" in props:
        textFile.write("parameter CVvec(i,jl);\n")
        textFile.write("CVvec(CV(i),'0') = itt('CV',i);\n")
        textFile.write("CVvec(CV(i),j) = -t(j)*(t(j)-1);\n\n")

    if "CP" in props:
        textFile.write("parameter CPvecA(i,jl), CPvecB(i,jl), CPvecCV(i,jl);\n\n")

        textFile.write("CPvecA(CP(i),'0')= 1;\n")
        textFile.write("CPvecA(CP(i),j)$(l(j) =0) = 2*d(j) + d(j)*(d(j)-1);\n")
        textFile.write(
            "CPvecA(CP(i),j)$(l(j) <> 0) = 2*(d(j)-l(j)*power(delta('CP',i),l(j))) + (d(j)-l(j)*power(delta('CP',i),l(j)))*(d(j) - 1 - l(j)*power(delta('CP',i),l(j))) - power(l(j),2)*power(delta('CP',i),l(j));\n\n"
        )

        textFile.write("CPvecB(CP(i),'0') = 1;\n")
        textFile.write("CPvecB(CP(i),j)$(l(j) = 0) = d(j)*(1-t(j));\n")
        textFile.write(
            "CPvecB(CP(i),j)$(l(j) <> 0) = (d(j) - l(j)*power(delta('CP',i),l(j)))*(1-t(j));\n\n"
        )

        textFile.write("CPvecCV(CP(i),'0') = 0;\n")
        textFile.write("CPvecCV(CP(i),j) = -t(j)*(t(j)-1);\n")

    if "SND" in props:

        textFile.write("parameter SNDvecA(i,jl), SNDvecB(i,jl), SNDvecCV(i,jl);\n\n")

        textFile.write("SNDvecA(SND(i),'0')= 1;\n")
        textFile.write("SNDvecA(SND(i),j)$(l(j) =0) = 2*d(j) + d(j)*(d(j)-1);\n")
        textFile.write(
            "SNDvecA(SND(i),j)$(l(j) <> 0) = 2*(d(j)-l(j)*power(delta('SND',i),l(j))) + (d(j)-l(j)*power(delta('SND',i),l(j)))*(d(j) - 1 - l(j)*power(delta('SND',i),l(j))) - power(l(j),2)*power(delta('SND',i),l(j));\n\n"
        )

        textFile.write("SNDvecB(SND(i),'0') = 1;\n")
        textFile.write("SNDvecB(SND(i),j)$(l(j) = 0) = d(j)*(1-t(j));\n")
        textFile.write(
            "SNDvecB(SND(i),j)$(l(j) <> 0) = (d(j) - l(j)*power(delta('SND',i),l(j)))*(1-t(j));\n\n"
        )

        textFile.write("SNDvecCV(SND(i),'0') = -itt('SND',i);\n")
        textFile.write("SNDvecCV(SND(i),j) = -t(j)*(t(j)-1);\n\n")

    Z = (float(critP) * 1000 / Rm / float(critT) / float(critD)) - 1.00
    textFile.write("parameter CRITvec(jl) ,CRITvec1(jl), CRITvec2(jl);\n\n")
    textFile.write("CRITvec('0') = %f;\n" % Z)
    textFile.write("CRITvec(j)$(l(j) = 0) = -d(j);\n")
    textFile.write(
        "CRITvec(j)$(l(j) <> 0) = -d(j) + l(j)*power(delta('Crit','1'),l(j));\n\n"
    )

    textFile.write("CRITvec1('0')= 1;\n")
    textFile.write("CRITvec1(j)$(l(j) =0) = d(j)*(d(j)+1);\n")
    textFile.write(
        "CRITvec1(j)$(l(j) <> 0) = d(j)*(d(j)+1) - 2*(d(j)+1)*l(j)*power(delta('Crit','1'),l(j)) - l(j)*(l(j)+1)*power(delta('Crit','1'),l(j)) + power(l(j),2)*power(delta('Crit','1'),(2*l(j)));\n\n"
    )

    textFile.write("CRITvec2('0') = 0;\n")
    textFile.write(
        "CRITvec2(j)$(l(j)=0) = 2*d(j)*(2*d(j)-1) + d(j)*(d(j)-1)*(d(j)-2);\n"
    )
    textFile.write(
        "CRITvec2(j)$(l(j)<>0) = 2*(d(j)-l(j)) + 4*((d(j)-l(j))*(d(j)-1-l(j)) - power(l(j),2)) + (d(j)*(d(j)-1)*(d(j)-2) + (-2*l(j) + 6*d(j)*l(j) - 3 * power(d(j),2)*l(j) - 3*d(j)*power(l(j),2) + 3*power(l(j),2) - power(l(j),3)) + (3*d(j)*power(l(j),2) - 3*power(l(j),2) + 3*power(l(j),3))-power(l(j),3));\n\n"
    )

    textFile.write("parameter ITvec(i,j);\n")
    textFile.write("ITvec(i,j)$(l(j) = 0) = d(j)*(d(j)-1)*(d(j)-2);\n")
    textFile.write(
        "ITvec(i,j)$(l(j) <> 0 ) =d(j)*(d(j)-1)*(d(j)-2) + power(delta('PVT',i),l(j))*(-2*l(j)+6*d(j)*l(j)-3*(d(j)**2)*l(j)-3*d(j)*(l(j)**2)+3*(l(j)**2) - (l(j)**3)) + power(delta('PVT',i),2*l(j))*(3*d(j)*(l(j)**2) - 3*(l(j)**2) + 3*(l(j)**3)) - (l(j)**3)*power(delta('PVT',i),3*l(j));\n\n"
    )


def writeBasisFunctions():
    """
    Write the values of the basis function terms
    """
    textFile.write("\n\n********** BASIS FUNCTIONS\n\n")

    textFile.write("parameter base(p,i,jl);\n\n")

    textFile.write("base(p,i,'0') = 1;\n")
    textFile.write(
        "base(p, i,j)$(data(p,i)=1 and l(j) <> 0) = (delta(p,i)**d(j))* (tau(p,i)**t(j)) * exp(- (delta(p,i)**l(j)));\n"
    )
    textFile.write(
        "base(p,i,j)$(data(p,i)=1 and l(j) = 0) = (delta(p,i)**d(j))* (tau(p,i)**t(j));\n\n"
    )


def writeEquationsAndVariablesB(props):
    """
    Writes multiple thermodynamic parameter equations and constants.
    :param props: Array containing the available properties.
    :type props: array
    """
    global constraints
    # properties = len(props);

    textFile.write("\n\n********  EQUATIONS and VARIABLES\n\n")

    textFile.write("binary variable y(j);\n")
    textFile.write("variables beta(jl), sse;\n")

    eqtnList = ""
    posVars = ""
    Vars = "gammaCrit, gammaCrit1, gammaCrit2"
    if "PVT" in props:
        eqtnList = "%s gammaAPVT," % eqtnList
        Vars = "%s, gammaPVT(i)" % Vars
        if DataImport.isothermIndex != 0:
            eqtnList = "%s ITgamma," % eqtnList
            Vars = "%s,  gammaIT(i)" % Vars
    if "CV" in props:
        eqtnList = "%s gammaACV," % eqtnList
        if posVars == "":
            posVars = "%s gammaCV(i)" % posVars
        else:
            posVars = "%s, gammaCV(i)" % posVars
    if "CP" in props:
        eqtnList = "%s gammaACP, gammaBCP, gammaCCP, " % eqtnList
        Vars = "%s, gammaCPB(i), gammaCPC(i)" % Vars
        if posVars == "":
            posVars = "%s gammaCPA(i)" % posVars
        else:
            posVars = "%s, gammaCPA(i)" % posVars
    if "SND" in props:
        eqtnList = "%s gammaASND, gammaBSND, gammaCSND, omegaSND, psiSND" % eqtnList
        Vars = "%s, gammaSNDA(i), gammaSNDB(i), SNDomega(i)" % Vars
        if posVars == "":
            posVars = "%s gammaSNDC(i), SNDpsi(i)" % posVars
        else:
            posVars = "%s, gammaSNDC(i), SNDpsi(i)" % posVars

    textFile.write("variables %s;\n" % Vars)
    textFile.write("positive variables %s;\n" % posVars)
    textFile.write("beta.lo(j) = betalo(j); beta.up(j) = betaup(j);\n")
    textFile.write(
        "equations eq1, eq2, eq3, ITgamma, Critgamma, Critgamma1, Critgamma2,  \n  %s \n obj;\n"
        % eqtnList
    )
    # eq5

    textFile.write("beta.FX('0')=1;\n")
    textFile.write("eq1(j).. betalo(j)*y(j) =l= beta(j);\n")
    textFile.write("eq2(j).. beta(j) =l= betaup(j)*y(j);\n")
    textFile.write("eq3.. sum(j, y(j)) =e= %d;\n\n" % num_terms)

    # Write property equations.
    if "PVT" in props:
        textFile.write(
            "gammaAPVT(PVT(i)).. gammaPVT(i) =e= sum(j,beta(j)*PVTvec(i,j)*base('PVT',i,j));\n\n"
        )
        if DataImport.isothermIndex != 0:
            textFile.write(
                "ITgamma(IT(i)).. gammaIT(i) =e= sum(j, beta(j)*ITvec(i,j)*base('PVT',i,j));\n"
            )
    if "CV" in props:
        textFile.write(
            "gammaACV(CV(i)).. gammaCV(i) =e= sum(jl, beta(jl)*CVvec(i,jl)*base('CV',i,jl));\n\n"
        )
    if "CP" in props:
        textFile.write(
            "gammaACP(CP(i)).. gammaCPA(i) =e= sum(j,beta(j)*CPvecCV(i,j)*base('CP',i,j));\n"
        )
        textFile.write(
            "gammaBCP(CP(i)).. gammaCPB(i) =e= sum(jl,  beta(jl)*CPvecB(i,jl)*base('CP',i,jl));\n"
        )
        textFile.write(
            "gammaCCP(CP(i)).. gammaCPC(i) =e= sum(jl,beta(jl)*CPvecA(i,jl)*base('CP',i,jl));\n\n"
        )
        # textFile.write("omegaCP(CP(i)).. CPomega(i) =e= gammaCPA(i)*gammaCPC(i);\n")
        # textFile.write("psiCP(CP(i)).. CPpsi(i) =e= gammaCPB(i)*gammaCPB(i);\n\n")
    if "SND" in props:
        textFile.write(
            "gammaASND(SND(i)).. gammaSNDA(i) =e= sum(jl,beta(jl)*SNDvecA(i,jl)*base('SND',i,jl));\n"
        )
        textFile.write(
            "gammaBSND(SND(i)).. gammaSNDB(i) =e= sum(jl,beta(jl)*SNDvecB(i,jl)*base('SND',i,jl));\n"
        )
        textFile.write(
            "gammaCSND(SND(i)).. gammaSNDC(i) =e= sum(jl,beta(jl)*SNDvecCV(i,jl)*base('SND',i,jl));\n"
        )
        textFile.write(
            "omegaSND(SND(i)).. SNDomega(i) =e= gammaSNDA(i)*gammaSNDC(i);\n"
        )
        textFile.write(
            "psiSND(SND(i))..   SNDpsi(i) =e= gammaSNDB(i)*gammaSNDB(i);\n\n"
        )

    textFile.write("Critgamma..  gammaCrit =e= sum(jl, beta(jl)*Critvec(jl));\n")
    textFile.write("Critgamma1.. gammaCrit1 =e= sum(jl, beta(jl)*Critvec1(jl));\n")
    textFile.write("Critgamma2.. gammaCrit2 =e= sum(jl, beta(jl)*Critvec2(jl));\n\n")


def writeObjectivesB(props):
    """
    Writes the objective equations for the different combination of properites and constraints
    """
    propsString = ""
    for p in props:
        propsString = "%s%s" % (propsString, p)

    objString = "obj%s.. sse =e= " % propsString

    textFile.write("EQUATION obj%s;\n" % propsString)

    if "PVT" in props:
        objString = (
            "%s + sum(PVT(i),power(z('PVT',i) - gammaPVT(i),2)/variance('PVT',i))"
            % objString
        )
    if "CV" in props:
        objString = (
            "%s + sum(CV(i), power( z('CV',i) - gammaCV(i),2)/variance('CV',i))"
            % objString
        )
    if "CP" in props:
        objString = (
            "%s + sum(CP(i),power(z('CP',i) - gammaCPA(i) - power(gammaCPB(i),2)/gammaCPC(i),2)/variance('CP',i))"
            % objString
        )
    if "SND" in props:
        objString = (
            "%s + sum(SND(i),power(z('SND',i)*gammaSNDC(i) - SNDomega(i) - SNDpsi(i) ,2)/variance('SND',i))"
            % objString
        )
    if DataImport.isothermIndex != 0:
        objString = "%s + 100*sum(IT(i),power(gammaIT(i),2)) " % objString

    textFile.write(
        "%s  + 5* power( gammaCrit,2) + 5 * power(gammaCrit1,2)+ 5* power(gammaCrit2,2);\n\n"
        % objString
    )

    textFile.write("\n\n")


def writeBoundsB(props):
    """
    Writes bounds on the variables
    """
    textFile.write("\n\n****** Property Bounds\n\n")

    # textFile.write("y.fx(jld) = 0;\n")
    textFile.write("beta.up(j) = 5;\n")
    textFile.write("beta.lo(j) = -5;\n\n")

    if "PVT" in props:
        PVTVals = [DataManipulation.PVT(x)[2] for x in DataImport.PVTValues]
        minZ = np.min(PVTVals)
        maxZ = np.max(PVTVals)
        print("Min", minZ, "Max", maxZ)

        textFile.write("gammaPVT.up(PVT)$(z('PVT',PVT)>0) = %f;\n" % (maxZ * 10))
        textFile.write("gammaPVT.lo(PVT)$(z('PVT',PVT)>0) = %f;\n\n" % (-1))
        # textFile.write("gammaPVT.up(PVT)$(z('PVT',PVT)<0) = 0;\n");
        # textFile.write("gammaPVT.lo(PVT)$(z('PVT',PVT)<0) = -1;\n\n");

    if "CV" in props:
        textFile.write("*gammaCV.up(i) = z('CV',i)*1.05;\n")
        textFile.write("*gammaCV.lo(i) = z('CV',i)*0.95;\n\n")

    if "CP" in props:
        textFile.write("*gammaCPA.lo(i) = 0;\n")
        textFile.write("*gammaCPA.up(i) = 20;\n")
        textFile.write("*gammaCPB.lo(i) = 0;\n")
        textFile.write("*gammaCPB.up(i) = 20;\n")
        textFile.write("*gammaCPC.lo(i) = 1E-16;\n")
        textFile.write("*gammaCPC.up(i) = 40;\n\n")

    if "SND" in props:
        textFile.write("*gammaSNDA.up(i) = 50;\n")
        textFile.write("*gammaSNDA.lo(i) = 0;\n")
        textFile.write("*gammaSNDB.up(i) = 20;\n")
        textFile.write("*gammaSNDB.lo(i) = 0;\n")
        textFile.write("*gammaSNDC.up(i) = 20;\n")
        textFile.write("*gammaSNDC.lo(i) = 0;\n\n")

        textFile.write("*SNDomega.lo(i) = 0;\n")
        textFile.write("*SNDomega.up(i)  = 120;\n")
        textFile.write("*SNDpsi.lo(i) = 0;\n")
        textFile.write("*SNDpsi.up(i) = 100;\n\n")


def writeModelB(reslim, props):
    """
    Writes the GAMS model definition
    """
    propsString = ""
    for p in props:
        propsString = "%s%s" % (propsString, p)

    textFile.write("option reslim = %i;\n" % reslim)
    modelName = "sse_%s" % propsString
    # modelString = "model sse_%s /"%propsString;

    eqtnList = "eq1, eq2, eq3, Critgamma, Critgamma1, Critgamma2, obj%s" % propsString
    if "PVT" in props:
        eqtnList = "%s, gammaAPVT" % eqtnList
    if "CV" in props:
        eqtnList = "%s, gammaACV" % eqtnList
    if "CP" in props:
        eqtnList = "%s, gammaACP, gammaBCP, gammaCCP" % eqtnList
    if "SND" in props:
        eqtnList = "%s, gammaASND, gammaBSND, gammaCSND" % eqtnList
    if DataImport.isothermIndex != 0:
        eqtnList = "%s, ITgamma " % eqtnList

    modelString = "model %s / %s/;\n" % (modelName, eqtnList)
    solveString = "solve %s minimizing sse using minlp;\n" % (modelName)
    displayString = "display beta.l, y.l;\n\n"

    textFile.write(modelString)
    textFile.write(solveString)
    textFile.write(displayString)


def writeCalculateIntermediates(props):
    """
    Calculate intermediate values of properties not in the regression
    """

    # Write property equations.
    if "PVT" in props:
        textFile.write(
            "gammaPVT.l(i) = sum(j,beta.l(j)*PVTvec(i,j)*base('PVT',i,j));\n\n"
        )
        if DataImport.isothermIndex != 0:
            textFile.write(
                "gammaIT.l(i) = sum(j, beta.l(j)*ITvec(i,j)*base('PVT',i,j));\n"
            )
    if "CV" in props:
        textFile.write(
            "gammaCV.l(i) = sum(jl, beta.l(jl)*CVvec(i,jl)*base('CV',i,jl));\n\n"
        )
    if "CP" in props:
        textFile.write(
            "gammaCPA.l(i) = sum(j,beta.l(j)*CPvecCV(i,j)*base('CP',i,j));\n"
        )
        textFile.write(
            "gammaCPB.l(i) = sum(jl,  beta.l(jl)*CPvecB(i,jl)*base('CP',i,jl));\n"
        )
        textFile.write(
            "gammaCPC.l(i) = sum(jl,beta.l(jl)*CPvecA(i,jl)*base('CP',i,jl));\n\n"
        )
        # textFile.write("omegaCP(CP(i)).. CPomega(i) =e= gammaCPA(i)*gammaCPC(i);\n")
        # textFile.write("psiCP(CP(i)).. CPpsi(i) =e= gammaCPB(i)*gammaCPB(i);\n\n")
    if "SND" in props:
        textFile.write(
            "gammaSNDA.l(i) = sum(jl,beta.l(jl)*SNDvecA(i,jl)*base('SND',i,jl));\n"
        )
        textFile.write(
            "gammaSNDB.l(i) = sum(jl,beta.l(jl)*SNDvecB(i,jl)*base('SND',i,jl));\n"
        )
        textFile.write(
            "gammaSNDC.l(i) = sum(jl,beta.l(jl)*SNDvecCV(i,jl)*base('SND',i,jl));\n"
        )
        textFile.write("SNDomega.l(i) = gammaSNDA.l(i)*gammaSNDC.l(i);\n")
        textFile.write("SNDpsi.l(i) = gammaSNDB.l(i)*gammaSNDB.l(i);\n\n")

    # textFile.write("Critgamma..  gammaCrit =e= sum(jl, beta(jl)*Critvec(jl));\n")
    # textFile.write("Critgamma1.. gammaCrit1 =e= sum(jl, beta(jl)*Critvec1(jl));\n")
    # textFile.write("Critgamma2.. gammaCrit2 =e= sum(jl, beta(jl)*Critvec2(jl));\n\n")


def writeModelPostEvaluations(props):
    """
    Calculates thermodynamic values after the regression
    """

    textFile.write("\n\n****************Post Evaluation\n\n")

    textFile.write("display gammaCrit.l, gammaCrit1.l, gammaCrit2.l;\n")

    if DataImport.isothermIndex != 0:
        textFile.write("display gammaIT.l;\n")

    if "PVT" in props:
        textFile.write("display gammaPVT.l;\n")
    if "CV" in props:
        textFile.write("display gammaCV.l;\n")

    if "CP" in props:
        textFile.write("display gammaCPA.l, gammaCPB.l, gammaCPC.l;\n")

        textFile.write("parameter mingammaCPA, mingammaCPB, mingammaCPC;\n")
        textFile.write("mingammaCPA = smin(CP(i), gammaCPA.l(i));\n")
        textFile.write("mingammaCPB = smin(CP(i), gammaCPB.l(i));\n")
        textFile.write("mingammaCPC = smin(CP(i), gammaCPC.l(i));\n")
        textFile.write("display mingammaCPA, mingammaCPB, mingammaCPC;\n\n")

        textFile.write("parameter maxgammaCPA, maxgammaCPB, maxgammaCPC;\n")
        textFile.write("maxgammaCPA = smax(CP(i), gammaCPA.l(i));\n")
        textFile.write("maxgammaCPB = smax(CP(i), gammaCPB.l(i));\n")
        textFile.write("maxgammaCPC = smax(CP(i), gammaCPC.l(i));\n")
        # textFile.write("maxCPomega = smax(CP(i), CPomega.l(i));\n");
        # textFile.write("maxCPpsi = smax(CP(i), CPpsi.l(i));\n");
        textFile.write("display maxgammaCPA, maxgammaCPB, maxgammaCPC;\n\n")

    if "SND" in props:
        textFile.write("display gammaSNDA.l, gammaSNDB.l, gammaSNDC.l;\n")

        textFile.write(
            "parameter mingammaSNDA, mingammaSNDB, mingammaSNDC, minSNDomega, minSNDpsi;\n"
        )
        textFile.write("mingammaSNDA = smin(SND(i), gammaSNDA.l(i));\n")
        textFile.write("mingammaSNDB = smin(SND(i), gammaSNDB.l(i));\n")
        textFile.write("mingammaSNDC = smin(SND(i), gammaSNDC.l(i));\n")
        textFile.write("minSNDomega = smin(SND(i), SNDomega.l(i));\n")
        textFile.write("minSNDpsi = smin(SND(i), SNDpsi.l(i));\n")
        textFile.write(
            "display mingammaSNDA, mingammaSNDB, mingammaSNDC, minSNDomega, minSNDpsi;\n\n"
        )

        textFile.write(
            "parameter maxgammaSNDA, maxgammaSNDB, maxgammaSNDC, maxSNDomega, maxSNDpsi;\n"
        )
        textFile.write("maxgammaSNDA = smax(SND(i), gammaSNDA.l(i));\n")
        textFile.write("maxgammaSNDB = smax(SND(i), gammaSNDB.l(i));\n")
        textFile.write("maxgammaSNDC = smax(SND(i), gammaSNDC.l(i));\n")
        textFile.write("maxSNDomega = smax(SND(i), SNDomega.l(i));\n")
        textFile.write("maxSNDpsi = smax(SND(i), SNDpsi.l(i));\n")
        textFile.write(
            "display maxgammaSNDA, maxgammaSNDB, maxgammaSNDC, maxSNDomega, maxSNDpsi;\n\n"
        )

    errString = ""
    TotErrString = ""
    for p in props:
        errString = "%s err%s," % (errString, p)
        TotErrString = "%s + err%s" % (TotErrString, p)
    errString = "%s totErr;\n" % errString
    textFile.write("parameter %s" % errString)

    if "PVT" in props:
        textFile.write(
            "errPVT = sum(PVT,power(z('PVT',PVT) - gammaPVT.l(PVT),2)/variance('PVT',PVT));\n"
        )
    if "CV" in props:
        textFile.write(
            "errCV =  sum(CV(i), power( z('CV',i) - gammaCV.l(i),2)/variance('CV',i));\n"
        )
    if "CP" in props:
        textFile.write(
            "errCP =  sum(CP(i),power(z('CP',i) - gammaCPA.l(i) - power(gammaCPB.l(i),2)/gammaCPC.l(i),2)/variance('CP',i));\n"
        )
    if "SND" in props:
        textFile.write(
            "errSND=  sum(SND(i),power(z('SND',i)*gammaSNDC.l(i) - SNDomega.l(i) - SNDpsi.l(i) ,2)/variance('SND',i));\n"
        )
    if DataImport.isothermIndex != 0:
        TotErrString = "%s + 100*sum(IT(i),power(gammaIT.l(i),2)) " % TotErrString

    textFile.write(
        "totErr = %s  + 5* power( gammaCrit.l,2) + 5 * power(gammaCrit1.l,2)+ 5* power(gammaCrit2.l,2);\n"
        % TotErrString
    )

    textFile.write("display %s" % errString)
