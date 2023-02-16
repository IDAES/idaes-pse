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
"""HELMET Plotting capabilities"""

import scipy.stats as stats
import numpy as np
import math

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib

from . import (
    DataImport,
    DataManipulation,
    BasisFunctions,
    parseGAMS,
    AncillaryEquations,
)

font = {"size": 12}
matplotlib.rc("font", **font)
not_parsed = True

CVValuesn, SNDValues = [], []
critT, critD, critP, acc, R, M, Rm = [0, 0, 0, 0, 0, 0, 0]
molecule = ""
triple = 0

global props
props = ["PVT", "CV", "CP", "SND"]
markers = [
    ".",
    ",",
    "o",
    "v",
    "^",
    "<",
    ">",
    "1",
    "2",
    "3",
    "4",
    "8",
    "s",
    "p",
    "P",
    "*",
    "x",
    "X",
    "D",
    "d",
    "|",
    "_",
]


def molData(fluidData, Dmolecule, RVal):
    """
    Shared data about the molecule and ideas gas constant R
    """
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    molecule = Dmolecule
    (critT, critP, critD, M, triple, acc) = fluidData
    R = RVal
    Rm = RVal


# View Ancillary Regressions
def viewAnc():
    """Plots all the ancillary equations for saturated density and vapor pressure"""
    plotPV()
    plotDV()
    plotDL()


def plotPV():
    """Plot vapor pressure"""
    global critT, critP, critD, M, triple, acc, R, Rm, S

    DataImport.PV(molecule)
    Values = DataImport.Values
    DataToWrite = []
    T = []
    D = []
    for x in Values:
        D.append(float(x[0]))
        T.append(float(x[1]))
        manVal = DataManipulation.PV(x)
        DataToWrite.append(manVal)

    Dcalc = []

    for i, t in zip(DataToWrite, T):
        PVf = AncillaryEquations.getPV()
        PV = np.exp(PVf.f(i[0], i[1])) * critP
        Dcalc.append(PV)

    fig = plt.figure()
    fig.suptitle("Vapor Pressure Data", fontsize=20, fontweight="bold")
    ax = fig.add_subplot(111)
    ax.scatter(T, D, s=10, edgecolor="r", facecolor="w", marker="o", label="data")
    ax.scatter(T, Dcalc[:], s=10, c="g", marker="o", label="calculated")
    ax.set_ylabel("Pressure (MPa)", fontsize=20)
    ax.set_xlabel("Temperature (K)", fontsize=20)
    ax.legend(loc="upper right", fontsize=20)
    plt.show()


def plotDV():
    """Plot saturated vapor density"""
    global critT, critP, critD, M, triple, acc, R, Rm, S

    DataImport.DV(molecule)
    Values = DataImport.DVValues
    DataToWrite = []
    T = []
    D = []
    for x in Values:
        D.append(float(x[0]))
        T.append(float(x[1]))

        manVal = DataManipulation.DV(x)
        DataToWrite.append(manVal)

    Dcalc = []

    DV = AncillaryEquations.getDV()
    for i, t in zip(DataToWrite, T):
        Dc = float(critD)
        Ts = 1 - float(t) / float(critT)
        Dcalc.append(np.exp(DV.f(Ts)) * Dc)

    fig = plt.figure()
    fig.suptitle("Saturated Vapor Density", fontsize=20, fontweight="bold")
    ax = fig.add_subplot(111)
    ax.scatter(T, D, s=10, edgecolor="r", facecolor="w", marker="o", label="data")
    ax.scatter(T, Dcalc[:], s=10, c="g", marker="o", label="calculated")
    ax.set_ylabel("Density (dM)", fontsize=20)
    ax.set_xlabel("Temperature (K)", fontsize=20)
    ax.legend(loc="upper right", fontsize=20)
    plt.show()


def plotDL():
    """Plot saturate liquid density"""
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    DataImport.DL(molecule)
    Values = DataImport.DLValues

    DataToWrite = []
    T = []
    D = []
    for x in Values:
        D.append(float(x[0]))
        T.append(float(x[1]))
        manVal = DataManipulation.DL(x)
        DataToWrite.append(manVal)

    Dcalc = []

    for i, t in zip(DataToWrite, T):
        Dc = critD
        DL = AncillaryEquations.getDL()
        Dcalc.append((DL.f(float(i[0])) + 1) * (Dc))

    fig = plt.figure()
    fig.suptitle("Saturated Liquid Density", fontsize=20, fontweight="bold")
    ax = fig.add_subplot(111)
    ax.scatter(T, D, s=10, edgecolor="r", facecolor="w", marker="o", label="data")
    ax.scatter(T, Dcalc[:], s=10, c="g", marker="o", label="calculated")
    ax.set_ylabel("Density (dM)", fontsize=20)
    ax.set_xlabel("Temperature (K)", fontsize=20)
    ax.legend(loc="upper right", fontsize=20)
    plt.show()


# View Data
def viewData():
    """
    View imported data
    """
    global props

    if "PVT" in props:
        plotPVT()
    if "CV" in props:
        plotCV()
    if "CP" in props:
        plotCP()
    if "SND" in props:
        plotSND()


def plotPVT():
    """
    Plot Pressure-Volume-Temperature data
    """
    DataImport.PVT(molecule)
    PVTValues = DataImport.PVTValues

    PVTnp = np.asarray(PVTValues)
    T = PVTnp[:, 2]
    D = PVTnp[:, 1]
    P = PVTnp[:, 0]

    fig = plt.figure()
    fig.suptitle("Pressure Data", fontsize=20, fontweight="bold")

    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(T, D, P, s=20, c=D, marker="o", label="data")

    # Plot Labels
    ax.set_zlabel("Pressure", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.set_xlabel("Temperature (K)", fontsize=12)
    ax.legend(loc="upper right", fontsize=15)
    plt.show()


def plotCV():
    """
    Plot isochoric heat capacity
    """
    DataImport.CV(molecule)
    CVValues = DataImport.CVValues
    DT = []
    for x in CVValues:
        Tau = float(critT) / float(x[1])
        Delta = float(x[0]) / float(critD)
        CVValuesn.append((Delta, Tau, float(x[2])))
        DT.append((Delta, Tau))
    CVnp = np.asarray(CVValuesn)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(
        CVnp[:, 1], CVnp[:, 0], CVnp[:, 2], s=20, c="b", marker="o", label="NIST data"
    )
    fig.suptitle("Isochoric Heat Capacity Data", fontsize=20, fontweight="bold")

    ax.legend(loc="upper right", fontsize=15)
    ax.set_zlabel("Isochoric Heat Capacity", fontsize=12)
    ax.set_ylabel("Reduced Density", fontsize=12)
    ax.set_xlabel("Inverse Reduced Temperature", fontsize=12)
    plt.show()


def plotCP():
    """
    Plot isobaric heat capacity
    """
    DataImport.CP(molecule)
    CPValues = DataImport.CPValues
    CPValuesn = []
    DT = []
    for x in CPValues:
        Theta = float(critT) / float(x[1])
        De = float(x[0]) / float(critD)
        z = float(x[2])
        CPValuesn.append((De, Theta, z, float(x[0])))
        DT.append((De, Theta))

    CPnp = np.asarray(CPValuesn)
    T = CPnp[:, 1]
    D = CPnp[:, 0]
    CP = CPnp[:, 2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(T, D, CP, s=20, c=T, marker="o", label="NIST data")
    fig.suptitle("Isobaric Heat Capacity Data", fontsize=20, fontweight="bold")
    ax.legend(loc="upper right")
    ax.set_zlabel("$C_p$", fontsize=12)
    ax.set_ylabel(r"$\delta$", fontsize=12)
    ax.set_xlabel("$ \\tau $", fontsize=12)

    plt.show()


def plotSND():
    """
    Plot speed of sound data
    """
    DataImport.SND(molecule)
    Values = DataImport.SNDValues
    WS = []

    for x in Values:
        Theta = float(critT) / float(x[1])
        Delta = float(x[0]) / float(critD)
        z = float(x[2])
        WS.append([Theta, Delta, z])

    SNDnp = np.asarray(WS)
    T = SNDnp[:, 0]
    D = SNDnp[:, 1]
    w = SNDnp[:, 2]

    w = (w**2) / R / 1000 * M / T * float(critT)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(T, D, w, s=25, c="r", marker="o", label="NIST data")
    ax.legend(loc="upper right", fontsize=15)
    fig.suptitle("Speed of Sound Data", fontsize=20, fontweight="bold")
    ax.set_zlabel("Dimensionless Speed of Sound", fontsize=12)
    ax.set_ylabel(r"$\delta$", fontsize=15)
    ax.set_xlabel("$\\tau$", fontsize=15)

    plt.show()


# View Regressed Equation
def sseCombo(lstFile=None, plot=False, report=False, surface=cm.coolwarm):
    """Plot regressed equation and data. Calculates statistical anlaysis metrics"""
    Y = []
    Beta = []

    parseGAMS.parser(lstFile)
    Y = parseGAMS.indexes[0]
    Beta = [float(b) for b in parseGAMS.betas[0]]

    BasisFunctions.formCustomBasis()
    parseGAMS.writeEquation(Y, Beta)

    saveFig = False
    showGraph = plot

    objZ, objCV, objCP, objSND = 0, 0, 0, 0

    if "PVT" in props:
        objZ = ssePVT(Y, Beta, saveFig, showGraph, report)
    if "CV" in props:
        objCV = sseCV(Y, Beta, saveFig, showGraph, report)
    if "CP" in props:
        objCP = sseCP(Y, Beta, saveFig, showGraph, report)
    if "SND" in props:
        objSND = sseSND(Y, Beta, saveFig, showGraph, report)

    if showGraph:
        HelmetSurface(Y, Beta, showGraph, surface)


def ssePVT(PVT1=[], PVT1Vals=[], saveFig=False, show=True, report=False):
    """Plots and metrics for Pressure-volume-temperature data"""
    Y, Beta = PVT1, PVT1Vals
    BasisFunctions.formCustomBasis()

    DT = []
    PlaceHolder = []
    Z = []

    DataImport.PVT(molecule)
    PVTValues = DataImport.PVTValues

    DataToWrite = []
    for x in PVTValues:
        manVal = DataManipulation.PVT(x)
        DataToWrite.append(manVal)

    for x in DataToWrite:
        Pval = BasisFunctions.drdRes(x[0], x[1], Y, Beta)
        PlaceHolder.append(Pval)
        Z.append(x[2])

    PVTnp = np.asarray(PVTValues)
    P = []
    DT = []
    for x in PVTValues:
        P.append(float(x[0]))
        T = float(critT) / float(x[2])
        D = float(x[1]) / float(critD)
        DT.append((D, T))

    anp = np.asarray(DT)
    D = anp[:, 0]
    T = anp[:, 1]
    Zvar = np.var(Z)
    Ze = [(y - x) for x, y in zip(Z, PlaceHolder)]
    Z = [(y - x) ** 2 for x, y in zip(Z, PlaceHolder)]
    PlaceHolder = [1 + x for x in PlaceHolder]
    PlaceHolder = [
        y * R / 1000 * float(x[2]) * float(x[1]) for x, y in zip(PVTValues, PlaceHolder)
    ]

    Ta = PVTnp[:, 2]
    Pa = np.array(P)

    delPcalc = (Pa - PlaceHolder) / Pa
    AAD = np.mean(delPcalc)

    critVal = BasisFunctions.drdRes(1, 1, Y, Beta)
    d2rdv = BasisFunctions.d2rdRes(1, 1, Y, Beta)
    d3rdv = BasisFunctions.d3rdRes(1, 1, Y, Beta)
    if report:
        print("Critical Values", critVal)
    CritPressureCalc = (critVal + 1) * R / 1000 * float(critT) * float(critD)
    if report:
        print("Critical Pressures", float(critP), CritPressureCalc)

    CritFirstDeriv = 1 + 2 * critVal + d2rdv
    CritSecondDeriv = 2 * critVal + 4 * d2rdv + d3rdv

    if report:
        print("First Deriv", 0, CritFirstDeriv)
        print("Second Deriv", 0, CritSecondDeriv)
    errors = [(m - x) ** 2 for m, x in zip(Pa, PlaceHolder)]

    sseCalc = sum(errors)
    avgY = sum(Pa) / len(Pa)
    sseTotal = sum((Pa - avgY) ** 2)
    R2calc = 1 - sseCalc / sseTotal
    meanRes = np.mean(Pa - PlaceHolder)
    residuals = [m - x for m, x in zip(Pa, PlaceHolder)]

    S2 = sum((residuals - meanRes) ** 2) / (len(residuals) - 1)
    S20 = sum(np.array(residuals) ** 2) / (len(residuals) - 1)
    tscore = meanRes / np.sqrt(S2 + S20)

    if report:
        print("Pressure ----------------------")
        print("SSE for Fit %f" % sseCalc)
        print("Average residual %f" % np.mean(Pa - PlaceHolder))
        print("R2 calculated is %f" % R2calc)
        RMSE = np.sqrt(np.mean(errors))
        print("RMSE calculated is %f" % RMSE)
        print("AAD is %f, %f, %f" % (min(delPcalc), AAD, max(delPcalc)))
        print("t value %f, df %f" % (tscore, len(residuals)))
        print(stats.ttest_1samp(residuals, 0))

    if show or saveFig:
        fig2 = plt.figure()
        fig2.suptitle("%s Pressure Volume Temperature Regression" % molecule)

        # Pressure Parity Plot
        ax2 = fig2.add_subplot(131)
        # ax2.set_title("%s PVT Regression" % molecule)
        ax2.set_xlabel("Observed pressure (MPa)")
        ax2.set_ylabel("Calculated pressure (MPa)")

        ax2.plot(
            [0, max(max(Pa), max(PlaceHolder))],
            [0, max(max(Pa), max(PlaceHolder))],
            c="b",
            label="$R^2 = %.3f$" % R2calc,
        )
        ax2.set_xlim([0, 120])
        ax2.set_ylim([0, 120])
        ax2.scatter(
            Pa,
            PlaceHolder,
            c="w",
            s=10,
            edgecolor="k",
            label="Developed equation\n (multivariable)",
        )
        ax2.legend(loc="upper left")

        # Pressure Data Plot
        ax = fig2.add_subplot(132)
        ax.scatter(
            Ta, Pa, s=25, c="r", marker="o", edgecolor="k", label="%s data" % molecule
        )
        ax.scatter(
            Ta,
            PlaceHolder,
            s=20,
            c="b",
            marker="o",
            edgecolor="k",
            label="Developed equation",
        )
        # ax.set_title('Data Fit');
        # ax1.set_title('Residuals')
        ax.set_ylabel("Pressure (MPa)")
        ax.set_xlabel("Temperature (K)")
        ax.legend(loc="upper right")

        #  Residuals Plot
        ax4 = fig2.add_subplot(133)
        ax4.scatter(T, Pa - PlaceHolder, c=residuals, edgecolor="k")
        ax4.set_xlabel("Temperature (K)")
        ax4.set_ylabel("Residuals")

        # 3d Plot
        fig3 = plt.figure()
        fig3.suptitle("%s Pressure Volume Temperature Regression" % molecule)
        ax3 = fig3.add_subplot(111, projection="3d")
        ax3.scatter(D, T, Pa, edgecolor="k", label="Data")
        ax3.scatter(D, T, PlaceHolder, edgecolor="k", label="Regression")
        ax3.set_xlabel("Reduced Density")
        ax3.set_ylabel("Inverse Reduced Temperature")
        ax3.set_zlabel("Pressure (MPa)")

        # figAAD = plt.figure()
        # axAAD = figAAD.add_subplot(111)
        # PVTAAD = [ (x-y)/x for x, y in zip(Pa,PlaceHolder)]
        # SAAD = [ x[0] for x in zip(Sources,PVTAAD) if abs(x[1]) >1]
        # TAAD = [ x[0] for x in zip(T,PVTAAD) if abs(x[1]) >1]
        # PVTAAD = [ x[1] for x in zip(T,PVTAAD) if abs(x[1]) >1]

        # df = pd.DataFrame(dict(x=TAAD, y=PVTAAD, label=SAAD))
        # groups = df.groupby('label')

        # i = 0
        # for name, group in groups:
        #     axAAD.scatter(group.x, group.y, label=name, marker=markers[i])
        #     i+=1
        # # axAAD.scatter(TAAD, PVTAAD, label=SAAD)
        # axAAD.legend(loc="right")

        if saveFig:
            fig2.savefig("%sPics/%sPVT.eps" % (molecule, molecule))
            fig3.savefig("%sPics/%sPVTR2.eps" % (molecule, molecule))
        if show:
            plt.show()

    if report:
        print("VARZ: ", Zvar)
    # l-1 norm
    sseZ1 = [abs(x) for x in Ze]
    # l-2 norm
    sseZ2 = [x**2 for x in Ze]
    # print "SSEZ2: ", np.sum(sseZ1)
    return [np.sum(sseZ1) / Zvar, np.sum(sseZ2) / Zvar, np.sum(sseZ1), np.sum(sseZ2)]


def sseCV(Y=[], Beta=[], saveFig=False, show=True, report=False):
    """Plots and metrics for isochoric heat capacity"""
    DT = []
    CVValuesn = []
    PlaceHolder = []

    global DataToWrite
    CVValues = []
    try:
        DataImport.CV(molecule)
        CVValues = DataImport.CVValues
    except Exception:
        print("Error plotting CV")
        return
    else:
        pass

    CVValues = DataImport.CVValues
    DataToWrite = []
    for x in CVValues:
        manVal = DataManipulation.CV(x)
        DataToWrite.append(manVal)
    BasisFunctions.formCustomBasis()

    for x in DataToWrite:
        val = BasisFunctions.rTTRes(x[0], x[1], Y, Beta)
        itt_val = BasisFunctions.iTT(x[0], x[1])
        CVR = val + itt_val
        PlaceHolder.append(CVR)

    for x in CVValues:
        Tau = float(critT) / float(x[1])
        Delta = float(x[0]) / float(critD)
        CVValuesn.append((Delta, Tau, float(x[2])))
        DT.append((Delta, Tau))

    PlaceHolder = [-z * R for T, z in zip(DT, PlaceHolder)]
    CVnp = np.asarray(CVValuesn)
    T = CVnp[:, 1]
    CV = CVnp[:, 2]

    if show or saveFig:
        fig = plt.figure()
        fig.suptitle("%s Isochoric Heat Capacity Regression" % molecule)
        ax = fig.add_subplot(121)
        ax.scatter(
            T, CV, s=10, c="b", edgecolors="k", marker="s", label="%s data" % molecule
        )
        ax.scatter(
            T,
            PlaceHolder,
            s=10,
            c="r",
            marker="o",
            edgecolors="k",
            label="Developed equation",
        )
        ax.legend(loc="upper left")
        ax.set_ylabel("Isochoric Heat Capacity")
        ax.set_xlabel("Temperature (K)")

        delCVcalc = (CV - PlaceHolder) / (CV)
        AAD = np.mean(delCVcalc)
        sseCalc = sum((CV - PlaceHolder) ** 2)
        avgY = sum(CV) / len(CV)
        sseTotal = sum((CV - avgY) ** 2)
        R2calc = 1 - sseCalc / sseTotal
        errors = (CV - PlaceHolder) ** 2

        meanRes = np.mean(CV - PlaceHolder)
        residuals = [m - x for m, x in zip(CV, PlaceHolder)]
        S2 = sum((residuals - meanRes) ** 2) / (len(residuals) - 1)
        S20 = sum(np.array(residuals) ** 2) / (len(residuals) - 1)
        t = meanRes / np.sqrt(S2 + S20)

        if report:
            print("Isochoric Heat Capacity --------------------------")
            print("SSE for Fit %f" % sseCalc)
            print("Average residual %f" % np.mean(CV - PlaceHolder))
            RMSE = np.sqrt(np.mean(errors))
            print("RMSE calculated is %f" % RMSE)
            # print "SSE total %f" % sseTotal;
            print("R2 calculated is %f" % R2calc)
            print("AAD is %f, %f, %f" % (min(delCVcalc), AAD, max(delCVcalc)))
            print("t value %f, df %i" % (t, len(residuals)))
            print(stats.ttest_1samp(residuals, 0))

        # Parity Plot
        # fig2 = plt.figure()
        ax2 = fig.add_subplot(122)
        ax2.set_xlabel(r"$C_V/(J/(mol\cdot K))$ observed")
        ax2.set_ylabel(r"$C_V/(J/(mol\cdot K))$ calculated")

        maxCV = math.ceil(max(max(CV), max(PlaceHolder)))
        ax2.plot([0, maxCV], [0, maxCV], label="$R^2 =  %.3f$" % R2calc, color="b")
        # ax2.plot([0,550], [0,550], label='$R^2 =  %.3f$' %R2calc, color='b')
        ax2.scatter(
            CV,
            PlaceHolder,
            edgecolors="k",
            facecolors="none",
            label="Developed equation",
        )
        # ax2.set_xlim([0,550])
        # ax2.set_ylim([0,550])
        ax2.legend(loc="upper left")

        # fig3 = plt.figure()
        # ax3 = fig3.add_subplot(111, projection='3d');
        # ax3.scatter(D, T, CV, edgecolor = 'k', label ='Data')
        # ax3.scatter(D,T, PlaceHolder, edgecolor='k', label='Regression')
        # ax3.set_xlabel('Delta');
        # ax3.set_ylabel('Tau');
        # ax3.set_zlabel('Pressure')

        if saveFig:
            fig.savefig("%sPics/%sCV.eps" % (molecule, molecule))
            # fig2.savefig('%sPics/%sCVR2.eps' %(molecule,molecule))
        if show:
            plt.show()

    CVR = [x / R for x in CV]
    PR = [x / R for x in PlaceHolder]
    # l-1 norm
    sseCV1 = [abs(x - y) for x, y in zip(CVR, PR)]
    # l-2 norm
    sseCV2 = [(x - y) ** 2 for x, y in zip(CVR, PR)]
    if report:
        print("VARCV: ", np.var(CVR))
    return [sum(sseCV1) / np.var(CVR), sum(sseCV2) / np.var(CVR)]


def sseCP(CP1=[], CP1Vals=[], saveFig=False, show=True, report=False):
    """Plots and calculates metrics for isobaric heat capacity"""
    global molecule
    Y = CP1
    Beta = CP1Vals

    BasisFunctions.formCustomBasis()
    DataImport.CP(molecule)

    CPValues = DataImport.CPValues
    CPValuesn = []
    DT = []

    rD = []
    rDD = []
    rTD = []
    PlaceHolder = []

    for x in CPValues:
        Theta = float(critT) / float(x[1])
        De = float(x[0]) / float(critD)
        z = float(x[2])
        CPValuesn.append((De, Theta, z, float(x[0])))
        DT.append((De, Theta))

    for x in DT:
        val = BasisFunctions.rTTRes(x[0], x[1], Y, Beta)
        itt_val = BasisFunctions.iTT(x[0], x[1])
        CPR = val + itt_val
        PlaceHolder.append(CPR)
        val = BasisFunctions.drdRes(x[0], x[1], Y, Beta)
        rD.append(val)
        rDD.append(BasisFunctions.d2rdRes(x[0], x[1], Y, Beta))
        rTD.append(BasisFunctions.dtrdtRes(x[0], x[1], Y, Beta))

    CPnp = np.asarray(CPValuesn)
    T = CPnp[:, 1]
    D = CPnp[:, 0]
    CP = CPnp[:, 2]

    dPD = [1 + 2 * y + z for x, y, z in zip(DT, rD, rDD)]  # 10 to 10^-2
    dPT = [1 + y - z for x, y, z in zip(DT, rD, rTD)]  # 10^-1

    CPcalc = [-x + (y**2) / z for x, y, z in zip(PlaceHolder, dPT, dPD)]
    CPcalc = [x * R for x in CPcalc]

    errors = (CP - CPcalc) ** 2
    delCPcalc = (CP - CPcalc) / CP
    AAD = np.mean(delCPcalc)
    sseCalc = sum((CP - CPcalc) ** 2)
    avgY = sum(CP) / len(CP)
    sseTotal = sum((CP - avgY) ** 2)
    R2calc = 1 - sseCalc / sseTotal

    meanRes = np.mean(CP - CPcalc)
    residuals = [m - x for m, x in zip(CP, CPcalc)]
    S2 = sum((residuals - meanRes) ** 2) / (len(residuals) - 1)
    S20 = sum(np.array(residuals) ** 2) / (len(residuals) - 1)
    t = meanRes / np.sqrt(S2 + S20)

    if report:
        print("Isobaric Heat Capacity ---------------------------")
        print("SSE for Fit %f" % sseCalc)
        print("Average residual %f" % np.mean(CP - CPcalc))
        RMSE = np.sqrt(np.mean(errors))
        print("RMSE calculated is %f" % RMSE)
        # print "SSE total %f" % sseTotal;
        print("R2 calculated is %f" % R2calc)
        print("AAD is %f, %f, %f" % (min(delCPcalc), AAD, max(delCPcalc)))
        print("t value %f, df %i" % (t, len(residuals)))
        print(stats.ttest_1samp(residuals, 0))

    if show or saveFig:

        fig2 = plt.figure()
        fig2.suptitle("%s Isobaric Heat Capacity Regression" % molecule)

        # Parity Plot
        ax2 = fig2.add_subplot(131)
        ax2.set_xlabel(r"$C_P/(J/(mol\cdot K))$ observed")
        ax2.set_ylabel(r"$C_P/(J/(mol\cdot K))$ calculated")
        minCP = math.floor(min(min(CP), min(CPcalc)))
        maxCP = math.ceil(max(max(CP), max(CPcalc)))
        ax2.plot(
            [minCP, maxCP], [minCP, maxCP], label="$R^2 = %.3f$" % R2calc, color="b"
        )
        ax2.set_xlim([minCP, maxCP])
        ax2.set_ylim([minCP, maxCP])
        ax2.scatter(
            CP,
            CPcalc,
            edgecolors="k",
            facecolors="none",
            label="Developed equation\n (multivariable)",
        )
        ax2.legend(loc="upper right")

        # Scatter Plot
        ax = fig2.add_subplot(132)
        ax.scatter(
            T,
            CP,
            s=15,
            c="r",
            marker="s",
            edgecolors="k",
            label="%s data (NIST)" % molecule,
        )
        ax.scatter(
            T,
            CPcalc,
            s=15,
            c="b",
            marker="o",
            edgecolors="k",
            label="Developed equation",
        )
        ax.set_ylabel(r"$C_P/(J/(mol\cdot K))$")
        ax.set_xlabel("Inverse Reduced Temperature")
        ax.legend(loc="upper right")
        ax.set_xlim((0.8, 3.4))
        ax.set_ylim((0, 900))

        ax4 = fig2.add_subplot(133)

        ax4.scatter(
            D,
            CP,
            s=15,
            c="r",
            marker="s",
            edgecolors="k",
            label="%s data (NIST)" % molecule,
        )
        ax4.scatter(
            D,
            CPcalc,
            s=15,
            c="b",
            marker="o",
            edgecolors="k",
            label="Developed equation",
        )
        ax4.set_ylabel(r"$C_P/(J/(mol\cdot K))$")
        ax4.set_xlabel("Reduced Density")
        ax4.legend(loc="upper right")
        # ax4.set_title("%s Isobaric Heat Capacity Regression" % molecule)
        ax.legend(loc="upper right")

        fig3 = plt.figure()
        ax3 = fig3.add_subplot(111, projection="3d")
        ax3.scatter(D, T, CP, edgecolor="k", label="Data")
        ax3.scatter(D, T, CPcalc, edgecolor="k", label="Regression")
        ax3.set_xlabel("Reduced Density")
        ax3.set_ylabel("Inverse Reduced Temperature")
        ax3.set_zlabel("Pressure")
        ax3.set_title("%s Isobaric Heat Capacity Data" % molecule)
        ax3.legend()

        if saveFig:
            fig2.savefig("%sPics/%sCP.eps" % (molecule, molecule))
        if show:
            plt.show()

    CPR = [x / R for x in CP]
    # l -1 norm
    sseCP1 = [abs(x / R - y / R) for x, y in zip(CP, CPcalc)]
    # l-2  norm
    sseCP2 = [(x / R - y / R) ** 2 for x, y in zip(CP, CPcalc)]

    if report:
        print("VARCP: ", np.var(CPR))
    return [sum(sseCP1) / np.var(CPR), sum(sseCP2) / np.var(CPR)]


def sseSND(SND1=[], SND1Vals=[], saveFig=False, show=True, report=False):
    """Plots and calculates metrics for speed of sound"""
    Y = SND1
    Beta = SND1Vals

    BasisFunctions.formCustomBasis()
    DataImport.SND(molecule)
    Values = DataImport.SNDValues

    DT = []
    WS = []

    rD = []
    rDD = []
    rTD = []
    rTT = []
    PlaceHolder = []
    SNDValues = []

    for x in Values:
        Theta = float(critT) / float(x[1])
        D = float(x[0]) / float(critD)
        z = float(x[2])
        SNDValues.append((D, float(x[1]), z, float(x[0])))
        DT.append((D, Theta))
        WS.append(z)

    SNDnp = np.asarray(SNDValues)
    T = SNDnp[:, 1]
    D = SNDnp[:, 0]
    w = SNDnp[:, 2]

    for x in DT:
        val = BasisFunctions.rTTRes(x[0], x[1], Y, Beta)
        itt_val = BasisFunctions.iTT(x[0], x[1])
        # was not negative, but .iTT returns - itt
        CPR = val + itt_val
        PlaceHolder.append(CPR)
        val = BasisFunctions.drdRes(x[0], x[1], Y, Beta)
        rD.append(val)
        rDD.append(BasisFunctions.d2rdRes(x[0], x[1], Y, Beta))
        rTD.append(BasisFunctions.dtrdtRes(x[0], x[1], Y, Beta))
        rTT.append(BasisFunctions.rTTRes(x[0], x[1], Y, Beta))

    dPD = [1 + 2 * x + y for x, y in zip(rD, rDD)]  # 10 to 10^-2
    dPT = [1 + x - y for x, y in zip(rD, rTD)]  # 10^-1
    Dwcalc = [x - (y**2) / z for x, y, z in zip(dPD, dPT, PlaceHolder)]
    Dwcalc = [x for x in Dwcalc]

    w = (w**2) / R / 1000 * M / T * float(critT)
    WS = [(x**2) / R / 1000 * M * y[1] / float(critT) for x, y in zip(WS, DT)]
    We = WS
    Wes = Dwcalc

    WSn = np.asarray(WS)
    Dw = np.asarray(Dwcalc)

    delWcalc = [(d - w) / d for d, w in zip(WS, Dwcalc)]
    AAD = np.mean(delWcalc)
    sseCalc = sum((WSn - Dw) ** 2)
    avgY = sum(WSn) / len(WSn)
    sseTotal = sum((WSn - avgY) ** 2)
    R2calc = 1 - sseCalc / sseTotal
    errors = [(d - w) ** 2 for d, w in zip(WS, Dwcalc)]
    residual = [d - w for d, w in zip(WS, Dwcalc)]

    meanRes = np.mean(residual)
    residuals = [m - x for m, x in zip(WS, Dwcalc)]
    S2 = sum((residuals - meanRes) ** 2) / (len(residuals) - 1)
    S20 = sum(np.array(residuals) ** 2) / (len(residuals) - 1)
    t = meanRes / np.sqrt(S2 + S20)

    if report:
        print("Speed of Sound ----------------------------")
        print("SSE for Fit %f" % sseCalc)
        print("Average residual %f" % np.mean(residual))
        RMSE = np.sqrt(np.mean(errors))
        print("RMSE calculated is %f" % RMSE)
        # print "SSE total %f" % sseTotal;
        print("R2 calculated is %f" % R2calc)
        print("AAD is %f, %f, %f" % (min(delWcalc), AAD, max(delWcalc)))
        print("t value %f, df %i" % (t, len(residuals)))
        print(stats.ttest_1samp(residuals, 0))

    if show or saveFig:
        fig2 = plt.figure()
        fig2.suptitle("%s Speed of Sound Regression" % molecule)

        # Parity Plot
        ax2 = fig2.add_subplot(121)
        ax2.set_xlabel("Speed of sound observed")
        ax2.set_ylabel("Speed of sound calculated")
        minSND = math.floor(min(min(WS), min(Dwcalc)))
        maxSND = math.ceil(max(max(WS), max(Dwcalc)))
        ax2.plot(
            [minSND, maxSND], [minSND, maxSND], label="$R^2 = %.3f$" % R2calc, color="b"
        )
        ax2.set_xlim([minSND, maxSND])
        ax2.set_ylim([minSND, maxSND])
        ax2.scatter(
            WS,
            Dwcalc,
            facecolors="none",
            edgecolors="k",
            label="Developed equation\n (multivariable)",
        )
        ax2.legend(loc="upper right")

        ax = fig2.add_subplot(122)

        ax.scatter(
            T, WS, s=25, c="r", marker="s", edgecolors="k", label="%s data" % molecule
        )
        ax.scatter(
            T,
            Dwcalc,
            s=20,
            c="b",
            marker="o",
            edgecolors="k",
            label="Developed equation",
        )
        ax.set_xlabel("Inverse reduced temperature")
        ax.set_ylabel("Speed of sound")
        ax.legend(loc="upper right")
        ax.set_xlim(min(T) * 0.8, max(T) * 1.1)
        ax.set_ylim(min(WS) * 0.8, max(WS) * 1.1)

        fig3 = plt.figure()
        fig3.suptitle("%s Speed of Sound Regression" % molecule)
        ax3 = fig3.add_subplot(111, projection="3d")
        ax3.scatter(D, T, WS, edgecolor="k", label="Data")
        ax3.scatter(D, T, Dwcalc, edgecolor="k", label="Regression")
        ax3.set_xlabel("Reduced Density")
        ax3.set_ylabel("Inverse Reduced Temperature")
        ax3.set_zlabel("Speed of sound")

        if saveFig:
            fig2.savefig("%sPics/%sSND.eps" % (molecule, molecule))
            fig3.savefig("%sPics/%sSNDR2.eps" % (molecule, molecule))
        if show:
            plt.show()

    if report:
        print("VAR", np.var(We))
    # l-1 norm
    sseW1 = [abs(x - y) for x, y in zip(Wes, We)]
    # l-2 norm
    sseW2 = [(x - y) ** 2 for x, y in zip(Wes, We)]
    return [sum(sseW1) / np.var(We), sum(sseW2) / np.var(We)]


def HelmetSurface(Y=[], Beta=[], show=True, surface=cm.coolwarm):
    """Plots Helmholtz Surface"""
    global critT, critP, critD, M, triple, acc, R, Rm, molecule
    BasisFunctions.formCustomBasis()
    # print(Y, Beta)

    fig = plt.figure()
    fig.suptitle("%s Helmholtz Energy Surface" % molecule)
    ax = fig.gca(projection="3d")

    D = np.arange(0.1, 2, 0.005)
    T = np.arange(0.1, 2, 0.005)
    D, T = np.meshgrid(D, T)

    PlaceHolder = np.zeros((len(D), len(T)))

    for i in range(len(D)):
        for j in range(len(D)):
            # print i,j, D[i,j], T[i,j]
            val = BasisFunctions.arBY(D[i, j], T[i, j], Y, Beta)
            itt_val = BasisFunctions.idealBY(D[i, j], T[i, j], Y, Beta)
            Helm = val + itt_val
            PlaceHolder[i, j] = Helm

    ax.plot_surface(D, T, PlaceHolder, cmap=surface, linewidth=0, antialiased=False)

    ax.xaxis.labelpad = 5
    ax.yaxis.labelpad = 5
    ax.zaxis.labelpad = 5

    ax.set_xlabel("Reduced density")
    ax.set_ylabel("Inverse reduced temperature")
    ax.set_zlabel("Helmholtz energy")

    ax.set_xticks([0.25, 0.75, 1.25, 1.75])
    ax.set_yticks([0.25, 0.75, 1.25, 1.75])
    # ax.set_zticks([-2, 0, 2])

    plt.show()
