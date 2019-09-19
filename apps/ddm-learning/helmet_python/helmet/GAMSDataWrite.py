""" .. module:: GAMSDataWrite
       :platform: Unix

    .. moduleauthor:: Marissa Engle <mengle@andrew.cmu.edu>

"""

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# import BasisFunctions
# import DataImport
# import helmet
from helmet import BasisFunctions, DataImport


Combination = False
R = 8.314472


def Crit(textFile, DataToWrite, Combination=False):
    BasisFunctions.drd(1, 1)
    drd_vals = BasisFunctions.drd_vals
    # print drd_vals;
    # textFile.write("crit0.. %f=E= %f*(1+ sum(j,beta(j)*crit('1',j,'drd')));\n" %(critP, val2))
    # textFile.write("crit1.. 0 =E= %f * (1+2* sum(j,beta(j)*crit('1',j,'drd')) + sum(j,beta(j)*crit('1',j,'d2rd')));\n" %val)
    # textFile.write("crit2.. 0  =E= %f * (2*sum(j,beta(j)*crit('1',j,'drd')) + 4*sum(j,beta(j)*crit('1',j,'d2rd')) + sum(j,beta(j)*crit('1',j,'d3rd')));\n" %val1)

    i = 1
    j = 1
    for y in drd_vals:
        if Combination:
            textFile.write(
                "crit('%s', '%d', '%d', '%s') = %f ;\n" % ("PVT", i, j, "drd", y)
            )
        else:
            textFile.write("crit('%d', '%d', '%s') = %f ;\n" % (i, j, "drd", y))
        j += 1
    i = 1
    j = 1
    BasisFunctions.d2rd(1, 1)
    d2rd_vals = BasisFunctions.d2rd_vals
    for y in d2rd_vals:
        if Combination:
            textFile.write(
                "crit('%s', '%d', '%d', '%s') = %f ;\n" % ("PVT", i, j, "d2rd", y)
            )
        else:
            textFile.write("crit('%d', '%d', '%s') = %f ;\n" % (i, j, "d2rd", y))
        j += 1
    i = 1
    j = 1
    BasisFunctions.d3rd(1, 1)
    d3rd_vals = BasisFunctions.d3rd_vals
    for y in d3rd_vals:
        if Combination:
            textFile.write(
                "crit('%s', '%d', '%d', '%s') = %f ;\n" % ("PVT", i, j, "d3rd", y)
            )
        else:
            textFile.write("crit('%d', '%d', '%s') = %f ;\n" % (i, j, "d3rd", y))
        j += 1
    i = 1


def InSat(textFile, DataToWrite, Combination=False):

    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.d3rd(x[0], x[1])
        BasisFunctions.d4rd(x[0], x[1])
        BasisFunctions.d5rd(x[0], x[1])
        d3rd_vals = BasisFunctions.d3rd_vals
        d4rd_vals = BasisFunctions.d4rd_vals
        d5rd_vals = BasisFunctions.d5rd_vals
        for x, y, z in zip(d3rd_vals, d4rd_vals, d5rd_vals):
            val = 12 * x + 8 * y + z
            if Combination:
                textFile.write(
                    "InSat('%s', '%d', '%d', '%s') = %f ;\n"
                    % ("PVT", i, j, "d3rd", val)
                )
            else:
                textFile.write("InSat('%d', '%d', '%s') = %f ;\n" % (i, j, "d3rd", val))
            j += 1
        i += 1


def Tangent(textFile, DataToWrite, Combination=False):
    # BasisFunctions.ideal(deltas,taus);
    # BasisFunctions.ar(deltas,taus);
    # BasisFunctions.drd(deltas,taus);
    # print DataToWrite
    # for x in DataToWrite:
    #   print x[0],x[1]
    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.drd(x[0], x[1])
        drd_vals = BasisFunctions.drd_vals
        if i % 10 > 5:
            for y in drd_vals:
                if Combination:
                    textFile.write(
                        "tan('%s', '%d', '%d', '%s') = %f ;\n" % ("PVT", i, j, "drd", y)
                    )
                else:
                    textFile.write("tan('%d', '%d', '%s') = %f ;\n" % (i, j, "drd", y))
                j += 1
        i += 1
    # i = 1
    # for x in DataToWrite:
    #   j=1
    #   BasisFunctions.ar(x[0],x[1]);
    #   ar_vals = BasisFunctions.ar_vals;
    #   for y in ar_vals:
    #       val = y/(x[0]**3)
    #       if Combination:
    #           textFile.write("tan('%s', '%d', '%d', '%s') = %f ;\n" %('PVT', i, j, 'ar', val))
    #       else:
    #           textFile.write("tan('%d', '%d', '%s') = %f ;\n" %(i, j, 'ar', val))
    #       j+=1
    # for x in DataToWrite:
    #   j=1
    #   BasisFunctions.ideal(x[0],x[1]);
    #   ideal_vals = BasisFunctions.ideal_vals;
    #   for y in d3rd_vals:
    #       val = y/(x[0]**3)
    #       if Combination:
    #           textFile.write("isoT('%s', '%d', '%d', '%s') = %f ;\n" %('PVT', i, j, 'd3rd', val))
    #       else:
    #           textFile.write("isoT('%d', '%d', '%s') = %f ;\n" %(i, j, 'd3rd', val))
    #       j+=1


# @staticmethod # trying to speed up
def PVT2(textFile, DataToWrite, Combination=False):
    """ Imports P-V-T data into the GAMS document, to be written in :func:'GamsWrite'

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """
    print("PVT: ")
    # PVTindexes = DataImport.PVTindexes
    IsothermIndexes = DataImport.isothermIndex
    # print Combination
    i = 1
    for x in DataToWrite:
        j = 1
        if Combination:
            textFile.write("zPVT('%d') = %f ;\n" % (i, x[2]))
        else:
            textFile.write("zPVT('%d') = %f ;\n" % (i, x[2]))
        i += 1
    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.drd(x[0], x[1])
        drd_vals = BasisFunctions.drd_vals
        for y in drd_vals:
            if Combination:
                textFile.write("xPVT('%d', '%d') = %f ;\n" % (i, j, y))
            else:
                textFile.write("xPVT('%d', '%d') = %f ;\n" % (i, j, y))
            j += 1
        i += 1
    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.d3rd(x[0], x[1])
        d3rd_vals = BasisFunctions.d3rd_vals
        for y in d3rd_vals:
            val = y / (x[0] ** 3)
            if i <= IsothermIndexes:
                if Combination:
                    textFile.write("isoT('%d', '%d') = %f ;\n" % (i, j, val))
                else:
                    textFile.write("isoT('%d', '%d') = %f ;\n" % (i, j, val))
            j += 1
        i += 1
    Crit(textFile, DataToWrite, Combination)


# @staticmethod
def PVT(textFile, DataToWrite, Combination=False):
    """ Imports P-V-T data into the GAMS document, to be written in :func:'GamsWrite'

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """
    print("PVT: ")
    # PVTindexes = DataImport.PVTindexes
    IsothermIndexes = DataImport.isothermIndex
    # print Combination
    i = 1
    for x in DataToWrite:
        j = 1
        if Combination:
            textFile.write("z('%s', '%d') = %f ;\n" % ("PVT", i, x[2]))
        else:
            textFile.write("z('%d') = %f ;\n" % (i, x[2]))
        i += 1
    i = 1
    drd_max = 0
    drd_min = 100
    for x in DataToWrite:
        j = 1
        BasisFunctions.drd(x[0], x[1])
        drd_vals = BasisFunctions.drd_vals

        if np.amax(drd_vals) > drd_max:
            drd_max = 5 * np.amax(drd_vals)
        if np.amin(drd_vals) < drd_min:
            drd_min = 5 * np.amin(drd_vals)

        for y in drd_vals:
            # if y < 100 and y > 10:
            #   print x[0], x[1], j, y
            if Combination:
                textFile.write(
                    "xijk('%s', '%d', '%d', '%s') = %f ;\n" % ("PVT", i, j, "drd", y)
                )
            else:
                textFile.write("xijk('%d', '%d', '%s') = %f ;\n" % (i, j, "drd", y))
            j += 1
        i += 1
    i = 1

    # print " PVT DRD Max : ", drd_max;
    # print " PVT DRD Min : ", drd_min;

    for x in DataToWrite:
        j = 1
        BasisFunctions.d3rd(x[0], x[1])
        # BasisFunctions.d4rd(x[0],x[1]);
        # BasisFunctions.d5rd(x[0],x[1]);
        d3rd_vals = BasisFunctions.d3rd_vals
        # d4rd_vals = BasisFunctions.d4rd_vals;
        for y in d3rd_vals:
            val = y / (x[0] ** 3)
            if i <= IsothermIndexes:  # if i <= PVTindexes[0]:
                if Combination:
                    textFile.write(
                        "isoT('%s', '%d', '%d', '%s') = %f ;\n"
                        % ("PVT", i, j, "d3rd", val)
                    )
                else:
                    textFile.write(
                        "isoT('%d', '%d', '%s') = %f ;\n" % (i, j, "d3rd", val)
                    )
            j += 1
        i += 1
    Crit(textFile, DataToWrite, Combination)
    InSat(textFile, DataImport.InSatValues, Combination)

    x, y, z = [], [], []
    for xs in DataToWrite:
        x.append(xs[0])
        y.append(xs[1])
        z.append(xs[2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, edgecolor="k", c="none", marker="o")
    ax.set_xlabel("delta")
    ax.set_ylabel("tau")
    ax.set_zlabel("z")
    plt.show()


# @staticmethod
def newPVT(textFile, DataToWrite, Combination=False):
    """ Imports P-V-T data into the GAMS document, to be written in :func:'GamsWrite'

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """
    print("PVT: ")
    PVTindexes = DataImport.PVTindexes
    # print Combination
    i = 1
    for x in DataToWrite:
        j = 1
        if Combination:
            textFile.write("z('%s', '%d') = %f ;\n" % ("PVT", i, x[2]))
        else:
            textFile.write("z('%d') = %f ;\n" % (i, x[2]))
        i += 1
    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.drd(x[0], x[1])
        drd_vals = BasisFunctions.drd_vals
        for y in drd_vals:
            if Combination:
                textFile.write(
                    "xijk('%s', '%d', '%d', '%s') = %f ;\n" % ("PVT", i, j, "drd", y)
                )
            else:
                textFile.write("xijk('%d', '%d', '%s') = %f ;\n" % (i, j, "drd", y))
            j += 1
        i += 1
    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.d3rd(x[0], x[1])
        d3rd_vals = BasisFunctions.d3rd_vals
        for y in d3rd_vals:
            val = y / (x[0] ** 3)
            if i <= PVTindexes[0]:
                if Combination:
                    textFile.write(
                        "isoT('%s', '%d', '%d', '%s') = %f ;\n"
                        % ("PVT", i, j, "d3rd", val)
                    )
                else:
                    textFile.write(
                        "isoT('%d', '%d', '%s') = %f ;\n" % (i, j, "d3rd", val)
                    )
            j += 1
        i += 1
    Crit(textFile, DataToWrite, Combination)


# @staticmethod
def CV(
    textFile, DataToWrite, Combination=False
):  # ERROR: including itt_val in every term
    """ Imports Isochoric Heat Capacity (CV) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    .. todo:: Check itt_val in every term addition.
    """

    i = 1
    maxCVI = 0
    minCVI = 23
    for x in DataToWrite:
        BasisFunctions.iTT(x[0], x[1])
        itt_val = -BasisFunctions.itt_val
        if itt_val > maxCVI:
            maxCVI = itt_val
        if itt_val < minCVI:
            minCVI = itt_val

        val = x[2] - itt_val
        # print itt_val, x[2], val
        if Combination:
            # print Combination
            textFile.write("z('%s', '%d') = %f ;\n" % ("CV", i, val))
            # print val
        else:
            textFile.write("z('%d') = %f ;\n" % (i, val))
        i += 1
    i = 1
    maxCVR = 0
    minCVR = 0
    for x in DataToWrite:
        j = 1
        BasisFunctions.rTT(x[0], x[1])
        rtt_vals = BasisFunctions.rtt_vals

        if max(rtt_vals) > maxCVR:
            maxCVR = max(rtt_vals)
        if min(rtt_vals) < minCVR:
            minCVR = min(rtt_vals)

        for y in rtt_vals:
            val = -(y)
            # print y, val
            if Combination:
                # textFile.write("xijk('%s', '%d', '%d', '%s') = %f ;\n" %('CV', i, j, 'itt', val))
                textFile.write(
                    "xijk('%s', '%d', '%d', '%s') = %f ;\n" % ("CV", i, j, "rtt", val)
                )
            else:
                # textFile.write("xijk('%d', '%d', '%s') = %f ;\n" %(i, j, 'itt', val))
                textFile.write("xijk('%d', '%d', '%s') = %f ;\n" % (i, j, "rtt", val))
            j += 1
        i += 1

    # print "Max Ideal CV ", maxCVI;
    # print "Min Ideal CV ", minCVI;
    # print "Max Res CV ", maxCVR;
    # print "Min Res CV ", minCVR;

    x, y, z = [], [], []
    for xs in DataToWrite:
        x.append(xs[0])
        y.append(xs[1])
        z.append(xs[2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, edgecolor="k", c="none", marker="o")
    ax.set_xlabel("delta")
    ax.set_ylabel("tau")
    ax.set_zlabel("CV")
    plt.show()


# @staticmethod
def CP(
    textFile, DataToWrite, Combination=False
):  # ERROR including itt_val in every term
    """ Imports Isobaric Heat Capacity(CP) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    .. todo:: Check formulation.
    """

    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.iTT(x[0], x[1])
        itt_val = -BasisFunctions.itt_val
        val = x[2] - itt_val
        # if i ==1:
        #   print x[2]
        #   print itt_val
        # if i ==1:
        #   print x[2]
        if Combination:
            textFile.write("z('%s', '%d') = %.13f ;\n" % ("CP", i, val))
            # print val
        else:
            textFile.write("z('%d') = %.13f ;\n" % (i, val))
        i += 1
    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.rTT(x[0], x[1])
        rtt_vals = BasisFunctions.rtt_vals
        BasisFunctions.drd(x[0], x[1])
        drd_vals = BasisFunctions.drd_vals
        BasisFunctions.dtrdt(x[0], x[1])
        dtrdt_vals = BasisFunctions.dtrdt_vals
        BasisFunctions.d2rd(x[0], x[1])
        d2rd_vals = BasisFunctions.d2rd_vals

        maxDenom = 0
        DenomCheck = [2 * x + y for x, y in zip(drd_vals, d2rd_vals)]
        upMax = [5 * x for x in DenomCheck]
        lowMax = [-5 * x for x in DenomCheck]
        denomMax = np.append(lowMax, upMax)
        denomMax = np.sort(denomMax)
        denomMaxRev = denomMax[::-1]
        if sum(denomMaxRev[:12] > maxDenom):
            maxDenom = sum(denomMaxRev[:12])

        maxNumer = 0
        NumerCheck = [x - y for x, y in zip(drd_vals, dtrdt_vals)]
        upMax = [5 * x for x in NumerCheck]
        lowMax = [-5 * x for x in NumerCheck]
        numerMax = np.append(lowMax, upMax)
        numerMax = np.sort(numerMax)
        numerMaxRev = numerMax[::-1]
        # print denomMaxRev;
        # print
        if sum(numerMaxRev[:12] > maxNumer):
            maxNumer = sum(numerMaxRev[:12])

        for rtt, drd, d2rd, dtrdt in zip(rtt_vals, drd_vals, d2rd_vals, dtrdt_vals):
            # val = -(rtt + itt_val)
            # print itt_val
            numerator = drd - dtrdt
            # print numerator
            denominator = 2 * drd + d2rd

            if Combination:
                # textFile.write("xijk('%s', '%d', '%d', '%s') = %f ;\n" %('CP', i, j, 'itt', itt_val))
                textFile.write(
                    "xijk('%s', '%d', '%d', '%s') = %.13f ;\n"
                    % ("CP", i, j, "rtt", -rtt)
                )
                textFile.write(
                    "xijk('%s','%d', '%d', '%s') = %.13f ;\n"
                    % ("CP", i, j, "num", numerator)
                )
                textFile.write(
                    "xijk('%s','%d', '%d', '%s') = %.13f ;\n"
                    % ("CP", i, j, "denom", denominator)
                )
            else:
                # textFile.write("xijk('%d', '%d', '%s') = %f ;\n" %(i, j, 'itt', itt_val))
                textFile.write(
                    "xijk('%d', '%d', '%s') = %.13f ;\n" % (i, j, "rtt", -rtt)
                )
                textFile.write(
                    "xijk('%d', '%d', '%s') = %.13f ;\n" % (i, j, "num", numerator)
                )
                textFile.write(
                    "xijk('%d', '%d', '%s') = %.13f ;\n" % (i, j, "denom", denominator)
                )
            j += 1
        i += 1
    # print "Maximum CP Denominator:", maxDenom
    # print "Maximum CP Numerator:", (1+maxNumer)**2

    x, y, z = [], [], []
    for xs in DataToWrite:
        x.append(xs[0])
        y.append(xs[1])
        z.append(xs[2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, edgecolor="k", c="none", marker="o")
    ax.set_xlabel("delta")
    ax.set_ylabel("tau")
    ax.set_zlabel("Isobaric heat capacity")
    ax.set_title("Molecule isobaric heat capacity data")
    plt.show()


# @staticmethod
def SND(textFile, DataToWrite, Combination=False):

    """ Imports Speed of Sound (SND) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    .. todo:: Check Density and formulation
    """
    # global DataToWrite
    i = 1
    for x in DataToWrite:
        j = 1
        if Combination:
            textFile.write("z('%s', '%d') = %.13f ;\n" % ("SND", i, x[2]))
        else:
            textFile.write("z('%d') = %.13f ;\n" % (i, x[2]))
        i += 1
    i = 1
    for x in DataToWrite:
        j = 1
        BasisFunctions.rTT(x[0], x[1])
        rtt_vals = BasisFunctions.rtt_vals
        BasisFunctions.iTT(x[0], x[1])
        itt_val = BasisFunctions.itt_val
        BasisFunctions.drd(x[0], x[1])
        drd_vals = BasisFunctions.drd_vals
        BasisFunctions.dtrdt(x[0], x[1])
        dtrdt_vals = BasisFunctions.dtrdt_vals
        BasisFunctions.d2rd(x[0], x[1])
        d2rd_vals = BasisFunctions.d2rd_vals

        maxDenom = 0
        DenomCheck = [2 * x + y for x, y in zip(drd_vals, d2rd_vals)]
        upMax = [5 * x for x in DenomCheck]
        lowMax = [-5 * x for x in DenomCheck]
        denomMax = np.append(lowMax, upMax)
        denomMax = np.sort(denomMax)
        denomMaxRev = denomMax[::-1]
        if sum(denomMaxRev[:12] > maxDenom):
            maxDenom = sum(denomMaxRev[:12])

        maxNumer = 0
        NumerCheck = [x - y for x, y in zip(drd_vals, dtrdt_vals)]
        upMax = [5 * x for x in NumerCheck]
        lowMax = [-5 * x for x in NumerCheck]
        numerMax = np.append(lowMax, upMax)
        numerMax = np.sort(numerMax)
        numerMaxRev = numerMax[::-1]
        # print denomMaxRev;
        # print
        if sum(numerMaxRev[:12] > maxNumer):
            maxNumer = sum(numerMaxRev[:12])

        for rtt, drd, d2rd, dtrdt in zip(rtt_vals, drd_vals, d2rd_vals, dtrdt_vals):
            numerator = drd - dtrdt
            denominator = 2 * drd + d2rd
            if Combination:
                if j == 1:
                    textFile.write(
                        "xijk('%s','%d', '%d', '%s') = %.13f ;\n"
                        % ("SND", i, j, "itt", itt_val)
                    )
                textFile.write(
                    "xijk('%s','%d', '%d', '%s') = %.13f ;\n"
                    % ("SND", i, j, "rtt", rtt)
                )
                textFile.write(
                    "xijk('%s','%d', '%d', '%s') = %.13f ;\n"
                    % ("SND", i, j, "num", numerator)
                )
                textFile.write(
                    "xijk('%s','%d', '%d', '%s') = %.13f ;\n"
                    % ("SND", i, j, "denom", denominator)
                )
            else:
                if j == 1:
                    textFile.write(
                        "xijk('%d', '%d', '%s') = %.13f ;\n" % (i, j, "itt", itt_val)
                    )
                textFile.write(
                    "xijk('%d', '%d', '%s') = %.13f ;\n" % (i, j, "rtt", rtt)
                )
                textFile.write(
                    "xijk('%d', '%d', '%s') = %.13f ;\n" % (i, j, "num", numerator)
                )
                textFile.write(
                    "xijk('%d', '%d', '%s') = %.13f ;\n" % (i, j, "denom", denominator)
                )
            j += 1
        i += 1
    # print "Maximum SND Denominator:", maxDenom
    # print "Maximum SND Numerator:", (1+maxNumer)**2

    x, y, z = [], [], []
    for xs in DataToWrite:
        x.append(xs[0])
        y.append(xs[1])
        z.append(xs[2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, edgecolor="k", c="none", marker="o")
    ax.set_xlabel("delta")
    ax.set_ylabel("tau")
    ax.set_zlabel("Speed of Sound")
    plt.show()


# @staticmethod
def PVTdt(textFile, DataToWrite, Combination=False, PlotData=False):
    """ Imports P-V-T data into the GAMS document, to be written in :func:'GamsWrite'

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """
    # print "PVT: "
    # PVTindexes = DataImport.PVTindexes
    IsothermIndexes = DataImport.isothermIndex
    # print Combination
    i = 1
    for x in DataToWrite:
        # j = 1
        if Combination:
            textFile.write("z('%s', '%d') = %f ;\n" % ("PVT", i, x[2]))
        else:
            textFile.write("z('%d') = %f ;\n" % (i, x[2]))
        i += 1
    i = 1
    # drd_max = 0
    # drd_min = 100
    for x in DataToWrite:
        if Combination:
            textFile.write("delta('%s','%d') = %f;\n" % ("PVT", i, x[0]))
            textFile.write("tau('%s','%d') = %f;\n" % ("PVT", i, x[1]))
        else:
            textFile.write("delta('%d') = %f;\n" % (i, x[0]))
            textFile.write("tau('%d') = %f;\n" % (i, x[1]))

        i += 1
    i = 1

    # print " PVT DRD Max : ", drd_max;
    # print " PVT DRD Min : ", drd_min;

    # for x in DataToWrite:
    #   j=1
    #   BasisFunctions.d3rd(x[0],x[1]);
    #   # BasisFunctions.d4rd(x[0],x[1]);
    #   # BasisFunctions.d5rd(x[0],x[1]);
    #   d3rd_vals = BasisFunctions.d3rd_vals;
    #   # d4rd_vals = BasisFunctions.d4rd_vals;
    #   for y in d3rd_vals:
    #       val = y/(x[0]**3)
    #       if i <= IsothermIndexes: #if i <= PVTindexes[0]:
    #           if Combination:
    #               textFile.write("isoT('%s', '%d', '%d', '%s') = %f ;\n" %('PVT', i, j, 'd3rd', val))
    #           else:
    #               textFile.write("isoT('%d', '%d', '%s') = %f ;\n" %(i, j, 'd3rd', val))
    #       j+=1
    #   i+=1;
    # Crit(textFile, DataToWrite, Combination)
    InSat(textFile, DataImport.InSatValues,Combination)

    if PlotData:
        x, y, z = [], [], []
        for xs in DataToWrite:
            x.append(xs[0])
            y.append(xs[1])
            z.append(xs[2])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x, y, z, edgecolor="k", c=y, marker="o", s=25)
        ax.xaxis.set_ticks(np.arange(min(x), max(x), (max(x) - min(x)) / 4))
        ax.yaxis.set_ticks(np.arange(min(y), max(y), (max(y) - min(y)) / 5))
        ax.zaxis.set_ticks(np.arange(min(z), max(z), (max(z) - min(z)) / 5))
        ax.set_xlabel("delta", labelpad=20, fontsize=20)
        ax.set_ylabel("tau", labelpad=20, fontsize=20)
        ax.set_zlabel("z", labelpad=10, fontsize=20)
        plt.show()


# @staticmethod
def CVdt(textFile, DataToWrite, Combination=False, PlotData=False):  # ERROR: including itt_val in every term
    """ Imports Isochoric Heat Capacity (CV) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    .. todo:: Check itt_val in every term addition.
    """

    i = 1
    # maxCVI = 0;
    # minCVI = 23;
    for x in DataToWrite:
        BasisFunctions.iTT(x[0], x[1])
        itt_val = -BasisFunctions.itt_val
        # if itt_val >maxCVI:
        #   maxCVI = itt_val;
        # if itt_val<minCVI:
        #   minCVI = itt_val;

        # val = x[2]; - itt_val;
        val = x[2]
        # print val
        # print itt_val, x[2], val
        if Combination:
            # print Combination
            textFile.write("z('%s', '%d') = %f ;\n" % ("CV", i, val))
            textFile.write("itt('%s', '%d') = %f ;\n" % ("CV", i, itt_val))
            textFile.write("delta('%s', '%d') = %.13f ;\n" % ("CV", i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % ("CV", i, x[1]))
            # print val
        else:
            textFile.write("z('%d') = %f ;\n" % (i, val))
        i += 1
    # i = 1
    # maxCVR = 0;
    # minCVR = 0;
    # for x in DataToWrite:
    #   j = 1
    #   BasisFunctions.rTT(x[0],x[1]);
    #   rtt_vals = BasisFunctions.rtt_vals;

    #   # if max(rtt_vals)>maxCVR:
    #   #   maxCVR = max(rtt_vals);
    #   # if min(rtt_vals)<minCVR:
    #   #   minCVR = min(rtt_vals);

    #   for y in rtt_vals:
    #       val = -(y)
    #       #print y, val
    #       if Combination:
    #           textFile.write("z('%s', '%d') = %f ;\n" %('CV', i, val))
    #           textFile.write("itt('%s', '%d') = %f ;\n" %('CV', i, itt_val))
    #       else:
    #           textFile.write("xijk('%d', '%d', '%s') = %f ;\n" %(i, j, 'itt', val))
    #           textFile.write("xijk('%d', '%d', '%s') = %f ;\n" %(i, j, 'rtt', val))
    #       j+=1
    #   i += 1

    # print "Max Ideal CV ", maxCVI;
    # print "Min Ideal CV ", minCVI;
    # print "Max Res CV ", maxCVR;
    # print "Min Res CV ", minCVR;

    if PlotData:
        x, y, z = [], [], []
        for xs in DataToWrite:
            x.append(xs[0])
            y.append(xs[1])
            z.append(xs[2])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x, y, z, edgecolor="k", c="none", marker="o")
        ax.set_xlabel("delta")
        ax.set_ylabel("tau")
        ax.set_zlabel("CV")
        plt.show()


# @staticmethod
def CPdt(textFile, DataToWrite, Combination=False, PlotData=False):  # ERROR including itt_val in every term
    """ Imports Isobaric Heat Capacity(CP) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    .. todo:: Check formulation.
    """

    i = 1
    for x in DataToWrite:
        # j = 1
        BasisFunctions.iTT(x[0], x[1])
        itt_val = -BasisFunctions.itt_val
        val = x[2] - itt_val
        # if i ==1:
        #   print x[2]
        #   print itt_val
        # if i ==1:
        #   print x[2]
        if Combination:
            textFile.write("z('%s', '%d') = %.13f ;\n" % ("CP", i, val))
            # print val
        else:
            textFile.write("z('%d') = %.13f ;\n" % (i, val))
        i += 1
    i = 1
    for x in DataToWrite:

        if Combination:
            textFile.write("delta('%s', '%d') = %.13f ;\n" % ("CP", i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % ("CP", i, x[1]))
        else:
            textFile.write("delta( '%d') = %.13f ;\n" % (i, x[0]))
            textFile.write("tau('%d') = %.13f ;\n" % (i, x[1]))

        i += 1
    # print "Maximum CP Denominator:", maxDenom
    # print "Maximum CP Numerator:", (1+maxNumer)**2

    if PlotData:
        x, y, z = [], [], []
        for xs in DataToWrite:
            x.append(xs[0])
            y.append(xs[1])
            z.append(xs[2])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x, y, z, edgecolor="k", c=z, marker="o", s=30)
        ax.xaxis.set_ticks(np.arange(min(x), max(x), (max(x) - min(x)) / 4))
        ax.yaxis.set_ticks(np.arange(min(y), max(y), (max(y) - min(y)) / 5))
        ax.zaxis.set_ticks(np.arange(min(z), max(z), (max(z) - min(z)) / 5))
        ax.set_xlabel("delta", labelpad=20, fontsize=20)
        ax.set_ylabel("tau", labelpad=20, fontsize=20)
        ax.set_zlabel("Isobaric heat capacity", labelpad=10, fontsize=20)

        # ax.set_title('Molecule isobaric heat capacity data')
        plt.show()


# @staticmethod
def CPdtMix(
    textFile, DataToWrite, Combination=False
):  # ERROR including itt_val in every term
    """ Imports Isobaric Heat Capacity(CP) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    .. todo:: Check formulation.
    """
    # delta,tau, x1, x2, Z, rTT, num, denom)
    i = 1
    for x in DataToWrite:
        if Combination:
            textFile.write("z('%s', '%d') = %.13f ;\n" % ("CP", i, x[4]))
            # print val
        else:
            textFile.write("z('%d') = %.13f ;\n" % (i, x[4]))
        i += 1
    i = 1
    for x in DataToWrite:

        if Combination:
            textFile.write("delta('%s', '%d') = %.13f ;\n" % ("CP", i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % ("CP", i, x[1]))
            textFile.write("x1('%s','%d') = %.13f ;\n" % ("CP", i, x[2]))
            textFile.write("x2('%s','%d') = %.13f ;\n" % ("CP", i, x[3]))
            textFile.write("rTT('%s','%d') = %.13f ;\n" % ("CP", i, x[5]))
            textFile.write("num('%s','%d') = %.13f ;\n" % ("CP", i, x[6]))
            textFile.write("denom('%s','%d') = %.13f ;\n" % ("CP", i, x[7]))
        else:
            textFile.write("delta( '%d') = %.13f ;\n" % (i, x[0]))
            textFile.write("tau('%d') = %.13f ;\n" % (i, x[1]))

        i += 1
    # print "Maximum CP Denominator:", maxDenom
    # print "Maximum CP Numerator:", (1+maxNumer)**2

    x, y, z = [], [], []
    for xs in DataToWrite:
        x.append(xs[0])
        y.append(xs[1])
        z.append(xs[4])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, edgecolor="k", c=z, marker="o", s=100)
    ax.xaxis.set_ticks(np.arange(min(x), max(x), (max(x) - min(x)) / 4))
    ax.yaxis.set_ticks(np.arange(min(y), max(y), (max(y) - min(y)) / 5))
    ax.zaxis.set_ticks(np.arange(min(z), max(z), (max(z) - min(z)) / 5))
    ax.set_xlabel("delta", labelpad=20, fontsize=20)
    ax.set_ylabel("tau", labelpad=20, fontsize=20)
    ax.set_zlabel("Isobaric heat capacity", labelpad=10, fontsize=20)
    # ax.set_title('Molecule isobaric heat capacity data')
    plt.show()


# @staticmethod
def SNDdt(textFile, DataToWrite, Combination=False, PlotData=False):

    """ Imports Speed of Sound (SND) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    .. todo:: Check Density and formulation
    """
    # global DataToWrite
    i = 1
    for x in DataToWrite:
        # j = 1
        if Combination:
            textFile.write("z('%s', '%d') = %.13f ;\n" % ("SND", i, x[2]))
        else:
            textFile.write("z('%d') = %.13f ;\n" % (i, x[2]))
        i += 1
    i = 1
    for x in DataToWrite:
        # j = 1
        BasisFunctions.iTT(x[0], x[1])
        itt_val = BasisFunctions.itt_val
        # print itt_val;#, BasisFunctions.iTT(x[0],x[1])# positive version
        BasisFunctions.drd(x[0], x[1])
        if Combination:
            textFile.write("itt('%s','%d') = %.13f ;\n" % ("SND", i, itt_val))
            textFile.write("delta('%s', '%d') = %.13f ;\n" % ("SND", i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % ("SND", i, x[1]))
        else:
            textFile.write("itt('%s','%d') = %.13f ;\n" % (i, itt_val))
            textFile.write("delta('%s', '%d') = %.13f ;\n" % (i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % (i, x[1]))

        # BasisFunctions.rTT(x[0],x[1]);
        # rtt_vals = BasisFunctions.rtt_vals;
        # BasisFunctions.iTT(x[0], x[1]);
        # itt_val = BasisFunctions.itt_val;
        # BasisFunctions.drd(x[0], x[1]);
        # drd_vals = BasisFunctions.drd_vals;
        # BasisFunctions.dtrdt(x[0], x[1]);
        # dtrdt_vals = BasisFunctions.dtrdt_vals;
        # BasisFunctions.d2rd(x[0], x[1]);
        # d2rd_vals = BasisFunctions.d2rd_vals;

        # maxDenom = 0;
        # DenomCheck = [2*x + y for x,y in zip(drd_vals, d2rd_vals)];
        # upMax = [5*x for x in DenomCheck];
        # lowMax = [-5*x for x in DenomCheck];
        # denomMax = np.append(lowMax, upMax)
        # denomMax = np.sort(denomMax);
        # denomMaxRev = denomMax[::-1]
        # if sum(denomMaxRev[:12] > maxDenom):
        #   maxDenom = sum(denomMaxRev[:12])

        # maxNumer = 0;
        # NumerCheck = [x - y for x,y in zip(drd_vals, dtrdt_vals)];
        # upMax = [5*x for x in NumerCheck];
        # lowMax = [-5*x for x in NumerCheck];
        # numerMax = np.append(lowMax, upMax)
        # numerMax = np.sort(numerMax);
        # numerMaxRev = numerMax[::-1]
        # #print denomMaxRev;
        # #print
        # if sum(numerMaxRev[:12] > maxNumer):
        #   maxNumer = sum(numerMaxRev[:12])

        # for rtt, drd, d2rd, dtrdt in zip(rtt_vals, drd_vals, d2rd_vals, dtrdt_vals):
        #   numerator = drd - dtrdt;
        #   denominator = 2*drd + d2rd;
        #   if Combination:
        #       if j ==1:
        #           textFile.write("xijk('%s','%d', '%d', '%s') = %.13f ;\n" %('SND', i, j, 'itt', itt_val))
        #       textFile.write("xijk('%s','%d', '%d', '%s') = %.13f ;\n" %('SND', i, j, 'rtt', rtt))
        #       textFile.write("xijk('%s','%d', '%d', '%s') = %.13f ;\n" %('SND', i, j, 'num', numerator))
        #       textFile.write("xijk('%s','%d', '%d', '%s') = %.13f ;\n" %('SND',i, j, 'denom', denominator))
        #   else:
        #       if j ==1:
        #           textFile.write("xijk('%d', '%d', '%s') = %.13f ;\n" %(i, j, 'itt', itt_val))
        #       textFile.write("xijk('%d', '%d', '%s') = %.13f ;\n" %(i, j, 'rtt', rtt))
        #       textFile.write("xijk('%d', '%d', '%s') = %.13f ;\n" %(i, j, 'num', numerator))
        #       textFile.write("xijk('%d', '%d', '%s') = %.13f ;\n" %(i, j, 'denom', denominator))
        #   j+=1
        i += 1
    # print "Maximum SND Denominator:", maxDenom
    # print "Maximum SND Numerator:", (1+maxNumer)**2

    if PlotData:
        x, y, z = [], [], []
        for xs in DataToWrite:
            x.append(xs[0])
            y.append(xs[1])
            z.append(xs[2])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x, y, z, edgecolor="k", c=z, marker="o", s=30)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x, y, z, edgecolor="k", c=z, marker="o")
        ax.xaxis.set_ticks(np.arange(min(x), max(x), (max(x) - min(x)) / 4))
        ax.yaxis.set_ticks(np.arange(min(y), max(y), (max(y) - min(y)) / 5))
        ax.zaxis.set_ticks(np.arange(min(z), max(z), (max(z) - min(z)) / 5))
        ax.set_xlabel("delta", labelpad=20, fontsize=20)
        ax.set_ylabel("tau", labelpad=20, fontsize=20)
        ax.set_zlabel("Speed of sound", labelpad=10, fontsize=20)

        plt.show()


def writeExp(textFile, Combination=False):
    coeffs = BasisFunctions.coeffs

    i = 1
    for d, t, l, m in coeffs:
        textFile.write("d('%d') = %i ;\n " % (i, d))
        textFile.write("t('%d') = %.6f ;\n " % (i, t))
        textFile.write("l('%d') = %i ;\n " % (i, l))
        i = i + 1
        # textFile.write("m('%d') = %.6f ;\n "%(i,d))


# def writeExpFreeT(textFile, Combination=False):
#     i = 1
#     # d_exp = range(0, 6)
#     # l_exp = range(0, 4)
#     for d in range(0, 6):
#         for l in range(0, 4):
#             t = 1
#             textFile.write("d('%d') = %i ;\n " % (i, d))
#             textFile.write("t('%d') = %.6f ;\n " % (i, t))
#             textFile.write("l('%d') = %i ;\n " % (i, l))
#             i = i + 1
#         # textFile.write("m('%d') = %.6f ;\n "%(i,d))


# @staticmethod
def VapVLEdt(textFile, DataToWrite, Combination=False):
    """ Imports P-V-T data into the GAMS document, to be written in :func:'GamsWrite'

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """
    # print "Vap: "
    # PVTindexes = DataImport.PVTindexes
    # IsothermIndexes = DataImport.isothermIndex
    # print Combination
    i = 1
    for x in DataToWrite:
        # j = 1
        if Combination:
            textFile.write("z('%s', '%d') = %f ;\n" % ("Vap", i, x[2]))
        else:
            textFile.write("z('%d') = %f ;\n" % (i, x[2]))
        i += 1
    i = 1
    # drd_max = 0
    # drd_min = 100
    for x in DataToWrite:
        if Combination:
            textFile.write("delta('%s','%d') = %f;\n" % ("Vap", i, x[0]))
            textFile.write("tau('%s','%d') = %f;\n" % ("Vap", i, x[1]))
        else:
            textFile.write("delta('%d') = %f;\n" % (i, x[0]))
            textFile.write("tau('%d') = %f;\n" % (i, x[1]))

        i += 1
    i = 1

    # print " PVT DRD Max : ", drd_max;
    # print " PVT DRD Min : ", drd_min;

    # for x in DataToWrite:
    #   j=1
    #   BasisFunctions.d3rd(x[0],x[1]);
    #   # BasisFunctions.d4rd(x[0],x[1]);
    #   # BasisFunctions.d5rd(x[0],x[1]);
    #   d3rd_vals = BasisFunctions.d3rd_vals;
    #   # d4rd_vals = BasisFunctions.d4rd_vals;
    #   for y in d3rd_vals:
    #       val = y/(x[0]**3)
    #       if i <= IsothermIndexes: #if i <= PVTindexes[0]:
    #           if Combination:
    #               textFile.write("isoT('%s', '%d', '%d', '%s') = %f ;\n" %('PVT', i, j, 'd3rd', val))
    #           else:
    #               textFile.write("isoT('%d', '%d', '%s') = %f ;\n" %(i, j, 'd3rd', val))
    #       j+=1
    #   i+=1;
    # Crit(textFile, DataToWrite, Combination)
    # InSat(textFile, DataImport.InSatValues,Combination)

    x, y, z = [], [], []
    for xs in DataToWrite:
        x.append(xs[0])
        y.append(xs[1])
        z.append(xs[2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, edgecolor="k", c="none", marker="o")
    ax.set_xlabel("delta")
    ax.set_ylabel("tau")
    ax.set_zlabel("z")
    plt.show()

    # @staticmethod


def LiqVLEdt(textFile, DataToWrite, Combination=False):
    """ Imports P-V-T data into the GAMS document, to be written in :func:'GamsWrite'

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """
    # print "Liq: "
    # PVTindexes = DataImport.PVTindexes
    # IsothermIndexes = DataImport.isothermIndex
    # print Combination
    i = 1
    for x in DataToWrite:
        # j = 1
        if Combination:
            textFile.write("z('%s', '%d') = %f ;\n" % ("Liq", i, x[2]))
        else:
            textFile.write("z('%d') = %f ;\n" % (i, x[2]))
        i += 1
    i = 1
    # drd_max = 0
    # drd_min = 100
    for x in DataToWrite:
        if Combination:
            textFile.write("delta('%s','%d') = %f;\n" % ("Liq", i, x[0]))
            textFile.write("tau('%s','%d') = %f;\n" % ("Liq", i, x[1]))
        else:
            textFile.write("delta('%d') = %f;\n" % (i, x[0]))
            textFile.write("tau('%d') = %f;\n" % (i, x[1]))

        i += 1
    i = 1

    # print " PVT DRD Max : ", drd_max;
    # print " PVT DRD Min : ", drd_min;

    # for x in DataToWrite:
    #   j=1
    #   BasisFunctions.d3rd(x[0],x[1]);
    #   # BasisFunctions.d4rd(x[0],x[1]);
    #   # BasisFunctions.d5rd(x[0],x[1]);
    #   d3rd_vals = BasisFunctions.d3rd_vals;
    #   # d4rd_vals = BasisFunctions.d4rd_vals;
    #   for y in d3rd_vals:
    #       val = y/(x[0]**3)
    #       if i <= IsothermIndexes: #if i <= PVTindexes[0]:
    #           if Combination:
    #               textFile.write("isoT('%s', '%d', '%d', '%s') = %f ;\n" %('PVT', i, j, 'd3rd', val))
    #           else:
    #               textFile.write("isoT('%d', '%d', '%s') = %f ;\n" %(i, j, 'd3rd', val))
    #       j+=1
    #   i+=1;
    # Crit(textFile, DataToWrite, Combination)
    # InSat(textFile, DataImport.InSatValues,Combination)

    x, y, z = [], [], []
    for xs in DataToWrite:
        x.append(xs[0])
        y.append(xs[1])
        z.append(xs[2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x, y, z, edgecolor="k", c="none", marker="o")
    ax.set_xlabel("delta")
    ax.set_ylabel("tau")
    ax.set_zlabel("z")
    plt.show()
