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
Writer of the data into the GAMS file
"""

import numpy as np
import matplotlib.pyplot as plt
from . import BasisFunctions, DataImport

Combination = False
R = 8.314472


def writeExp(textFile, Combination=False):
    """
    Writes into the GDX file the basis function parameters
    """
    coeffs = BasisFunctions.coeffs

    i = 1
    for d, t, l, m in coeffs:
        textFile.write("d('%d') = %i ;\n " % (i, d))
        textFile.write("t('%d') = %.6f ;\n " % (i, t))
        textFile.write("l('%d') = %i ;\n " % (i, l))
        i = i + 1


def PVTdt(textFile, DataToWrite, Combination=False, PlotData=False):
    """Imports P-V-T data into the GAMS document, to be written in :func:'GamsWrite'

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """

    IsothermIndexes = DataImport.isothermIndex

    i = 1
    for x in DataToWrite:
        if Combination:
            textFile.write("z('%s', '%d') = %f ;\n" % ("PVT", i, x[2]))
        else:
            textFile.write("z('%d') = %f ;\n" % (i, x[2]))
        i += 1
    i = 1

    for x in DataToWrite:
        if Combination:
            textFile.write("delta('%s','%d') = %f;\n" % ("PVT", i, x[0]))
            textFile.write("tau('%s','%d') = %f;\n" % ("PVT", i, x[1]))
        else:
            textFile.write("delta('%d') = %f;\n" % (i, x[0]))
            textFile.write("tau('%d') = %f;\n" % (i, x[1]))

        i += 1
    i = 1

    InSat(textFile, DataImport.InSatValues, Combination)

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


def CVdt(textFile, DataToWrite, Combination=False, PlotData=False):
    """Imports Isochoric Heat Capacity (CV) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.
    """

    i = 1

    for x in DataToWrite:
        BasisFunctions.iTT(x[0], x[1])
        itt_val = -BasisFunctions.itt_val
        val = x[2]
        if Combination:
            textFile.write("z('%s', '%d') = %f ;\n" % ("CV", i, val))
            textFile.write("itt('%s', '%d') = %f ;\n" % ("CV", i, itt_val))
            textFile.write("delta('%s', '%d') = %.13f ;\n" % ("CV", i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % ("CV", i, x[1]))
        else:
            textFile.write("z('%d') = %f ;\n" % (i, val))
        i += 1

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


def CPdt(textFile, DataToWrite, Combination=False, PlotData=False):
    """Imports Isobaric Heat Capacity(CP) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    """

    i = 1
    for x in DataToWrite:
        BasisFunctions.iTT(x[0], x[1])
        itt_val = -BasisFunctions.itt_val
        val = x[2] - itt_val
        if Combination:
            textFile.write("z('%s', '%d') = %.13f ;\n" % ("CP", i, val))
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


def SNDdt(textFile, DataToWrite, Combination=False, PlotData=False):
    """Imports Speed of Sound (SND) data into the GAMS document

    :param textFile: Gams File written to.
    :type textFile: str.
    :param DataToWrite: Data prepared for the GAMS file.
    :type DataToWrite: array.
    :returns: void.

    """
    i = 1
    for x in DataToWrite:
        if Combination:
            textFile.write("z('%s', '%d') = %.13f ;\n" % ("SND", i, x[2]))
        else:
            textFile.write("z('%d') = %.13f ;\n" % (i, x[2]))
        i += 1
    i = 1
    for x in DataToWrite:
        BasisFunctions.iTT(x[0], x[1])
        itt_val = BasisFunctions.itt_val
        BasisFunctions.drd(x[0], x[1])
        if Combination:
            textFile.write("itt('%s','%d') = %.13f ;\n" % ("SND", i, itt_val))
            textFile.write("delta('%s', '%d') = %.13f ;\n" % ("SND", i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % ("SND", i, x[1]))
        else:
            pass  # required by pylint so that the block directive applies
            # PYLINT-TODO-FIX use correct number of args in format strings
            # pylint: disable=too-few-format-args
            textFile.write("itt('%s','%d') = %.13f ;\n" % (i, itt_val))
            textFile.write("delta('%s', '%d') = %.13f ;\n" % (i, x[0]))
            textFile.write("tau('%s','%d') = %.13f ;\n" % (i, x[1]))
        i += 1

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


def Crit(textFile, DataToWrite, Combination=False):
    """Import Critical data points"""
    BasisFunctions.drd(1, 1)
    drd_vals = BasisFunctions.drd_vals

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
    """Import saturation density values"""
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
