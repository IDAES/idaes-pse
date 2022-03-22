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
import numpy as np

from ..geometry import RectPrism
from ..atom import Atom


def readPointsAndAtomsFromCFG(filename):
    Pts = []
    Atoms = []
    GS = 1.0
    Vx = np.array([0, 0, 0], dtype=float)
    Vy = np.array([0, 0, 0], dtype=float)
    Vz = np.array([0, 0, 0], dtype=float)
    blnAtomToAdd = False  # Used to track
    Elem = None  # Used to store most recent element in the new format
    with open(filename, "r") as infile:
        for line in infile:
            splitLine = line.split()
            if len(splitLine) == 0:
                continue
            blnAtomToAdd = False
            s1 = s2 = s3 = None
            Tagname = splitLine[0]
            if Tagname[0] == "#":
                continue
            elif Tagname == "A":
                GS = GS * float(splitLine[2])
            elif Tagname == "H0(1,1)":
                Vx[0] = GS * float(splitLine[2])
            elif Tagname == "H0(1,2)":
                Vx[1] = GS * float(splitLine[2])
            elif Tagname == "H0(1,3)":
                Vx[2] = GS * float(splitLine[2])
            elif Tagname == "H0(2,1)":
                Vy[0] = GS * float(splitLine[2])
            elif Tagname == "H0(2,2)":
                Vy[1] = GS * float(splitLine[2])
            elif Tagname == "H0(2,3)":
                Vy[2] = GS * float(splitLine[2])
            elif Tagname == "H0(3,1)":
                Vz[0] = GS * float(splitLine[2])
            elif Tagname == "H0(3,2)":
                Vz[1] = GS * float(splitLine[2])
            elif Tagname == "H0(3,3)":
                Vz[2] = GS * float(splitLine[2])
            elif Tagname[0].isdigit() and len(splitLine) == 1:
                # Start of new format, line for mol weight
                line = next(infile)
                splitLine = line.split()
                Elem = Atom(splitLine[0])
            elif Tagname[0].isdigit() and splitLine[1][0].isdigit():
                # New format atom
                blnAtomToAdd = True
                s1 = float(splitLine[0])
                s2 = float(splitLine[1])
                s3 = float(splitLine[2])
            elif Tagname[0].isdigit() and splitLine[1][0].isalpha():
                # Old format atom
                blnAtomToAdd = True
                Elem = Atom(splitLine[1])
                s1 = float(splitLine[2])
                s2 = float(splitLine[3])
                s3 = float(splitLine[4])
            else:
                # Other entry, not usefule
                pass
            if blnAtomToAdd:
                assert Elem is not None
                x = s1 * Vx[0] + s2 * Vy[0] + s3 * Vz[0]
                y = s1 * Vx[1] + s2 * Vy[1] + s3 * Vz[1]
                z = s1 * Vx[2] + s2 * Vy[2] + s3 * Vz[2]
                Pts.append(np.array([x, y, z], dtype=float))
                Atoms.append(Elem)
    return Pts, Atoms


def writeDesignToCFG(
    D, filename, GS=None, BBox=None, AuxPropMap=None, blnGroupByType=True
):
    if BBox is None:
        BBox = RectPrism.fromPointsBBox(D.Canvas.Points)
        BBox.scale(2.0)
    with open(filename, "w") as outfile:
        outfile.write("Number of particles = {}\n".format(D.NonVoidCount))
        outfile.write(
            "A = {} Angstrom (basic length-scale)\n".format(1.0 if GS is None else GS)
        )
        outfile.write("H0(1,1) = {} A\n".format(BBox.Vx[0]))
        outfile.write("H0(1,2) = {} A\n".format(BBox.Vx[1]))
        outfile.write("H0(1,3) = {} A\n".format(BBox.Vx[2]))
        outfile.write("H0(2,1) = {} A\n".format(BBox.Vy[0]))
        outfile.write("H0(2,2) = {} A\n".format(BBox.Vy[1]))
        outfile.write("H0(2,3) = {} A\n".format(BBox.Vy[2]))
        outfile.write("H0(3,1) = {} A\n".format(BBox.Vz[0]))
        outfile.write("H0(3,2) = {} A\n".format(BBox.Vz[1]))
        outfile.write("H0(3,3) = {} A\n".format(BBox.Vz[2]))
        outfile.write(".NO_VELOCITY.\n")
        if AuxPropMap is not None:
            outfile.write("entry_count = {}\n".format(len(AuxPropMap) + 3))
            for i, AuxProp in enumerate(AuxPropMap):
                PropName = AuxProp[0]
                PropUnit = AuxProp[1]
                outfile.write("auxiliary[{}] = {} [{}]".format(i, PropName, PropUnit))
            outfile.write("\n")
        else:
            outfile.write("entry_count = 3\n")

        if blnGroupByType:
            for Elem in D.NonVoidElems:
                outfile.write("{}\n".format(Elem.Mass))
                outfile.write("{}\n".format(Elem.Symbol))
                for i in range(len(D)):
                    if D.Contents[i] == Elem:
                        P = BBox.getFractionalCoords(D.Canvas.Points[i])
                        outfile.write("{} {} {}".format(P[0], P[1], P[2]))
                        if AuxPropMap is not None:
                            for AuxProp in AuxPropMap:
                                if i in AuxPropMap[AuxProp]:
                                    outfile.write(" {}".format(AuxPropMap[AuxProp][i]))
                        outfile.write("\n")
        else:
            for i in range(len(D)):
                if not (D.Contents[i] is None or D.Contents[i] == Atom()):
                    P = BBox.getFractionalCoords(D.Canvas.Points[i])
                    Elem = D.Contents[i]
                    outfile.write(
                        "{} {} {} {} {}".format(
                            Elem.Mass, Elem.Symbol, P[0], P[1], P[2]
                        )
                    )
                    if AuxPropMap is not None:
                        for AuxProp in AuxPropMap:
                            if i in AuxPropMap[AuxProp]:
                                outfile.write(" {}".format(AuxPropMap[AuxProp][i]))
                    outfile.write("\n")
