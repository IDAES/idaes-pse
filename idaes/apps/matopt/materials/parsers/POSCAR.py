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

from ..geometry import Parallelepiped, RectPrism
from ..atom import Atom


def readPointsAndAtomsFromPOSCAR(filename, ImpliedElems=None):
    Pts = []
    Atoms = []
    BBox = Parallelepiped.fromPOSCAR(filename)
    blnAtomToAdd = False  # Used to track
    Elem = None  # Used to store most recent element in the new format
    with open(filename, "r") as infile:
        CommentLine = infile.readline()
        GSLine = infile.readline().split()
        GS = float(GSLine[0])
        VxLine = infile.readline()
        VyLine = infile.readline()
        VzLine = infile.readline()
        ElementNamesOrCountsLine = infile.readline().split()
        if ElementNamesOrCountsLine[0].isdigit():
            if ImpliedElems is not None:
                Elems = ImpliedElems
                ElemCountsLine = ElementNamesOrCountsLine
            else:
                raise ValueError(
                    "Tried to read a POSCAR without explicit elements provided"
                )
        else:
            assert ElementNamesOrCountsLine[0].isalpha()
            ElemsLine = ElementNamesOrCountsLine
            Elems = [Atom(Elem) for Elem in ElemsLine]
            ElemCountsLine = infile.readline().split()
        ElemCounts = [int(ElemCount) for ElemCount in ElemCountsLine]
        CartesianOrDirectLine = infile.readline().split()
        CorDFlag = CartesianOrDirectLine[0][0].lower()
        blnIsCartesian = CorDFlag == "c" or CorDFlag == "k"
        for t, Elem in enumerate(Elems):
            for i in range(ElemCounts[t]):
                XYZLine = infile.readline().split()
                if blnIsCartesian:
                    P = GS * np.array(
                        [float(XYZLine[0]), float(XYZLine[1]), float(XYZLine[2])],
                        dtype=float,
                    )
                else:
                    # NOTE: The global scaling GS has already been taken in the BBox
                    #       so no need to multiply it here
                    Pfrac = np.array(
                        [float(XYZLine[0]), float(XYZLine[1]), float(XYZLine[2])],
                        dtype=float,
                    )
                    P = np.inner(np.array([BBox.Vx, BBox.Vy, BBox.Vz]).T, Pfrac)
                Pts.append(P)
                Atoms.append(Elem)
    return Pts, Atoms


def writeDesignToPOSCAR(
    D, filename, CommentLine=None, GS=None, BBox=None, Elems=None, blnUseDirect=True
):
    if GS is None:
        GS = 1.0
    if BBox is None:
        BBox = RectPrism.fromPointsBBox(D.Canvas.Points)
        BBox.scale(2.0)
    if Elems is None:
        Elems = D.NonVoidElems
    with open(filename, "w") as outfile:
        outfile.write("{}\n".format(CommentLine if CommentLine is not None else ""))
        outfile.write("{}\n".format(GS))
        outfile.write(
            "{} {} {}\n".format(BBox.Vx[0] / GS, BBox.Vx[1] / GS, BBox.Vx[2] / GS)
        )
        outfile.write(
            "{} {} {}\n".format(BBox.Vy[0] / GS, BBox.Vy[1] / GS, BBox.Vy[2] / GS)
        )
        outfile.write(
            "{} {} {}\n".format(BBox.Vz[0] / GS, BBox.Vz[1] / GS, BBox.Vz[2] / GS)
        )
        outfile.write("{}\n".format(" ".join([Elem.Symbol for Elem in Elems])))
        outfile.write(
            "{}\n".format(" ".join([str(D.Contents.count(Elem)) for Elem in Elems]))
        )
        outfile.write("{}\n".format("Direct" if blnUseDirect else "Cartesian"))
        for Elem in Elems:
            for i in range(len(D)):
                P = (
                    D.Canvas.Points[i]
                    if blnUseDirect
                    else BBox.getFractionalCoords(D.Canvas.Points[i])
                )
                if D.Contents[i] == Elem:
                    outfile.write("{} {} {}\n".format(P[0] / GS, P[1] / GS, P[2] / GS))
