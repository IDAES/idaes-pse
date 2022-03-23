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

from ..atom import Atom


def isLineAtomRecord(line):
    return line[0:4] == "ATOM"


def readPointsFromPDB(filename):
    Points = []
    with open(filename, "r") as infile:
        for line in infile:
            if isLineAtomRecord(line):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                Points.append(np.array([x, y, z], dtype=float))
    return Points


def readAtomsFromPDB(filename):
    Atoms = []
    with open(filename, "r") as infile:
        for line in infile:
            if isLineAtomRecord(line):
                Atoms.append(Atom((line[12:16]).strip()))
    return Atoms


def readPointsAndAtomsFromPDB(filename):
    return readPointsFromPDB(filename), readAtomsFromPDB(filename)


def writeDesignToPDB(D, filename):
    with open(filename, "w") as outfile:
        for i in range(len(D)):
            if not (D.Contents[i] is None or D.Contents[i] == Atom()):
                outfile.write(
                    "ATOM  {:>5d} {:<4s}{:14}{:>8.3f}{:>8.3f}{:>8.3f}{:26}\n".format(
                        i,
                        D.Contents[i].Symbol,
                        "",
                        D.Canvas.Points[i][0],
                        D.Canvas.Points[i][1],
                        D.Canvas.Points[i][2],
                        "",
                    )
                )
