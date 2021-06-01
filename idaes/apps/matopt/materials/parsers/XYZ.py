###############################################################################
# ** Copyright Notice **
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
###############################################################################
import numpy as np

from ..atom import Atom


def readPointsFromXYZ(filename):
    Points = []
    with open(filename, 'r') as infile:
        nAtoms = int(infile.readline().split()[0])
        _ = next(infile)  # skip the comment line
        for _ in range(nAtoms):
            line = next(infile).split()
            x = float(line[1])
            y = float(line[2])
            z = float(line[3])
            Points.append(np.array([x, y, z], dtype=float))
    return Points


def readAtomsFromXYZ(filename):
    Atoms = []
    with open(filename, 'r') as infile:
        nAtoms = int(infile.readline().split()[0])
        _ = next(infile)  # skip the comment line
        for _ in range(nAtoms):
            line = next(infile).split()
            Atoms.append(Atom(line[0]))
    return Atoms


def readPointsAndAtomsFromXYZ(filename):
    return readPointsFromXYZ(filename), readAtomsFromXYZ(filename)


def writeDesignToXYZ(D, filename, comment_line=None):
    with open(filename, 'w') as outfile:
        outfile.write('{:d}\n'.format(D.NonVoidCount))  # number of atoms
        outfile.write('{}\n'.format(comment_line if comment_line is not None else ''))
        for i in range(len(D)):
            if not (D.Contents[i] is None or D.Contents[i] == Atom()):
                outfile.write('{} {:.8f} {:.8f} {:.8f}\n'.format(D.Contents[i].Symbol,
                                                                 D.Canvas.Points[i][0],
                                                                 D.Canvas.Points[i][1],
                                                                 D.Canvas.Points[i][2]))
