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
