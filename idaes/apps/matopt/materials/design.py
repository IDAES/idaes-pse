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
import os
from copy import deepcopy

from .atom import Atom
from .canvas import Canvas
from .parsers.PDB import readPointsAndAtomsFromPDB, writeDesignToPDB
from .parsers.XYZ import readPointsAndAtomsFromXYZ, writeDesignToXYZ
from .parsers.CFG import readPointsAndAtomsFromCFG, writeDesignToCFG
from .parsers.POSCAR import readPointsAndAtomsFromPOSCAR, writeDesignToPOSCAR


class Design(object):
    """A class used to represent material designs.

    This class combines a ``Canvas`` objects and a list of contents.
    It assigns an element (possibly None) to each point in the ``Canvas``.
    This generally works for any type of content, but it is intended
    to work with ``Atom`` objects and can be used to generate CFG, PDB, POSCAR,
    and XYZ files.
    """

    # === STANDARD CONSTRUCTOR
    def __init__(self, Canvas_=None, Contents=None):
        """Initialize a Design object from standard data."""
        if Canvas_ is not None and Contents is None:
            Contents = [None] * len(Canvas_)
        elif Canvas_ is not None and isinstance(Contents, Atom):
            Contents = [Contents] * len(Canvas_)
        elif Canvas_ is None:
            Canvas_ = Canvas()
            Contents = []
        self._Canvas = Canvas_
        self._Contents = Contents

    # === CONSTRUCTOR - From PDB
    @classmethod
    def fromPDB(cls, filename, DefaultNN=0):
        """Make a Design by reading from PDB file.

        Args:
            filename(str): PDB file to read from.
            DefaultNN(int, optional): Optional, the default number of nearest
                neighbors to initialize the Canvas with.

        Returns:
            (Design) A new Design.

        """
        Pts, Atoms = readPointsAndAtomsFromPDB(filename)
        return cls(Canvas(Points=Pts, DefaultNN=DefaultNN), Atoms)

    # === CONSTRUCTOR - From XYZ
    @classmethod
    def fromXYZ(cls, filename, DefaultNN=0):
        """Make a Design by reading from XYZ file.

        Args:
            filename(str): XYZ file to read from.
            DefaultNN(int, optional): Optional, the default number of nearest
                neighbors to initialize the Canvas with.

        Returns:
            (Design) A new Design.

        """
        Pts, Atoms = readPointsAndAtomsFromXYZ(filename)
        return cls(Canvas(Points=Pts, DefaultNN=DefaultNN), Atoms)

    # === CONSTRUCTOR - From CFG
    @classmethod
    def fromCFG(cls, filename, DefaultNN=0):
        """Makes a Design by reading from CFG file.

        Args:
            filename(str): CFG file to read from.
            DefaultNN(int, optional): Optional, the default number of nearest
                neighbors to initialize the Canvas with.

        Returns:
            (Design) A new Design.

        """
        Pts, Atoms = readPointsAndAtomsFromCFG(filename)
        return cls(Canvas(Points=Pts, DefaultNN=DefaultNN), Atoms)

    # === CONSTRUCTOR - From POSCAR/CONTCAR
    @classmethod
    def fromPOSCAR(cls, filename, DefaultNN=0):
        """Make a Design by reading from POSCAR file.

        Args:
            filename(str): POSCAR file to read from.
            DefaultNN(int, optional): Optional, the default number of nearest
                neighbors to initialize the Canvas with.

        Returns:
            (Design) A new Design.

        """
        Pts, Atoms = readPointsAndAtomsFromPOSCAR(filename)
        return cls(Canvas(Points=Pts, DefaultNN=DefaultNN), Atoms)

    fromCONTCAR = fromPOSCAR
    """Makes a Design by reading from CONTCAR file.
       See documentation for fromPOSCAR. """

    # === MANIPULATION METHODS
    def setContent(self, i, Elem):
        """Set content at a particular location.

        Args:
            i(int): Index to set.
            Elem(Atom/Any): Element to place in Contents.

        Returns:
            None.

        """
        self._Contents[i] = Elem

    def setContents(self, Elem):
        """Set content across all locations.

        Args:
            Elem(Atom/Any): Element to place in each index of Contents.

        Returns:
            None.

        """
        for i in range(len(self.Contents)):
            self.setContent(i, Elem)

    def transform(self, TransF):
        """Transform the this Design according to a functor.

        Args:
            TransF(TransformFunc): Transformation to apply to Canvas.

        Returns:
            None.

        """
        self._Canvas.transform(TransF)

    def getTransformed(self, TransF):
        """Copy and transform this Design.

        Args:
            TransF(TransformFunc): Transformation to apply to Canvas.

        Returns:
            Design) Deepcopy of Design after transformation is applied.

        """
        result = deepcopy(self)
        result.transform(TransF)
        return result

    def add(self, P, Elem):
        """Add point and element to this Design.

        Args:
            P(numpy.ndarray): Point to add.
            Elem(Atom/Any): Content to add.

        Returns:

        """
        self._Canvas.addLocation(P)
        self._Contents.append(Elem)

    def addOther(self, other, blnAssertNotAlreadyInDesign=True):
        """Add another Design to this one.

        Args:
            other(Design): Design to append.
            blnAssertNotAlreadyInDesign(bool, optional): Optional, flag to enable
                assertion that all locations were new and unique.
                (Default value = True)

        Returns:
            None.

        """
        self._Canvas.addOther(other.Canvas, blnAssertNotAlreadyInDesign)
        self._Contents.extend(other.Contents)

    # === PROPERTY EVALUATION METHODS
    def __len__(self):
        """Get the size of this Canvas for this Design."""
        return len(self.Canvas)

    def __eq__(self, other):
        """Compare strict equality of two Designs."""
        return self.Contents == other.Contents and self.Canvas == other.Canvas

    def isEquivalentTo(self, other, blnPreserveIndexing=False, blnIgnoreVoid=True):
        """Compare equivilancy of two Designs.

        Args:
            blnIgnoreVoid:
            other(Design): other Design to compare against.
            blnPreserveIndexing(bool, optional): Optional, flag to determine if
                index order is considered. (Default value = False)
            blnIgnoreVoid(bool, optional): Optional, flag to determine if void
                (i.e., None or Atom() contents) should be considered.
                For example, if all solid atoms are equivalent, but there
                are void locations that do not match. (Default value = True)

        Returns:
            bool) True if Designs are considered equivalent.

        """
        if not blnIgnoreVoid and len(self) != len(other):
            return False
        if blnPreserveIndexing:
            return self == other
        for i, P in enumerate(self.Canvas.Points):
            if not other.Canvas.hasPoint(P):
                return False
            else:
                j = other.Canvas.getPointIndex(P)
                if self.Contents[i] != other.Contents[j]:
                    return False
        return True

    @property
    def NonVoidCount(self):
        """Count number of contents that are not considered void."""
        return sum(Elem is not None and Elem != Atom() for Elem in self.Contents)

    @property
    def NonVoidElems(self):
        """Get the set of non-void contents."""
        result = set(self.Contents)
        result.discard(None)
        result.discard(Atom())
        return result

    # === BASIC QUERY METHODS
    @property
    def Canvas(self):
        """Get the Canvas for this Design."""
        return self._Canvas

    @property
    def Contents(self):
        """Get the Contents for this Design."""
        return self._Contents

    # === REPORTING METHODS
    def toPDB(self, filename):
        """Write a Design to PDB file.

        Args:
            filename(str): PDB file to write to.

        Returns:
            None.

        """
        writeDesignToPDB(self, filename)

    def toXYZ(self, filename):
        """Write a Design to XYZ file.

        Args:
            filename(str): XYZ file to write to.

        Returns:
            None.

        """
        writeDesignToXYZ(self, filename)

    def toCFG(self, filename, GS=None, BBox=None, AuxPropMap=None, blnGroupByType=True):
        """Write a Design to CFG file.

        Args:
            blnGroupByType:
            AuxPropMap:
            BBox:
            filename(str): CFG file to write to.
            GS(float, optional): Optional, Global scaling to write to file
                (see CFG file format). (Default value = None)
            BBox(Parallelepiped, optional): Optional, Bounding box to write to file
                (see CFG file format). If not provided, calculates a
                rectangular prism 2x the necesary size to encompass points.
                (Default value = None)
            AuxPropMap(dict<tuple<str, optional): Optional, Auxilliary
                property map. Example: {('Energy','eV'):[0.0, 1.0, ... ]}
                (Default value = None)
            blnGroupByType(bool, optional): Optional, flag to group atoms by element.
                (see CFG file format). (Default value = True)

        Returns:
            None.

        """
        writeDesignToCFG(
            self,
            filename,
            GS=GS,
            BBox=BBox,
            AuxPropMap=AuxPropMap,
            blnGroupByType=blnGroupByType,
        )

    def toPOSCAR(
        self,
        filename,
        CommentLine=None,
        GS=None,
        BBox=None,
        Elems=None,
        blnUseDirect=True,
    ):
        """Write a Design to CFG file.

        Args:
            blnUseDirect:
            Elems:
            BBox:
            filename(str): CFG file to write to.
            CommentLine(str, optional): Optional, line to write at top of file. (Default value = None)
            GS(float, optional): Optional, Global scaling to write to file
                (see POSCAR file format). (Default value = None)
            BBox(Parallelepiped, optional): Optional, Bounding box to write to file
                (see POSCAR file format). If not provided, calculates a
                rectangular prism 2x the necesary size to encompass points.
                (Default value = None)
            Elems(list<Atom>, optional): Optional, order of elements to write to file.
                Only important because they are sometimes implicitly defined
                in other VASP files and may need to be in a definite order.
                (Default value = None)
            blnUseDirect(bool, optional): Optional, flag to switch between direct and
                cartesian flag of POSCAR file. (Default value = True)

        Returns:
            None.

        """
        writeDesignToPOSCAR(
            self,
            filename,
            CommentLine=CommentLine,
            GS=GS,
            BBox=BBox,
            Elems=Elems,
            blnUseDirect=blnUseDirect,
        )


def loadFromPDBs(filenames, folder=None):
    """Load a list of Designs from PDB files.

    Args:
        filenames (list<str>): List of files to read.
        folder (str): Optional, folder to prepend to filenames.
            (Default value = None)

    Returns:
        (list<Design>): List of Designs created.

    """
    result = []
    for filename in filenames:
        if folder is not None:
            filename = os.path.join(folder, filename)
        result.append(Design.fromPDB(filename))
    return result


def loadFromXYZs(filenames, folder=None):
    """Load a list of Designs from XYZ files.

    Args:
        filenames (list<str>): List of files to read.
        folder (str): Optional, folder to prepend to filenames.
            (Default value = None)

    Returns:
        (list<Design>): List of Designs created.

    """
    result = []
    for filename in filenames:
        if folder is not None:
            filename = os.path.join(folder, filename)
        result.append(Design.fromXYZ(filename))
    return result


def loadFromCFGs(filenames, folder=None):
    """Load a list of Designs from CFG files.

    Args:
        filenames (list<str>): List of files to read.
        folder (str): Optional, folder to prepend to filenames.
            (Default value = None)

    Returns:
        (list<Design>): List of Designs created.

    """
    result = []
    for filename in filenames:
        if folder is not None:
            filename = os.path.join(folder, filename)
        result.append(Design.fromCFG(filename))
    return result
