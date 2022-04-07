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
from copy import deepcopy

from ..util.util import myArrayEq, myPointsEq, ListHasPoint
from .parsers.PDB import readPointsAndAtomsFromPDB
from .parsers.XYZ import readPointsAndAtomsFromXYZ
from .parsers.CFG import readPointsAndAtomsFromCFG
from .geometry import RectPrism


class Canvas(object):
    """A class for combining geometric points and neighbors.

    This class contains a list of Cartesian points coupled with a graph of nodes for sites and arcs
    for bonds. A ``Canvas`` object establishes a mapping from the abstract, mathematical modeling of
    materials as graphs to the geometry of the material lattice. The list of points and neighbor
    connections necessary to create a ``Canvas`` object can be obtained from the combination of
    ``Lattice``, ``Shape``, and ``Tiling`` objects.
    """

    DBL_TOL = 1e-5

    # === STANDARD CONSTRUCTOR
    def __init__(self, Points=None, NeighborhoodIndexes=None, DefaultNN=0):
        if Points is None and NeighborhoodIndexes is None:
            Points = []
            NeighborhoodIndexes = []
        elif Points is None:
            Points = [None] * len(NeighborhoodIndexes)
        elif NeighborhoodIndexes is None:
            NeighborhoodIndexes = [[None] * DefaultNN for _ in range(len(Points))]
        self._Points = Points
        self._NeighborhoodIndexes = NeighborhoodIndexes
        self.__DefaultNN = DefaultNN
        assert self.isConsistentWithDesign()

    # === CONSTRUCTOR - From PDB File
    @classmethod
    def fromPDB(cls, filename, Lat=None, DefaultNN=0):
        """Make Canvas by reading from PDB file.

        Args:
            filename(str): Name of PDB file to read.
            Lat (Lattice, optional): A lattice to define neighbor connections
            DefaultNN(int, optional): Number of neighbors to allocate
                in neighborhood. (Default value = 0)

        Returns:
            Canvas: A new Canvas object.

        """
        Pts, _ = readPointsAndAtomsFromPDB(filename)
        result = cls(Points=Pts, DefaultNN=DefaultNN)
        if Lat is not None:
            result.setNeighborsFromFunc(Lat.getNeighbors)
        return result

    @classmethod
    def fromXYZ(cls, filename, Lat=None, DefaultNN=0):
        """Make Canvas by reading from XYZ file.

        Args:
            filename(str): Name of XYZ file to read.
            Lat (Lattice, optional): A lattice to define neighbor connections
            DefaultNN(int, optional): Number of neighbors to allocate
                in neighborhood. (Default value = 0)

        Returns:
            Canvas: A new Canvas object.

        """
        Pts, _ = readPointsAndAtomsFromXYZ(filename)
        result = cls(Points=Pts, DefaultNN=DefaultNN)
        if Lat is not None:
            result.setNeighborsFromFunc(Lat.getNeighbors)
        return result

    @classmethod
    def fromCFG(cls, filename, Lat=None, DefaultNN=0):
        """Make Canvas by reading from CFG file.

        Args:
            filename(str): Name of CFG file to read.
            Lat (Lattice, optional): A lattice to define neighbor connections
            DefaultNN(int, optional): Number of neighbors to allocate
                in neighborhood. (Default value = 0)

        Returns:
            Canvas: A new Canvas object.

        """
        Pts, _ = readPointsAndAtomsFromCFG(filename)
        result = cls(Points=Pts, DefaultNN=DefaultNN)
        if Lat is not None:
            result.setNeighborsFromFunc(Lat.getNeighbors)
        return result

    # === CONSTRUCTOR - From Lattice and Shape
    @classmethod
    def fromLatticeAndShape(
        cls, Lat, S, Seed=np.array([0, 0, 0], dtype=float), DefaultNN=0
    ):
        """Make Canvas by iterating over Lattice points that fit in Shape.

        This constructor starts from a seed location and repeatedly adds
        neighbors until there are no more neighbors that lie inside the
        provided Shape object. It is potentially very slow and should only
        be considered if other methods are unavailable.
        Prefer Canvas.fromLatticeAndShapeScan.

        Args:
            DefaultNN:
            Lat(Lattice): Lattice to provide getNeighbors function.
            S(Shape): Shape to iterate over.
            Seed(numpy.ndarray, optional): Location to begin adding neighbors from.
                Should be on the Lattice. (Default value =
                np.array([0,0,0],dtype=float))
            DefaultNN(int, optional): Number of neighbors to allocate
                in neighborhood. (Default value = 0)

        Returns:
            Canvas: A new Canvas object.

        """
        result = cls(DefaultNN=DefaultNN)
        Stack = [Seed]
        while len(Stack) > 0:
            P = Stack.pop()
            if not result.hasPoint(P) and P in S:
                PNs = Lat.getNeighbors(P)
                result.addLocation(P, len(PNs))
                # NOTE: At first, we checked if P needed to be added
                #       (i.e., if it was not alreay in the Stack)
                #       but doing so was significantly slower than
                #       just extending the Stack without checks and
                #       just checking if P was already in the result
                Stack.extend(PNs)
        result.setNeighborsFromFunc(Lat.getNeighbors)
        return result

    @classmethod
    def fromLatticeAndShapeScan(cls, Lat, argPolyhedron, DefaultNN=0):
        """Make Canvas by iterating over Lattice points that fit in Shape.

        This constructor takes advantage of methods in Lattice to
        efficiently scan over sites. This requires the Lattice.Scan
        method to produce a generator object.

        Args:
            DefaultNN:
            Lat(Lattice): Lattice with has defined Scan method.
            argPolyhedron(Polyhedron: Polyhedron): Shape to iterate over.
                NOTE: A Polyhedron is required because there is a
                complicated step in finding bounds for the Shape
                in the reference lattice space that is not
                generally valid for all Shapes.
            DefaultNN(int, optional): Number of neighbors to allocate
                in neighborhood. (Default value = 0)

        Returns:
            Canvas: A new Canvas object.

        """
        result = cls(DefaultNN=DefaultNN)
        BBox = RectPrism.fromPointsBBox(argPolyhedron.getBounds())
        for P in Lat.Scan(BBox):
            if P in argPolyhedron:
                result.addLocation(P, len(Lat.getNeighbors(P)))
        result.setNeighborsFromFunc(Lat.getNeighbors)
        return result

    @classmethod
    def fromLatticeAndTiling(
        cls, Lat, T, Seed=np.array([0, 0, 0], dtype=float), DefaultNN=0
    ):
        """Make Canvas by iterating over Lattice points that fit in Tiling.

        See documentation for fromLatticeAndShape.
        This constructor additionally makes the resulting Canvas periodic.

        Args:
            DefaultNN:
            Lat(Lattice): Lattice to provide getNeighbors function.
            T(Tiling): Tiling which provides a Shape to iterate over.
            Seed(numpy.ndarray, optional): Location to begin adding neighbors from.
                Should be on the Lattice. (Default value = np.array([0,0,0],dtype=float))
            DefaultNN(int, optional): Number of neighbors to allocate
                in neighborhood. (Default value = 0)

        Returns:
            Canvas: A new Canvas object.

        """
        result = cls.fromLatticeAndShape(
            Lat, T.TileShape, Seed=Seed, DefaultNN=DefaultNN
        )
        result.makePeriodic(T, Lat.getNeighbors)
        return result

    @classmethod
    def fromLatticeAndTilingScan(cls, Lat, T, DefaultNN=0):
        """Make Canvas by iterating over Lattice points that fit in Tiling.

        See documentation for fromLatticeAndTiling.
        This constructor additionally makes the resulting Canvas periodic.

        Args:
            Lat(Lattice): Lattice with has defined Scan method.
            T(Tiling): Tiling that provides a Polyhedron shape.
            DefaultNN(int, optional): Number of neighbors to allocate
                in neighborhood. (Default value = 0)

        Returns:
            Canvas: A new Canvas object.

        """
        result = cls.fromLatticeAndShapeScan(Lat, T.TileShape, DefaultNN=DefaultNN)
        result.makePeriodic(T, Lat.getNeighbors)
        return result

    # === ASSERTION OF CLASS DESIGN
    def isConsistentWithDesign(self):
        """Determine if object is consistent with class assumptions."""
        if len(self.Points) != len(self.NeighborhoodIndexes):
            return False
        return True

    # === MANIPULATION METHODS
    def addLocation(self, P, NNeighbors=None):
        """Add new location to Points.

        Args:
            P(numpy.ndarray): Point to add.
            NNeighbors(int, optional): Number of neighbors to allocate
                in a neighborhood. If None, use the instance default
                self.DefaultNN (Default value = None)

        Returns:
            None.

        """
        assert not self.hasPoint(P)
        self._Points.append(P)
        self._NeighborhoodIndexes.append([None] * (NNeighbors or self.__DefaultNN))
        assert self.isConsistentWithDesign()

    def setNeighbors(self, P1, P2, l=None):
        """Set a (directed) neighbor connection between two points.

        Args:
            P1(numpy.ndarray): Canvas Point to set neighbor for.
            P2(numpy.ndarray): Canvas Point to set as neighbor.
            l(int, optional): Index of neighbor to set. For example, if
                l=3, then P2 is set to the fourth neighbor of P1.
                If None, appends to the neighborhood.
                (Default value = None)

        Returns:
            None.

        """
        assert self.hasPoint(P1)
        assert self.hasPoint(P2)
        i = self.getPointIndex(P1)
        j = self.getPointIndex(P2)
        self.setNeighborsIJ(i, j, l=l)
        assert self.isConsistentWithDesign()

    def setNeighborsIJ(self, i, j, l=None):
        """Set a (directed) neighbor connection between two indices.

        Args:
            i(int): Index of location to set neighbor for.
            j(int): Index of location to set as neighbor.
            l(int, optional): Index of neighbor to set. For example, if
                l=3, then j is set to the fourth neighbor of i.
                If None, appends to the neighborhood.
                (Default value = None)

        Returns:
            None.

        """
        if l is None:
            self._NeighborhoodIndexes[i].append(j)
        else:
            self._NeighborhoodIndexes[i][l] = j
        assert self.isConsistentWithDesign()

    def setNeighborLofI(self, PN, l, i, blnSetNoneOtherwise=True):
        """Set a (directed) neighbor connection between a Point and index.

        Args:
            blnSetNoneOtherwise:
            i:
            PN(numpy.ndarray): Canvas Point to set as neighbor.
            l(int): Index of neighbor to set. For example, if
                l=3, then PN is set as the fourth neighbor of i.
            i(int): Index of location to set neighbor for.
            blnSetNoneOtherwise(bool, optional): Flag to control behavior
                if PN is not found in Canvas. (Default value = True)

        Returns:
            None.

        """
        assert i < len(self._NeighborhoodIndexes)
        assert l < len(self._NeighborhoodIndexes[i])
        if self.hasPoint(PN):
            self._NeighborhoodIndexes[i][l] = self.getPointIndex(PN)
        elif blnSetNoneOtherwise:
            self._NeighborhoodIndexes[i][l] = None
        assert self.isConsistentWithDesign()

    def setNeighborsOfI(self, PNs, i):
        """Set a list of points as neighbors to an index.

        Args:
            PNs(list<numpy.ndarray>): Points to set as neighbors of i.
            i(int): Index to set neighbors for.

        Returns:
            None.

        """
        self._NeighborhoodIndexes[i] = [None] * len(PNs)
        for l, P in enumerate(PNs):
            self.setNeighborLofI(P, l, i)
        assert self.isConsistentWithDesign()

    def setNeighborsFromFunc(self, NeighborsFunc):
        """Set neighbors across the Canvas from a functor.

        Args:
            NeighborsFunc(function): Function that takes as input a
                point (numpy.ndarray) and returns a list of points.

        Returns:
            None.

        """
        for i, P in enumerate(self.Points):
            PNs = NeighborsFunc(P)
            self.setNeighborsOfI(PNs, i)
        assert self.isConsistentWithDesign()

    def getNeighborsFromFunc(self, NeighborsFunc, layer=1):
        """Get neighbors across the Canvas from a functor.

        Args:
            NeighborsFunc(function): Function that takes as input a point
                (numpy.ndarray) and returns a list of points.
            layer(int): Target layer of neighbors.

        Returns:
            list<list<int>> Indexes matrix

        """
        result = [[] for _ in range(len(self.Points))]
        for i, P in enumerate(self.Points):
            PNs = NeighborsFunc(P, layer)
            result[i] = [None] * len(PNs)
            for j, PN in enumerate(PNs):
                if self.hasPoint(PN):
                    result[i][j] = self.getPointIndex(PN)
        return result

    def getNeighborsFromFuncAndTiling(self, NeighborsFunc, argTiling, layer=1):
        """Get neighbors across the Canvas from a functor and Tiling.

        Args:
            NeighborsFunc(function): Function that takes as input a point
                (numpy.ndarray) and returns a list of points.
            argTiling(Tiling): Tiling object that specifies the periodicity
                of the Canvas.
            layer(int): Target layer of neighbors.

        Returns:
            list<list<int>> Indexes matrix
        """
        result = [[] for _ in range(len(self.Points))]
        for i, P in enumerate(self.Points):
            PNs = NeighborsFunc(P, layer)
            result[i] = [None] * len(PNs)
            for j, PN in enumerate(PNs):
                if self.hasPoint(PN):
                    result[i][j] = self.getPointIndex(PN)
                else:
                    for TilingDirection in argTiling.TilingDirections:
                        PtoTry = PN + TilingDirection
                        if self.hasPoint(PtoTry):
                            result[i][j] = self.getPointIndex(PtoTry)
        return result

    def makePeriodic(self, argTiling, NeighborsFunc):
        """Make connections periodic accross the edges of Tiling.

        Args:
            argTiling(Tiling: Tiling): Tiling to provide TilingDirections.
            NeighborsFunc(function): Function that takes as input a
                point (numpy.ndarray) and returns a list of points.

        Returns:
            None.

        """
        for i, P in enumerate(self.Points):
            LatNeighbors = NeighborsFunc(P)
            for l, Index in enumerate(self.NeighborhoodIndexes[i]):
                if Index is None:
                    assert not self.hasPoint(
                        LatNeighbors[l]
                    )  # else, Canvas constructed incorrectly
                    for TilingDirection in argTiling.TilingDirections:
                        PtoTry = LatNeighbors[l] + TilingDirection
                        if self.hasPoint(PtoTry):
                            self._NeighborhoodIndexes[i][l] = self.getPointIndex(PtoTry)
                            break

    def addShells(self, n, NeighborsFunc):
        """Add locations in n-shells around current Points.

        Args:
            n(int): Number of shells to add.
            NeighborsFunc(function): Function that takes as input a
                point (numpy.ndarray) and returns a list of points.

        Returns:
            None.

        """
        for _ in range(n):
            self.addShell(NeighborsFunc)

    def addShell(self, NeighborsFunc):
        """Add locations in a shell around current Points.

        Args:
            NeighborsFunc(function): Function that takes as input a
                point (numpy.ndarray) and returns a list of points.

        Returns:
            None.

        """
        Shell = self.getShell(NeighborsFunc)
        for P in Shell:
            self.addLocation(P)
        self.setNeighborsFromFunc(NeighborsFunc)

    def transform(self, TransF):
        """Transform the Points in Canvas accroding to a functor.

        Args:
            TransF(TransformFunc): Transformation to apply to Points.

        Returns:
            None.

        """
        for P in self._Points:
            TransF.transform(P)

    def getTransformed(self, TransF):
        """Copy and transform this Canvas.

        Args:
            TransF(TransformFunc): Transformation to apply to Points.

        Returns:
            None.

        """
        result = deepcopy(self)
        result.transform(TransF)
        return result

    def addOther(self, other, blnAssertNotAlreadyInCanvas=True):
        """Add other Canvas to this one.

        Args:
            other(Canvas): Canvas to append.
            blnAssertNotAlreadyInCanvas(bool, optional): Flag to enable
                assertion that all locations were new and unique.
                (Default value = True)

        Returns:
            None.

        """
        for P in other.Points:
            if blnAssertNotAlreadyInCanvas:
                assert not self.hasPoint(P)
            self.addLocation(P)

    # === PROPERTY EVALUATION METHODS
    def __len__(self):
        """Get the number of Points for this Canvas."""
        return len(self.Points)

    def __eq__(self, other):
        """Compare strict equality of Canvas data."""
        return (
            myPointsEq(self.Points, other.Points, Canvas.DBL_TOL)
            and self.NeighborhoodIndexes == other.NeighborhoodIndexes
        )

    def hasPoint(self, P):
        """Identify if point is in Canvas.

        Args:
            P(numpy.ndarray): Point to identify in Canvas.

        Returns:
            bool) True if Points has P.

        """
        for Q in self.Points:
            # NOTE: There are several ways to test point membership.
            #       This is optimized for speed.
            if myArrayEq(P, Q, Canvas.DBL_TOL):
                return True
        return False

    def getPointIndex(self, P):
        """Identify the index of a point in the Canvas.

        Args:
            P(numpy.ndarray): Point in the Canvas.

        Returns:
            int) Index of P in Points.

        """
        for i, Q in enumerate(self.Points):
            # NOTE: There are several ways to test point membership.
            #       This is optimized for speed.
            if myArrayEq(P, Q, Canvas.DBL_TOL):
                return i
        return None

    def getNeighbors(self, P):
        """Identify set of neighbors to a point in Canvas.

        Args:
            P(numpy.ndarray): Point to get neighbors for.

        Returns:
            list<int>: Neighborhood of P.

        """
        return self.NeighborhoodIndexes[self.getPointIndex(P)]

    def getNeighborLofI(self, l, i):
        """Identify neighbor for specific index and neighbor order.

        Args:
            l(int): Order of neighbor in neighborhood
            i(int): Index to consider the neighborhood of.

        Returns:
            int: Index for neighbor l of i.

        """
        assert i < len(self.NeighborhoodIndexes)
        assert l < len(self.NeighborhoodIndexes[i])
        return self.NeighborhoodIndexes[i][l]

    def getShell(self, NeighborsFunc):
        """Identify set of neighboring points not in Canvas.

        Args:
            NeighborsFunc(function): Function that takes as input a
                point (numpy.ndarray) and returns a list of points.

        Returns:
            list<numpy.ndarray>: Set of points to consider as
            neighboring shell.

        """
        result = []
        for i, P in enumerate(self.Points):
            Neighs = NeighborsFunc(P)
            for Neigh in Neighs:
                if not self.hasPoint(Neigh) and not ListHasPoint(
                    result, Neigh, Canvas.DBL_TOL
                ):
                    result.append(Neigh)
        return result

    def getNeighborhoodIndexes(self, Lat, layer=1, T=None):
        """Wrapper of functions that returns neighbors across the Canvas

        Args:
            Lat(Lattice): Lattice of the Canvas
            layer(int): optional, target layer of neighbors
            T(Tiling): optional, Tiling object that specifies the periodicity

        Returns:
            list<list<int>> Indexes matrix
        """
        if layer == 1:
            return self.NeighborhoodIndexes
        elif T is None:
            return self.getNeighborsFromFunc(Lat.getNeighbors, layer)
        else:
            return self.getNeighborsFromFuncAndTiling(Lat.getNeighbors, T, layer)

    # === BASIC QUERY METHODS
    @property
    def Points(self):
        """Get Points in Canvas."""
        return self._Points

    @property
    def NeighborhoodIndexes(self):
        """Get description of neighborhoods in Canvas."""
        return self._NeighborhoodIndexes

    # === REPORTING METHODS
    def printPoints(self):
        """Pretty-print Points."""
        for i, P in enumerate(self.Points):
            print("{}: {}".format(i, P))

    def printNeighborhoodIndexes(self, layer=1):
        """Pretty-print NeighborhoodIndexes."""
        for i, Neighborhood in enumerate(self.getNeighborhoodIndexes(layer)):
            print("{}: {}".format(i, Neighborhood))
